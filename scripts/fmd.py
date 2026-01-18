"""
Compute and cache FMD maps using Sionna RT + real CIR data.
Supports partial Rx updates via a bitmask.

Cache is a json with the format:
{
  "(tx_x, tx_y, tx_z)": {
    "(rx_x, rx_y, rx_z)": { "fmd": ..., "rssi": ... },
    ...
  },
  ...
}
"""

import time, gc, json
from typing import Tuple, Optional, Sequence, Dict, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

from sionna.rt import load_scene, PlanarArray, Receiver, Transmitter, PathSolver

#Constants
C = 3e8
EPS = 1e-18
K_DEFAULT = 100

#Helpers
def add_tx_rx(scene, receiver_positions, transmitter_positions):
    scene.rx_array = PlanarArray(
        num_rows=1, num_cols=1,
        vertical_spacing=0.5, horizontal_spacing=0.5,
        pattern="dipole", polarization="V"
    )
    rx_nodes = []
    for i, pos in enumerate(receiver_positions):
        rx = Receiver(name=f"node_{i}", position=pos, display_radius=0.1)
        scene.add(rx)
        rx_nodes.append(rx)

    scene.tx_array = PlanarArray(
        num_rows=1, num_cols=1,
        vertical_spacing=0.5, horizontal_spacing=0.5,
        pattern="tr38901", polarization="V"
    )
    tx_nodes = []
    for i, pos in enumerate(transmitter_positions):
        tx = Transmitter(name=f"tx_{i}", position=pos, display_radius=0.5, look_at=[0,0,0])
        scene.add(tx)
        tx_nodes.append(tx)
    return rx_nodes, tx_nodes


def calculate_paths(scene):
    p_solver = PathSolver()
    paths_refraction = p_solver(scene=scene,
                                max_depth=10, los=True,
                                specular_reflection=False,
                                diffuse_reflection=False,
                                refraction=True,
                                synthetic_array=False, seed=41)
    paths_reflection = p_solver(scene=scene,
                                max_depth=10, los=True,
                                specular_reflection=True,
                                diffuse_reflection=True,
                                refraction=False,
                                synthetic_array=False, seed=41)
    return paths_reflection, paths_refraction


def calculate_taps(paths_obj, bandwidth=500e6, sampling_frequency=1e9,
                   l_min=0, l_max=99, normalize_delays=False):
    a_taps = paths_obj.taps(
        bandwidth=bandwidth, l_min=l_min, l_max=l_max,
        sampling_frequency=sampling_frequency,
        normalize_delays=normalize_delays, out_type="numpy"
    )
    a_taps_s = np.squeeze(a_taps)
    if a_taps_s.ndim == 1:
        a_taps_s = a_taps_s[np.newaxis, :]
    taps = np.arange(l_min, l_max + 1) / sampling_frequency
    taps_ns = taps * 1e9
    return a_taps_s, taps_ns

def build_sparse_with_neighbors_iq(df, frac=0.1, radius=0.2, z_fixed=1.2, K=100):
    i_cols = [f"I_{i}" for i in range(K)]
    q_cols = [f"Q_{i}" for i in range(K)]
    sampled = df.sample(frac=frac, random_state=42)

    xs_all = df['x_position'].to_numpy(float)
    ys_all = df['y_positions'].to_numpy(float)
    I_all = np.nan_to_num(df[i_cols].to_numpy(float))
    Q_all = np.nan_to_num(df[q_cols].to_numpy(float))

    results = []
    for _, center in sampled.iterrows():
        cx, cy = float(center['x_position']), float(center['y_positions'])
        dx = xs_all - cx
        dy = ys_all - cy
        dist = np.sqrt(dx ** 2 + dy ** 2)
        idx = np.where(dist <= radius)[0]
        if len(idx) == 0:
            continue
        cirs_complex = I_all[idx, :] + 1j * Q_all[idx, :]
        neigh_pos = np.stack([ys_all[idx], np.full(len(idx), z_fixed), xs_all[idx]], axis=1)
        pos = (cy, z_fixed, cx)  # (y,z,x)
        results.append((pos, cirs_complex, neigh_pos))
    return results

def generate_grid_receivers(xlim, ylim, z=1.5, spacing=0.1):
    xs = np.arange(xlim[0], xlim[1], spacing) + spacing/2
    ys = np.arange(ylim[0], ylim[1], spacing) + spacing/2
    return np.array([(x, y, z) for x in xs for y in ys])

def lookup_rx_power(coord, receivers_df):
    xq, yq = coord[0], coord[1]
    xs = receivers_df["x_position"].to_numpy(float)
    ys = receivers_df["y_positions"].to_numpy(float)
    dists = np.sqrt((xs - xq) ** 2 + (ys - yq) ** 2)
    return float(receivers_df["RX_Power_dBm"].iat[int(np.argmin(dists))])

def compute_real_mu_sigma(sparse_iq, receivers_df, l0=25, l1=None):
    coords, mu_all, var_all = [], [], []
    for entry in sparse_iq:
        coord = np.array(entry[0], float)  # (cy,z,cx)
        cir = np.array(entry[1])
        if cir.ndim == 1: cir = cir[np.newaxis, :]
        roi = cir[:, l0:] if l1 is None else cir[:, l0:l1]
        rx_dbm = lookup_rx_power((coord[2], coord[0]), receivers_df)
        rx_mw = 10 ** (rx_dbm / 10)
        E = np.sum(np.abs(roi) ** 2, axis=1, keepdims=True)
        scale = rx_mw / np.maximum(E, 1e-12)
        roi_mw = np.abs(roi) ** 2 * scale
        roi_dbm = 10 * np.log10(np.maximum(roi_mw, 1e-18))
        mu_all.append(np.mean(roi_dbm, axis=0))
        var_all.append(np.var(roi_dbm, axis=0))
        coords.append(np.array([coord[2], coord[0], coord[1]], float))  # (x,y,z)
    return np.array(coords), np.array(mu_all), np.array(var_all)

def inverse_distance_weights(target, sources, power=2, eps=1e-9):
    dists = np.linalg.norm(sources - target, axis=1) + eps
    w = 1.0 / (dists ** power)
    return w / np.sum(w)

def fuse_mu_sigma(target_coords, sim_dbm, real_coords, mu_real, var_real, power=2):
    N, T = sim_dbm.shape
    mu_out, var_out = np.zeros((N, T)), np.zeros((N, T))
    for i, coord in enumerate(target_coords):
        w = inverse_distance_weights(coord, real_coords, power=power)
        mu_r = np.sum(w[:, None] * mu_real, axis=0)
        var_r = np.sum(w[:, None] * var_real, axis=0)
        mu_out[i] = 0.5 * sim_dbm[i] + 0.5 * mu_r
        var_out[i] = var_r
    return mu_out, var_out

def compute_fmd(mu_taps, sigma2_taps, tx_coords, rx_coords, Ts, gamma_dbm=-90):
    N, T = mu_taps.shape
    fmd = np.zeros(N)
    tx = np.array(tx_coords, float)
    for n in range(N):
        mu, sigma = mu_taps[n], np.sqrt(np.maximum(sigma2_taps[n], 1e-12))
        q = norm.sf((gamma_dbm - mu) / sigma)
        probs, prod = np.zeros(T), 1.0
        for j in range(T):
            probs[j] = prod * q[j]
            prod *= (1 - q[j])
        E_X = np.sum(np.arange(T) * probs)
        d = np.linalg.norm(tx - rx_coords[n])
        fmd[n] = d / (C * Ts) - E_X
    return fmd

def plot_fmd(rx_coords, fmd_vals, Ts, unit="meters"):
    x, y = rx_coords[:, 0], rx_coords[:, 1]
    z = fmd_vals * (C * Ts) if unit == "meters" else fmd_vals
    label = "FMD (m)" if unit == "meters" else "FMD (taps)"
    plt.figure(figsize=(8, 6))
    sc = plt.scatter(x, y, c=z, cmap="coolwarm", s=80, edgecolors="k")
    plt.colorbar(sc, label=label)
    plt.xlabel("x [m]"); plt.ylabel("y [m]")
    plt.title("FMD map")
    plt.axis("equal"); plt.grid(True)
    plt.show()

class FMDCalculator:
    def __init__(self, arena_size, spacing, real_csv_path, scene_path,
                 rx_height=1.5, roi_l0=25, roi_l1=None,
                 sampling_frequency=1e9, sim_frequency=6.5e9):
        t0 = time.time()
        self.xlim, self.ylim = (0, arena_size[0]), (0, arena_size[1])
        self.spacing, self.rx_height = spacing, rx_height
        self.real_csv_path, self.scene_path = real_csv_path, scene_path
        self.roi_l0, self.roi_l1 = roi_l0, roi_l1
        self.fs, self.Ts, self.sim_frequency = sampling_frequency, 1/sampling_frequency, sim_frequency

        self.receivers = pd.read_csv(real_csv_path)
        self.rx_positions = generate_grid_receivers(self.xlim, self.ylim, z=rx_height, spacing=spacing)
        self.grid_shape = (
            len(np.arange(*self.ylim, spacing)),
            len(np.arange(*self.xlim, spacing))
        )
        self._cache: Dict[str, Dict[str, Any]] = {}
        self.setup_time = time.time() - t0
        self.scene = None

    def update_fmd(self, tx_coords, gamma_dbm=-90,
                   frac_sample=0.1, radius_sample=0.2, K=K_DEFAULT,
                   bitmask: Optional[np.ndarray] = None):
        """
        Compute/update FMD for a given Tx.
        If bitmask is provided (shape=grid_shape), only Rx where mask=1 are updated.
        """
        timers = {}
        t0 = time.time()
        tx_tuple, tx_key = tuple(map(float, tx_coords)), str(tuple(map(float, tx_coords)))

        # -------------------------------
        # Select receivers
        # -------------------------------
        t_sel = time.time()
        if bitmask is None:
            sel_rx_positions = self.rx_positions
        else:
            ny, nx = self.grid_shape
            if bitmask.shape != (ny, nx):
                raise ValueError(f"bitmask must have shape {self.grid_shape}, got {bitmask.shape}")
            sel_rx_positions = []
            for rx in self.rx_positions:
                x, y = rx[0], rx[1]
                x_idx = int(round((x - self.xlim[0] - self.spacing/2) / self.spacing))
                y_idx = int(round((y - self.ylim[0] - self.spacing/2) / self.spacing))
                if 0 <= x_idx < nx and 0 <= y_idx < ny and bitmask[y_idx, x_idx] == 1:
                    sel_rx_positions.append(rx)
            sel_rx_positions = np.array(sel_rx_positions)
        timers["receiver_selection"] = time.time() - t_sel
        if len(sel_rx_positions) == 0:
            raise ValueError("No receivers selected (check bitmask)")

        # Scene + paths
        t1 = time.time()
        self.scene = load_scene(self.scene_path)
        self.scene.frequency = self.sim_frequency
        add_tx_rx(self.scene, sel_rx_positions, [tx_tuple])
        paths_reflection, paths_refraction = calculate_paths(self.scene)
        timers["scene_and_paths"] = time.time() - t1

        # Tap extraction
        t2 = time.time()
        a_tap, _ = calculate_taps(paths_reflection, sampling_frequency=self.fs)
        a_tap_refract, _ = calculate_taps(paths_refraction, sampling_frequency=self.fs)
        a_comb = np.zeros_like(a_tap, complex)
        both = (a_tap != 0) & (a_tap_refract != 0)
        a_comb[both] = 0.5 * (a_tap[both] + a_tap_refract[both])
        a_comb[(a_tap != 0) & ~both] = a_tap[(a_tap != 0) & ~both]
        a_comb[(a_tap_refract != 0) & ~both] = a_tap_refract[(a_tap_refract != 0) & ~both]
        timers["tap_processing"] = time.time() - t2

        # Sparse real stats
        t3 = time.time()
        sparse_iq = build_sparse_with_neighbors_iq(
            self.receivers, frac=frac_sample, radius=radius_sample,
            z_fixed=self.rx_height, K=K
        )
        real_coords, mu_real, var_real = compute_real_mu_sigma(
            sparse_iq, self.receivers, l0=self.roi_l0, l1=self.roi_l1
        )
        timers["real_data_processing"] = time.time() - t3

        # Sim → dBm
        t4 = time.time()
        roi_len = mu_real.shape[1] if mu_real.size > 0 else a_comb.shape[1]
        sim_cut = np.array([row[:roi_len] for row in a_comb])
        rx_dbm_sim = np.array([lookup_rx_power(rp[:2], self.receivers) for rp in sel_rx_positions])
        rx_mw_sim = 10 ** (rx_dbm_sim / 10)
        E_sim = np.sum(np.abs(sim_cut) ** 2, axis=1)
        scale = rx_mw_sim / np.maximum(E_sim, EPS)
        per_tap_mw_sim = (np.abs(sim_cut) ** 2) * scale[:, None]
        per_tap_dbm_sim = 10 * np.log10(np.maximum(per_tap_mw_sim, EPS))
        timers["sim_to_dbm"] = time.time() - t4

        # Fusion
        t5 = time.time()
        if real_coords.shape[0] > 0:
            mu_taps, sigma2_taps = fuse_mu_sigma(sel_rx_positions, per_tap_dbm_sim,
                                                 real_coords, mu_real, var_real)
        else:
            mu_taps, sigma2_taps = per_tap_dbm_sim, np.full_like(per_tap_dbm_sim, 1e-6)
        timers["fusion"] = time.time() - t5

        # FMD
        t6 = time.time()
        fmd_vals = compute_fmd(mu_taps, sigma2_taps, tx_tuple, sel_rx_positions, self.Ts, gamma_dbm)
        timers["fmd_computation"] = time.time() - t6

        # Cache update
        t7 = time.time()
        if tx_key not in self._cache:
            self._cache[tx_key] = {}
        for rx, fmd_val, rssi_val in zip(sel_rx_positions, fmd_vals, rx_dbm_sim):
            rx_key = str(tuple(np.round(rx, 6)))
            self._cache[tx_key][rx_key] = {"fmd": float(fmd_val), "rssi": float(rssi_val)}
        timers["cache_update"] = time.time() - t7

        timers["total"] = time.time() - t0
        gc.collect()
        return fmd_vals, timers

    def get_fmd_map(self, tx_coords, unit="meters"):
        """Return FMD values as a 2D numpy array [ny, nx]."""
        tx_key = str(tuple(map(float, tx_coords)))
        if tx_key not in self._cache:
            raise ValueError("Tx not found in cache, run update_fmd first")

        fmd_map = np.full(self.grid_shape, np.nan)
        for rx_key, vals in self._cache[tx_key].items():
            rx = np.array(eval(rx_key))
            x_idx = int(round((rx[0] - self.xlim[0] - self.spacing/2) / self.spacing))
            y_idx = int(round((rx[1] - self.ylim[0] - self.spacing/2) / self.spacing))
            fmd_val = vals["fmd"]
            if unit == "meters":
                fmd_val *= (C * self.Ts)
            fmd_map[y_idx, x_idx] = fmd_val
        return fmd_map

    @property
    def cache(self): return self._cache

    def export_cache_to_json(self, path, pretty=True):
        with open(path, "w") as f:
            json.dump(self._cache, f, indent=2 if pretty else None)

    def import_cache_from_json(self, path):
        with open(path, "r") as f:
            data = json.load(f)
        for tx_key, rx_dict in data.items():
            if tx_key not in self._cache:
                self._cache[tx_key] = {}
            self._cache[tx_key].update(rx_dict)

# -------------------------------
# Example usage
# -------------------------------
if __name__ == "__main__":
    np.random.seed(0)

    fmd = FMDCalculator(
        arena_size=(15, 9), spacing=0.5,
        real_csv_path="../data/raw/rx_samples.csv",
        scene_path="../data/raw/robolab_sionna_materials/roboscene_with_materials.xml"
    )
    print(f"Initial setup time: {fmd.setup_time:.3f} s")

    tx1 = (7, 1.2, 3)
    ny, nx = fmd.grid_shape

    for frac in [0.2, 0.5, 0.8]:
        mask = (np.random.rand(ny, nx) < frac).astype(int)
        fmd_vals, timers = fmd.update_fmd(tx1, bitmask=mask)
        print(f"\n--- Tx {tx1}, mask={int(frac*100)}% ---")
        for k, v in timers.items():
            print(f"{k:20s}: {v:.3f} s")

    # Heatmap
    fmd_map = fmd.get_fmd_map(tx1, unit="meters")
    plt.figure(figsize=(8,6))
    plt.imshow(fmd_map, origin="lower", cmap="plasma",
               extent=[fmd.xlim[0], fmd.xlim[1], fmd.ylim[0], fmd.ylim[1]])
    plt.colorbar(label="FMD [m]")
    plt.xlabel("x [m]"); plt.ylabel("y [m]")
    plt.title("FMD Heatmap (plasma)")
    plt.show()