"""
Argos Optimizer
Authors: Arko Datta and Ayon Chakraborty, SeNSE Lab, IIT Madras

This helper module takes in 
(i) the scene and its ∆FMD cache, 
(ii) candidate anchor positions, 
(iii) the anchor budget, and 
(iv) a prior location map of the object(s) to localize, 
and returns the anchor subset fitting the budget, which reduces the localization error. 
The module runs the optimization algorithm, incrementally selecting anchors to minimise
the mean ∆FMD weighted by the prior and the geometric dilution of precision (GDoP). 
It returns the recommended anchor positions along with the map of localization errors.
"""
import numpy as np
class ArgosOptimization:
    def __init__(self,
                 fmd_calc,
                 candidate_txs,
                 prior_map,
                 rssi_threshold=-90.0,
                 min_anchors=3):
        self.fmd_calc = fmd_calc #fmd calculator instance -> has the fmd maps

        # keep both world and internal coords
        self.candidate_txs_world = [tuple(map(float, tx)) for tx in candidate_txs] # (x,z,y)
        self.candidate_txs = [(tx[0], tx[1], tx[2]) for tx in self.candidate_txs_world]  # (x,y,z)

        self.prior_map = np.array(prior_map, dtype=float)
        self.rssi_threshold = float(rssi_threshold)
        self.min_anchors = int(min_anchors)
        # grid sizes
        self.ny, self.nx = self.prior_map.shape
        self.Nrx = self.ny * self.nx #number of receivers
        self.Nc = len(self.candidate_txs) #number of candidate transmitters

        # flatten prior
        self.prior_flat = self.prior_map.flatten().astype(float)
        self.total_prior = np.sum(self.prior_flat)

        # RX coordinates (grid centers)
        xs = np.linspace(self.fmd_calc.xlim[0] + 0.5*self.fmd_calc.spacing,
                         self.fmd_calc.xlim[1] - 0.5*self.fmd_calc.spacing, self.nx)
        ys = np.linspace(self.fmd_calc.ylim[0] + 0.5*self.fmd_calc.spacing,
                         self.fmd_calc.ylim[1] - 0.5*self.fmd_calc.spacing, self.ny)
        xx, yy = np.meshgrid(xs, ys)
        z = getattr(self.fmd_calc, 'rx_height', 0.0)
        self.rx_coords = np.column_stack([xx.ravel(), yy.ravel(), np.full(xx.size, z)])

        # Precompute FMD and RSSI [Nc, Nrx]
        self.fmd_flat = np.full((self.Nc, self.Nrx), np.nan, dtype=float)
        self.rssi_flat = np.full((self.Nc, self.Nrx), -np.inf, dtype=float)
        for ci, tx in enumerate(self.candidate_txs):
            tx_key = str(tuple(map(float, tx)))
            rx_dict = self.fmd_calc.cache.get(tx_key, {})
            for ri in range(self.Nrx):
                rx = tuple(np.round(self.rx_coords[ri], 6))
                rx_key = str(rx)
                if rx_key in rx_dict:
                    vals = rx_dict[rx_key]
                    self.fmd_flat[ci, ri] = float(vals.get("fmd", np.nan))
                    self.rssi_flat[ci, ri] = float(vals.get("rssi", -np.inf))

    def _score_candidate(self, c_idx, kappa, selected):
        mask = self.rssi_flat[c_idx] >= self.rssi_threshold #Indicates the cells which have ben covered by the selected anchor
        if not np.any(mask):
            print("UNCOVERED")
            return np.inf, dict(
                score=np.inf,
                avg_fmd=np.inf,
                delta_prior_cov=0.0,
                cover_total=0,
                cover_active=0,
                penalty=0.0,
            )

        before_trilat = (kappa >= self.min_anchors) #Indicate the cells which have been covered by minimum anchors before selection
        after_trilat  = (kappa + mask.astype(int) >= self.min_anchors) #Indicates the cells which have been covered by minimum anchors after selection

        #Gets the sum of prior weights being serviced and covered by minimum number of anchors. Our goal is to ensure that the parts of the prior map having higher weight must be covered.
        prior_before = np.sum(self.prior_flat[before_trilat]) #Before being covered by selected anchor
        prior_after  = np.sum(self.prior_flat[after_trilat]) #After being covered by selected anchor
        delta_prior_cov = (prior_after - prior_before) / (self.total_prior + 1e-12) #Gets the weight of cells gained after being covered
        
        # FMD quality for whole map
        #fmd_vals = self.fmd_flat[c_idx, mask]
        fmd_vals_mask = (self.fmd_flat[c_idx, mask] * self.prior_flat[mask]) > 0
        fmd_vals = (self.fmd_flat[c_idx, mask] * self.prior_flat[mask])[fmd_vals_mask]

        #fmd_vals = np.clip(fmd_vals * self.range_scale, self.sigma_floor, None)
        #avg_fmd = np.average(fmd_vals, weights=(self.prior_flat[mask])) if np.sum(self.prior_flat[mask]) > 0 else 1e9 #Avg fmd is set to infinity if the anchor cannot reach that cell.
        avg_fmd = np.percentile(fmd_vals[np.isfinite(fmd_vals)], 75)
        avg_fmd = 1/avg_fmd
        
        # Multilateration penalty
        under_mask = (kappa < self.min_anchors) & mask #Cells that are not yet covered
        multi_pen = np.sum(self.prior_flat[under_mask]) / (self.total_prior + 1e-12)

        #Collinearity penalty: collinear anchors incur higher GDOP
        pts = np.array(self.candidate_txs_world)[:, [0, 2]]  # (x,z) only
        col_pen = 0
        p = pts[c_idx]
        
        for i in range(len(selected)):
            for j in range(i+1, len(selected)):
                a = pts[selected[i]]
                b = pts[selected[j]]
                area2 = abs((b[0]-a[0])*(p[1]-a[1]) - (b[1]-a[1])*(p[0]-a[0]))
                col_pen += (area2 < 1e-6)

        metrics = dict(
            avg_fmd=avg_fmd,
            delta_prior_cov=(1 - delta_prior_cov),
            cover_total=int(mask.sum()),
            cover_active=int(((kappa < self.min_anchors) & mask).sum()),
            penalty=multi_pen,
            collinearity_penalty=col_pen
        )
        return metrics

    def run(self, budget):
        kappa = np.zeros(self.Nrx, dtype=int)
        selected = []
        history = []

        for step in range(1, budget + 1):
            cand_metrics = {}
            for c in range(self.Nc):
                if c in selected: 
                    continue
                metrics = self._score_candidate(c, kappa, selected)
                if metrics is not None:
                    cand_metrics[c] = metrics
        
            if not cand_metrics:
                break
        
            keys = ["delta_prior_cov", "avg_fmd", "penalty","collinearity_penalty"]
            norm = {}
        
            for k in keys:
                vals = np.array([cand_metrics[c][k] for c in cand_metrics])
                lo, hi = np.nanmin(vals), np.nanmax(vals)
                norm[k] = {c: (cand_metrics[c][k] - lo) / (hi - lo + 1e-12)
                           for c in cand_metrics}
                
            scores = {
                c: norm["delta_prior_cov"][c]
                   + norm["avg_fmd"][c]
                   + norm["penalty"][c]
                   + norm["collinearity_penalty"][c]
                for c in cand_metrics
            }

            best_idx = min(scores, key=scores.get)
            best_metrics = cand_metrics[best_idx]
        
            selected.append(best_idx)
            cover_mask = (self.rssi_flat[best_idx] >= self.rssi_threshold)
            kappa = kappa + cover_mask.astype(int)
        
            history.append(dict(
                step=step,
                chosen_idx=best_idx,
                metrics=best_metrics,
                score=scores[best_idx]
            ))
        selected_tx = [self.candidate_txs_world[i] for i in selected]
        return dict(
            selected=selected_tx,
            selected_indices=selected,
            history=history,
            fmd_map = self.compute_median_fmd_map(selected_tx)
        )

    def compute_median_fmd_map(self, selected_tx):
        return np.median(np.stack([self.fmd_calc.get_fmd_map(txcoord) for txcoord in selected_tx]),axis=0)