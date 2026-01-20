"""
Miscellaneous Utilities
Authors: Arko Datta and Ayon Chakraborty, SeNSE Lab, IIT Madras

This module has miscellaneous utilities for rendering the outputs of the user-facing applications
"""

import numpy as np
from sionna.rt import RadioMap
import mitsuba as mi
import drjit as dr
import matplotlib.pyplot as plt
import matplotlib as mpl

def make_custom_radiomap(scene, rss_matrix_dBm, center, size):
    """
    Create a RadioMap from a custom matrix of RSS values in dBm.
    Works even if no transmitters exist in the scene.
    """
    m, n = rss_matrix_dBm.shape
    cell_size = (size[0]/n, size[1]/m)

    rm = RadioMap(
        scene=scene,
        center=mi.Point3f(*center),
        orientation=mi.Point3f((np.pi)/2, 0.0, 0.0),
        size=mi.Point2f(*size),
        cell_size=mi.Point2f(*cell_size)
    )

    # If scene has no TX, fake one with power=1W
    if dr.width(rm._tx_powers) == 0:
        rm._tx_powers = mi.Float([1.0])   # one TX, 1 Watt
        rm._tx_positions = mi.Point3f([0.0], [0.0], [0.0])

    # Convert dBm -> Watts
    rss_W = 10**((rss_matrix_dBm - 30)/10.0)

    # Use first TX power
    tx_power_W = float(rm._tx_powers[0])
    path_gain = rss_W / tx_power_W

    # Inject into pathgain map
    rm._pathgain_map = mi.TensorXf(path_gain[None, :, :])

    return rm

def render_custom_radiomap(custom_map, scene, center, size, invert=False):
    if invert==True:
        for i in range(custom_map.shape[0]):
            for j in range(custom_map.shape[1]):
                custom_map[i,j] = np.max(custom_map)-custom_map[i,j]
    rss_matrix_dBm = custom_map * (-40) + (1 - custom_map) * (-100)
    rss_matrix_dBm = rss_matrix_dBm/np.max(rss_matrix_dBm)
    rm = make_custom_radiomap(
        scene,
        rss_matrix_dBm=rss_matrix_dBm,
        center=center,
        size=size
    )
    return rm

def render_colorbar(values, title="", invert=False):
    colormap = "viridis" if invert==False else "viridis_r"
    sm = mpl.cm.ScalarMappable(
        cmap=colormap,
        norm=mpl.colors.Normalize(vmin=np.min(values), vmax=np.max(values))
    )
    sm.set_array([])
    
    fig, ax = plt.subplots(figsize=(12, 1))
    ax.set_title(title)
    plt.colorbar(sm, cax=ax, orientation="horizontal")
    plt.show()