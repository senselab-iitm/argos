#Specify the paths to jointly collected visual + RF data
cfg_mesh_directory="../data/visual_priors/factory/"
cfg_cir_path="../data/rf_priors/rx_samples.csv" #Please refer to the provided sheet for examples
cfg_prior_map = "../data/prior_map.npy" #The prior map for mobile node(s) expressed as a 2D array spanning 0-1
cfg_fmd_cache_name = "../data/fmd_cache.json" #The generated fmd cache

#In meters, if .blend file is available. Can be obtained from Blender.
cfg_scene_length = 14
cfg_scene_breadth = 9

#Transmitter grid specifications
cfg_tx_grid_length = 4
cfg_tx_grid_breadth = 3

#Other aesthetic options
cfg_bounding_box_name="elm__28"
cfg_set_bounding_box_transparent=True
cfg_floor_name = "elm__6"
cfg_set_floor_transparent=True
cfg_clip_at = 1
cfg_clip_plane_orientation = (0,-1.9,0)
cfg_color_dict = {
    "metal":(0.60, 0.65, 0.70),
    "wood":(0.76, 0.60, 0.42),
    "marble":(0.98, 0.95, 0.75),
    "concrete":(0.75, 0.75, 0.75)
}