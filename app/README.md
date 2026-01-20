<h1>Application Specific Notebooks for Argos</h1>
The app folder has 4 tools in form of IPYTHON NOTEBOOKS to facilitate clear demonstration.
<ol>
<li>create_fmd_map.ipynb: provides a step-by-step method to initialize an ∆FMD cache<br>
You can define the grid of anchors from which the selection is to be done.<br>

![anchor_grid](../readme_assets/anchor_grid.gif)
<br>
You can also view ∆FMD for each transmitter after it has been computed<br>

![fmd_map](../readme_assets/single_anchor_fmd.gif)
<br>

</li><br>

<li>argos_optimize.ipynb: helps the user to perform anchor selection given the ∆FMD cache. It takes in a prior trajectory probability map of mobile nodes and shows the anchors to be deployed in order to achieve superior quality localization given a budget.<br>

![fmd_map](../readme_assets/single_anchor_fmd.gif)
<br>

</li><br>

<li>iter_objects.ipynb: a utility which allows the user to iterate through the objects in the scene and view their identifiers for debugging purposes.

![fmd_map](../readme_assets/item_iter.gif)
<br>

</li><br>

<li>benchmark_updates.ipynb: benchmarks percentage of updates for the hardware Argos is executing on. Generates a benchmark.csv file in the app folder by updating FMD cache for single Tx with a certain percentage receivers. The meaning of different columns are as follows:
	<ol>
		<li>receiver_selection: Time spent selecting receiver grid points based on the random bitmask.</li><br>
		<li>scene_and_paths: Time spent loading the scene and computing propagation paths.</li><br>
		<li>tap_processing: Time spent extracting and combining multipath channel taps.</li><br>
		<li>real_data_processing: Time spent sampling receivers and computing real-data statistics (mean and variance).</li><br>
		<li>sim_to_dbm: Time spent converting simulated channel energy into received power in dBm.</li><br>
		<li>fusion: Time spent fusing simulated and real signal statistics.</li><br>
		<li>fmd_computation: Time spent computing the Fisher Metric Distance values.</li><br>
		<li>cache_update: Time spent storing computed FMD and RSSI values in the cache.</li><br>
		<li>total: Total time for the full FMD update pipeline.	</li><br>
	</ol>
</ol> 
Users are provided with configurations in the form of a configs.py file so that user-centric configurations are not tightly coupled with
raw code.
