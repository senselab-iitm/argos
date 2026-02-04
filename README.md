<h1>Argos: Leveraging Visual Priors for Scalable Wireless Navigation in Dynamic Environments</h1>

![banner](./readme_assets/artifact_workflow.png)

This artefact demonstrates the anchor selection for a material-aware scalable digital twin of a factory. This repository provides a segmented digital twin of a factory, where users can define the prior path of a movable tag, and the number of anchors required to localize it throughout its journey. Our algorithm provides the optimal subset of anchors and the FMD (which affects the ranging error) map of the deployment. This repository has three folders:
<ol>
	<li>data: collected test data from a factory environment. A material aware 3D model of the factory is provided in the subdirectory visual_priors, covering 14 m × 9 m. UWB Channel Impulse Responses (CIRs), recorded from this environment using QORVO DECAWAVE DW-1000 UWB tags, is provided in the subdirectory rf_priors.</li>
	<li>scripts: scripts for generating ∆FMD caches and optimizers for anchor placement. We have also provided other visualization helpers  and utilities in this directory for visualizing resultant 3D models.</li>
	<li>apps: directory contains user-facing scripts for: (a)generating the ∆FMD cache given an RF-aware 3D model and prior UWB transmitter grid configurations, and (b)finding the optimal anchor subset given an anchor budget and the ∆FMD cache. It also contains helper scripts for the ease-of-use for the user.</li>
</ol>

<h2>Installation and setup</h2>

The code supports **Python 3.10 - 3.12** and has been developed and tested with **Python 3.10** on **Ubuntu 22.04 LTS**; using the same setup is recommended for reproducibility. The current implementation is compatible **only with Sionna-RT v1.0.1** and is not expected to work with other versions.

### Setup Instructions

1. Install python3.10, python3-venv and llvm (for running Sionna-RT in CPU)
   ```bash
   sudo apt install python3.10
   sudo apt install python3-venv
   sudo apt install llvm
   ```
2. Install git and clone the repository
   ```bash
   sudo apt install git
   git clone https://github.com/senselab-iitm/argos.git
   ```
3. Change directory to project root
   ```bash
   cd argos
   ```   
4. Inside the project root, create a virtual environment  
   ```bash
   python3 -m venv .venv
   ```
5. Enable this virtual environment<br>
	```bash
	source .venv/bin/activate
	```
6. In this virtual environment install necessary packages as specified in `requirements.txt`<br>
	```bash
	pip3 install -r requirements.txt
	```
7. Then make the virtual environment accessible to jupyter lab<br>
	```bash
	python -m ipykernel install --user --name=.venv
	```
8. From this virtual environment, start jupyter lab<br>
	```bash
	jupyter lab
	```

8. Select the .venv kernel in jupyter lab by clicking the `kernel->change kernel menu`. If the `kernel` menu appears disabled/unclickable, open an ipynb notebook from the `app` directory

<h1>Application Specific Notebooks for Argos</h1>
The `app` folder has 4 tools in form of IPYTHON NOTEBOOKS to facilitate clear demonstration.
<ol>
<li>create_fmd_map.ipynb: provides a step-by-step method to initialize an ∆FMD cache<br>
You can define the grid of anchors from which the selection is to be done.<br>

![anchor_grid](readme_assets/anchor_grid.gif)
<br>
You can also view ∆FMD for each transmitter after it has been computed<br>

![fmd_map](readme_assets/single_anchor_fmd.gif)
<br>

</li><br>

<li>select_anchors.ipynb: helps the user to perform anchor selection given the ∆FMD cache. It takes in a prior trajectory probability map of mobile nodes and shows the anchors to be deployed in order to achieve superior quality localization given a budget.<br>

![fmd_map](readme_assets/anchors_selected.gif)
<br>

</li><br>

<li>iter_objects.ipynb: a utility which allows the user to iterate through the objects in the scene and view their identifiers for debugging purposes.

![fmd_map](readme_assets/item_iter.gif)
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

### User Configurations
Users are provided with configurations in the form of a configs.py file so that user-centric configurations are not tightly coupled with raw code.