<h1>Argos: Leveraging Visual Priors for Scalable Wireless Navigation in Dynamic Environments</h1>
This artefact demonstrates the anchor selection for a material-aware scalable digital twin of a factory. This repository provides a segmented digital twin of a factory, where users can define the prior path of a movable tag, and the number of anchors required to localize it throughout its journey. Our algorithm provides the optimal subset of anchors and the FMD (which affects the ranging error) map of the deployment. This repository has three folders:
<ol>
	<li>data: collected test data from a factory environment. A material aware 3D model of the factory is provided in the subdirectory visual_priors, covering 14 m × 9 m. UWB Channel Impulse Responses (CIRs), recorded from this environment using QORVO DECAWAVE DW-1000 UWB tags, is provided in the subdirectory rf_priors.</li>
	<li>scripts: scripts for generating ∆FMD caches and optimizers for anchor placement. We have also provided other visualization helpers  and utilities in this directory for visualizing resultant 3D models.</li>
	<li>apps: directory contains user-facing scripts for: (a)generating the ∆FMD cache given an RF-aware 3D model and prior UWB transmitter grid configurations, and (b)finding the optimal anchor subset given an anchor budget and the ∆FMD cache. It also contains helper scripts for the ease-of-use for the user.</li>
</ol>

<h2>Installation</h2>
Kindly install the required packages from requirements.txt. Our current version works with Sionna-RT 1.0.1, as is specified and will not work with other versions.<br>
1. Inside the project root, create a virtual environment<br>

`python3 -m venv .venv` <br>
2. Enable this virtual environment<br>
`source .venv/bin/activate`<br>
3. In this virtual environment install necessary packages<br>
`pip3 install -r requirements.txt`<br>
4. Then make the virtual environment accessible to jupyter lab<br>
`ipykernel install --user name=.venv`<br>
5. From this virtual environment, start jupyter lab<br>
`jupyter lab`<br>
Select the .venv kernel in jupyter lab by clicking the `kernel->change kernel menu`
