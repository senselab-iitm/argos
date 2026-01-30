<h1>Argos: Leveraging Visual Priors for Scalable Wireless Navigation in Dynamic Environments</h1>

![banner](./readme_assets/artifact_workflow.png)

This artefact demonstrates the anchor selection for a material-aware scalable digital twin of a factory. This repository provides a segmented digital twin of a factory, where users can define the prior path of a movable tag, and the number of anchors required to localize it throughout its journey. Our algorithm provides the optimal subset of anchors and the FMD (which affects the ranging error) map of the deployment. This repository has three folders:
<ol>
	<li>data: collected test data from a factory environment. A material aware 3D model of the factory is provided in the subdirectory visual_priors, covering 14 m × 9 m. UWB Channel Impulse Responses (CIRs), recorded from this environment using QORVO DECAWAVE DW-1000 UWB tags, is provided in the subdirectory rf_priors.</li>
	<li>scripts: scripts for generating ∆FMD caches and optimizers for anchor placement. We have also provided other visualization helpers  and utilities in this directory for visualizing resultant 3D models.</li>
	<li>apps: directory contains user-facing scripts for: (a)generating the ∆FMD cache given an RF-aware 3D model and prior UWB transmitter grid configurations, and (b)finding the optimal anchor subset given an anchor budget and the ∆FMD cache. It also contains helper scripts for the ease-of-use for the user.</li>
</ol>

<h2>Installation and setup</h2>

The code has been developed and tested with **Python 3.10** on **Ubuntu 22.04 LTS**; using the same setup is recommended for reproducibility. The current implementation is compatible **only with Sionna-RT v1.0.1** and is not expected to work with other versions. Installation instructions for Sionna-RT v1.0.1 are available in the [official documentation](https://github.com/NVlabs/sionna/blob/v1.0.1/README.md), which is recommended to follow in order to avoid dependency issues.

### Sionna-RT Installation

> Sionna RT has the same requirements as [Mitsuba 3](https://github.com/mitsuba-renderer/mitsuba3) and we refer to its
> [installation guide](https://mitsuba.readthedocs.io/en/stable/) for further
> information. To run Sionna RT on CPU, [LLVM](https://llvm.org) is required by
> [Dr.Jit](https://drjit.readthedocs.io/en/stable/). Please check the
> [installation instructions for the LLVM backend](https://drjit.readthedocs.io/en/latest/what.html#backends).
> The source code of Sionna RT is located in a separate
> [GitHub repository](https://github.com/NVlabs/sionna-rt).

_Source: [Sionna-RT Official README](https://github.com/NVlabs/sionna/blob/v1.0.1/README.md)_

### Setup Instructions

1. Inside the project root, create a virtual environment  
   ```bash
   python3 -m venv .venv
   ```
2. Enable this virtual environment<br>
```bash
source .venv/bin/activate
```

3. In this virtual environment install necessary packages as specified in `requirements.txt`<br>
```bash
pip3 install -r requirements.txt
```

5. Then make the virtual environment accessible to jupyter lab<br>
```bash
python -m ipykernel install --user --name=.venv
```

7. From this virtual environment, start jupyter lab<br>
```bash
jupyter lab
```

8. Select the .venv kernel in jupyter lab by clicking the `kernel->change kernel menu`
