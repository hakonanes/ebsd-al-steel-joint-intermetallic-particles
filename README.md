# Analysis of electron backscatter diffraction (EBSD) patterns of five phases from an Al-steel joint 

This repository contains Jupyter notebooks and a MATLAB script necessary to reproduce the EBSD results in the conference paper "" which was recently submitted to X.

The content in this repository is licensed under the GPLv3+, since many of the softwares used have this license.

## Contents

The Jupyter notebook and Python files are numbered according to the steps taken in the EBSD analysis:
1. `ebsd1_preprocess.ipynb`: Increase the signal-to-noise ratio of patterns by background subtraction and averaging. Generate indexing-independent views of EBSD datasets (mean intensity map, image quality map, and average neighbour dot product map), and calibrate the detector-sample geometry via projection center (PC) optimization with the [PyEBSDIndex Python package](https://github.com/USNavalResearchLaboratory/PyEBSDIndex) (cubic phases only!). An average PC is used in dictionary indexing.
2. `ebsd2_dictionary_indexing.py`: Obtain crystal orientations from the EBSD patterns via dictionary indexing (DI) as implemented in `kikuchipy`. Requires master patterns of each phase, generated with EMsoft.
3. `ebsd3_orientation_refinement.py`: Refine crystal orientations obtained from DI.
4. `ebsd4_postprocess_indexing_results.ipynb`: Create multi-phase crystal map from the single-phase maps obtained from DI.
5. `ebsd5_plot_simulations.ipynb`: Plot geometrical and dynamical simulations for each phase.

Python packages used in the notebooks and scripts are listed in `requirements.txt` and can be installed into a virtual or conda environment:

```bash
pip install -r requirements.txt
```

Analysis of the EBSD indexing results is done with the script `orientation_analysis.m`. MATLAB packages used are [MTEX](https://mtex-toolbox.github.io/) and [export_fig](https://se.mathworks.com/matlabcentral/fileexchange/23629-export_fig).