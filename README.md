# Risk analysis of natural hazards to power grids in Southeast and East Asia
This repository provides the code to perform the risk analysis of tropical cyclones and coastal floods to power grids in Southeast and East Asia on an asset level.

It also provides Jupyter Notebooks to reproduce the figures and tables in Ye et al. (in review).

## Data requirements
- Tropical cyclone wind speed data: [Synthetic Tropical cyclOne geneRation Model (STORM, version 3)](https://data.4tu.nl/articles/dataset/STORM_climate_change_tropical_cyclone_wind_speed_return_periods/14510817)
- Coastal flood data: [Aqueduct Floods Hazard Maps (version 2)](http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2/index.html)
- Electricity infrastructure data:
  - OpenStreetMap (OSM). Download the latest .osm.pbf file, then edit the osmconf.ini file to extract the specific assets.
  - Government power grid maps. This compiled dataset of power plants, substations, and power lines is created based on power grid maps of each country in the study area from government reports and organizations (e.g., World Bank, World Resources Institute). It is open-access at https://zenodo.org/deposit/7550620#.
- Vulnerability curves: The detailed information and description of how to process these curves from various sources are available at https://zenodo.org/deposit/7550620# and Ye et al..

## Python requirements
The recommended option is to use a miniconda environment to work in for this project, relying on conda to handle some of the trickier library dependencies.
  <# Create a conda environment for the project and install packages>
  <conda env create -f environment.yml>
  <conda activate py310>

## How to cite:
If you use this repository in your work, please cite the corresponding paper:
Ye et al.
