# Risk analysis of natural hazards to power grids in Southeast and East Asia
This repository provides the code to perform the risk analysis of tropical cyclones and coastal floods to power grids in Southeast and East Asia on an asset level.

It also provides Jupyter Notebooks to reproduce the figures and tables in [Ye et al (2024)](Ye, M., Ward, P.J., Bloemendaal, N. et al. Risk of Tropical Cyclones and Floods to Power Grids in Southeast and East Asia. Int J Disaster Risk Sci 15, 494–507 (2024). https://doi.org/10.1007/s13753-024-00573-7).

## Data requirements
- Tropical cyclone wind speed data: [Synthetic Tropical cyclOne geneRation Model (STORM, version 3)](https://data.4tu.nl/articles/dataset/STORM_climate_change_tropical_cyclone_wind_speed_return_periods/14510817)
- Coastal flood data: [Aqueduct Floods Hazard Maps (version 2)](http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2/index.html)
- Electricity infrastructure data:
  - OpenStreetMap (OSM). Download the latest .osm.pbf file, then edit the osmconf.ini file to extract the specific assets.
  - Government power grid maps. This compiled dataset of power plants, substations, and power lines is created based on power grid maps of each country in the study area from government reports and organizations (e.g., World Bank, World Resources Institute). It is open-access at [Electricity Infrastructure and Vulnerability Database for Power Gird Risk Assessment](https://zenodo.org/records/7550620).
- Vulnerability curves: The detailed information and description of how to process these curves from various sources are available at [Electricity Infrastructure and Vulnerability Database for Power Gird Risk Assessment](https://zenodo.org/records/7550620).

## Python requirements
The recommended option is to use a miniconda environment to work in for this project, relying on conda to handle some of the trickier library dependencies.
```
# Create a conda environment for the project and install packages
conda env create -f environment.yml
conda activate py310
```

## How to cite:
If you use this repository in your work, please cite the corresponding paper:<br>
Ye, M., Ward, P.J., Bloemendaal, N. et al. Risk of Tropical Cyclones and Floods to Power Grids in Southeast and East Asia. International Journal of Disaster Risk Science 15, 494–507 (2024). https://doi.org/10.1007/s13753-024-00573-7

If you would like to refer to the [Electricity Infrastructure and Vulnerability Database for Power Gird Risk Assessment](https://zenodo.org/records/7550620), please cite:<br>
Ye, M., Ward, P., Bloemendaal, N., Nirandjan, S., & Koks, E. (2023). Electricity Infrastructure and Vulnerability Database for Power Gird Risk Assessment [Data set]. Zenodo. https://doi.org/10.5281/zenodo.7550620
