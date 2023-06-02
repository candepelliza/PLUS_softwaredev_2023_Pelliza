# Relative Elevation Model: Mendoza River

This workflow generates a Relative Elevation Model (REM) for a river watershed, starting from a DEM file, making possible to visualize more clearly the river changes and waterflow along time. 

The script is developed following the [tutorial by Matt Forrest](https://github.com/mbforr/youtube-examples/tree/main/relative-elevation-model), and applied to the case of Mendoza River the main river at the city with the same name in Argentina. The DEM files were obtained from the ALOS-1 satellite with a resolution of 30m, which gave a much lower resolution result than the original example.

The libraries utilized for the workflow are pathlib, IPython, numpy, pandas, geopandas, osmnx, shapely, scipy, xarray, xrspatial, rioxarray, matloplib, datadasher and geojson. The environment file is available for download.

