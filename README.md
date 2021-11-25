# Prerequisites
This processing chain works on Python 2.7 and works mainly on GDAL and OTB.
Before running this script your have ton install OTB following this [link](https://www.orfeo-toolbox.org/ "link").
We recommand you to use Anaconda or Miniconda in order to install multiple python environments.

## Frequent errors or warnings
- "ERROR 4: Unable to open EPSG support file gcs.csv.  Try setting the GDAL_DATA environment variable to point to the directory containing EPSG csv files." This error may occurs while using Rasterio. You just need to unset GDAL_DATA (Source : https://rasterio.readthedocs.io/en/latest/faq.html)

## Notes
- Automatic S1 download doesn't work. Prefer using [s1tiling](https://s1-tiling.pages.orfeo-toolbox.org/s1tiling/latest/install.html)
- Reference data can be downloaded following this [link](https://www.geograndest.fr/portail/fr/projets/occupation-du-sol "link"). If you're working on multiple departments, you have to convert gpkg to shp and merge them if needed

## TODO
- Accelerate treatments to create the dataset in less than 2 or 3 days
- Propose some function to extract every patch available per tile, in a time interval ...
