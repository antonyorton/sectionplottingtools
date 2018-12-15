# sectionplottingtools
Python functions to visualise boreholes and DEMs in section
Functions are provided to do the following:
  * Load ESRI .asc files to pandas dataframes
  * Load borehole coordinate data and lithology (geology) to pandas
  * Load shapefiles to shapely objects
  * Plot borehole data and DEM surface elevations on cross sections
  
Please read individual function docstrings for more information

DEMs can be obtained at http://elevation.fsdf.org.au/ for Australia

Borehole data should be contained in two csv files:
  1. containing the word 'Hole' in the file name and column names 'borehole, x, y, top_rl, EOH_depth, dip, dip_direction'
    provides the coordinate information for each hole
  2. containing the word 'Geology' in the file name and column names 'borehole, fromDepth, material'
    specifies the top depth of each geological unit
