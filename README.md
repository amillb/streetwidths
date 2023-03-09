# The width and value of residential streets

Code to replicate the analysis in:
Millard-Ball, Adam. 2021. "The width and value of residential streets." *Journal of the American Planning Association*, in press.

For more information visit [https://streetwidths.its.ucla.edu](https://streetwidths.its.ucla.edu)

## Running this project

This guide provides, hopefully, the easiest way to set up your environment to
run this project. It assumes you have the following installed and are familiar
with using them:
- [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
- [PostGIS](https://postgis.net/install/)
- [Java (64-bit)](https://www.java.com/en/download/help/download_options.html)

Setting up the conda environment:

1. After cloning this repository, create a new conda environment using the
   provided environment.yml file (`conda env create -f environment.yml`).
2. Clone this repository anywhere and keep note of the directory it is in:
   [https://github.com/amillb/pgMapMatch](https://github.com/amillb/pgMapMatch)
3. Within pgMapMatch, duplicate config_template.py and name the copy config.py.
   The contents of config.py do not have to be changed.
4. Activate the new conda environment (`conda activate streetwidths`), and add
   pgMapMatch to the environment's PATH (`conda develop
   /path/to/the/directory/containing/pgMapMatch`, not pgMapMatch itself).

Setting up the PostgreSQL database:

1. Using pgAdmin or psql, create a database named streetwidths.
2. Connect to the database and perform the following:
   - Create a schema named main and another named rawdata.
   - Run the SQL `CREATE EXTENSION postgis; CREATE EXTENSION postgis_raster;`.
3. Open pg_hba.conf and write the following line at the top of the file:
```
host    streetwidths    all             ::1/128                 trust
```
- (pg_hba.conf can be located by running the SQL `SHOW data_directory;`.)
4. Make sure your PostgreSQL binaries directory is added to your PATH
   environment variable.

Setting up the data and osm2po:

1. Locate the README's inside the following directories, and download (and
   extract) the files from the sources listed by the README's into the same
   directories:
   - Data
   - Data/Tracts
   - Data/Parcel_data
   - Data/Roads
   - Data/Other
   - osm2po-5
2. For each county, you'll need to know the EPSG code for the State Plane (in meters). Add it to the dictionary called `projections` at the start of `calculate_street_value.py`. (If anyone wants to share a lookup file, please do and we can automate this step.)
3. After downloading all the data, make sure the directories resemble this
   structure:
- Data
  - FBUY
    - data
      - FBUY.tif
      - ...
    - ...
  - Other
    - Land-Prices_DLOS_2020_26Oct.xlsx
    - tl_2019_us_uac10.shp
    - ...
  - Parcel_data
    - (your_counties)
      - (parcel_boundaries).shp
      - ...
    - ...
  - Roads
    - (your_states)-latest.osm.pbf
    - ...
  - Tracts
    - cb_2019_us_tract_500k.shp
    - ...
- osm2po-5
  - osm2po-core-5.2.43-signed.jar
  - ...

Note that the parcels data should have voids for the street rights of way, rather than representing streets with their own parcel(s). If you don't have avoids, you will probably get an out of memory error. (Feel free to adapt the code to make it robust to this situation.)

Additional Windows setup:

1. Make sure Powershell is installed and added to your PATH. (This is usually
   already done in a normal Windows installation.)
2. Install the GDAL package from [OSGeo4W](https://trac.osgeo.org/osgeo4w/).
   (The ogr2ogr that comes with conda and PostGIS usually do not work or have
   the necessary drivers.)

Calling calculate_street_value.py:

- At this point, the project can be run: make sure the conda environment is
  activated (`conda activate streetwidths`), and simply run
  calculate_street_value.py using python (navigate to this project's root folder
  and call `python calculate_street_value.py`).
- If you wish to use a different PostgreSQL user, change the `postgres_user`
  variable near the top of calculate_street_value.py accordingly, and make sure
  the user has the same privileges as postgres.
- calculate_street_value.py has been set up to automatically detect what
  counties are present in the Parcel_data directory, but if you want to manually
  choose what counties to run, remove the code under `if __name__ ==
  '__main__':` and replace it with a call to `run_counties()`, passing in a list
  of county names matching those in Parcel_data (i.e.
  `run_counties(['san_francisco_ca','tarrant_tx'])`).
- If you want to recreate the JAPA analysis, you can use `runAll()`. However, you will need the parcel shapefiles for all 20 counties, and the underlying census data too. Feel free to contact the author for assistance. 

## Where are my results?
They will be in the `main` schema as a postgres table: `osm_metrics_XXXXX`, where XXXXX is the county fips code. You can also see some of the intermediate tables for the county.

See the [data dictionary](https://github.com/amillb/streetwidths/blob/master/data_dictionary.csv) for details of how to interpret each column.

## Credits

This study was made possible through funding received by the University of California Institute of Transportation Studies Statewide Transportation Research Program, funded by the Road Repair and Accountability Act of 2017 (SB 1). The authors would like to thank the State of California for its support of university-based research, and especially for the funding provided in support of this project.

Thanks to Theodore Lau for programming assistance and documentation. 
