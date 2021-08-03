# The width and value of residential streets

Code to replicate the analysis in:
Millard-Ball, Adam. 2021. "The width and value of residential streets." *Journal of the American Planning Association*, in press.

For more information visit [https://streetwidths.its.ucla.edu](https://streetwidths.its.ucla.edu)

Requires Postgres/PostGIS and Python. 

To run (abbreviated guide - please contact me for assistance):
- set up your directory structure as in definePaths()
- change your Postgres credentials in the `tools.dbConnection()` calls 
- add the parcel data, OSM roads data, HISDAC, land price data, and tract boundaries to the relevant directories
- call `run_counties(list_of_counties)` to create the output files with street width and value

This study was made possible through funding received by the University of California Institute of Transportation Studies Statewide Transportation Research Program, funded by the Road Repair and Accountability Act of 2017 (SB 1). The authors would like to thank the State of California for its support of university-based research, and especially for the funding provided in support of this project.
