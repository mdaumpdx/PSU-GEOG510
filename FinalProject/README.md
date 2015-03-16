This folder contains the final project source code for georeferencing 
Rapid Bio-Assessment survey data. 

It contains two python scripts, define_RBA_dist_adj_factors.py and
georef_RBA_survey_dat.py, and a utility module, RBA_georef_util.py, 
containing code common to both scripts.


Steps for use with RBA survey data:

1. Add X, Y coordinates for sync points to survey data.

2. Calculate adjustment factors:
Run define_RBA_dist_adj_factors, yielding stream distance info 
in a csv file.

3. Review the adjustment factors, and possibly modify 
sync points or manually change adjustment factors.

4. Georeference Survey Data:
Run georef_RBA_survey_data, yielding point feature
class for surveyed data, completely populated. 

5. Review locations visually.



