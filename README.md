# HBM-restoring

This project includes codes to do restoring for salinity and temperature for HBM ocean model. The restoring is a common practice to force ocean models to maintain realistic conditions. See Guideline.md for mode details.

The restoring tools is built as an HBM tool.
The strategy is to run as a separate program in a gpc node (avoid queuing), which means a program to change the restart file. After testing, it could be built in HBM at the later stage.
The restoring frequency will be flexible in the tool. In the testing run, we will do monthly restoring, which means the relaxation time is a month to the climatological fields to avoid spikes in the results
