###
## Description for the example of analysis for CARDAMOM

Author: T. L. Smallman (t.l.smallman@ed.ac.uk)
Date: 08/11/2021

Contents:
This example contains an example control file, time series meteorological inputs, time series observational constraints and 
observational information for the initial state for a Fluxnet site (FI-Hyy).
The analysis covers a 16 year period (1999-2014) at a monthly time step. Obervational constrains available in the inputs are LAI, 
NEE and Evaporation.
NOTE: that while the example here uses DALEC_CDEA_LU_FIRES (also known as DALEC2, C1 or M1) the model being used can be changed 
within the control script. 

To Run:
Using R within a linux environment this example can be operated by iterating through the different stages (-1->4) found in the 
control file ("control_CDEA_FI-Hyy.r"). First set the working directory to the location the source code has been storged in. 
Next set each stage in turn within the control file, save the file and then run the file from R.

> source("control_CDEA_FI-Hyy.r")

On First running of CARDAMOM you will need to specify the location of your output directory (which must have already been created).
You will also be asked to provide the address of a remote cluster to work on. The R code is currently only set up to use 
Edinburgh University's high performance computer (HPC) Eddie (eddie.ecdf.ed.ac.uk). A new user can run on your local machine 
"out of the box" while running on a remote server is possible, it will require system specific modifications to the R code.
Even if you will using your local machine and not Eddie (or other HPC) you can just fill in Eddie's address to fill the question.
The actual choice of running a specific job on the local machine or remote cluster is made using the control script during stage -1.

When running a job using "site_specific" input files (i.e. not drawing from global databases) the time step and time period covered by 
the input files must exactly match that expected by the source code and specified in the control file.


