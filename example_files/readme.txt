###
## Description for the example of analysis for CARDAMOM

Author: T. L. Smallman
Date: 04/03/2021

Contents:
This example contains an example control file, meteorological inputs and observational constraints for a Fluxnet site (FI-Hyy).
The analysis covers a 16 year period (1999-2014) at a monthly time step. Obervational constrains available in the inputs are LAI, NEE and Evaporation.
NOTE: that the example here uses DALEC_CDEA_LU_FIRES (also known as DALEC2, C1 or M1) simulates the carbon cycle only.

To Run:
Using R within a linux environment this example can be operated by iterating through the different stages (-1->4) found in the control file ("control_CDEA_FI-Hyy.r").
First set the working directory to the location the source code has been storged in. Next set each stage in turn within the control file, save the file and then run the file from R.

> source("control_CDEA_FI-Hyy.r")


