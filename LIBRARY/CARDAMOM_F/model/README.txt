
This file aims to provide a brief summary of the various models found within this directory.
Information on these models is important as not all of them have been kept up to day with current developments and therefore may not run.
However, the reason for many of these models continuing to exist is due to a desire to keep a record of abandoned development lines which may, in the future, aid in developing better models.

Actively maintained models
  1) ACM (WORKING MODEL)
     This model contains the ACM-GPP-ET source code for use in parameterising ACM-GPP-ET not DALEC.
     The model is currently coded and tested to assimilate consistent time series of GPP, Etrans, Esoil, Ewet.
     The current form is not set to deal with "real" data, for example from Eddy covariance.
     Be aware of the need to keep this code in step with changes to the code as applied in the ACM-GPP-ET repository and the BUCKET models below.

  2) DALEC_GSI_DFOL_CWD_FR (WORKING MODEL)
     The active development version of the C-cycle model.
     This model uses a growing season index related model for phenology.
     The revised ACM-GPP-ET models is present here, but has been modified to assume field capacity at all times.
     Therefore, the model responds to changes in atmospheric demand to evaporation but not the development of long term drought.

  3) DALEC_GSI_BUCKET (WORKING MODEL)
     As DALEC_GSI_DFOL_CWD_FR except that the full ACM-GPP-ETv1 along with water cycle is coupled.

  4) DALEC_BUCKET (IN DEVELOPMENT)
     The active development of the C+H2O cycles model.
     Full implementation of ACM-GPP-ET.
     Implementation of maintenance respiration model linked to retrieved tissue C:N based on Reich et al (2008).
     The phenology model is a simplification of the GSI approach, where the trajectory of GPP increment is used to determine growth and scenesence.
     This system will form an intermediate step between the C-cycle only system and the C+H2O coupled system with mechanistic phenology.

  5) DALECN_BUCKET (IN DEVELOPMENT)
     The active development of the C+H2O cycles model.
     Full implementation of ACM-GPP-ET.
     Implementation of maintenance respiration model linked to retrieved tissue C:N based on Reich et al (2008).
     Based on revised ACM-GPP-ET, a maintenance respiration and leaf age model the phenology model is being revised to make growth and mortality choices based on marginal return calculation.
     This system is still in development, as a result there is a simpler intermediate model between the DALEC_GSI_DFOL_CWD_FR and DALEC_BUCKET.

  6) DALEC_CDEA_LU_FIRES (WORKING MODEL)
     Application of the DALEC-CDEA model which has biomass removal and fire sub-models imposed.
     This model is the closed to that used in Bloom et al., 2016
     Note that in all other actively developed models the biomass removal and fire sub-models are also implemented.

Preserved models
   1) ACM_TESSEL
   2) AT_DALEC
   3) DALEC_CDEA
   4) DALEC_CDEA_FR
   5) DALEC_CDEA_LU_FIRES_HBV
      Does not run
      Does exists in C also which does run
   7) DALEC_GSI_DBio_FR
   8) DALEC_GSI_DFOL_FR
   9) DALEC_GSI_DFOL_FROOT_FR
  10) DALEC_GSI_DFO_LABILE_FR
  11) DALEC_GSI_FR
  12) DALECN_CYCLE_GSI_BUCKET
  13) DALECN_GSI_BUCKET
  14) DALECN_GSI_DFOL_LABILE_FR
  15) DALECN_GSI_FR
