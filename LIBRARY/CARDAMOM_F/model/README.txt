Author: T. Luke Smallman
Last updated: 22/06/2022

This file aims to provide a brief summary of the various models found within this, the Fortran code, directory.
Information on these models is important as not all of them have been kept up to day with current developments
and therefore may not run. However, the reason for many of these models continuing to exist is due to a desire
to keep a record of abandoned development lines which may, in the future, aid in developing better models.

For detailed information and appropriate references see the relevant source codes.

Primary developers / maintainers:
Mathew Williams (MW)
Anthony Alexis Bloom (AAB)
Jean-François Exbrayat (JFE)
Thomas Luke Smallman (TLS)
David Thomas Milowdowski (DTM)
Andrew Revill (AR)

Maintained models
  1) ACM (TLS/MW: MODEL; Forests2020 / NCEO)
     This model contains the ACM-GPP-ET source code for use in parameterising ACM-GPP-ET not DALEC.
     The model is currently coded and tested to assimilate consistent time series of GPP, Etrans, Esoil, Ewet.
     The current form is not set to deal with "real" data, for example from Eddy covariance.
     Be aware of the need to keep this code in step with changes to the code as applied in the ACM-GPP-ET
     repository and the BUCKET models below.
     NOTE: the implemented ACM-GPP-ET code into various DALEC models (aka ACM2) has diverged
     due to more robust representations afforded by being coupled into a C-cycle representation

  2) DALEC_CDEA_LU_FIRES (TLS/JFE/AAB/MW: WORKING MODEL; NCEO LTSS)
     Application of the DALEC-CDEA model which has biomass removal and fire sub-models imposed.
     These sub-models are developed substantially since their original implementations.
     This model is the closest to that used in Bloom et al., 2016
     Note that in all other actively developed models the biomass removal and fire sub-models are also implemented.

  3) DALEC_CDEA_ACM2 (TLS/DTM: WORKING MODEL; NCEO LTSS)
     Application of the DALEC-CDEA model which has biomass removal and fire sub-models imposed.
     These sub-models are developed substantially since their original implementations.
     This model is the closest to that used in Bloom et al., 2016 in terms of the C-cycle.
     The ACM-GPP-ET (aka ACM2) model has been implemented assuming soil moisture is not limiting.
     Note that in all other actively developed models the biomass removal and fire sub-models are also implemented.

  4) DALEC_CDEA_ACM2_BUCKET (TLS: WORKING MODEL; NCEO LTSS)
     Application of the DALEC-CDEA model which has biomass removal and fire sub-models imposed.
     These sub-models are developed substantially since their original implementations.
     This model is the closest to that used in Bloom et al., 2016 in terms of the C-cycle.
     The ACM-GPP-ET (aka ACM2) model has been implemented in full.

  5) DALEC_CDEA_ACM2_BUCKET_RmRg (TLS: WORKING MODEL; NCEO LTSS)
     Application of the DALEC-CDEA model which has biomass removal and fire sub-models imposed.
     These sub-models are developed substantially since their original implementations.
     Autotrophic respiration has been partitioned between growth and maintenance. Growth is assumed to
     by 28 % of NPP while maintenance continues to be a retrieved parameter.
     The ACM-GPP-ET (aka ACM2) model has been implemented in full.

  6) DALEC_CDEA_ACM2_BUCKET_LAB (AAB/TLS: WORKING MODEL; NCEO LTSS)
     Application of the DALEC-CDEA model which has biomass removal and fire sub-models imposed.
     These sub-models are developed substantially since their original implementations.
     The labile pool in this model version is allowed to have an residence time of greater than 1 year.
     The ACM-GPP-ET (aka ACM2) model has been implemented in full.

  7) DALEC_CDEA_ACM2_BUCKET_LAB_wMRT (AAB/TLS/MW: WORKING MODEL; NCEO LTSS)
     Application of the DALEC-CDEA model which has biomass removal and fire sub-models imposed.
     These sub-models are developed substantially since their original implementations.
     The labile pool in this model version is allowed to have an residence time of greater than 1 year.
     Wood turnover follows a density dependence function (Michaelis-Menten) reducing losses at low biomass.
     The ACM-GPP-ET (aka ACM2) model has been implemented in full.

  8) DALEC_CDEA_ACM2_BUCKET_wMRT (TLS/MW: WORKING MODEL; NCEO LTSS)
     Application of the DALEC-CDEA model which has biomass removal and fire sub-models imposed.
     These sub-models are developed substantially since their original implementations.
     Wood turnover follows a density dependence function (Michaelis-Menten) reducing losses at low biomass.
     The ACM-GPP-ET (aka ACM2) model has been implemented in full.

  9) DALEC_CDEA_ACM2_BUCKET_RmRg_CWD (AAB/TLS: WORKING MODEL; NCEO LTSS)
     Application of the DALEC-CDEA model which has biomass removal and fire sub-models imposed.
     These sub-models are developed substantially since their original implementations.
     Autotrophic respiration has been partitioned between growth and maintenance. Growth is assumed to
     by 28 % of NPP while maintenance continues to be a retrieved parameter.
     A wood litter pool (indicated by CWD) has been added in addition to the existing litter and som pools.
     The ACM-GPP-ET (aka ACM2) model has been implemented in full.

 10) DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT (TLS: WORKING MODEL; NCEO LTSS)
     Application of the DALEC-CDEA model which has biomass removal and fire sub-models imposed.
     These sub-models are developed substantially since their original implementations.
     Autotrophic respiration has been partitioned between growth and maintenance. Growth is assumed to
     by 28 % of NPP while maintenance continues to be a retrieved parameter.
     Wood turnover follows a density dependence function (Michaelis-Menten) reducing losses at low biomass.
     A wood litter pool (indicated by CWD) has been added in addition to the existing litter and som pools.
     The ACM-GPP-ET (aka ACM2) model has been implemented in full.

 11) DALEC_CDEA_ACM2_BUCKET_RmHeskel_Rg_CWD_wMRT (TLS: WORKING MODEL; NCEO LTSS)
     Application of the DALEC-CDEA model which has biomass removal and fire sub-models imposed.
     These sub-models are developed substantially since their original implementations.
     Autotrophic respiration has been partitioned between growth and maintenance. Growth is assumed to
     by 28 % of NPP. Canopy maintenance is assumed to follow Mary Heskel's maintenance respiration model,
     while wood and root maintenance are assumed to be calibrated fixed parameters.
     A wood litter pool (indicated by CWD) has been added in addition to the existing litter and som pools.
     The ACM-GPP-ET (aka ACM2) model has been implemented in full.

 12) DALEC_GSI_DFOL_CWD_FR (TLS: WORKING MODEL; GREENHOUSE)
     The active development version of the C-cycle model.
     This model uses a growing season index related model for phenology.
     The trend in GSI is used to determine whether canopy growth or mortality are being considered.
     For canopy growth there is a second order condition based on whether estimated GPP return is above a critical threshold.
     The revised ACM-GPP-ET models is present here, but has been modified to assume field capacity at all times.
     Therefore, the model responds to changes in atmospheric demand to evaporation but not the development of long term drought.

 13) DALEC_GSI_BUCKET (TLS: WORKING MODEL; Forests2020 / NCEO)
     As DALEC_GSI_DFOL_CWD_FR except that the full ACM-GPP-ETv1 along with water cycle is coupled.
     Except, root growth is linked to temperature and wood growth is linked to temperature and water stress.
     The trend in GSI is used to determine whether canopy growth or mortality are being considered.
     For canopy growth there is a second order condition based on whether estimated GPP return is above a critical threshold.

 14) DALEC (TLS: NCEO)
     The revised ACM-GPP-ET models is present here, but has been modified to assume field capacity at all times.
     Implementation of maintenance respiration model linked to retrieved tissue C:N based on Reich et al (2008).
     The phenology model is a simplification of the GSI approach, where the trajectory of GPP increment is used to determine growth and scenesence.
     This system will form an intermediate step between the C-cycle only system and the C+H2O coupled system with mechanistic phenology.

 15) DALEC_BUCKET (TLS: NCEO)
     The active development of the C+H2O cycles model.
     Full implementation of ACM-GPP-ET.
     Implementation of maintenance respiration model linked to retrieved tissue C:N based on Reich et al (2008).
     The phenology model is a simplification of the GSI approach, where the trajectory of GPP increment is used to determine growth and scenesence.
     This system will form an intermediate step between the C-cycle only system and the C+H2O coupled system with mechanistic phenology.

Preserved models (ASSUME NON-FUNCTIONING), project / purpose of code development in parenthesise
   1) ACM_TESSEL (TLS: GREENHOUSE)
      The DALEC_CDEA model coupled to a special re-calibration of the ACM-GPP model to emulate photosynthetic responses of the CTESSEL model (ECMWF land surface scheme)
   2) AT_DALEC (TLS: GREENHOUSE)
      The DALEC_CDEA model coupled to random forest based emulation photosynthetic responses of the CTESSEL model (ECMWF land surface scheme)
   3) DALEC_CDEA (AAB: GEOCARBON / NCEO)
      Fortran version of the DALEC-CDEA model as found in Bloom & Williams (2015)
      EDCs also as found in Bloom & Williams (2015)
   4) DALEC_CDEA_FR (TLS: Innovate UK, formally TSB)
      Fortran version of the DALEC-CDEA model as found in Bloom & Williams (2015)
      EDCs modifie to weaken the steady state assumption an experiment into simulating aggrading forest ecosystems (_FR; Forest Rotation)
   5) DALEC_CDEA_LU_FIRES_HBV (JFE: NCEO LTSS)
      Does not run
      Does exists in C also which does run
   7) DALEC_GSI_DBio_FR (TLS: GREENHOUSE)
      Coupling of C-only DALEC to the microbial decomposition model (Xenakis & Williams 2014).
      The microbial model required a large number of parameters for which good priors and inter-relations is difficult to quantify.
      Preserved to give basic model structure of microbial model should one be included in future, particularly if linked to N cycle model instead.
   8) DALEC_GSI_DFOL_FR (TLS: GREENHOUSE)
      Development of the GSI model where growth or scenesence are based on whether GSI is increasing or decreasing beyond a parameterised threshold.
      This model built on to create the DALEC_GSI_DFOL_CWD_FR into which a coarse woody debris pool was added
   9) DALEC_GSI_DFOL_FROOT_FR (TLS: GREENHOUSE)
      Development to split the root pool into 2 responsible for uptake and transport roles.
      This is a C-cycle model only so the physiological roles were un-constrained.
      However the basic influence of low CN short lived uptake roots vs long lived high CN transport roots should be reinvestigated in due course.
  10) DALEC_GSI_DFOL_LABILE_FR (TLS: GREENHOUSE)
      C-cycle exploration of model structure to have single labile pool supplying Rm, Rg, fol, wood, roots.
      Due to lack of water or N cycles creating demand for tissues this was uncontrainable.
      This model structure was used as the template for the DALECN_BUCKET model which has a water cycle to impose some requirement on allocation.
      Moreover, the Reich et al (2008) was added to represent Rm in later models which provides division between Rg and Rm demand.
  11) DALEC_GSI_FR (JFE/TLS: GREENHOUSE)
      First implentation of the GSI phenology sub-model. Built directory on DALEC-CDEA replacing the fixed phenology model.
      This implentation assumes concurrent leaf growth and loss based on the GSI / 1-GSI.
  12) DALECN_CYCLE_GSI_BUCKET (TLS: Forests2020)
      Attempt to include a full N cycle based on ACONITE.
  13) DALECN_GSI_BUCKET (TLS: Forests2020)
  14) DALECN_GSI_DFOL_LABILE_FR (TLS: Forests2020)
      Stepping stone from DALEC_GSI_DFOL_CWD_FR, with an implementation of the Reich et al., (2008) maintenance respiratoin model.
      Used as part basis for the DALECN_BUCKET model.
  15) DALECN_GSI_FR (TLS: Forests2020)
  5) DALECN_BUCKET (TLS: HIATUS DEVELOPMENT; Forests2020 / NCEO)
     The active development of the C+H2O cycles model.
     Full implementation of ACM-GPP-ET.
     Implementation of maintenance respiration model linked to retrieved tissue C:N based on Reich et al (2008).
     Based on revised ACM-GPP-ET, a maintenance respiration and leaf age model the phenology model is being revised to make growth and mortality choices based on marginal return calculation.
     This system is still in development, as a result there is a simpler intermediate model between the DALEC_GSI_DFOL_CWD_FR and DALEC_BUCKET.
