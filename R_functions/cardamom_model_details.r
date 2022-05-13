
###
## Function which contains basic information about the models
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

cardamom_model_details <-function(modelname,specific_pft,ctessel_pft) {

  if (modelname == "ACM") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(2,dim=c(length(ctessel_pft)))
    nopars=array(20,dim=c(length(ctessel_pft)))
    nofluxes=array(4,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="ACM",nopools=nopools,nofluxes=nofluxes,nomet=19+4,nopars=nopars)
  } else if (modelname == "DALEC_CROP") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(35,dim=c(length(ctessel_pft)))
    nofluxes=array(21,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_CROP",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_CROP_BUCKET") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(9,dim=c(length(ctessel_pft)))
    nopars=array(38,dim=c(length(ctessel_pft)))
    nofluxes=array(21,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_CROP_BUCKET",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_CDEA_LU_FIRES") {
    # Information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(6,dim=c(length(ctessel_pft)))
    nopars=array(28,dim=c(length(ctessel_pft)))
    nofluxes=array(39,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_CDEA_LU_FIRES",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_CDEA_ACM2") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(6,dim=c(length(ctessel_pft)))
    nopars=array(28,dim=c(length(ctessel_pft)))
    nofluxes=array(39,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_CDEA_ACM2",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_CDEA_ACM2_BUCKET") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(32,dim=c(length(ctessel_pft)))
    nofluxes=array(40,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_CDEA_ACM2_BUCKET",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_CDEA_ACM2_BUCKET_wMRT") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(33,dim=c(length(ctessel_pft)))
    nofluxes=array(40,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_CDEA_ACM2_BUCKET_wMRT",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_CDEA_ACM2_BUCKET_LAB") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(33,dim=c(length(ctessel_pft)))
    nofluxes=array(40,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_CDEA_ACM2_BUCKET_LAB",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_CDEA_ACM2_BUCKET_LAB_wMRT") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(34,dim=c(length(ctessel_pft)))
    nofluxes=array(40,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_CDEA_ACM2_BUCKET_LAB_wMRT",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_CDEA_ACM2_BUCKET_RmRg") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(32,dim=c(length(ctessel_pft)))
    nofluxes=array(40,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_CDEA_ACM2_BUCKET_RmRg",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_CDEA_ACM2_BUCKET_RmRg_CWD") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(35,dim=c(length(ctessel_pft)))
    nofluxes=array(44,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_CDEA_ACM2_BUCKET_RmRg_CWD",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(36,dim=c(length(ctessel_pft)))
    nofluxes=array(44,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_CDEA_ACM2_BUCKET_RmHeskel_Rg_CWD_wMRT") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(38,dim=c(length(ctessel_pft)))
    nofluxes=array(44,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_CDEA_ACM2_BUCKET_RmHeskel_Rg_CWD_wMRT",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_CDEA_no_lit_root") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(4,dim=c(length(ctessel_pft)))
    nopars=array(17,dim=c(length(ctessel_pft)))
    nofluxes=array(28,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_CDEA_no_lit_root",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_EVERGREEN") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(5,dim=c(length(ctessel_pft)))
    nopars=array(17,dim=c(length(ctessel_pft)))
    nofluxes=array(28,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_EVERGREEN",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_EVERGREEN_no_lit_root") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(3,dim=c(length(ctessel_pft)))
    nopars=array(11,dim=c(length(ctessel_pft)))
    nofluxes=array(28,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_EVERGREEN_no_lit_root",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(43,dim=c(length(ctessel_pft)))
    nofluxes=array(25,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_BUCKET") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(43,dim=c(length(ctessel_pft)))
    nofluxes=array(25,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_BUCKET",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_G5") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(46,dim=c(length(ctessel_pft)))
    nofluxes=array(25,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_G5",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_G6") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(46,dim=c(length(ctessel_pft)))
    nofluxes=array(25,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_G6",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_BUCKET_CanAGE") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(48,dim=c(length(ctessel_pft)))
    nofluxes=array(25,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_BUCKET_CanAGE",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_GSI_BUCKET") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(46,dim=c(length(ctessel_pft)))
    nofluxes=array(25,dim=c(length(ctessel_pft)))
    if (specific_pft == "pft_specific") {nopars[which(ctessel_pft == 1)]=38 ; nofluxes[which(ctessel_pft == 1)]=21 ; nopools[which(ctessel_pft == 1)]=9}
    cardamom_model_details=list(name="DALEC_GSI_BUCKET",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_GSI_DFOL_CWD_FR") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(37,dim=c(length(ctessel_pft)))
    nofluxes=array(25,dim=c(length(ctessel_pft)))
    if (specific_pft == "pft_specific") {nopars[which(ctessel_pft == 1)]=35 ; nofluxes[which(ctessel_pft == 1)]=21 ; nopools[which(ctessel_pft == 1)]=8}
    cardamom_model_details=list(name="DALEC_GSI_DFOL_CWD_FR",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_1005") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(37,dim=c(length(ctessel_pft)))
    nofluxes=array(32,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_1005",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_1005a") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(37,dim=c(length(ctessel_pft)))
    nofluxes=array(32,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_1005a",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } # If modelname == "..."

} # end function cardamom_model_details

## Use byte compile
cardamom_model_details<-cmpfun(cardamom_model_details)
