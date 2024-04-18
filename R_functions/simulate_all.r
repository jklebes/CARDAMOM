
###
## Function to generate mean state variable information by running the parameters and model choice
###

# Function to actually restrict values of a double precision numeric into single precision allowed range
# Note: modified from "readBrukerFlexData" library
double2single <-function(x) {

    # Ensure the incoming variable is a double precision variable
    stopifnot(is.double(x))
    # Create a virtual connection which we will use to write / read the single precision variable version
    virtualCon = raw()
    # Write out to virtual connection in 4L (i.e. 4 byte single precision)
    virtualCon = writeBin(object = x, con = virtualCon, size = 4L) # 4L defines single precision
    # Read back into the new variable. what = double() as R does not have a format for single precision
    y = readBin(con = virtualCon, what = double(), size = 4L, n = length(x))
    # Return back to user
    return(y)

} # end function double2single

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

simulate_all<- function (site,PROJECT,model_name,met,pars,lat,pft,parameter_type,exepath,soil_info) {

  output_dim=17 ; aNPP_dim = 3 ; MTT_dim = 5 ; SS_dim = 5 ; fire_dim = 6
  noyears = length(c(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)))

  # restructure pars
  if (length(dim(pars)) > 2) {
      pars_in = array(0,dim=c(dim(pars)[1],dim(pars)[2]*dim(pars)[3]))
      nos_iter = dim(pars)[2]*dim(pars)[3]
      for (n in seq(1,dim(pars)[1])){
           pars_in[n,] = pars[n,,]
      }
  } else if (length(dim(pars)) == 2) {
      nos_iter = dim(pars)[2]
      pars_in = pars
  } else {
      nos_iter = 1
      pars_in = array(pars, dim=c(length(pars),1))
  }

  # loop through combinations
  if (model_name == "ACM") {
      # load the function code from the dalec shared object.
      # NOTE: that the name of the shared object is hardcoded to dalec.so
      # while the function call within will be specific to the actual model being called
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      output_dim = 11
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "racm",output_dim=as.integer(output_dim),met=as.double(t(met)),pars=as.double(pars_in)
                          ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                          ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                          ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                          ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                          ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                          ,soil_frac_clay=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                          ,soil_frac_sand=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3))) )
      output=tmp$out_var
      output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # construct output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      ET_kgH2Om2day = output[,,3],soilET_kgH2Om2day = output[,,4],
                      Rtot_MPasm2mmol = output[,,5],wetcanET_kgH2Om2day = output[,,6],
                      gs_demand_supply = output[,,7], gs_total_canopy = output[,,8],
                      APAR_MJm2day = output[,,9], gb_total_canopy = output[,,10],
                      CiCa = output[,,11])
  } else if (model_name == "DALEC.C3.M1.#") {
# THIS CODE AND THE R INTERFACE NEED UPDATING TO MAKE OPTIMAL USE OF THE VARIABLES HERE
# The C-only crop model had not yet been integrated
      output_dim=28
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalec14",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
                                      ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                      ,fire_dim=as.integer(fire_dim)
                                      ,met=as.double(t(met)),pars=as.double(pars_in)
                                      ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                      ,out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                      ,out_var3=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                      ,out_var4=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                      ,out_var5=as.double(array(0,dim=c(nos_iter,MTT_dim,noyears)))
                                      ,out_var6=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                      ,out_var7=as.double(array(0,dim=c(nos_iter,fire_dim,noyears)))
                                      ,out_var8=as.double(array(0,dim=c(nos_iter,fire_dim,noyears)))
                                      ,out_var9=as.double(array(0,dim=c(nos_iter,fire_dim,noyears)))
                                      ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site])
                                      ,nomet=as.integer(dim(met)[2]),nofluxes=as.integer(PROJECT$model$nofluxes[site])
                                      ,nopools=as.integer(PROJECT$model$nopools[site]),pft=as.integer(pft)
                                      ,pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                      ,nos_years=as.integer(noyears)
                                      ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                                      ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                                      ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3)))
                                      ,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
      output = tmp$out_var   ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2    ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3     ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var4 ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5    ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6  ; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      FIREemiss = tmp$out_var7; FIREemiss = array(FIREemiss, dim=c(nos_iter,fire_dim,noyears))
      FIRElit = tmp$out_var8  ; FIRElit = array(FIRElit, dim=c(nos_iter,fire_dim,noyears))
      outflux_nat = tmp$out_var9  ; outflux_nat = array(outflux_nat, dim=c(nos_iter,fire_dim,noyears))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object (14,15,16,17 unused)
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], litter_gCm2 = output[,,10],
                      labile_gCm2 = output[,,11], foliage_gCm2 = output[,,12],
                      harvest_gCm2day = output[,,13], ET_kgH2Om2day = output[,,18],
                      SurfWater_kgH2Om2 = output[,,19], wSWP_MPa = output[,,20],
                      woodlitter_gCm2 = output[,,21], auto_gCm2 = output[,,22],
                      fire_gCm2day = output[,,23],
                      gs_demand_supply = output[,,24], gs_total_canopy = output[,,25],
                      APAR_MJm2day = output[,,26], gb_total_canopy = output[,,27],
                      CiCa = output[,,28],
                      aNPP = aNPP, MTT = MTT, SS_gCm2 = SS_gCm2, aMTT = aMTT, natMTT = MTTnat,
                      FIREemiss_gCm2yr = FIREemiss, FIRElit_gCm2yr = FIRElit,
                      NAToutflux_gCm2yr = outflux_nat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
      # Final tidy
      rm(output,MTT_gCm2,SS_gCm2)
  } else if (model_name == "DALEC.A3.C3.H2.M1.#") {
      output_dim = 58 ; MTT_dim = 8 ; SS_dim = 8
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      #crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
      crop_type = 1 # Winter Wheat
      wd_old = getwd() ; setwd(PROJECT$exepath)
      tmp=.Fortran( "rdalec15",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                             
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                             ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2]))
                             ,pathlength=as.integer(crop_type))
                             #,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc() ; setwd(wd_old)
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output[,,3],
                      mean_annual_rhet_litter_gCm2day = output[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      harvest_gCm2day = output[,,5],
                      mean_harvest_gCm2day = output_mean[,5],
                      mean_annual_harvest_gCm2day = output_annual[,,5],
                      extracted_residue_gCm2day = output[,,6],
                      mean_extracted_residue_gCm2day = output_mean[,6],
                      mean_annual_extracted_residue_gCm2day = output_annual[,,6],
                      fire_gCm2day = array(0, dim(output)[1:2]), # kept because fire is a common management routine outside of europe
                      mean_fire_gCm2day = array(0, dim(output_mean)[1]), # kept because fire is a common management routine outside of europe
                      mean_annual_fire_gCm2day = array(0, dim(output_annual)[1:2]), # kept because fire is a common management routine outside of europe
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      alloc_autotrophic_gCm2day = output[,,12],
                      mean_alloc_autotrophic_gCm2day = output_mean[,12],
                      mean_annual_alloc_autotrophic_gCm2day = output_annual[,,12],
                      alloc_StorageOrgan_gCm2day = output[,,13],
                      mean_alloc_StorageOrgan_gCm2day = output_mean[,13],
                      mean_annual_alloc_StorageOrgan_gCm2day = output_annual[,,13],
                      foliage_to_litter_gCm2day = output[,,14],
                      mean_foliage_to_litter_gCm2day = output_mean[,14],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,14],
                      roots_to_litter_gCm2day = output[,,15],
                      mean_roots_to_litter_gCm2day = output_mean[,15],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,15],
                      wood_to_litter_gCm2day = output[,,16],
                      mean_wood_to_litter_gCm2day = output_mean[,16],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,16],
                      litter_to_som_gCm2day = output[,,17],
                      mean_litter_to_som_gCm2day = output_mean[,17],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,17],
                      rauto_maintenance_gCm2day = output[,,18],
                      mean_rauto_maintenance_gCm2day = output_mean[,18],
                      mean_annual_rauto_maintenance_gCm2day = output_annual[,,18],
                      rauto_labile_to_foliage_gCm2day = output[,,19],
                      mean_rauto_labile_to_foliage_gCm2day = output_mean[,19],
                      mean_annual_rauto_labile_to_foliage_gCm2day = output_annual[,,19],
                      rauto_npp_to_labile_gCm2day = output[,,20],
                      mean_rauto_npp_to_labile_gCm2day = output_mean[,20],
                      mean_annual_rauto_npp_to_labile_gCm2day = output_annual[,,20],
                      rauto_foliage_to_litter_gCm2day = output[,,21],
                      mean_rauto_foliage_to_litter_gCm2day = output_mean[,21],
                      mean_annual_rauto_foliage_to_litter_gCm2day = output_annual[,,21],
                      rauto_wood_to_litter_gCm2day = output[,,22],
                      mean_rauto_wood_to_litter_gCm2day = output_mean[,22],
                      mean_annual_rauto_wood_to_litter_gCm2day = output_annual[,,22],
                      HARVESTextracted_foliage_gCm2day = output[,,23],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,23],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,23],
                      HARVESTextracted_wood_gCm2day = output[,,24],
                      mean_HARVESTextracted_wood_gCm2day = output_annual[,,24],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,24],
                      HARVESTextracted_DeadFoliage_gCm2day = output[,,25],
                      mean_HARVESTextracted_DeadFoliage_gCm2day = output_mean[,25],
                      mean_annual_HARVESTextracted_DeadFoliage_gCm2day = output_annual[,,25],
                      HARVESTextracted_labile_gCm2day = output[,,26],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,26],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,26],
                      HARVESTlitter_foliage_gCm2day = output[,,27],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,27],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,27],
                      HARVESTlitter_wood_gCm2day = output[,,28],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,28],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,28],
                      HARVESTlitter_DeadFoliage_gCm2day = output[,,29],
                      mean_HARVESTlitter_DeadFoliage_gCm2day = output_mean[,29],
                      mean_annual_HARVESTlitter_DeadFoliage_gCm2day = output_annual[,,29],
                      HARVESTlitter_autotrophic_gCm2day = output[,,30],
                      mean_HARVESTlitter_autotrophic_gCm2day = output_mean[,30],
                      mean_annual_HARVESTlitter_autotrophic_gCm2day = output_annual[,,30],
                      HARVESTlitter_labile_gCm2day = output[,,31],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,31],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,31],
                      PLOUGHlitter_roots_gCm2day = output[,,32],
                      mean_PLOUGHlitter_roots_gCm2day = output_mean[,32],
                      mean_annual_PLOUGHlitter_roots_gCm2day = output_annual[,,32],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,33],
                      mean_labile_gCm2 = output_mean[,33],
                      mean_annual_labile_gCm2 = output_annual[,,33],
                      foliage_gCm2 = output[,,34],
                      mean_foliage_gCm2 = output_mean[,34],
                      mean_annual_foliage_gCm2 = output_annual[,,34],
                      roots_gCm2 = output[,,35],
                      mean_roots_gCm2 = output_mean[,35],
                      mean_annual_roots_gCm2 = output_annual[,,35],
                      wood_gCm2 = output[,,36],
                      mean_wood_gCm2 = output_mean[,36],
                      mean_annual_wood_gCm2 = output_annual[,,36],
                      litter_gCm2 = output[,,37],
                      mean_litter_gCm2 = output_mean[,37],
                      mean_annual_litter_gCm2 = output_annual[,,37],
                      som_gCm2 = output[,,38],
                      mean_som_gCm2 = output_mean[,38],
                      mean_annual_som_gCm2 = output_annual[,,38],
                      autotrophic_gCm2 = output[,,39],
                      mean_autotrophic_gCm2 = output_mean[,39],
                      mean_annual_autotrophic_gCm2 = output_annual[,,39],
                      StorageOrgan_gCm2 = output[,,40],
                      mean_StorageOrgan_gCm2 = output_mean[,40],
                      mean_annual_StorageOrgan_gCm2 = output_annual[,,40],
                      DeadFoliage_gCm2 = output[,,41],
                      mean_DeadFoliage_gCm2 = output_mean[,41],
                      mean_annual_DeadFoliage_gCm2 = output_annual[,,41],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,42],
                      mean_ET_kgH2Om2day = output_mean[,42],
                      mean_annual_ET_kgH2Om2day = output_annual[,,42],
                      Etrans_kgH2Om2day = output[,,43],
                      mean_Etrans_kgH2Om2day = output_mean[,43],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,43],
                      Esoil_kgH2Om2day = output[,,44],
                      mean_Esoil_kgH2Om2day = output_mean[,44],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,44],
                      Ewetcanopy_kgH2Om2day = output[,,45],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,45],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,45],
                      runoff_kgH2Om2day = output[,,46],
                      mean_runoff_kgH2Om2day = output_mean[,46],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,46],
                      underflow_kgH2Om2day = output[,,47],
                      mean_underflow_kgH2Om2day = output_mean[,47],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,47],
                      SurfWater_kgH2Om2 = output[,,48],
                      mean_SurfWater_kgH2Om2 = output_mean[,48],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,48],
                      wSWP_MPa = output[,,49],
                      mean_wSWP_MPa = output_mean[,49],
                      mean_annual_wSWP_MPa = output_annual[,,49],
                      snow_kgH2Om2 = output[,,50],
                      mean_snow_kgH2Om2 = output_mean[,50],
                      mean_annual_snow_kgH2Om2 = output_annual[,,50],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,51],
                      mean_lai_m2m2 = output_mean[,51],
                      mean_annual_lai_m2m2 = output_annual[,,51],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,52],
                      mean_gs_demand_supply_ratio = output_mean[,52],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,52],
                      gs_mmolH2Om2s = output[,,53],
                      mean_gs_mmolH2Om2s = output_mean[,53],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,53],
                      APAR_MJm2day = output[,,54],
                      mean_APAR_MJm2day = output_mean[,54],
                      mean_annual_APAR_MJm2day = output_annual[,,54],
                      gb_mmolH2Om2s = output[,,55],
                      mean_gb_mmolH2Om2s = output_mean[,55],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,55],
                      CiCa = output[,,56],
                      mean_CiCa = output_mean[,56],
                      mean_annual_CiCa = output_annual[,,56],
                      # Misc
                      RootDepth_m = output[,,57],
                      mean_RootDepth_m = output_mean[,57],
                      mean_annual_RootDepth_m = output_annual[,,57],
                      DevelopmentStage = output[,,58],
                      mean_DevelopmentStage = output_mean[,58],
                      mean_annual_DevelopmentStage = output_annual[,,58],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      MTT_autotrophic_years = MTT_years[,7],
                      MTT_DeadFoliage_years = MTT_years[,8],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6],
                      SS_autotrophic_gCm2 = SS_gCm2[,7],
                      SS_DeadFoliage_gCm2 = SS_gCm2[,8])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day +
                           states_all$alloc_autotrophic_gCm2day +
                           states_all$alloc_StorageOrgan_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean),
                           apply(states_all$alloc_autotrophic_gCm2day,1,mean),
                           apply(states_all$alloc_StorageOrgan_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      states_all$NPP_autotrophic_fraction = NPP_fraction[,4]
      states_all$NPP_StorageOrgan_fraction = NPP_fraction[,5]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P4.R2.#") {
    output_dim = 62 ; MTT_dim = 7 ; SS_dim = 7
    # Load the required dalec shared object
    dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
    tmp=.Fortran( "rdalec11",output_dim=as.integer(output_dim)
                            ,aNPP_dim=as.integer(aNPP_dim)
                            ,MTT_dim=as.integer(MTT_dim)
                            ,SS_dim = as.integer(SS_dim)
                            ,fire_dim=as.integer(fire_dim)
                            ,met=as.double(t(met))
                            ,pars=as.double(pars_in)
                            ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                            ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                            ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                            ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                            ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                            
                            ,lat=as.double(lat)
                            ,nopars=as.integer(PROJECT$model$nopars[site])
                            ,nomet=as.integer(dim(met)[2])
                            ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                            ,nopools=as.integer(PROJECT$model$nopools[site])
                            ,pft=as.integer(pft)
                            ,nodays=as.integer(dim(met)[1])
                            ,nos_years=as.integer(noyears)
                            ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                            ,nos_iter=as.integer(nos_iter)
                            ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                            ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3))))
    # Extract the different output variables
    output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
    output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    mean_gpp_gCm2day = output_mean[,1],
                    mean_annual_gpp_gCm2day = output_annual[,,1],
                    rauto_gCm2day = output[,,2],
                    mean_rauto_gCm2day = output_mean[,2],
                    mean_annual_rauto_gCm2day = output_annual[,,2],
                    rg_foliage_gCm2day = output[,,3],
                    mean_rg_foliage_gCm2day = output_mean[,3],
                    mean_annual_rg_foliage_gCm2day = output_annual[,,3],
                    rhet_litter_gCm2day = output[,,4],
                    mean_rhet_litter_gCm2day = output_mean[,4],
                    mean_annual_rhet_litter_gCm2day = output_annual[,,4],
                    rhet_som_gCm2day = output[,,5],
                    mean_rhet_som_gCm2day = output_mean[,5],
                    mean_annual_rhet_som_gCm2day = output_annual[,,5],
                    rhet_woodlitter_gCm2day = output[,,6],
                    mean_rhet_woodlitter_gCm2day = output_mean[,6],
                    mean_annual_rhet_woodlitter_gCm2day = output_annual[,,6],
                    fire_gCm2day = output[,,7],
                    mean_fire_gCm2day = output_mean[,7],
                    mean_annual_fire_gCm2day = output_annual[,,7],
                    harvest_gCm2day = output[,,8],
                    mean_harvest_gCm2day = output_mean[,8],
                    mean_annual_harvest_gCm2day = output_annual[,,8],
                    # Internal fluxes
                    alloc_labile_gCm2day = output[,,9],
                    mean_alloc_labile_gCm2day = output_mean[,9],
                    mean_annual_alloc_labile_gCm2day = output_annual[,,9],
                    alloc_roots_gCm2day = output[,,10],
                    mean_alloc_roots_gCm2day = output_mean[,10],
                    mean_annual_alloc_roots_gCm2day = output_annual[,,10],
                    alloc_wood_gCm2day = output[,,11],
                    mean_alloc_wood_gCm2day = output_mean[,11],
                    mean_annual_alloc_wood_gCm2day = output_annual[,,11],
                    labile_to_foliage_gCm2day = output[,,12],
                    mean_labile_to_foliage_gCm2day = output_mean[,12],
                    mean_annual_labile_to_foliage_gCm2day = output_annual[,,12],
                    foliage_to_litter_gCm2day = output[,,13],
                    mean_foliage_to_litter_gCm2day = output_mean[,13],
                    mean_annual_foliage_to_litter_gCm2day = output_annual[,,13],
                    roots_to_litter_gCm2day = output[,,14],
                    mean_roots_to_litter_gCm2day = output_mean[,14],
                    mean_annual_roots_to_litter_gCm2day = output_annual[,,14],
                    wood_to_litter_gCm2day = output[,,15],
                    mean_wood_to_litter_gCm2day = output_mean[,15],
                    mean_annual_wood_to_litter_gCm2day = output_annual[,,15],
                    litter_to_som_gCm2day = output[,,16],
                    mean_litter_to_som_gCm2day = output_mean[,16],
                    mean_annual_litter_to_som_gCm2day = output_annual[,,16],
                    woodlitter_to_som_gCm2day = output[,,17],
                    mean_woodlitter_to_som_gCm2day = output_mean[,17],
                    mean_annual_woodlitter_to_som_gCm2day = output_annual[,,17],
                    # Disturbance fluxes
                    FIREemiss_labile_gCm2day = output[,,18],
                    mean_FIREemiss_labile_gCm2day = output_mean[,18],
                    mean_annual_FIREemiss_labile_gCm2day = output_annual[,,18],
                    FIRElitter_labile_gCm2day = output[,,19],
                    mean_FIRElitter_labile_gCm2day = output_mean[,19],
                    mean_annual_FIRElitter_labile_gCm2day = output_annual[,,19],
                    FIREemiss_foliage_gCm2day = output[,,20],
                    mean_FIREemiss_foliage_gCm2day = output_mean[,20],
                    mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,20],
                    FIRElitter_foliage_gCm2day = output[,,21],
                    mean_FIRElitter_foliage_gCm2day = output_mean[,21],
                    mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,21],
                    FIREemiss_roots_gCm2day = output[,,22],
                    mean_FIREemiss_roots_gCm2day = output_mean[,22],
                    mean_annual_FIREemiss_roots_gCm2day = output_annual[,,22],
                    FIRElitter_roots_gCm2day = output[,,23],
                    mean_FIRElitter_roots_gCm2day = output_mean[,23],
                    mean_annual_FIRElitter_roots_gCm2day = output_annual[,,23],
                    FIREemiss_wood_gCm2day = output[,,24],
                    mean_FIREemiss_wood_gCm2day = output_mean[,24],
                    mean_annual_FIREemiss_wood_gCm2day = output_annual[,,24],
                    FIRElitter_wood_gCm2day = output[,,25],
                    mean_FIRElitter_wood_gCm2day = output_mean[,25],
                    mean_annual_FIRElitter_wood_gCm2day = output_annual[,,25],
                    FIREemiss_litter_gCm2day = output[,,26],
                    mean_FIREemiss_litter_gCm2day = output_mean[,26],
                    mean_annual_FIREemiss_litter_gCm2day = output_annual[,,26],
                    FIRElitter_litter_gCm2day = output[,,27],
                    mean_FIRElitter_litter_gCm2day = output_mean[,27],
                    mean_annual_FIRElitter_litter_gCm2day = output_annual[,,27],
                    FIREemiss_woodlitter_gCm2day = output[,,28],
                    mean_FIREemiss_woodlitter_gCm2day = output_mean[,28],
                    mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,28],
                    FIRElitter_woodlitter_gCm2day = output[,,29],
                    mean_FIRElitter_woodlitter_gCm2day = output_mean[,29],
                    mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,29],
                    FIREemiss_som_gCm2day = output[,,30],
                    mean_FIREemiss_som_gCm2day = output_mean[,30],
                    mean_annual_FIREemiss_som_gCm2day = output_annual[,,30],
                    HARVESTextracted_labile_gCm2day = output[,,31],
                    mean_HARVESTextracted_labile_gCm2day = output_mean[,31],
                    mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,31],
                    HARVESTextracted_foliage_gCm2day = output[,,32],
                    mean_HARVESTextracted_foliage_gCm2day = output_mean[,32],
                    mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,32],
                    HARVESTextracted_roots_gCm2day = output[,,33],
                    mean_HARVESTextracted_roots_gCm2day = output_mean[,33],
                    mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,33],
                    HARVESTextracted_wood_gCm2day = output[,,34],
                    mean_HARVESTextracted_wood_gCm2day = output_mean[,34],
                    mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,34],
                    HARVESTextracted_litter_gCm2day = output[,,35],
                    mean_HARVESTextracted_litter_gCm2day = output_mean[,35],
                    mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,35],
                    HARVESTextracted_woodlitter_gCm2day = output[,,36],
                    mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,36],
                    mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,36],
                    HARVESTextracted_som_gCm2day = output[,,37],
                    mean_HARVESTextracted_som_gCm2day = output_mean[,37],
                    mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,37],
                    HARVESTlitter_labile_gCm2day = output[,,38],
                    mean_HARVESTlitter_labile_gCm2day = output_mean[,38],
                    mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,38],
                    HARVESTlitter_foliage_gCm2day = output[,,39],
                    mean_HARVESTlitter_foliage_gCm2day = output_mean[,39],
                    mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,39],
                    HARVESTlitter_roots_gCm2day = output[,,40],
                    mean_HARVESTlitter_roots_gCm2day = output_mean[,40],
                    mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,40],
                    HARVESTlitter_wood_gCm2day = output[,,41],
                    mean_HARVESTlitter_wood_gCm2day = output_mean[,41],
                    mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,41],
                    # C pools (gC/m2)
                    labile_gCm2 = output[,,42],
                    mean_labile_gCm2 = output_mean[,42],
                    mean_annual_labile_gCm2 = output_annual[,,42],
                    foliage_gCm2 = output[,,43],
                    mean_foliage_gCm2 = output_mean[,43],
                    mean_annual_foliage_gCm2 = output_annual[,,43],
                    roots_gCm2 = output[,,44],
                    mean_roots_gCm2 = output_mean[,44],
                    mean_annual_roots_gCm2 = output_annual[,,44],
                    wood_gCm2 = output[,,45],
                    mean_wood_gCm2 = output_mean[,45],
                    mean_annual_wood_gCm2 = output_annual[,,45],
                    litter_gCm2 = output[,,46],
                    mean_litter_gCm2 = output_mean[,46],
                    mean_annual_litter_gCm2 = output_annual[,,46],
                    woodlitter_gCm2 = output[,,47],
                    mean_woodlitter_gCm2 = output_mean[,47],
                    mean_annual_woodlitter_gCm2 = output_annual[,,47],
                    som_gCm2 = output[,,48],
                    mean_som_gCm2 = output_mean[,48],
                    mean_annual_som_gCm2 = output_annual[,,48],
                    # Water cycle related
                    ET_kgH2Om2day = output[,,49],
                    mean_ET_kgH2Om2day = output_mean[,49],
                    mean_annual_ET_kgH2Om2day = output_annual[,,49],
                    SurfWater_kgH2Om2 = output[,,50],
                    mean_SurfWater_kgH2Om2 = output_mean[,50],
                    mean_annual_SurfWater_kgH2Om2 = output_annual[,,50],
                    wSWP_MPa = output[,,51],
                    mean_wSWP_MPa = output_mean[,51],
                    mean_annual_wSWP_MPa = output_annual[,,51],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,52],
                    mean_lai_m2m2 = output_mean[,52],
                    mean_annual_lai_m2m2 = output_annual[,,52],
                    gsi = output[,,53],
                    mean_gsi = output_mean[,53],
                    mean_annual_gsi = output_annual[,,53],
                    gsi_itemp = output[,,54],
                    mean_gsi_itemp = output_mean[,54],
                    mean_annual_gsi_itemp = output_annual[,,54],
                    gsi_iphoto = output[,,55],
                    mean_gsi_iphoto = output_mean[,55],
                    mean_annual_gsi_iphoto = output_annual[,,55],
                    gsi_ivpd = output[,,56],
                    mean_gsi_ivpd = output_mean[,56],
                    mean_annual_gsi_ivpd = output_annual[,,56],
                    # Photosynthesis / C~water coupling related
                    gs_demand_supply_ratio = output[,,57],
                    mean_gs_demand_supply_ratio = output_mean[,57],
                    mean_annual_gs_demand_supply_ratio = output_annual[,,57],
                    gs_mmolH2Om2s = output[,,58],
                    mean_gs_mmolH2Om2s = output_mean[,58],
                    mean_annual_gs_mmolH2Om2s = output_annual[,,58],
                    APAR_MJm2day = output[,,59],
                    mean_APAR_MJm2day = output_mean[,59],
                    mean_annual_APAR_MJm2day = output_annual[,,59],
                    gb_mmolH2Om2s = output[,,60],
                    mean_gb_mmolH2Om2s = output_mean[,60],
                    mean_annual_gb_mmolH2Om2s = output_annual[,,60],
                    CiCa = output[,,61],
                    mean_CiCa = output_mean[,61],
                    mean_annual_CiCa = output_annual[,,61],
                    # Misc
                    RootDepth_m = output[,,62],
                    mean_RootDepth_m = output_mean[,62],
                    mean_annual_RootDepth_m = output_annual[,,62],
                    ## Aggregated variables
                    # Mean Transit times
                    MTT_labile_years = MTT_years[,1],
                    MTT_foliage_years = MTT_years[,2],
                    MTT_roots_years = MTT_years[,3],
                    MTT_wood_years = MTT_years[,4],
                    MTT_litter_years = MTT_years[,5],
                    MTT_woodlitter_years = MTT_years[,6],
                    MTT_som_years = MTT_years[,7],
                    # Steady state estimates
                    SS_labile_gCm2 = SS_gCm2[,1],
                    SS_foliage_gCm2 = SS_gCm2[,2],
                    SS_roots_gCm2 = SS_gCm2[,3],
                    SS_wood_gCm2 = SS_gCm2[,4],
                    SS_litter_gCm2 = SS_gCm2[,5],
                    SS_woodlitter_gCm2 = SS_gCm2[,6],
                    SS_som_gCm2 = SS_gCm2[,7])
    # Determine the NPP fraction of expressed NPP
    # i.e. actual growth not GPP-Ra
    NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                         states_all$alloc_roots_gCm2day +
                         states_all$alloc_wood_gCm2day,1,mean)
    NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day,1,mean),
                         apply(states_all$alloc_roots_gCm2day,1,mean),
                         apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
    states_all$NPP_foliage_fraction = NPP_fraction[,1]
    states_all$NPP_roots_fraction = NPP_fraction[,2]
    states_all$NPP_wood_fraction = NPP_fraction[,3]
    # Tidy up variables
    rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P7.R2.#") {
      output_dim = 59 ; MTT_dim = 7 ; SS_dim = 7
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec23",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                              
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      rhet_woodlitter_gCm2day = output[,,5],
                      mean_rhet_woodlitter_gCm2day = output_mean[,5],
                      mean_annual_rhet_woodlitter_gCm2day = output_annual[,,5],
                      fire_gCm2day = output[,,6],
                      mean_fire_gCm2day = output_mean[,6],
                      mean_annual_fire_gCm2day = output_annual[,,6],
                      harvest_gCm2day = output[,,7],
                      mean_harvest_gCm2day = output_mean[,7],
                      mean_annual_harvest_gCm2day = output_annual[,,7],
                      # Internal fluxes
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      woodlitter_to_som_gCm2day = output[,,16],
                      mean_woodlitter_to_som_gCm2day = output_mean[,16],
                      mean_annual_woodlitter_to_som_gCm2day = output_annual[,,16],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,17],
                      mean_FIREemiss_labile_gCm2day = output_mean[,17],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,17],
                      FIRElitter_labile_gCm2day = output[,,18],
                      mean_FIRElitter_labile_gCm2day = output_mean[,18],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,18],
                      FIREemiss_foliage_gCm2day = output[,,19],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,19],
                      FIRElitter_foliage_gCm2day = output[,,20],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,20],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,20],
                      FIREemiss_roots_gCm2day = output[,,21],
                      mean_FIREemiss_roots_gCm2day = output_mean[,21],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,21],
                      FIRElitter_roots_gCm2day = output[,,22],
                      mean_FIRElitter_roots_gCm2day = output_mean[,22],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,22],
                      FIREemiss_wood_gCm2day = output[,,23],
                      mean_FIREemiss_wood_gCm2day = output_mean[,23],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,23],
                      FIRElitter_wood_gCm2day = output[,,24],
                      mean_FIRElitter_wood_gCm2day = output_mean[,24],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,24],
                      FIREemiss_litter_gCm2day = output[,,25],
                      mean_FIREemiss_litter_gCm2day = output_mean[,25],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,25],
                      FIRElitter_litter_gCm2day = output[,,26],
                      mean_FIRElitter_litter_gCm2day = output_mean[,26],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,26],
                      FIREemiss_woodlitter_gCm2day = output[,,27],
                      mean_FIREemiss_woodlitter_gCm2day = output_mean[,27],
                      mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,27],
                      FIRElitter_woodlitter_gCm2day = output[,,28],
                      mean_FIRElitter_woodlitter_gCm2day = output_mean[,28],
                      mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,28],
                      FIREemiss_som_gCm2day = output[,,29],
                      mean_FIREemiss_som_gCm2day = output_mean[,29],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,29],
                      HARVESTextracted_labile_gCm2day = output[,,30],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,30],
                      HARVESTextracted_foliage_gCm2day = output[,,31],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,31],
                      HARVESTextracted_roots_gCm2day = output[,,32],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,32],
                      HARVESTextracted_wood_gCm2day = output[,,33],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,33],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,33],
                      HARVESTextracted_litter_gCm2day = output[,,34],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,34],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,34],
                      HARVESTextracted_woodlitter_gCm2day = output[,,35],
                      mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,35],
                      mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,35],
                      HARVESTextracted_som_gCm2day = output[,,36],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,36],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,36],
                      HARVESTlitter_labile_gCm2day = output[,,37],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,37],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,37],
                      HARVESTlitter_foliage_gCm2day = output[,,38],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,38],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,38],
                      HARVESTlitter_roots_gCm2day = output[,,39],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,39],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,39],
                      HARVESTlitter_wood_gCm2day = output[,,40],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,40],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,40],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,41],
                      mean_labile_gCm2 = output_mean[,41],
                      mean_annual_labile_gCm2 = output_annual[,,41],
                      foliage_gCm2 = output[,,42],
                      mean_foliage_gCm2 = output_mean[,42],
                      mean_annual_foliage_gCm2 = output_annual[,,42],
                      roots_gCm2 = output[,,43],
                      mean_roots_gCm2 = output_mean[,43],
                      mean_annual_roots_gCm2 = output_annual[,,43],
                      wood_gCm2 = output[,,44],
                      mean_wood_gCm2 = output_mean[,44],
                      mean_annual_wood_gCm2 = output_annual[,,44],
                      litter_gCm2 = output[,,45],
                      mean_litter_gCm2 = output_mean[,45],
                      mean_annual_litter_gCm2 = output_annual[,,45],
                      woodlitter_gCm2 = output[,,46],
                      mean_woodlitter_gCm2 = output_mean[,46],
                      mean_annual_woodlitter_gCm2 = output_annual[,,46],
                      som_gCm2 = output[,,47],
                      mean_som_gCm2 = output_mean[,47],
                      mean_annual_som_gCm2 = output_annual[,,47],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,48],
                      mean_ET_kgH2Om2day = output_mean[,48],
                      mean_annual_ET_kgH2Om2day = output_annual[,,48],
                      SurfWater_kgH2Om2 = output[,,49],
                      mean_SurfWater_kgH2Om2 = output_mean[,49],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,49],
                      wSWP_MPa = output[,,50],
                      mean_wSWP_MPa = output_mean[,50],
                      mean_annual_wSWP_MPa = output_annual[,,50],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,51],
                      mean_lai_m2m2 = output_mean[,51],
                      mean_annual_lai_m2m2 = output_annual[,,51],
                      cgi = output[,,52],
                      mean_cgi = output_mean[,52],
                      mean_annual_cgi = output_annual[,,52],
                      ncce_gCm2day = output[,,53],
                      mean_ncce_gCm2day = output_mean[,53],
                      mean_annual_ncce_gCm2day = output_annual[,,53],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,54],
                      mean_gs_demand_supply_ratio = output_mean[,54],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,54],
                      gs_mmolH2Om2s = output[,,55],
                      mean_gs_mmolH2Om2s = output_mean[,55],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,55],
                      APAR_MJm2day = output[,,56],
                      mean_APAR_MJm2day = output_mean[,56],
                      mean_annual_APAR_MJm2day = output_annual[,,56],
                      gb_mmolH2Om2s = output[,,57],
                      mean_gb_mmolH2Om2s = output_mean[,57],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,57],
                      CiCa = output[,,58],
                      mean_CiCa = output_mean[,58],
                      mean_annual_CiCa = output_annual[,,58],
                      # Misc
                      RootDepth_m = output[,,59],
                      mean_RootDepth_m = output_mean[,59],
                      mean_annual_RootDepth_m = output_annual[,,59],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_woodlitter_years = MTT_years[,6],
                      MTT_som_years = MTT_years[,7],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_woodlitter_gCm2 = SS_gCm2[,6],
                      SS_som_gCm2 = SS_gCm2[,7])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P8.R2.#") {
      output_dim = 60 ; MTT_dim = 7 ; SS_dim = 7
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec24",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                              
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      rhet_woodlitter_gCm2day = output[,,5],
                      mean_rhet_woodlitter_gCm2day = output_mean[,5],
                      mean_annual_rhet_woodlitter_gCm2day = output_annual[,,5],
                      fire_gCm2day = output[,,6],
                      mean_fire_gCm2day = output_mean[,6],
                      mean_annual_fire_gCm2day = output_annual[,,6],
                      harvest_gCm2day = output[,,7],
                      mean_harvest_gCm2day = output_mean[,7],
                      mean_annual_harvest_gCm2day = output_annual[,,7],
                      # Internal fluxes
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      woodlitter_to_som_gCm2day = output[,,16],
                      mean_woodlitter_to_som_gCm2day = output_mean[,16],
                      mean_annual_woodlitter_to_som_gCm2day = output_annual[,,16],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,17],
                      mean_FIREemiss_labile_gCm2day = output_mean[,17],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,17],
                      FIRElitter_labile_gCm2day = output[,,18],
                      mean_FIRElitter_labile_gCm2day = output_mean[,18],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,18],
                      FIREemiss_foliage_gCm2day = output[,,19],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,19],
                      FIRElitter_foliage_gCm2day = output[,,20],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,20],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,20],
                      FIREemiss_roots_gCm2day = output[,,21],
                      mean_FIREemiss_roots_gCm2day = output_mean[,21],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,21],
                      FIRElitter_roots_gCm2day = output[,,22],
                      mean_FIRElitter_roots_gCm2day = output_mean[,22],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,22],
                      FIREemiss_wood_gCm2day = output[,,23],
                      mean_FIREemiss_wood_gCm2day = output_mean[,23],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,23],
                      FIRElitter_wood_gCm2day = output[,,24],
                      mean_FIRElitter_wood_gCm2day = output_mean[,24],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,24],
                      FIREemiss_litter_gCm2day = output[,,25],
                      mean_FIREemiss_litter_gCm2day = output_mean[,25],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,25],
                      FIRElitter_litter_gCm2day = output[,,26],
                      mean_FIRElitter_litter_gCm2day = output_mean[,26],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,26],
                      FIREemiss_woodlitter_gCm2day = output[,,27],
                      mean_FIREemiss_woodlitter_gCm2day = output_mean[,27],
                      mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,27],
                      FIRElitter_woodlitter_gCm2day = output[,,28],
                      mean_FIRElitter_woodlitter_gCm2day = output_mean[,28],
                      mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,28],
                      FIREemiss_som_gCm2day = output[,,29],
                      mean_FIREemiss_som_gCm2day = output_mean[,29],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,29],
                      HARVESTextracted_labile_gCm2day = output[,,30],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,30],
                      HARVESTextracted_foliage_gCm2day = output[,,31],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,31],
                      HARVESTextracted_roots_gCm2day = output[,,32],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,32],
                      HARVESTextracted_wood_gCm2day = output[,,33],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,33],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,33],
                      HARVESTextracted_litter_gCm2day = output[,,34],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,34],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,34],
                      HARVESTextracted_woodlitter_gCm2day = output[,,35],
                      mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,35],
                      mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,35],
                      HARVESTextracted_som_gCm2day = output[,,36],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,36],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,36],
                      HARVESTlitter_labile_gCm2day = output[,,37],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,37],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,37],
                      HARVESTlitter_foliage_gCm2day = output[,,38],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,38],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,38],
                      HARVESTlitter_roots_gCm2day = output[,,39],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,39],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,39],
                      HARVESTlitter_wood_gCm2day = output[,,40],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,40],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,40],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,41],
                      mean_labile_gCm2 = output_mean[,41],
                      mean_annual_labile_gCm2 = output_annual[,,41],
                      foliage_gCm2 = output[,,42],
                      mean_foliage_gCm2 = output_mean[,42],
                      mean_annual_foliage_gCm2 = output_annual[,,42],
                      roots_gCm2 = output[,,43],
                      mean_roots_gCm2 = output_mean[,43],
                      mean_annual_roots_gCm2 = output_annual[,,43],
                      wood_gCm2 = output[,,44],
                      mean_wood_gCm2 = output_mean[,44],
                      mean_annual_wood_gCm2 = output_annual[,,44],
                      litter_gCm2 = output[,,45],
                      mean_litter_gCm2 = output_mean[,45],
                      mean_annual_litter_gCm2 = output_annual[,,45],
                      woodlitter_gCm2 = output[,,46],
                      mean_woodlitter_gCm2 = output_mean[,46],
                      mean_annual_woodlitter_gCm2 = output_annual[,,46],
                      som_gCm2 = output[,,47],
                      mean_som_gCm2 = output_mean[,47],
                      mean_annual_som_gCm2 = output_annual[,,47],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,48],
                      mean_ET_kgH2Om2day = output_mean[,48],
                      mean_annual_ET_kgH2Om2day = output_annual[,,48],
                      SurfWater_kgH2Om2 = output[,,49],
                      mean_SurfWater_kgH2Om2 = output_mean[,49],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,49],
                      wSWP_MPa = output[,,50],
                      mean_wSWP_MPa = output_mean[,50],
                      mean_annual_wSWP_MPa = output_annual[,,50],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,51],
                      mean_lai_m2m2 = output_mean[,51],
                      mean_annual_lai_m2m2 = output_annual[,,51],
                      cgi = output[,,52],
                      mean_cgi = output_mean[,52],
                      mean_annual_cgi = output_annual[,,52],
                      cmi = output[,,53],
                      mean_cmi = output_mean[,53],
                      mean_annual_cmi = output_annual[,,53],
                      ncce_gCm2day = output[,,54],
                      mean_ncce_gCm2day = output_mean[,54],
                      mean_annual_ncce_gCm2day = output_annual[,,54],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,55],
                      mean_gs_demand_supply_ratio = output_mean[,55],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,55],
                      gs_mmolH2Om2s = output[,,56],
                      mean_gs_mmolH2Om2s = output_mean[,56],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,56],
                      APAR_MJm2day = output[,,57],
                      mean_APAR_MJm2day = output_mean[,57],
                      mean_annual_APAR_MJm2day = output_annual[,,57],
                      gb_mmolH2Om2s = output[,,58],
                      mean_gb_mmolH2Om2s = output[,,58],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,58],
                      CiCa = output[,,59],
                      mean_CiCa = output_mean[,59],
                      mean_annual_CiCa = output_annual[,,59],
                      # Misc
                      RootDepth_m = output[,,60],
                      mean_RootDepth_m = output_mean[,60],
                      mean_annual_RootDepth_m = output_annual[,,60],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_woodlitter_years = MTT_years[,6],
                      MTT_som_years = MTT_years[,7],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_woodlitter_gCm2 = SS_gCm2[,6],
                      SS_som_gCm2 = SS_gCm2[,7])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P9.R2.#") {
      output_dim = 60 ; MTT_dim = 7 ; SS_dim = 7
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec25",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                              
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      rhet_woodlitter_gCm2day = output[,,5],
                      mean_rhet_woodlitter_gCm2day = output_mean[,5],
                      mean_annual_rhet_woodlitter_gCm2day = output_annual[,,5],
                      fire_gCm2day = output[,,6],
                      mean_fire_gCm2day = output_mean[,6],
                      mean_annual_fire_gCm2day = output_annual[,,6],
                      harvest_gCm2day = output[,,7],
                      mean_harvest_gCm2day = output_mean[,7],
                      mean_annual_harvest_gCm2day = output_annual[,,7],
                      # Internal fluxes
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      woodlitter_to_som_gCm2day = output[,,16],
                      mean_woodlitter_to_som_gCm2day = output_mean[,16],
                      mean_annual_woodlitter_to_som_gCm2day = output_annual[,,16],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,17],
                      mean_FIREemiss_labile_gCm2day = output_mean[,17],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,17],
                      FIRElitter_labile_gCm2day = output[,,18],
                      mean_FIRElitter_labile_gCm2day = output_mean[,18],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,18],
                      FIREemiss_foliage_gCm2day = output[,,19],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,19],
                      FIRElitter_foliage_gCm2day = output[,,20],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,20],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,20],
                      FIREemiss_roots_gCm2day = output[,,21],
                      mean_FIREemiss_roots_gCm2day = output_mean[,21],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,21],
                      FIRElitter_roots_gCm2day = output[,,22],
                      mean_FIRElitter_roots_gCm2day = output_mean[,22],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,22],
                      FIREemiss_wood_gCm2day = output[,,23],
                      mean_FIREemiss_wood_gCm2day = output_mean[,23],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,23],
                      FIRElitter_wood_gCm2day = output[,,24],
                      mean_FIRElitter_wood_gCm2day = output_mean[,24],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,24],
                      FIREemiss_litter_gCm2day = output[,,25],
                      mean_FIREemiss_litter_gCm2day = output_mean[,25],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,25],
                      FIRElitter_litter_gCm2day = output[,,26],
                      mean_FIRElitter_litter_gCm2day = output_mean[,26],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,26],
                      FIREemiss_woodlitter_gCm2day = output[,,27],
                      mean_FIREemiss_woodlitter_gCm2day = output_mean[,27],
                      mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,27],
                      FIRElitter_woodlitter_gCm2day = output[,,28],
                      mean_FIRElitter_woodlitter_gCm2day = output_mean[,28],
                      mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,28],
                      FIREemiss_som_gCm2day = output[,,29],
                      mean_FIREemiss_som_gCm2day = output_mean[,29],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,29],
                      HARVESTextracted_labile_gCm2day = output[,,30],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,30],
                      HARVESTextracted_foliage_gCm2day = output[,,31],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,31],
                      HARVESTextracted_roots_gCm2day = output[,,32],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,32],
                      HARVESTextracted_wood_gCm2day = output[,,33],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,33],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,33],
                      HARVESTextracted_litter_gCm2day = output[,,34],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,34],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,34],
                      HARVESTextracted_woodlitter_gCm2day = output[,,35],
                      mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,35],
                      mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,35],
                      HARVESTextracted_som_gCm2day = output[,,36],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,36],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,36],
                      HARVESTlitter_labile_gCm2day = output[,,37],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,37],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,37],
                      HARVESTlitter_foliage_gCm2day = output[,,38],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,38],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,38],
                      HARVESTlitter_roots_gCm2day = output[,,39],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,39],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,39],
                      HARVESTlitter_wood_gCm2day = output[,,40],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,40],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,40],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,41],
                      mean_labile_gCm2 = output_mean[,41],
                      mean_annual_labile_gCm2 = output_annual[,,41],
                      foliage_gCm2 = output[,,42],
                      mean_foliage_gCm2 = output_mean[,42],
                      mean_annual_foliage_gCm2 = output_annual[,,42],
                      roots_gCm2 = output[,,43],
                      mean_roots_gCm2 = output_mean[,43],
                      mean_annual_roots_gCm2 = output_annual[,,43],
                      wood_gCm2 = output[,,44],
                      mean_wood_gCm2 = output_mean[,44],
                      mean_annual_wood_gCm2 = output_annual[,,44],
                      litter_gCm2 = output[,,45],
                      mean_litter_gCm2 = output_mean[,45],
                      mean_annual_litter_gCm2 = output_annual[,,45],
                      woodlitter_gCm2 = output[,,46],
                      mean_woodlitter_gCm2 = output_mean[,46],
                      mean_annual_woodlitter_gCm2 = output_annual[,,46],
                      som_gCm2 = output[,,47],
                      mean_som_gCm2 = output_mean[,47],
                      mean_annual_som_gCm2 = output_annual[,,47],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,48],
                      mean_ET_kgH2Om2day = output_mean[,48],
                      mean_annual_ET_kgH2Om2day = output_annual[,,48],
                      SurfWater_kgH2Om2 = output[,,49],
                      mean_SurfWater_kgH2Om2 = output_mean[,49],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,49],
                      wSWP_MPa = output[,,50],
                      mean_wSWP_MPa = output_mean[,50],
                      mean_annual_wSWP_MPa = output_annual[,,50],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,51],
                      mean_lai_m2m2 = output_mean[,51],
                      mean_annual_lai_m2m2 = output_annual[,,51],
                      cgi = output[,,52],
                      mean_cgi = output_mean[,52],
                      mean_annual_cgi = output_annual[,,52],
                      cmi = output[,,53],
                      mean_cmi = output_mean[,53],
                      mean_annual_cmi = output_annual[,,53],
                      ncce_gCm2day = output[,,54],
                      mean_ncce_gCm2day = output_mean[,54],
                      mean_annual_ncce_gCm2day = output_annual[,,54],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,55],
                      mean_gs_demand_supply_ratio = output_mean[,55],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,55],
                      gs_mmolH2Om2s = output[,,56],
                      mean_gs_mmolH2Om2s = output_mean[,56],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,56],
                      APAR_MJm2day = output[,,57],
                      mean_APAR_MJm2day = output_mean[,57],
                      mean_annual_APAR_MJm2day = output_annual[,,57],
                      gb_mmolH2Om2s = output[,,58],
                      mean_gb_mmolH2Om2s = output_mean[,58],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,58],
                      CiCa = output[,,59],
                      mean_CiCa = output_mean[,59],
                      mean_annual_CiCa = output_annual[,,59],
                      # Misc
                      RootDepth_m = output[,,60],
                      mean_RootDepth_m = output_mean[,60],
                      mean_annual_RootDepth_m = output_annual[,,60],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_woodlitter_years = MTT_years[,6],
                      MTT_som_years = MTT_years[,7],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_woodlitter_gCm2 = SS_gCm2[,6],
                      SS_som_gCm2 = SS_gCm2[,7])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P10.R2.#") {
# Interface needs updating
      output_dim=28
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalec26",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met)),pars=as.double(pars_in)
                              ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,MTT_dim,noyears)))
                              ,out_var6=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site])
                              ,nomet=as.integer(dim(met)[2]),nofluxes=as.integer(PROJECT$model$nofluxes[site])
                              ,nopools=as.integer(PROJECT$model$nopools[site]),pft=as.integer(pft)
                              ,pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                              ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3)))
                              ,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
      output = tmp$out_var   ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      aNPP = tmp$out_var2    ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3     ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var4 ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5    ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6  ; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], litter_gCm2 = output[,,10],
                      labile_gCm2 = output[,,11], foliage_gCm2 = output[,,12],
                      harvest_gCm2day = output[,,13], gsi = output[,,14],
                      gsi_itemp = output[,,15], gsi_iphoto = output[,,16],
                      gsi_ivpd = output[,,17], ET_kgH2Om2day = output[,,18],
                      SurfWater_kgH2Om2 = output[,,19], wSWP_MPa = output[,,20],
                      woodlitter_gCm2 = output[,,21], fire_gCm2day = output[,,23],
                      gs_demand_supply = output[,,24], gs_total_canopy = output[,,25],
                      APAR_MJm2day = output[,,26], gb_total_canopy = output[,,27],
                      CiCa = output[,,28],
                      aNPP = aNPP, MTT = MTT, SS_gCm2 = SS_gCm2, aMTT = aMTT, natMTT = MTTnat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
      # Final tidy
      rm(output,MTT_gCm2,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H1.P4.R2.#") {
    output_dim = 59 ; MTT_dim = 7 ; SS_dim = 7
    # Load the required dalec shared object
    dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
    tmp=.Fortran( "rdalec10",output_dim=as.integer(output_dim)
                            ,aNPP_dim=as.integer(aNPP_dim)
                            ,MTT_dim=as.integer(MTT_dim)
                            ,SS_dim = as.integer(SS_dim)
                            ,fire_dim=as.integer(fire_dim)
                            ,met=as.double(t(met))
                            ,pars=as.double(pars_in)
                            ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                            ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                            ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                            ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                            ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                           
                            ,lat=as.double(lat)
                            ,nopars=as.integer(PROJECT$model$nopars[site])
                            ,nomet=as.integer(dim(met)[2])
                            ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                            ,nopools=as.integer(PROJECT$model$nopools[site])
                            ,pft=as.integer(pft)
                            ,nodays=as.integer(dim(met)[1])
                            ,nos_years=as.integer(noyears)
                            ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                            ,nos_iter=as.integer(nos_iter))
    # Extract the different output variables
    output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
    output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    mean_gpp_gCm2day = output_mean[,1],
                    mean_annual_gpp_gCm2day = output_annual[,,1],
                    rauto_gCm2day = output[,,2],
                    mean_rauto_gCm2day = output_mean[,2],
                    mean_annual_rauto_gCm2day = output_annual[,,2],
                    rg_foliage_gCm2day = output[,,3],
                    mean_rg_foliage_gCm2day = output_mean[,3],
                    mean_annual_rg_foliage_gCm2day = output_annual[,,3],
                    rhet_litter_gCm2day = output[,,4],
                    mean_rhet_litter_gCm2day = output_mean[,4],
                    mean_annual_rhet_litter_gCm2day = output_annual[,,4],
                    rhet_som_gCm2day = output[,,5],
                    mean_rhet_som_gCm2day = output_mean[,5],
                    mean_annual_rhet_som_gCm2day = output_annual[,,5],
                    rhet_woodlitter_gCm2day = output[,,6],
                    mean_rhet_woodlitter_gCm2day = output_mean[,6],
                    mean_annual_rhet_woodlitter_gCm2day = output_annual[,,6],
                    fire_gCm2day = output[,,7],
                    mean_fire_gCm2day = output_mean[,7],
                    mean_annual_fire_gCm2day = output_annual[,,7],
                    harvest_gCm2day = output[,,8],
                    mean_harvest_gCm2day = output_mean[,8],
                    mean_annual_harvest_gCm2day = output_annual[,,8],
                    # Internal fluxes
                    alloc_labile_gCm2day = output[,,9],
                    mean_alloc_labile_gCm2day = output_mean[,9],
                    mean_annual_alloc_labile_gCm2day = output_annual[,,9],
                    alloc_roots_gCm2day = output[,,10],
                    mean_alloc_roots_gCm2day = output_mean[,10],
                    mean_annual_alloc_roots_gCm2day = output_annual[,,10],
                    alloc_wood_gCm2day = output[,,11],
                    mean_alloc_wood_gCm2day = output_mean[,11],
                    mean_annual_alloc_wood_gCm2day = output_annual[,,11],
                    labile_to_foliage_gCm2day = output[,,12],
                    mean_labile_to_foliage_gCm2day = output_mean[,12],
                    mean_annual_labile_to_foliage_gCm2day = output_annual[,,12],
                    foliage_to_litter_gCm2day = output[,,13],
                    mean_foliage_to_litter_gCm2day = output_mean[,13],
                    mean_annual_foliage_to_litter_gCm2day = output_annual[,,13],
                    roots_to_litter_gCm2day = output[,,14],
                    mean_roots_to_litter_gCm2day = output_mean[,14],
                    mean_annual_roots_to_litter_gCm2day = output_annual[,,14],
                    wood_to_litter_gCm2day = output[,,15],
                    mean_wood_to_litter_gCm2day = output_mean[,15],
                    mean_annual_wood_to_litter_gCm2day = output_annual[,,15],
                    litter_to_som_gCm2day = output[,,16],
                    mean_litter_to_som_gCm2day = output_mean[,16],
                    mean_annual_litter_to_som_gCm2day = output_annual[,,16],
                    woodlitter_to_som_gCm2day = output[,,17],
                    mean_woodlitter_to_som_gCm2day = output_mean[,17],
                    mean_annual_woodlitter_to_som_gCm2day = output_annual[,,17],
                    # Disturbance fluxes
                    FIREemiss_labile_gCm2day = output[,,18],
                    mean_FIREemiss_labile_gCm2day = output_mean[,18],
                    mean_annual_FIREemiss_labile_gCm2day = output_annual[,,18],
                    FIRElitter_labile_gCm2day = output[,,19],
                    mean_FIRElitter_labile_gCm2day = output_mean[,19],
                    mean_annual_FIRElitter_labile_gCm2day = output_annual[,,19],
                    FIREemiss_foliage_gCm2day = output[,,20],
                    mean_FIREemiss_foliage_gCm2day = output_mean[,20],
                    mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,20],
                    FIRElitter_foliage_gCm2day = output[,,21],
                    mean_FIRElitter_foliage_gCm2day = output_mean[,21],
                    mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,21],
                    FIREemiss_roots_gCm2day = output[,,22],
                    mean_FIREemiss_roots_gCm2day = output_mean[,22],
                    mean_annual_FIREemiss_roots_gCm2day = output_annual[,,22],
                    FIRElitter_roots_gCm2day = output[,,23],
                    mean_FIRElitter_roots_gCm2day = output_mean[,23],
                    mean_annual_FIRElitter_roots_gCm2day = output_annual[,,23],
                    FIREemiss_wood_gCm2day = output[,,24],
                    mean_FIREemiss_wood_gCm2day = output_mean[,24],
                    mean_annual_FIREemiss_wood_gCm2day = output_annual[,,24],
                    FIRElitter_wood_gCm2day = output[,,25],
                    mean_FIRElitter_wood_gCm2day = output_mean[,25],
                    mean_annual_FIRElitter_wood_gCm2day = output_annual[,,25],
                    FIREemiss_litter_gCm2day = output[,,26],
                    mean_FIREemiss_litter_gCm2day = output_mean[,26],
                    mean_annual_FIREemiss_litter_gCm2day = output_annual[,,26],
                    FIRElitter_litter_gCm2day = output[,,27],
                    mean_FIRElitter_litter_gCm2day = output_mean[,27],
                    mean_annual_FIRElitter_litter_gCm2day = output_annual[,,27],
                    FIREemiss_woodlitter_gCm2day = output[,,28],
                    mean_FIREemiss_woodlitter_gCm2day = output_mean[,28],
                    mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,28],
                    FIRElitter_woodlitter_gCm2day = output[,,29],
                    mean_FIRElitter_woodlitter_gCm2day = output_mean[,29],
                    mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,29],
                    FIREemiss_som_gCm2day = output[,,30],
                    mean_FIREemiss_som_gCm2day = output_mean[,30],
                    mean_annual_FIREemiss_som_gCm2day = output_annual[,,30],
                    HARVESTextracted_labile_gCm2day = output[,,31],
                    mean_HARVESTextracted_labile_gCm2day = output_mean[,31],
                    mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,31],
                    HARVESTextracted_foliage_gCm2day = output[,,32],
                    mean_HARVESTextracted_foliage_gCm2day = output_mean[,32],
                    mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,32],
                    HARVESTextracted_roots_gCm2day = output[,,33],
                    mean_HARVESTextracted_roots_gCm2day = output_mean[,33],
                    mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,33],
                    HARVESTextracted_wood_gCm2day = output[,,34],
                    mean_HARVESTextracted_wood_gCm2day = output_mean[,34],
                    mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,34],
                    HARVESTextracted_litter_gCm2day = output[,,35],
                    mean_HARVESTextracted_litter_gCm2day = output_mean[,35],
                    mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,35],
                    HARVESTextracted_woodlitter_gCm2day = output[,,36],
                    mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,36],
                    mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,36],
                    HARVESTextracted_som_gCm2day = output[,,37],
                    mean_HARVESTextracted_som_gCm2day = output_mean[,37],
                    mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,37],
                    HARVESTlitter_labile_gCm2day = output[,,38],
                    mean_HARVESTlitter_labile_gCm2day = output_mean[,38],
                    mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,38],
                    HARVESTlitter_foliage_gCm2day = output[,,39],
                    mean_HARVESTlitter_foliage_gCm2day = output_mean[,39],
                    mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,39],
                    HARVESTlitter_roots_gCm2day = output[,,40],
                    mean_HARVESTlitter_roots_gCm2day = output_mean[,40],
                    mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,40],
                    HARVESTlitter_wood_gCm2day = output[,,41],
                    mean_HARVESTlitter_wood_gCm2day = output_mean[,41],
                    mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,41],
                    # C pools (gC/m2)
                    labile_gCm2 = output[,,42],
                    mean_labile_gCm2 = output_mean[,42],
                    mean_annual_labile_gCm2 = output_annual[,,42],
                    foliage_gCm2 = output[,,43],
                    mean_foliage_gCm2 = output_mean[,43],
                    mean_annual_foliage_gCm2 = output_annual[,,43],
                    roots_gCm2 = output[,,44],
                    mean_roots_gCm2 = output_mean[,44],
                    mean_annual_roots_gCm2 = output_annual[,,44],
                    wood_gCm2 = output[,,45],
                    mean_wood_gCm2 = output_mean[,45],
                    mean_annual_wood_gCm2 = output_annual[,,45],
                    litter_gCm2 = output[,,46],
                    mean_litter_gCm2 = output_mean[,46],
                    mean_annual_litter_gCm2 = output_annual[,,46],
                    woodlitter_gCm2 = output[,,47],
                    mean_woodlitter_gCm2 = output_mean[,47],
                    mean_annual_woodlitter_gCm2 = output_annual[,,47],
                    som_gCm2 = output[,,48],
                    mean_som_gCm2 = output_mean[,48],
                    mean_annual_som_gCm2 = output_annual[,,48],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,49],
                    mean_lai_m2m2 = output_mean[,49],
                    mean_annual_lai_m2m2 = output_annual[,,49],
                    gsi = output[,,50],
                    mean_gsi = output_mean[,50],
                    mean_annual_gsi = output_annual[,,50],
                    gsi_itemp = output[,,51],
                    mean_gsi_itemp = output_mean[,51],
                    mean_annual_gsi_itemp = output_annual[,,51],
                    gsi_iphoto = output[,,52],
                    mean_gsi_iphoto = output_mean[,52],
                    mean_annual_gsi_iphoto = output_annual[,,52],
                    gsi_ivpd = output[,,53],
                    mean_gsi_ivpd = output_mean[,53],
                    mean_annual_gsi_ivpd = output_annual[,,53],
                    # Photosynthesis / C~water coupling related
                    gs_demand_supply_ratio = output[,,54],
                    mean_gs_demand_supply_ratio = output_mean[,54],
                    mean_annual_gs_demand_supply_ratio = output_annual[,,54],
                    gs_mmolH2Om2s = output[,,55],
                    mean_gs_mmolH2Om2s = output_mean[,55],
                    mean_annual_gs_mmolH2Om2s = output_annual[,,55],
                    APAR_MJm2day = output[,,56],
                    mean_APAR_MJm2day = output_mean[,56],
                    mean_annual_APAR_MJm2day = output_annual[,,56],
                    gb_mmolH2Om2s = output[,,57],
                    mean_gb_mmolH2Om2s = output_mean[,57],
                    mean_annual_gb_mmolH2Om2s = output_annual[,,57],
                    CiCa = output[,,58],
                    mean_CiCa = output_mean[,58],
                    mean_annual_CiCa = output_annual[,,58],
                    # Misc
                    RootDepth_m = output[,,59],
                    mean_RootDepth_m = output_mean[,59],
                    mean_annual_RootDepth_m = output_annual[,,59],
                    ## Aggregated variables
                    # Mean Transit times
                    MTT_labile_years = MTT_years[,1],
                    MTT_foliage_years = MTT_years[,2],
                    MTT_roots_years = MTT_years[,3],
                    MTT_wood_years = MTT_years[,4],
                    MTT_litter_years = MTT_years[,5],
                    MTT_woodlitter_years = MTT_years[,6],
                    MTT_som_years = MTT_years[,7],
                    # Steady state estimates
                    SS_labile_gCm2 = SS_gCm2[,1],
                    SS_foliage_gCm2 = SS_gCm2[,2],
                    SS_roots_gCm2 = SS_gCm2[,3],
                    SS_wood_gCm2 = SS_gCm2[,4],
                    SS_litter_gCm2 = SS_gCm2[,5],
                    SS_woodlitter_gCm2 = SS_gCm2[,6],
                    SS_som_gCm2 = SS_gCm2[,7])
    # Determine the NPP fraction of expressed NPP
    # i.e. actual growth not GPP-Ra
    NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                         states_all$alloc_roots_gCm2day +
                         states_all$alloc_wood_gCm2day,1,mean)
    NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day,1,mean),
                         apply(states_all$alloc_roots_gCm2day,1,mean),
                         apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
    states_all$NPP_foliage_fraction = NPP_fraction[,1]
    states_all$NPP_roots_fraction = NPP_fraction[,2]
    states_all$NPP_wood_fraction = NPP_fraction[,3]
    # Tidy up variables
    rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P3.R1.#") {
      output_dim = 62 ; MTT_dim = 7 ; SS_dim = 7
      # Load the required dalec shared object
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec9",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,fire_dim=as.integer(fire_dim)
                             ,met=as.double(t(met)),pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                            
                             ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site])
                             ,nomet=as.integer(dim(met)[2]),nofluxes=as.integer(PROJECT$model$nofluxes[site])
                             ,nopools=as.integer(PROJECT$model$nopools[site]),pft=as.integer(pft)
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                             ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3))))
      # Extract the different output variables
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rg_foliage_gCm2day = output[,,3],
                      mean_rg_foliage_gCm2day = output_mean[,3],
                      mean_annual_rg_foliage_gCm2day = output_annual[,,3],
                      rhet_litter_gCm2day = output[,,4],
                      mean_rhet_litter_gCm2day = output_mean[,4],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,4],
                      rhet_som_gCm2day = output[,,5],
                      mean_rhet_som_gCm2day = output_mean[,5],
                      mean_annual_rhet_som_gCm2day = output_annual[,,5],
                      rhet_woodlitter_gCm2day = output[,,6],
                      mean_rhet_woodlitter_gCm2day = output_mean[,6],
                      mean_annual_rhet_woodlitter_gCm2day = output_annual[,,6],
                      fire_gCm2day = output[,,7],
                      mean_fire_gCm2day = output_mean[,7],
                      mean_annual_fire_gCm2day = output_annual[,,7],
                      harvest_gCm2day = output[,,8],
                      mean_harvest_gCm2day = output_mean[,8],
                      mean_annual_harvest_gCm2day = output_annual[,,8],
                      # Internal fluxes
                      alloc_labile_gCm2day = output[,,9],
                      mean_alloc_labile_gCm2day = output_mean[,9],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,9],
                      alloc_roots_gCm2day = output[,,10],
                      mean_alloc_roots_gCm2day = output_mean[,10],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,10],
                      alloc_wood_gCm2day = output[,,11],
                      mean_alloc_wood_gCm2day = output_mean[,11],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,11],
                      labile_to_foliage_gCm2day = output[,,12],
                      mean_labile_to_foliage_gCm2day = output_mean[,12],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,12],
                      foliage_to_litter_gCm2day = output[,,13],
                      mean_foliage_to_litter_gCm2day = output_mean[,13],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,13],
                      roots_to_litter_gCm2day = output[,,14],
                      mean_roots_to_litter_gCm2day = output_mean[,14],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,14],
                      wood_to_litter_gCm2day = output[,,15],
                      mean_wood_to_litter_gCm2day = output_mean[,15],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,15],
                      litter_to_som_gCm2day = output[,,16],
                      mean_litter_to_som_gCm2day = output_mean[,16],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,16],
                      woodlitter_to_som_gCm2day = output[,,17],
                      mean_woodlitter_to_som_gCm2day = output_mean[,17],
                      mean_annual_woodlitter_to_som_gCm2day = output_annual[,,17],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,18],
                      mean_FIREemiss_labile_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,18],
                      FIRElitter_labile_gCm2day = output[,,19],
                      mean_FIRElitter_labile_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,19],
                      FIREemiss_foliage_gCm2day = output[,,20],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,20],
                      FIRElitter_foliage_gCm2day = output[,,21],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,21],
                      FIREemiss_roots_gCm2day = output[,,22],
                      mean_FIREemiss_roots_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,22],
                      FIRElitter_roots_gCm2day = output[,,23],
                      mean_FIRElitter_roots_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,23],
                      FIREemiss_wood_gCm2day = output[,,24],
                      mean_FIREemiss_wood_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,24],
                      FIRElitter_wood_gCm2day = output[,,25],
                      mean_FIRElitter_wood_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,25],
                      FIREemiss_litter_gCm2day = output[,,26],
                      mean_FIREemiss_litter_gCm2day = output_mean[,26],
                      Fmean_annual_IREemiss_litter_gCm2day = output_annual[,,26],
                      FIRElitter_litter_gCm2day = output[,,27],
                      mean_FIRElitter_litter_gCm2day = output_mean[,27],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,27],
                      FIREemiss_woodlitter_gCm2day = output[,,28],
                      mean_FIREemiss_woodlitter_gCm2day = output_mean[,28],
                      mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,28],
                      FIRElitter_woodlitter_gCm2day = output[,,29],
                      mean_FIRElitter_woodlitter_gCm2day = output_mean[,29],
                      mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,29],
                      FIREemiss_som_gCm2day = output[,,30],
                      mean_FIREemiss_som_gCm2day = output_mean[,30],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,30],
                      HARVESTextracted_labile_gCm2day = output[,,31],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,31],
                      HARVESTextracted_foliage_gCm2day = output[,,32],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,32],
                      HARVESTextracted_roots_gCm2day = output[,,33],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,33],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,33],
                      HARVESTextracted_wood_gCm2day = output[,,34],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,34],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,34],
                      HARVESTextracted_litter_gCm2day = output[,,35],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,35],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,35],
                      HARVESTextracted_woodlitter_gCm2day = output[,,36],
                      mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,36],
                      mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,36],
                      HARVESTextracted_som_gCm2day = output[,,37],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,37],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,37],
                      HARVESTlitter_labile_gCm2day = output[,,38],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,38],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,38],
                      HARVESTlitter_foliage_gCm2day = output[,,39],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,39],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,39],
                      HARVESTlitter_roots_gCm2day = output[,,40],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,40],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,40],
                      HARVESTlitter_wood_gCm2day = output[,,41],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,41],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,41],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,42],
                      mean_labile_gCm2 = output_mean[,42],
                      mean_annual_labile_gCm2 = output_annual[,,42],
                      foliage_gCm2 = output[,,43],
                      mean_foliage_gCm2 = output_mean[,43],
                      mean_annual_foliage_gCm2 = output_annual[,,43],
                      roots_gCm2 = output[,,44],
                      mean_roots_gCm2 = output_mean[,44],
                      mean_annual_roots_gCm2 = output_annual[,,44],
                      wood_gCm2 = output[,,45],
                      mean_wood_gCm2 = output_mean[,45],
                      mean_annual_wood_gCm2 = output_annual[,,45],
                      litter_gCm2 = output[,,46],
                      mean_litter_gCm2 = output_mean[,46],
                      mean_annual_litter_gCm2 = output_annual[,,46],
                      woodlitter_gCm2 = output[,,47],
                      mean_woodlitter_gCm2 = output_mean[,47],
                      mean_annual_woodlitter_gCm2 = output_annual[,,47],
                      som_gCm2 = output[,,48],
                      mean_som_gCm2 = output_mean[,48],
                      mean_annual_som_gCm2 = output_annual[,,48],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,49],
                      mean_ET_kgH2Om2day = output_mean[,49],
                      mean_annual_ET_kgH2Om2day = output_annual[,,49],
                      SurfWater_kgH2Om2 = output[,,50],
                      mean_SurfWater_kgH2Om2 = output_mean[,50],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,50],
                      wSWP_MPa = output[,,51],
                      mean_wSWP_MPa = output_mean[,51],
                      mean_annual_wSWP_MPa = output_annual[,,51],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,52],
                      mean_lai_m2m2 = output_mean[,52],
                      mean_annual_lai_m2m2 = output_annual[,,52],
                      gsi = output[,,53],
                      mean_gsi = output_mean[,53],
                      mean_annual_gsi = output_annual[,,53],
                      gsi_itemp = output[,,54],
                      mean_gsi_itemp = output_mean[,54],
                      mean_annual_gsi_itemp = output_annual[,,54],
                      gsi_iphoto = output[,,55],
                      mean_gsi_iphoto = output_mean[,55],
                      mean_annual_gsi_iphoto = output_annual[,,55],
                      gsi_ivpd = output[,,56],
                      mean_gsi_ivpd = output_mean[,56],
                      mean_annual_gsi_ivpd = output_annual[,,56],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,57],
                      mean_gs_demand_supply_ratio = output_mean[,57],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,57],
                      gs_mmolH2Om2s = output[,,58],
                      mean_gs_mmolH2Om2s = output_mean[,58],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,58],
                      APAR_MJm2day = output[,,59],
                      mean_APAR_MJm2day = output_mean[,59],
                      mean_annual_APAR_MJm2day = output_annual[,,59],
                      gb_mmolH2Om2s = output[,,60],
                      mean_gb_mmolH2Om2s = output_mean[,60],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,60],
                      CiCa = output[,,61],
                      mean_CiCa = output_mean[,61],
                      mean_annual_CiCa = output_annual[,,61],
                      # Misc
                      RootDepth_m = output[,,62],
                      mean_RootDepth_m = output_mean[,62],
                      mean_annual_RootDepth_m = output_annual[,,62],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_woodlitter_years = MTT_years[,6],
                      MTT_som_years = MTT_years[,7],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_woodlitter_gCm2 = SS_gCm2[,6],
                      SS_som_gCm2 = SS_gCm2[,7])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.C1.D1.F2.P1.#") {
      output_dim = 44 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec2",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                                 
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter) )
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,43],
                      mean_lai_m2m2 = output_mean[,43],
                      mean_annual_lai_m2m2 = output_annual[,,43],
                      # Photosynthesis / C~water coupling related
                      CiCa = output[,,44],
                      mean_CiCa = output_mean[,44],
                      mean_annual_CiCa = output_annual[,,44],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC_1005") {
      output_dim = 46 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec1005",output_dim=as.integer(output_dim)
                                ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                ,met=as.double(t(met))
                                ,pars=as.double(pars_in)
                                ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                                ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                       
                                ,lat=as.double(lat)
                                ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                ,nodays=as.integer(dim(met)[1])
                                ,nos_years=as.integer(noyears)
                                ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      # extract output variables
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output_mean[,11],
                      mean_labile_to_foliage_gCm2day = output[,,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      SurfWater_kgH2Om2 = output[,,44],
                      mean_SurfWater_kgH2Om2 = output_mean[,44],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,44],
                      # wSWP_MPa = output[,,45],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,45],
                      mean_lai_m2m2 = output_mean[,45],
                      mean_annual_lai_m2m2 = output_annual[,,45],
                      # Photosynthesis / C~water coupling related
                      CiCa = output[,,46],
                      mean_CiCa = output_mean[,46],
                      mean_annual_CiCa = output_annual[,,46],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC_1005a") {
      output_dim = 46 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      # if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalec1005a",output_dim=as.integer(output_dim)
                                 ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                 ,met=as.double(t(met))
                                 ,pars=as.double(pars_in)
                                 ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                 ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                 ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                 ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                                 ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                        
                                 ,lat=as.double(lat)
                                 ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                 ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                 ,nodays=as.integer(dim(met)[1])
                                 ,nos_years=as.integer(noyears)
                                 ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      # extract output variables
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output_mean[,11],
                      mean_labile_to_foliage_gCm2day = output[,,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      SurfWater_kgH2Om2 = output[,,44],
                      mean_SurfWater_kgH2Om2 = output_mean[,44],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,44],
                      # wSWP_MPa = output[,,45],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,45],
                      mean_lai_m2m2 = output_mean[,45],
                      mean_annual_lai_m2m2 = output_annual[,,45],
                      # Photosynthesis / C~water coupling related
                      CiCa = output[,,46],
                      mean_CiCa = output_mean[,46],
                      mean_annual_CiCa = output_annual[,,46],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # add newly calculated variables
      # states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Final tidy
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C1.D2.F2.H1.P1.#") {
      output_dim = 48 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec3",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                    
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                             ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3))))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output_mean[,11],
                      mean_labile_to_foliage_gCm2day = output[,,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter__annualgCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,43],
                      mean_lai_m2m2 = output_mean[,43],
                      mean_annual_lai_m2m2 = output_annual[,,43],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,44],
                      mean_gs_demand_supply_ratio = output_mean[,44],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,44],
                      gs_mmolH2Om2s = output[,,45],
                      mean_gs_mmolH2Om2s = output_mean[,45],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,45],
                      APAR_MJm2day = output[,,46],
                      mean_APAR_MJm2day = output_mean[,46],
                      mean_annual_APAR_MJm2day = output_annual[,,46],
                      gb_mmolH2Om2s = output[,,47],
                      mean_gb_mmolH2Om2s = output_mean[,47],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,47],
                      CiCa = output[,,48],
                      mean_CiCa = output_mean[,48],
                      mean_annual_CiCa = output_annual[,,48],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPPHARVESTextracted_roots_gCm2day = output[,,29],
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C1.D2.F2.H2.P1.#") {
      output_dim = 58 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec4",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                    
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                             ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2     ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3       ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      Etrans_kgH2Om2day = output[,,44],
                      mean_Etrans_kgH2Om2day = output_mean[,44],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,44],
                      Esoil_kgH2Om2day = output[,,45],
                      mean_Esoil_kgH2Om2day = output_mean[,45],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,45],
                      Ewetcanopy_kgH2Om2day = output[,,46],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,46],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,46],
                      runoff_kgH2Om2day = output[,,47],
                      mean_runoff_kgH2Om2day = output_mean[,47],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,47],
                      underflow_kgH2Om2day = output[,,48],
                      mean_underflow_kgH2Om2day = output_mean[,48],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,48],
                      SurfWater_kgH2Om2 = output[,,49],
                      mean_SurfWater_kgH2Om2 = output_mean[,49],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,49],
                      wSWP_MPa = output[,,50],
                      mean_wSWP_MPa = output_mean[,50],
                      mean_annual_wSWP_MPa = output_annual[,,50],
                      snow_kgH2Om2 = output[,,51],
                      mean_snow_kgH2Om2 = output_mean[,51],
                      mean_annual_snow_kgH2Om2 = output_annual[,,51],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,52],
                      mean_lai_m2m2 = output_mean[,52],
                      mean_annual_lai_m2m2 = output_annual[,,52],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,53],
                      mean_gs_demand_supply_ratio = output_mean[,53],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,53],
                      gs_mmolH2Om2s = output[,,54],
                      mean_gs_mmolH2Om2s = output_mean[,54],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,54],
                      APAR_MJm2day = output[,,55],
                      mean_APAR_MJm2day = output_mean[,55],
                      mean_annual_APAR_MJm2day = output_annual[,,55],
                      gb_mmolH2Om2s = output[,,56],
                      mean_gb_mmolH2Om2s = output_mean[,56],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,56],
                      CiCa = output[,,57],
                      mean_CiCa = output_mean[,57],
                      mean_annual_CiCa = output_annual[,,57],
                      # Misc
                      RootDepth_m = output[,,58],
                      mean_RootDepth_m = output_mean[,58],
                      mean_annual_RootDepth_m = output_annual[,,58],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C1.D2.F2.H3.P1.#") {
      output_dim = 59 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec29",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                    
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                             ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],  
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output[,,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mea_meann[,,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      Etrans_kgH2Om2day = output[,,44],
                      mean_Etrans_kgH2Om2day = output_mean[,44],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,44],
                      Esoil_kgH2Om2day = output[,,45],
                      mean_Esoil_kgH2Om2day = output_mean[,45],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,45],
                      Ewetcanopy_kgH2Om2day = output[,,46],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,46],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,46],
                      runoff_kgH2Om2day = output[,,47],
                      mean_runoff_kgH2Om2day = output_mean[,47],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,47],
                      underflow_kgH2Om2day = output[,,48],
                      mean_underflow_kgH2Om2day = output[,,48],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,48],
                      SurfWater_kgH2Om2 = output[,,49],
                      mean_SurfWater_kgH2Om2 = output_mean[,49],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,49],
                      wSWP_MPa = output[,,50],
                      mean_wSWP_MPa = output_mean[,50],
                      mean_annual_wSWP_MPa = output_annual[,,50],
                      snow_kgH2Om2 = output[,,51],
                      mean_snow_kgH2Om2 = output_mean[,51],
                      mean_annual_snow_kgH2Om2 = output_annual[,,51],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,52],
                      mean_lai_m2m2 = output_mean[,52],
                      mean_annual_lai_m2m2 = output_annual[,,52],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,53],
                      mean_gs_demand_supply_ratio = output_mean[,53],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,53],
                      gs_mmolH2Om2s = output[,,54],
                      mean_gs_mmolH2Om2s = output_mean[,54],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,54],
                      APAR_MJm2day = output[,,55],
                      mean_APAR_MJm2day = output_mean[,55],
                      mean_annual_APAR_MJm2day = output_annual[,,55],
                      gb_mmolH2Om2s = output[,,56],
                      mean_gb_mmolH2Om2s = output_mean[,56],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,56],
                      CiCa = output[,,57],
                      mean_CiCa = output_mean[,57],
                      mean_annual_CiCa = output_annual[,,57],
                      # Misc
                      RootDepth_m = output[,,58],
                      mean_RootDepth_m = output_mean[,58],
                      mean_annual_RootDepth_m = output_annual[,,58],
                      LWP_MPa = output[,,59],
                      mean_LWP_MPa = output_mean[,59],
                      mean_annual_LWP_MPa = output_annual[,,59],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A2.C1.D2.F2.H2.P1.#") {
      output_dim = 58 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec20",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                     
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],  
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output[,,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      Etrans_kgH2Om2day = output[,,44],
                      mean_Etrans_kgH2Om2day = output_mean[,44],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,44],
                      Esoil_kgH2Om2day = output[,,45],
                      mean_Esoil_kgH2Om2day = output_mean[,45],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,45],
                      Ewetcanopy_kgH2Om2day = output[,,46],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,46],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,46],
                      runoff_kgH2Om2day = output[,,47],
                      mean_runoff_kgH2Om2day = output_mean[,47],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,47],
                      underflow_kgH2Om2day = output[,,48],
                      mean_underflow_kgH2Om2day = output_mean[,48],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,48],
                      SurfWater_kgH2Om2 = output[,,49],
                      mean_SurfWater_kgH2Om2 = output_mean[,49],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,49],
                      wSWP_MPa = output[,,50],
                      mean_wSWP_MPa = output_mean[,50],
                      mean_annual_wSWP_MPa = output_annual[,,50],
                      snow_kgH2Om2 = output[,,51],
                      mean_snow_kgH2Om2 = output_mean[,51],
                      mean_annual_snow_kgH2Om2 = output_annual[,,51],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,52],
                      mean_lai_m2m2 = output_mean[,52],
                      mean_annual_lai_m2m2 = output_annual[,,52],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,53],
                      mean_gs_demand_supply_ratio = output_mean[,53],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,53],
                      gs_mmolH2Om2s = output[,,54],
                      mean_gs_mmolH2Om2s = output_mean[,54],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,54],
                      APAR_MJm2day = output[,,55],
                      mean_APAR_MJm2day = output_mean[,55],
                      mean_annual_APAR_MJm2day = output_annual[,,55],
                      gb_mmolH2Om2s = output[,,56],
                      mean_gb_mmolH2Om2s = output_mean[,56],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,56],
                      CiCa = output[,,57],
                      mean_CiCa = output_mean[,57],
                      mean_annual_CiCa = output_annual[,,57],
                      # Misc
                      RootDepth_m = output[,,58],
                      mean_RootDepth_m = output_mean[,58],
                      mean_annual_RootDepth_m = output_annual[,,58],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A4.C6.D2.F2.H2.P11.#") {
      output_dim = 58 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec31",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,nos_years = as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      MTT_years = tmp$out_var2     ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3       ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      Etrans_kgH2Om2day = output[,,44],
                      mean_Etrans_kgH2Om2day = output_mean[,44],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,44],
                      Esoil_kgH2Om2day = output[,,45],
                      mean_Esoil_kgH2Om2day = output_mean[,45],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,45],
                      Ewetcanopy_kgH2Om2day = output[,,46],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,46],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,46],
                      runoff_kgH2Om2day = output[,,47],
                      mean_runoff_kgH2Om2day = output_mean[,47],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,47],
                      underflow_kgH2Om2day = output[,,48],
                      mean_underflow_kgH2Om2day = output_mean[,48],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,48],
                      SurfWater_kgH2Om2 = output[,,49],
                      mean_SurfWater_kgH2Om2 = output_mean[,49],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,49],
                      wSWP_MPa = output[,,50],
                      mean_wSWP_MPa = output_mean[,50],
                      mean_annual_wSWP_MPa = output_annual[,,50],
                      snow_kgH2Om2 = output[,,51],
                      mean_snow_kgH2Om2 = output_mean[,51],
                      mean_annual_snow_kgH2Om2 = output_annual[,,51],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,52],
                      mean_lai_m2m2 = output_mean[,52],
                      mean_annual_lai_m2m2 = output_annual[,,52],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,53],
                      mean_gs_demand_supply_ratio = output_mean[,53],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,53],
                      gs_mmolH2Om2s = output[,,54],
                      mean_gs_mmolH2Om2s = output_mean[,54],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,54],
                      APAR_MJm2day = output[,,55],
                      mean_APAR_MJm2day = output_mean[,55],
                      mean_annual_APAR_MJm2day = output_annual[,,55],
                      gb_mmolH2Om2s = output[,,56],
                      mean_gb_mmolH2Om2s = output_mean[,56],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,56],
                      CiCa = output[,,57],
                      mean_CiCa = output_mean[,57],
                      mean_annual_CiCa = output_annual[,,57],
                      # Misc
                      RootDepth_m = output[,,58],
                      mean_RootDepth_m = output_mean[,58],
                      mean_annual_RootDepth_m = output_annual[,,58],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,output_mean,output_annual,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A3.C1.D2.F2.H2.P1.#") {
      output_dim = 58 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec30",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                    
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                             ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      Etrans_kgH2Om2day = output[,,44],
                      mean_Etrans_kgH2Om2day = output_mean[,44],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,44],
                      Esoil_kgH2Om2day = output[,,45],
                      mean_Esoil_kgH2Om2day = output_mean[,45],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,45],
                      Ewetcanopy_kgH2Om2day = output[,,46],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,46],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,46],
                      runoff_kgH2Om2day = output[,,47],
                      mean_runoff_kgH2Om2day = output_mean[,47],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,47],
                      underflow_kgH2Om2day = output[,,48],
                      mean_underflow_kgH2Om2day = output_mean[,48],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,48],
                      SurfWater_kgH2Om2 = output[,,49],
                      mean_SurfWater_kgH2Om2 = output_mean[,49],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,49],
                      wSWP_MPa = output[,,50],
                      mean_wSWP_MPa = output_mean[,50],
                      mean_annual_wSWP_MPa = output_annual[,,50],
                      snow_kgH2Om2 = output[,,51],
                      mean_snow_kgH2Om2 = output_mean[,51],
                      mean_annual_snow_kgH2Om2 = output_annual[,,51],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,52],
                      mean_lai_m2m2 = output_mean[,52],
                      mean_annual_lai_m2m2 = output_annual[,,52],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,53],
                      mean_gs_demand_supply_ratio = output_mean[,53],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,53],
                      gs_mmolH2Om2s = output[,,54],
                      mean_gs_mmolH2Om2s = output_mean[,54],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,54],
                      APAR_MJm2day = output[,,55],
                      mean_APAR_MJm2day = output_mean[,55],
                      mean_annual_APAR_MJm2day = output_annual[,,55],
                      gb_mmolH2Om2s = output[,,56],
                      mean_gb_mmolH2Om2s = output_mean[,56],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,56],
                      CiCa = output[,,57],
                      mean_CiCa = output_mean[,57],
                      mean_annual_CiCa = output_annual[,,57],
                      # Misc
                      RootDepth_m = output[,,58],
                      mean_RootDepth_m = output_mean[,58],
                      mean_annual_RootDepth_m = output_annual[,,58],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C1.D2.F2.H2.P2.#") {
      output_dim = 51 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec18",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                     
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                              ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3))))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      SurfWater_kgH2Om2 = output[,,44],
                      mean_SurfWater_kgH2Om2 = output_mean[,44],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,44],
                      wSWP_MPa = output[,,45],
                      mean_wSWP_MPa = output_mean[,45],
                      mean_annual_wSWP_MPa = output_annual[,,45],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,46],
                      mean_lai_m2m2 = output_mean[,46],
                      mean_annual_lai_m2m2 = output_annual[,,46],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,47],
                      mean_gs_demand_supply_ratio = output_mean[,47],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,47],
                      gs_mmolH2Om2s = output[,,48],
                      mean_gs_mmolH2Om2s = output_mean[,48],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,48],
                      APAR_MJm2day = output[,,49],
                      mean_APAR_MJm2day = output_mean[,49],
                      mean_annual_APAR_MJm2day = output_annual[,,49],
                      gb_mmolH2Om2s = output[,,50],
                      mean_gb_mmolH2Om2s = output_mean[,50],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,50],
                      CiCa = output[,,51],
                      mean_CiCa = output_mean[,51],
                      mean_annual_CiCa = output_annual[,,51],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C1.D2.F2.H2.P5.#") {
      output_dim = 52 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec21",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                     
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                              ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3))))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      SurfWater_kgH2Om2 = output[,,44],
                      mean_SurfWater_kgH2Om2 = output_mean[,44],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,44],
                      wSWP_MPa = output[,,45],
                      mean_wSWP_MPa = output_mean[,45],
                      mean_annual_wSWP_MPa = output_annual[,,45],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,46],
                      mean_lai_m2m2 = output_mean[,46],
                      mean_annual_lai_m2m2 = output_annual[,,46],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,47],
                      mean_gs_demand_supply_ratio = output_mean[,47],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,47],
                      gs_mmolH2Om2s = output[,,48],
                      mean_gs_mmolH2Om2s = output_mean[,48],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,48],
                      APAR_MJm2day = output[,,49],
                      mean_APAR_MJm2day = output_mean[,49],
                      mean_annual_APAR_MJm2day = output_annual[,,49],
                      gb_mmolH2Om2s = output[,,50],
                      mean_gb_mmolH2Om2s = output_mean[,50],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,50],
                      CiCa = output[,,51],
                      mean_CiCa = output_mean[,51],
                      mean_annual_CiCa = output_annual[,,51],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C1.D2.F2.H2.P6.#") {
      output_dim = 51 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec22",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                     
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                              ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3))))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      SurfWater_kgH2Om2 = output[,,44],
                      mean_SurfWater_kgH2Om2 = output_mean[,44],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,44],
                      wSWP_MPa = output[,,45],
                      mean_wSWP_MPa = output_mean[,45],
                      mean_annual_wSWP_MPa = output_annual[,,45],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,46],
                      mean_lai_m2m2 = output_mean[,46],
                      mean_annual_lai_m2m2 = output_annual[,,46],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,47],
                      mean_gs_demand_supply_ratio = output_mean[,47],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,47],
                      gs_mmolH2Om2s = output[,,48],
                      mean_gs_mmolH2Om2s = output_mean[,48],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,48],
                      APAR_MJm2day = output[,,49],
                      mean_APAR_MJm2day = output_mean[,49],
                      mean_annual_APAR_MJm2day = output_annual[,,49],
                      gb_mmolH2Om2s = output[,,50],
                      mean_gb_mmolH2Om2s = output_mean[,50],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,50],
                      CiCa = output[,,51],
                      mean_CiCa = output_mean[,51],
                      mean_annual_CiCa = output_annual[,,51],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C1.D2.F2.H2.P1.R1.#") {
      output_dim = 51 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec5",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                    
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                             ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3))))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      SurfWater_kgH2Om2 = output[,,44],
                      mean_SurfWater_kgH2Om2 = output_mean[,44],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,44],
                      wSWP_MPa = output[,,45],
                      mean_wSWP_MPa = output_mean[,45],
                      mean_annual_wSWP_MPa = output_annual[,,45],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,46],
                      mean_lai_m2m2 = output_mean[,46],
                      mean_annual_lai_m2m2 = output_annual[,,46],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,47],
                      mean_gs_demand_supply_ratio = output_mean[,47],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,47],
                      gs_mmolH2Om2s = output[,,48],
                      mean_gs_mmolH2Om2s = output_mean[,48],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,48],
                      APAR_MJm2day = output[,,49],
                      mean_APAR_MJm2day = output_mean[,49],
                      mean_annual_APAR_MJm2day = output_annual[,,49],
                      gb_mmolH2Om2s = output[,,50],
                      mean_gb_mmolH2Om2s = output_mean[,50],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,50],
                      CiCa = output[,,51],
                      mean_CiCa = output_mean[,51],
                      mean_annual_CiCa = output_annual[,,51],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                     states_all$alloc_foliage_gCm2day +
                     states_all$alloc_roots_gCm2day +
                     states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                     apply(states_all$alloc_roots_gCm2day,1,mean),
                     apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P1.R1.#") {
      output_dim = 57 ; MTT_dim = 7 ; SS_dim = 7
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec6",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                    
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                             ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      rhet_woodlitter_gCm2day = output[,,5],
                      mean_rhet_woodlitter_gCm2day = output_mean[,5],
                      mean_annual_rhet_woodlitter_gCm2day = output_annual[,,5],
                      fire_gCm2day = output[,,6],
                      mean_fire_gCm2day = output_mean[,6],
                      mean_annual_fire_gCm2day = output_annual[,,6],
                      harvest_gCm2day = output[,,7],
                      mean_harvest_gCm2day = output_mean[,7],
                      mean_annual_harvest_gCm2day = output_annual[,,7],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,8],
                      mean_alloc_foliage_gCm2day = output_mean[,8],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,8],
                      alloc_labile_gCm2day = output[,,9],
                      mean_alloc_labile_gCm2day = output_mean[,9],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,9],
                      alloc_roots_gCm2day = output[,,10],
                      mean_alloc_roots_gCm2day = output_mean[,10],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,10],
                      alloc_wood_gCm2day = output[,,11],
                      mean_alloc_wood_gCm2day = output_mean[,11],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,11],
                      labile_to_foliage_gCm2day = output[,,12],
                      mean_labile_to_foliage_gCm2day = output_mean[,12],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,12],
                      foliage_to_litter_gCm2day = output[,,13],
                      mean_foliage_to_litter_gCm2day = output_mean[,13],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,13],
                      roots_to_litter_gCm2day = output[,,14],
                      mean_roots_to_litter_gCm2day = output_mean[,14],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,14],
                      wood_to_litter_gCm2day = output[,,15],
                      mean_wood_to_litter_gCm2day = output_mean[,15],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,15],
                      litter_to_som_gCm2day = output[,,16],
                      mean_litter_to_som_gCm2day = output_mean[,16],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,16],
                      woodlitter_to_som_gCm2day = output[,,17],
                      mean_woodlitter_to_som_gCm2day = output_mean[,17],
                      mean_annual_woodlitter_to_som_gCm2day = output_annual[,,17],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,18],
                      mean_FIREemiss_labile_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,18],
                      FIRElitter_labile_gCm2day = output[,,19],
                      mean_FIRElitter_labile_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,19],
                      FIREemiss_foliage_gCm2day = output[,,20],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,20],
                      FIRElitter_foliage_gCm2day = output[,,21],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,21],
                      FIREemiss_roots_gCm2day = output[,,22],
                      mean_FIREemiss_roots_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,22],
                      FIRElitter_roots_gCm2day = output[,,23],
                      mean_FIRElitter_roots_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,23],
                      FIREemiss_wood_gCm2day = output[,,24],
                      mean_FIREemiss_wood_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,24],
                      FIRElitter_wood_gCm2day = output[,,25],
                      mean_FIRElitter_wood_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,25],
                      FIREemiss_litter_gCm2day = output[,,26],
                      mean_FIREemiss_litter_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,26],
                      FIRElitter_litter_gCm2day = output[,,27],
                      mean_FIRElitter_litter_gCm2day = output_mean[,27],
                      mean_annual_FIRElitter_litter_gCm2day = output[,,27],
                      FIREemiss_woodlitter_gCm2day = output[,,28],
                      mean_FIREemiss_woodlitter_gCm2day = output_mean[,28],
                      mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,28],
                      FIRElitter_woodlitter_gCm2day = output[,,29],
                      mean_FIRElitter_woodlitter_gCm2day = output_mean[,29],
                      mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,29],
                      FIREemiss_som_gCm2day = output[,,30],
                      mean_FIREemiss_som_gCm2day = output_mean[,30],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,30],
                      HARVESTextracted_labile_gCm2day = output[,,31],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,31],
                      HARVESTextracted_foliage_gCm2day = output[,,32],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,32],
                      HARVESTextracted_roots_gCm2day = output[,,33],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,33],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,33],
                      HARVESTextracted_wood_gCm2day = output[,,34],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,34],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,34],
                      HARVESTextracted_litter_gCm2day = output[,,35],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,35],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,35],
                      HARVESTextracted_woodlitter_gCm2day = output[,,36],
                      mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,36],
                      mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,36],
                      HARVESTextracted_som_gCm2day = output[,,37],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,37],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,37],
                      HARVESTlitter_labile_gCm2day = output[,,38],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,38],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,38],
                      HARVESTlitter_foliage_gCm2day = output[,,39],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,39],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,39],
                      HARVESTlitter_roots_gCm2day = output[,,40],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,40],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,40],
                      HARVESTlitter_wood_gCm2day = output[,,41],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,41],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,41],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,42],
                      mean_labile_gCm2 = output_mean[,42],
                      mean_annual_labile_gCm2 = output_annual[,,42],
                      foliage_gCm2 = output[,,43],
                      mean_foliage_gCm2 = output_mean[,43],
                      mean_annual_foliage_gCm2 = output_annual[,,43],
                      roots_gCm2 = output[,,44],
                      mean_roots_gCm2 = output_mean[,44],
                      mean_annual_roots_gCm2 = output_annual[,,44],
                      wood_gCm2 = output[,,45],
                      mean_wood_gCm2 = output_mean[,45],
                      mean_annual_wood_gCm2 = output_annual[,,45],
                      litter_gCm2 = output[,,46],
                      mean_litter_gCm2 = output_mean[,46],
                      mean_annual_litter_gCm2 = output_annual[,,46],
                      woodlitter_gCm2 = output[,,47],
                      mean_woodlitter_gCm2 = output_mean[,47],
                      mean_annual_woodlitter_gCm2 = output_annual[,,47],
                      som_gCm2 = output[,,48],
                      mean_som_gCm2 = output_mean[,48],
                      mean_annual_som_gCm2 = output_annual[,,48],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,49],
                      mean_ET_kgH2Om2day = output_mean[,49],
                      mean_annual_ET_kgH2Om2day = output_annual[,,49],
                      SurfWater_kgH2Om2 = output[,,50],
                      mean_SurfWater_kgH2Om2 = output_mean[,50],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,50],
                      wSWP_MPa = output[,,51],
                      mean_wSWP_MPa = output_mean[,51],
                      mean_annual_wSWP_MPa = output_annual[,,51],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,52],
                      mean_lai_m2m2 = output_mean[,52],
                      mean_annual_lai_m2m2 = output_annual[,,52],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,53],
                      mean_gs_demand_supply_ratio = output_mean[,53],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,53],
                      gs_mmolH2Om2s = output[,,54],
                      mean_gs_mmolH2Om2s = output_mean[,54],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,54],
                      APAR_MJm2day = output[,,55],
                      mean_APAR_MJm2day = output_mean[,55],
                      mean_annual_APAR_MJm2day = output_annual[,,55],
                      gb_mmolH2Om2s = output[,,56],
                      mean_gb_mmolH2Om2s = output_mean[,56],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,56],
                      CiCa = output[,,57],
                      mean_CiCa = output_mean[,57],
                      mean_annual_CiCa = output_annual[,,57],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_woodlitter_years = MTT_years[,6],
                      MTT_som_years = MTT_years[,7],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_woodlitter_gCm2 = SS_gCm2[,6],
                      SS_som_gCm2 = SS_gCm2[,7])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P2.R1.#") {
      output_dim = 57 ; MTT_dim = 7 ; SS_dim = 7
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec7",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                    
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                             ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))            
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      rhet_woodlitter_gCm2day = output[,,5],
                      mean_rhet_woodlitter_gCm2day = output_mean[,5],
                      mean_annual_rhet_woodlitter_gCm2day = output_annual[,,5],
                      fire_gCm2day = output[,,6],
                      mean_fire_gCm2day = output_mean[,6],
                      mean_annual_fire_gCm2day = output_annual[,,6],
                      harvest_gCm2day = output[,,7],
                      mean_harvest_gCm2day = output_mean[,7],
                      mean_annual_harvest_gCm2day = output_annual[,,7],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,8],
                      mean_alloc_foliage_gCm2day = output_mean[,8],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,8],
                      alloc_labile_gCm2day = output[,,9],
                      mean_alloc_labile_gCm2day = output_mean[,9],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,9],
                      alloc_roots_gCm2day = output[,,10],
                      mean_alloc_roots_gCm2day = output_mean[,10],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,10],
                      alloc_wood_gCm2day = output[,,11],
                      mean_alloc_wood_gCm2day = output_mean[,11],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,11],
                      labile_to_foliage_gCm2day = output[,,12],
                      mean_labile_to_foliage_gCm2day = output_mean[,12],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,12],
                      foliage_to_litter_gCm2day = output[,,13],
                      mean_foliage_to_litter_gCm2day = output_mean[,13],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,13],
                      roots_to_litter_gCm2day = output[,,14],
                      mean_roots_to_litter_gCm2day = output_mean[,14],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,14],
                      wood_to_litter_gCm2day = output[,,15],
                      mean_wood_to_litter_gCm2day = output_mean[,15],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,15],
                      litter_to_som_gCm2day = output[,,16],
                      mean_litter_to_som_gCm2day = output_mean[,16],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,16],
                      woodlitter_to_som_gCm2day = output[,,17],
                      mean_woodlitter_to_som_gCm2day = output_mean[,17],
                      mean_annual_woodlitter_to_som_gCm2day = output_annual[,,17],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,18],
                      mean_FIREemiss_labile_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,18],
                      FIRElitter_labile_gCm2day = output[,,19],
                      mean_FIRElitter_labile_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,19],
                      FIREemiss_foliage_gCm2day = output[,,20],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,20],
                      FIRElitter_foliage_gCm2day = output[,,21],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,21],
                      FIREemiss_roots_gCm2day = output[,,22],
                      mean_FIREemiss_roots_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,22],
                      FIRElitter_roots_gCm2day = output[,,23],
                      mean_FIRElitter_roots_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,23],
                      FIREemiss_wood_gCm2day = output[,,24],
                      mean_FIREemiss_wood_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,24],
                      FIRElitter_wood_gCm2day = output[,,25],
                      mean_FIRElitter_wood_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,25],
                      FIREemiss_litter_gCm2day = output[,,26],
                      mean_FIREemiss_litter_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,26],
                      FIRElitter_litter_gCm2day = output[,,27],
                      mean_FIRElitter_litter_gCm2day = output_mean[,27],
                      mean_annual_FIRElitter_litter_gCm2day = output[,,27],
                      FIREemiss_woodlitter_gCm2day = output[,,28],
                      mean_FIREemiss_woodlitter_gCm2day = output_mean[,28],
                      mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,28],
                      FIRElitter_woodlitter_gCm2day = output[,,29],
                      mean_FIRElitter_woodlitter_gCm2day = output_mean[,29],
                      mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,29],
                      FIREemiss_som_gCm2day = output[,,30],
                      mean_FIREemiss_som_gCm2day = output_mean[,30],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,30],
                      HARVESTextracted_labile_gCm2day = output[,,31],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,31],
                      HARVESTextracted_foliage_gCm2day = output[,,32],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,32],
                      HARVESTextracted_roots_gCm2day = output[,,33],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,33],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,33],
                      HARVESTextracted_wood_gCm2day = output[,,34],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,34],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,34],
                      HARVESTextracted_litter_gCm2day = output[,,35],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,35],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,35],
                      HARVESTextracted_woodlitter_gCm2day = output[,,36],
                      mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,36],
                      mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,36],
                      HARVESTextracted_som_gCm2day = output[,,37],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,37],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,37],
                      HARVESTlitter_labile_gCm2day = output[,,38],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,38],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,38],
                      HARVESTlitter_foliage_gCm2day = output[,,39],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,39],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,39],
                      HARVESTlitter_roots_gCm2day = output[,,40],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,40],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,40],
                      HARVESTlitter_wood_gCm2day = output[,,41],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,41],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,41],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,42],
                      mean_labile_gCm2 = output_mean[,42],
                      mean_annual_labile_gCm2 = output_annual[,,42],
                      foliage_gCm2 = output[,,43],
                      mean_foliage_gCm2 = output_mean[,43],
                      mean_annual_foliage_gCm2 = output_annual[,,43],
                      roots_gCm2 = output[,,44],
                      mean_roots_gCm2 = output_mean[,44],
                      mean_annual_roots_gCm2 = output_annual[,,44],
                      wood_gCm2 = output[,,45],
                      mean_wood_gCm2 = output_mean[,45],
                      mean_annual_wood_gCm2 = output_annual[,,45],
                      litter_gCm2 = output[,,46],
                      mean_litter_gCm2 = output_mean[,46],
                      mean_annual_litter_gCm2 = output_annual[,,46],
                      woodlitter_gCm2 = output[,,47],
                      mean_woodlitter_gCm2 = output_mean[,47],
                      mean_annual_woodlitter_gCm2 = output_annual[,,47],
                      som_gCm2 = output[,,48],
                      mean_som_gCm2 = output_mean[,48],
                      mean_annual_som_gCm2 = output_annual[,,48],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,49],
                      mean_ET_kgH2Om2day = output_mean[,49],
                      mean_annual_ET_kgH2Om2day = output_annual[,,49],
                      SurfWater_kgH2Om2 = output[,,50],
                      mean_SurfWater_kgH2Om2 = output_mean[,50],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,50],
                      wSWP_MPa = output[,,51],
                      mean_wSWP_MPa = output_mean[,51],
                      mean_annual_wSWP_MPa = output_annual[,,51],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,52],
                      mean_lai_m2m2 = output_mean[,52],
                      mean_annual_lai_m2m2 = output_annual[,,52],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,53],
                      mean_gs_demand_supply_ratio = output_mean[,53],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,53],
                      gs_mmolH2Om2s = output[,,54],
                      mean_gs_mmolH2Om2s = output_mean[,54],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,54],
                      APAR_MJm2day = output[,,55],
                      mean_APAR_MJm2day = output_mean[,55],
                      mean_annual_APAR_MJm2day = output_annual[,,55],
                      gb_mmolH2Om2s = output[,,56],
                      mean_gb_mmolH2Om2s = output_mean[,56],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,56],
                      CiCa = output[,,57],
                      mean_CiCa = output_mean[,57],
                      mean_annual_CiCa = output_annual[,,57],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_woodlitter_years = MTT_years[,6],
                      MTT_som_years = MTT_years[,7],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_woodlitter_gCm2 = SS_gCm2[,6],
                      SS_som_gCm2 = SS_gCm2[,7])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P2.R3.#") {
      output_dim = 57 ; MTT_dim = 7 ; SS_dim = 7
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec19",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                     
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))            
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      rhet_woodlitter_gCm2day = output[,,5],
                      mean_rhet_woodlitter_gCm2day = output_mean[,5],
                      mean_annual_rhet_woodlitter_gCm2day = output_annual[,,5],
                      fire_gCm2day = output[,,6],
                      mean_fire_gCm2day = output_mean[,6],
                      mean_annual_fire_gCm2day = output_annual[,,6],
                      harvest_gCm2day = output[,,7],
                      mean_harvest_gCm2day = output_mean[,7],
                      mean_annual_harvest_gCm2day = output_annual[,,7],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,8],
                      mean_alloc_foliage_gCm2day = output_mean[,8],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,8],
                      alloc_labile_gCm2day = output[,,9],
                      mean_alloc_labile_gCm2day = output_mean[,9],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,9],
                      alloc_roots_gCm2day = output[,,10],
                      mean_alloc_roots_gCm2day = output_mean[,10],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,10],
                      alloc_wood_gCm2day = output[,,11],
                      mean_alloc_wood_gCm2day = output_mean[,11],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,11],
                      labile_to_foliage_gCm2day = output[,,12],
                      mean_labile_to_foliage_gCm2day = output_mean[,12],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,12],
                      foliage_to_litter_gCm2day = output[,,13],
                      mean_foliage_to_litter_gCm2day = output_mean[,13],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,13],
                      roots_to_litter_gCm2day = output[,,14],
                      mean_roots_to_litter_gCm2day = output_mean[,14],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,14],
                      wood_to_litter_gCm2day = output[,,15],
                      mean_wood_to_litter_gCm2day = output_mean[,15],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,15],
                      litter_to_som_gCm2day = output[,,16],
                      mean_litter_to_som_gCm2day = output_mean[,16],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,16],
                      woodlitter_to_som_gCm2day = output[,,17],
                      mean_woodlitter_to_som_gCm2day = output_mean[,17],
                      mean_annual_woodlitter_to_som_gCm2day = output_annual[,,17],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,18],
                      mean_FIREemiss_labile_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,18],
                      FIRElitter_labile_gCm2day = output[,,19],
                      mean_FIRElitter_labile_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,19],
                      FIREemiss_foliage_gCm2day = output[,,20],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,20],
                      FIRElitter_foliage_gCm2day = output[,,21],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,21],
                      FIREemiss_roots_gCm2day = output[,,22],
                      mean_FIREemiss_roots_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,22],
                      FIRElitter_roots_gCm2day = output[,,23],
                      mean_FIRElitter_roots_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,23],
                      FIREemiss_wood_gCm2day = output[,,24],
                      mean_FIREemiss_wood_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,24],
                      FIRElitter_wood_gCm2day = output[,,25],
                      mean_FIRElitter_wood_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,25],
                      FIREemiss_litter_gCm2day = output[,,26],
                      mean_FIREemiss_litter_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,26],
                      FIRElitter_litter_gCm2day = output[,,27],
                      mean_FIRElitter_litter_gCm2day = output_mean[,27],
                      mean_annual_FIRElitter_litter_gCm2day = output[,,27],
                      FIREemiss_woodlitter_gCm2day = output[,,28],
                      mean_FIREemiss_woodlitter_gCm2day = output_mean[,28],
                      mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,28],
                      FIRElitter_woodlitter_gCm2day = output[,,29],
                      mean_FIRElitter_woodlitter_gCm2day = output_mean[,29],
                      mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,29],
                      FIREemiss_som_gCm2day = output[,,30],
                      mean_FIREemiss_som_gCm2day = output_mean[,30],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,30],
                      HARVESTextracted_labile_gCm2day = output[,,31],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,31],
                      HARVESTextracted_foliage_gCm2day = output[,,32],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,32],
                      HARVESTextracted_roots_gCm2day = output[,,33],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,33],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,33],
                      HARVESTextracted_wood_gCm2day = output[,,34],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,34],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,34],
                      HARVESTextracted_litter_gCm2day = output[,,35],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,35],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,35],
                      HARVESTextracted_woodlitter_gCm2day = output[,,36],
                      mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,36],
                      mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,36],
                      HARVESTextracted_som_gCm2day = output[,,37],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,37],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,37],
                      HARVESTlitter_labile_gCm2day = output[,,38],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,38],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,38],
                      HARVESTlitter_foliage_gCm2day = output[,,39],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,39],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,39],
                      HARVESTlitter_roots_gCm2day = output[,,40],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,40],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,40],
                      HARVESTlitter_wood_gCm2day = output[,,41],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,41],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,41],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,42],
                      mean_labile_gCm2 = output_mean[,42],
                      mean_annual_labile_gCm2 = output_annual[,,42],
                      foliage_gCm2 = output[,,43],
                      mean_foliage_gCm2 = output_mean[,43],
                      mean_annual_foliage_gCm2 = output_annual[,,43],
                      roots_gCm2 = output[,,44],
                      mean_roots_gCm2 = output_mean[,44],
                      mean_annual_roots_gCm2 = output_annual[,,44],
                      wood_gCm2 = output[,,45],
                      mean_wood_gCm2 = output_mean[,45],
                      mean_annual_wood_gCm2 = output_annual[,,45],
                      litter_gCm2 = output[,,46],
                      mean_litter_gCm2 = output_mean[,46],
                      mean_annual_litter_gCm2 = output_annual[,,46],
                      woodlitter_gCm2 = output[,,47],
                      mean_woodlitter_gCm2 = output_mean[,47],
                      mean_annual_woodlitter_gCm2 = output_annual[,,47],
                      som_gCm2 = output[,,48],
                      mean_som_gCm2 = output_mean[,48],
                      mean_annual_som_gCm2 = output_annual[,,48],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,49],
                      mean_ET_kgH2Om2day = output_mean[,49],
                      mean_annual_ET_kgH2Om2day = output_annual[,,49],
                      SurfWater_kgH2Om2 = output[,,50],
                      mean_SurfWater_kgH2Om2 = output_mean[,50],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,50],
                      wSWP_MPa = output[,,51],
                      mean_wSWP_MPa = output_mean[,51],
                      mean_annual_wSWP_MPa = output_annual[,,51],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,52],
                      mean_lai_m2m2 = output_mean[,52],
                      mean_annual_lai_m2m2 = output_annual[,,52],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,53],
                      mean_gs_demand_supply_ratio = output_mean[,53],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,53],
                      gs_mmolH2Om2s = output[,,54],
                      mean_gs_mmolH2Om2s = output_mean[,54],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,54],
                      APAR_MJm2day = output[,,55],
                      mean_APAR_MJm2day = output_mean[,55],
                      mean_annual_APAR_MJm2day = output_annual[,,55],
                      gb_mmolH2Om2s = output[,,56],
                      mean_gb_mmolH2Om2s = output_mean[,56],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,56],
                      CiCa = output[,,57],
                      mean_CiCa = output_mean[,57],
                      mean_annual_CiCa = output_annual[,,57],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_woodlitter_years = MTT_years[,6],
                      MTT_som_years = MTT_years[,7],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_woodlitter_gCm2 = SS_gCm2[,6],
                      SS_som_gCm2 = SS_gCm2[,7])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.C5.D1.F2.P1.#") {
    output_dim = 22 ; MTT_dim = 4 ; SS_dim = 4
    dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
    tmp=.Fortran( "rdalec13",output_dim=as.integer(output_dim)
                            ,MTT_dim=as.integer(MTT_dim)
                            ,SS_dim = as.integer(SS_dim)
                            ,met=as.double(t(met))
                            ,pars=as.double(pars_in)
                            ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                            ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                            ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                            ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                            ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                   
                            ,lat=as.double(lat)
                            ,nopars=as.integer(PROJECT$model$nopars[site])
                            ,nomet=as.integer(dim(met)[2])
                            ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                            ,nopools=as.integer(PROJECT$model$nopools[site])
                            ,nodays=as.integer(dim(met)[1])
                            ,nos_years=as.integer(noyears)
                            ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                            ,nos_iter=as.integer(nos_iter) )
    output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
    output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))          
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    mean_gpp_gCm2day = output_mean[,1],
                    mean_annual_gpp_gCm2day = output_annual[,,1],
                    rauto_gCm2day = output[,,2],
                    mean_rauto_gCm2day = output_mean[,2],
                    mean_annual_rauto_gCm2day = output_annual[,,2],
                    rhet_dom_gCm2day = output[,,3],
                    mean_rhet_dom_gCm2day = output_mean[,3],
                    mean_annual_rhet_dom_gCm2day = output_annual[,,3],
                    fire_gCm2day = output[,,4],
                    mean_fire_gCm2day = output_mean[,4],
                    mean_annual_fire_gCm2day = output_annual[,,4],
                    # Internal fluxes
                    alloc_foliage_gCm2day = output[,,5],
                    mean_alloc_foliage_gCm2day = output_mean[,5],
                    mean_annual_alloc_foliage_gCm2day = output_annual[,,5],
                    alloc_labile_gCm2day = output[,,6],
                    mean_alloc_labile_gCm2day = output_mean[,6],
                    mean_annual_alloc_labile_gCm2day = output_annual[,,6],
                    alloc_roots_wood_gCm2day = output[,,7],
                    mean_alloc_roots_wood_gCm2day = output_mean[,7],
                    mean_annual_alloc_roots_wood_gCm2day = output_annual[,,7],
                    labile_to_foliage_gCm2day = output[,,8],
                    mean_labile_to_foliage_gCm2day = output_mean[,8],
                    mean_annual_labile_to_foliage_gCm2day = output_annual[,,8],
                    foliage_to_litter_gCm2day = output[,,9],
                    mean_foliage_to_litter_gCm2day = output_mean[,9],
                    mean_annual_foliage_to_litter_gCm2day = output_annual[,,9],
                    roots_wood_to_litter_gCm2day = output[,,10],
                    mean_roots_wood_to_litter_gCm2day = output_mean[,10],
                    mean_annual_roots_wood_to_litter_gCm2day = output_annual[,,10],
                    # Disturbance fluxes
                    FIREemiss_labile_gCm2day = output[,,11],
                    mean_FIREemiss_labile_gCm2day = output_mean[,11],
                    mean_annual_FIREemiss_labile_gCm2day = output_annual[,,11],
                    FIRElitter_labile_gCm2day = output[,,12],
                    mean_FIRElitter_labile_gCm2day = output_mean[,12],
                    mean_annual_FIRElitter_labile_gCm2day = output_annual[,,12],
                    FIREemiss_foliage_gCm2day = output[,,13],
                    mean_FIREemiss_foliage_gCm2day = output_mean[,13],
                    mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,13],
                    FIRElitter_foliage_gCm2day = output[,,14],
                    mean_FIRElitter_foliage_gCm2day = output_mean[,14],
                    mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,14],
                    FIREemiss_roots_wood_gCm2day = output[,,15],
                    mean_FIREemiss_roots_wood_gCm2day = output_mean[,15],
                    mean_annual_FIREemiss_roots_wood_gCm2day = output_annual[,,15],
                    FIRElitter_roots_wood_gCm2day = output[,,16],
                    mean_FIRElitter_roots_wood_gCm2day = output_mean[,16],
                    mean_annual_FIRElitter_roots_wood_gCm2day = output_annual[,,16],
                    FIREemiss_dom_gCm2day = output[,,17],
                    mean_FIREemiss_dom_gCm2day = output_mean[,17],
                    mean_annual_FIREemiss_dom_gCm2day = output_annual[,,17],
                    # C pools (gC/m2)
                    labile_gCm2 = output[,,18],
                    mean_labile_gCm2 = output_mean[,18],
                    mean_annual_labile_gCm2 = output_annual[,,18],
                    foliage_gCm2 = output[,,19],
                    mean_foliage_gCm2 = output_mean[,19],
                    mean_annual_foliage_gCm2 = output_annual[,,19],
                    roots_wood_gCm2 = output[,,20],
                    mean_roots_wood_gCm2 = output_mean[,20],
                    mean_annual_roots_wood_gCm2 = output_annual[,,20],
                    dom_gCm2 = output[,,21],
                    mean_dom_gCm2 = output_mean[,21],
                    mean_annual_dom_gCm2 = output_annual[,,21],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,22],
                    mean_lai_m2m2 = output_mean[,22],
                    mean_annual_lai_m2m2 = output_annual[,,22],
                    ## Aggregated variables
                    # Mean Transit times
                    MTT_labile_years = MTT_years[,1],
                    MTT_foliage_years = MTT_years[,2],
                    MTT_roots_wood_years = MTT_years[,3],
                    MTT_dom_years = MTT_years[,4],
                    # Steady state estimates
                    SS_labile_gCm2 = SS_gCm2[,1],
                    SS_foliage_gCm2 = SS_gCm2[,2],
                    SS_roots_wood_gCm2 = SS_gCm2[,3],
                    SS_dom_gCm2 = SS_gCm2[,4])
    # Determine the NPP fraction of expressed NPP
    # i.e. actual growth not GPP-Ra
    NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                         states_all$alloc_foliage_gCm2day +
                         states_all$alloc_roots_wood_gCm2day,1,mean)
    NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                         apply(states_all$alloc_roots_wood_gCm2day,1,mean)) / NPP_fraction
    states_all$NPP_foliage_fraction = NPP_fraction[,1]
    states_all$NPP_roots_wood_fraction = NPP_fraction[,2]
    # Tidy up variables
    rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.D1.F2.#") {
      output_dim = 36 ; MTT_dim = 5 ; SS_dim = 5
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec1",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                    
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))            
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # Create the output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_roots_gCm2day = output[,,8],
                      mean_alloc_roots_gCm2day = output_mean[,8],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,8],
                      alloc_wood_gCm2day = output[,,9],
                      mean_alloc_wood_gCm2day = output_mean[,9],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,9],
                      foliage_to_litter_gCm2day = output[,,10],
                      mean_foliage_to_litter_gCm2day = output_mean[,10],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,10],
                      roots_to_litter_gCm2day = output[,,11],
                      mean_roots_to_litter_gCm2day = output_mean[,11],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,11],
                      wood_to_litter_gCm2day = output[,,12],
                      mean_wood_to_litter_gCm2day = output_mean[,12],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,12],
                      litter_to_som_gCm2day = output[,,13],
                      mean_litter_to_som_gCm2day = output_mean[,13],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,13],
                      # Disturbance fluxes
                      FIREemiss_foliage_gCm2day = output[,,14],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,14],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,14],
                      FIRElitter_foliage_gCm2day = output[,,15],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,15],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,15],
                      FIREemiss_roots_gCm2day = output[,,16],
                      mean_FIREemiss_roots_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,16],
                      FIRElitter_roots_gCm2day = output[,,17],
                      mean_FIRElitter_roots_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,17],
                      FIREemiss_wood_gCm2day = output[,,18],
                      mean_FIREemiss_wood_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,18],
                      FIRElitter_wood_gCm2day = output[,,19],
                      mean_FIRElitter_wood_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,19],
                      FIREemiss_litter_gCm2day = output[,,20],
                      mean_FIREemiss_litter_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,20],
                      FIRElitter_litter_gCm2day = output[,,21],
                      mean_FIRElitter_litter_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,21],
                      FIREemiss_som_gCm2day = output[,,22],
                      mean_FIREemiss_som_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,22],
                      HARVESTextracted_foliage_gCm2day = output[,,23],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,23],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,23],
                      HARVESTextracted_roots_gCm2day = output[,,24],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,24],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,24],
                      HARVESTextracted_wood_gCm2day = output[,,25],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,25],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,25],
                      HARVESTextracted_litter_gCm2day = output[,,26],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,26],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,26],
                      HARVESTextracted_som_gCm2day = output[,,27],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,27],
                      HARVESTlitter_foliage_gCm2day = output[,,28],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,28],
                      HARVESTlitter_roots_gCm2day = output[,,29],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,29],
                      HARVESTlitter_wood_gCm2day = output[,,30],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,30],
                      # C pools (gC/m2)
                      foliage_gCm2 = output[,,31],
                      mean_foliage_gCm2 = output_mean[,31],
                      mean_annual_foliage_gCm2 = output_annual[,,31],
                      roots_gCm2 = output[,,32],
                      mean_roots_gCm2 = output_mean[,32],
                      mean_annual_roots_gCm2 = output_annual[,,32],
                      wood_gCm2 = output[,,33],
                      mean_wood_gCm2 = output_mean[,33],
                      mean_annual_wood_gCm2 = output_annual[,,33],
                      litter_gCm2 = output[,,34],
                      mean_litter_gCm2 = output_mean[,34],
                      mean_annual_litter_gCm2 = output_annual[,,34],
                      som_gCm2 = output[,,35],
                      mean_som_gCm2 = output_mean[,35],
                      mean_annual_som_gCm2 = output_annual[,,35],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,36],
                      mean_lai_m2m2 = output_mean[,36],
                      mean_annual_lai_m2m2 = output_annual[,,36],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_foliage_years = MTT_years[,1],
                      MTT_roots_years = MTT_years[,2],
                      MTT_wood_years = MTT_years[,3],
                      MTT_litter_years = MTT_years[,4],
                      MTT_som_years = MTT_years[,5],
                      # Steady state estimates
                      SS_foliage_gCm2 = SS_gCm2[,1],
                      SS_roots_gCm2 = SS_gCm2[,2],
                      SS_wood_gCm2 = SS_gCm2[,3],
                      SS_litter_gCm2 = SS_gCm2[,4],
                      SS_som_gCm2 = SS_gCm2[,5])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Final tidy
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.C4.D1.F2.#") {
    output_dim = 23 ; MTT_dim = 3 ; SS_dim = 3
    dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
    tmp=.Fortran( "rdalec12",output_dim=as.integer(output_dim)
                            ,MTT_dim=as.integer(MTT_dim)
                            ,SS_dim = as.integer(SS_dim)
                            ,met=as.double(t(met))
                            ,pars=as.double(pars_in)
                            ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                            ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                            ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                            ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                            ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                   
                            ,lat=as.double(lat)
                            ,nopars=as.integer(PROJECT$model$nopars[site])
                            ,nomet=as.integer(dim(met)[2])
                            ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                            ,nopools=as.integer(PROJECT$model$nopools[site])
                            ,nodays=as.integer(dim(met)[1])
                            ,nos_years=as.integer(noyears)
                            ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                            ,nos_iter=as.integer(nos_iter) )
    output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
    output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))          
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    mean_gpp_gCm2day = output_mean[,1],
                    mean_annual_gpp_gCm2day = output_annual[,,1],
                    rauto_gCm2day = output[,,2],
                    mean_rauto_gCm2day = output_mean[,2],
                    mean_annual_rauto_gCm2day = output_annual[,,2],
                    rhet_dom_gCm2day = output[,,3],
                    mean_rhet_dom_gCm2day = output_mean[,3],
                    mean_annual_rhet_dom_gCm2day = output_annual[,,3],
                    fire_gCm2day = output[,,4],
                    mean_fire_gCm2day = output_mean[,4],
                    mean_annual_fire_gCm2day = output_annual[,,4],
                    harvest_gCm2day = output[,,5],
                    mean_harvest_gCm2day = output_mean[,5],
                    mean_annual_harvest_gCm2day = output_annual[,,5],
                    # Internal fluxes
                    alloc_foliage_gCm2day = output[,,6],
                    mean_alloc_foliage_gCm2day = output_mean[,6],
                    mean_annual_alloc_foliage_gCm2day = output_annual[,,6],
                    alloc_roots_wood_gCm2day = output[,,7],
                    mean_alloc_roots_wood_gCm2day = output_mean[,7],
                    mean_annual_alloc_roots_wood_gCm2day = output_annual[,,7],
                    foliage_to_litter_gCm2day = output[,,8],
                    mean_foliage_to_litter_gCm2day = output_mean[,8],
                    mean_annual_foliage_to_litter_gCm2day = output_annual[,,8],
                    roots_wood_to_litter_gCm2day = output[,,9],
                    mean_roots_wood_to_litter_gCm2day = output_mean[,9],
                    mean_annual_roots_wood_to_litter_gCm2day = output_annual[,,9],
                    # Disturbance fluxes
                    FIREemiss_foliage_gCm2day = output[,,10],
                    mean_FIREemiss_foliage_gCm2day = output_mean[,10],
                    mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,10],
                    FIRElitter_foliage_gCm2day = output[,,11],
                    mean_FIRElitter_foliage_gCm2day = output_mean[,11],
                    mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,11],
                    FIREemiss_roots_wood_gCm2day = output[,,12],
                    mean_FIREemiss_roots_wood_gCm2day = output_mean[,12],
                    mean_annual_FIREemiss_roots_wood_gCm2day = output_annual[,,12],
                    FIRElitter_roots_wood_gCm2day = output[,,13],
                    mean_FIRElitter_roots_wood_gCm2day = output_mean[,13],
                    mean_annual_FIRElitter_roots_wood_gCm2day = output_annual[,,13],
                    FIREemiss_dom_gCm2day = output[,,14],
                    mean_FIREemiss_dom_gCm2day = output_mean[,14],
                    mean_annual_FIREemiss_dom_gCm2day = output_annual[,,14],
                    HARVESTextracted_foliage_gCm2day = output[,,15],
                    mean_HARVESTextracted_foliage_gCm2day = output_mean[,15],
                    mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,15],
                    HARVESTextracted_roots_wood_gCm2day = output[,,16],
                    mean_HARVESTextracted_roots_wood_gCm2day = output_mean[,16],
                    mean_annual_HARVESTextracted_roots_wood_gCm2day = output_annual[,,16],
                    HARVESTextracted_som_gCm2day = output[,,17],
                    mean_HARVESTextracted_som_gCm2day = output_mean[,17],
                    mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,17],
                    HARVESTlitter_foliage_gCm2day = output[,,18],
                    mean_HARVESTlitter_foliage_gCm2day = output_mean[,18],
                    mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,18],
                    HARVESTlitter_roots_wood_gCm2day = output[,,19],
                    mean_HARVESTlitter_roots_wood_gCm2day = output_mean[,19],
                    mean_annual_HARVESTlitter_roots_wood_gCm2day = output_annual[,,19],
                    # C pools (gC/m2)
                    foliage_gCm2 = output[,,20],
                    mean_foliage_gCm2 = output[,,20],
                    mean_annual_foliage_gCm2 = output_annual[,,20],
                    roots_wood_gCm2 = output[,,21],
                    mean_roots_wood_gCm2 = output[,,21],
                    mean_annual_roots_wood_gCm2 = output_annual[,,21],
                    dom_gCm2 = output[,,22],
                    mean_dom_gCm2 = output[,,22],
                    mean_annual_dom_gCm2 = output_annual[,,22],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,23],
                    mean_lai_m2m2 = output[,,23],
                    mean_annual_lai_m2m2 = output_annual[,,23],
                    ## Aggregated variables
                    # Mean Transit times
                    MTT_foliage_years = MTT_years[,1],
                    MTT_roots_wood_years = MTT_years[,2],
                    MTT_dom_years = MTT_years[,3],
                    # Steady state estimates
                    SS_foliage_gCm2 = SS_gCm2[,1],
                    SS_roots_wood_gCm2 = SS_gCm2[,2],
                    SS_dom_gCm2 = SS_gCm2[,3])
    # Determine the NPP fraction of expressed NPP
    # i.e. actual growth not GPP-Ra
    NPP_fraction = apply(states_all$alloc_foliage_gCm2day +
                         states_all$alloc_roots_wood_gCm2day,1,mean)
    NPP_fraction = cbind(apply(states_all$alloc_foliage_gCm2day,1,mean),
                         apply(states_all$alloc_roots_wood_gCm2day,1,mean)) / NPP_fraction
    states_all$NPP_foliage_fraction = NPP_fraction[,1]
    states_all$NPP_roots_wood_fraction = NPP_fraction[,2]
    # Tidy up variables
    rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H1.P3.R1.#") {
    output_dim = 59 ; MTT_dim = 7 ; SS_dim = 7
    # Load the required dalec shared object
    dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
    tmp=.Fortran( "rdalec8",output_dim=as.integer(output_dim)
                           ,aNPP_dim=as.integer(aNPP_dim)
                           ,MTT_dim=as.integer(MTT_dim)
                           ,SS_dim = as.integer(SS_dim)
                           ,fire_dim=as.integer(fire_dim)
                           ,met=as.double(t(met))
                           ,pars=as.double(pars_in)
                           ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                           ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                           ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                           ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                           ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                  
                           ,lat=as.double(lat)
                           ,nopars=as.integer(PROJECT$model$nopars[site])
                           ,nomet=as.integer(dim(met)[2])
                           ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                           ,nopools=as.integer(PROJECT$model$nopools[site])
                           ,pft=as.integer(pft)
                           ,nodays=as.integer(dim(met)[1])
                           ,nos_years=as.integer(noyears)
                           ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                           ,nos_iter=as.integer(nos_iter))
    # Extract the different output variables
    output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
    output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))          
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    mean_gpp_gCm2day = output_mean[,1],
                    mean_annual_gpp_gCm2day = output_annual[,,1],
                    rauto_gCm2day = output[,,2],
                    mean_rauto_gCm2day = output_mean[,2],
                    mean_annual_rauto_gCm2day = output_annual[,,2],
                    rg_foliage_gCm2day = output[,,3],
                    mean_rg_foliage_gCm2day = output_mean[,3],
                    mean_annual_rg_foliage_gCm2day = output_annual[,,3],
                    rhet_litter_gCm2day = output[,,4],
                    mean_rhet_litter_gCm2day = output_mean[,4],
                    mean_annual_rhet_litter_gCm2day = output_annual[,,4],
                    rhet_som_gCm2day = output[,,5],
                    mean_rhet_som_gCm2day = output_mean[,5],
                    mean_annual_rhet_som_gCm2day = output_annual[,,5],
                    rhet_woodlitter_gCm2day = output[,,6],
                    mean_rhet_woodlitter_gCm2day = output_mean[,6],
                    mean_annual_rhet_woodlitter_gCm2day = output_annual[,,6],
                    fire_gCm2day = output[,,7],
                    mean_fire_gCm2day = output_mean[,7],
                    mean_annual_fire_gCm2day = output_annual[,,7],
                    harvest_gCm2day = output[,,8],
                    mean_harvest_gCm2day = output_mean[,8],
                    mean_annual_harvest_gCm2day = output_annual[,,8],
                    # Internal fluxes
                    alloc_labile_gCm2day = output[,,9],
                    mean_alloc_labile_gCm2day = output_mean[,9],
                    mean_annual_alloc_labile_gCm2day = output_annual[,,9],
                    alloc_roots_gCm2day = output[,,10],
                    mean_alloc_roots_gCm2day = output_mean[,10],
                    mean_annual_alloc_roots_gCm2day = output_annual[,,10],
                    alloc_wood_gCm2day = output[,,11],
                    mean_alloc_wood_gCm2day = output_mean[,11],
                    mean_annual_alloc_wood_gCm2day = output_annual[,,11],
                    labile_to_foliage_gCm2day = output[,,12],
                    mean_labile_to_foliage_gCm2day = output_mean[,12],
                    mean_annual_labile_to_foliage_gCm2day = output_annual[,,12],
                    foliage_to_litter_gCm2day = output[,,13],
                    mean_foliage_to_litter_gCm2day = output_mean[,13],
                    mean_annual_foliage_to_litter_gCm2day = output_annual[,,13],
                    roots_to_litter_gCm2day = output[,,14],
                    mean_roots_to_litter_gCm2day = output_mean[,14],
                    mean_annual_roots_to_litter_gCm2day = output_annual[,,14],
                    wood_to_litter_gCm2day = output[,,15],
                    mean_wood_to_litter_gCm2day = output_mean[,15],
                    mean_annual_wood_to_litter_gCm2day = output_annual[,,15],
                    litter_to_som_gCm2day = output[,,16],
                    mean_litter_to_som_gCm2day = output_mean[,16],
                    mean_annual_litter_to_som_gCm2day = output_annual[,,16],
                    woodlitter_to_som_gCm2day = output[,,17],
                    mean_woodlitter_to_som_gCm2day = output_mean[,17],
                    mean_annual_woodlitter_to_som_gCm2day = output_annual[,,17],
                    # Disturbance fluxes
                    FIREemiss_labile_gCm2day = output[,,18],
                    mean_FIREemiss_labile_gCm2day = output_mean[,18],
                    mean_annual_FIREemiss_labile_gCm2day = output_annual[,,18],
                    FIRElitter_labile_gCm2day = output[,,19],
                    mean_FIRElitter_labile_gCm2day = output_mean[,19],
                    mean_annual_FIRElitter_labile_gCm2day = output_annual[,,19],
                    FIREemiss_foliage_gCm2day = output[,,20],
                    mean_FIREemiss_foliage_gCm2day = output_mean[,20],
                    mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,20],
                    FIRElitter_foliage_gCm2day = output[,,21],
                    mean_FIRElitter_foliage_gCm2day = output_mean[,21],
                    mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,21],
                    FIREemiss_roots_gCm2day = output[,,22],
                    mean_FIREemiss_roots_gCm2day = output[,,22],
                    mean_annual_FIREemiss_roots_gCm2day = output_annual[,,22],
                    FIRElitter_roots_gCm2day = output[,,23],
                    mean_FIRElitter_roots_gCm2day = output_mean[,23],
                    mean_annual_FIRElitter_roots_gCm2day = output_annual[,,23],
                    FIREemiss_wood_gCm2day = output[,,24],
                    mean_FIREemiss_wood_gCm2day = output_mean[,24],
                    mean_annual_FIREemiss_wood_gCm2day = output_annual[,,24],
                    FIRElitter_wood_gCm2day = output[,,25],
                    mean_FIRElitter_wood_gCm2day = output_mean[,25],
                    mean_annual_FIRElitter_wood_gCm2day = output_annual[,,25],
                    FIREemiss_litter_gCm2day = output[,,26],
                    mean_FIREemiss_litter_gCm2day = output_mean[,26],
                    mean_annual_FIREemiss_litter_gCm2day = output_annual[,,26],
                    FIRElitter_litter_gCm2day = output[,,27],
                    mean_FIRElitter_litter_gCm2day = output[,,27],
                    mean_annual_FIRElitter_litter_gCm2day = output_annual[,,27],
                    FIREemiss_woodlitter_gCm2day = output[,,28],
                    mean_FIREemiss_woodlitter_gCm2day = output_mean[,28],
                    mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,28],
                    FIRElitter_woodlitter_gCm2day = output[,,29],
                    mean_FIRElitter_woodlitter_gCm2day = output_mean[,29],
                    mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,29],
                    FIREemiss_som_gCm2day = output[,,30],
                    mean_FIREemiss_som_gCm2day = output_mean[,30],
                    mean_annual_FIREemiss_som_gCm2day = output_annual[,,30],
                    HARVESTextracted_labile_gCm2day = output[,,31],
                    mean_HARVESTextracted_labile_gCm2day = output_mean[,31],
                    mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,31],
                    HARVESTextracted_foliage_gCm2day = output[,,32],
                    mean_HARVESTextracted_foliage_gCm2day = output_mean[,32],
                    mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,32],
                    HARVESTextracted_roots_gCm2day = output[,,33],
                    mean_HARVESTextracted_roots_gCm2day = output_mean[,33],
                    mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,33],
                    HARVESTextracted_wood_gCm2day = output[,,34],
                    mean_HARVESTextracted_wood_gCm2day = output_mean[,34],
                    mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,34],
                    HARVESTextracted_litter_gCm2day = output[,,35],
                    mean_HARVESTextracted_litter_gCm2day = output_mean[,35],
                    mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,35],
                    HARVESTextracted_woodlitter_gCm2day = output[,,36],
                    mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,36],
                    mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,36],
                    HARVESTextracted_som_gCm2day = output[,,37],
                    mean_HARVESTextracted_som_gCm2day = output_mean[,37],
                    mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,37],
                    HARVESTlitter_labile_gCm2day = output[,,38],
                    mean_HARVESTlitter_labile_gCm2day = output_mean[,38],
                    mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,38],
                    HARVESTlitter_foliage_gCm2day = output[,,39],
                    mean_HARVESTlitter_foliage_gCm2day = output_mean[,39],
                    mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,39],
                    HARVESTlitter_roots_gCm2day = output[,,40],
                    mean_HARVESTlitter_roots_gCm2day = output_mean[,40],
                    mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,40],
                    HARVESTlitter_wood_gCm2day = output[,,41],
                    mean_HARVESTlitter_wood_gCm2day = output_mean[,41],
                    mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,41],
                    # C pools (gC/m2)
                    labile_gCm2 = output[,,42],
                    mean_labile_gCm2 = output[,,42],
                    mean_annual_labile_gCm2 = output_annual[,,42],
                    foliage_gCm2 = output[,,43],
                    mean_foliage_gCm2 = output[,,43],
                    mean_annual_foliage_gCm2 = output_annual[,,43],
                    roots_gCm2 = output[,,44],
                    mean_roots_gCm2 = output[,,44],
                    mean_annual_roots_gCm2 = output_annual[,,44],
                    wood_gCm2 = output[,,45],
                    mean_wood_gCm2 = output[,,45],
                    mean_annual_wood_gCm2 = output_annual[,,45],
                    litter_gCm2 = output[,,46],
                    mean_litter_gCm2 = output[,,46],
                    mean_annual_litter_gCm2 = output_annual[,,46],
                    woodlitter_gCm2 = output[,,47],
                    mean_woodlitter_gCm2 = output[,,47],
                    mean_annual_woodlitter_gCm2 = output_annual[,,47],
                    som_gCm2 = output[,,48],
                    mean_som_gCm2 = output[,,48],
                    mean_annual_som_gCm2 = output_annual[,,48],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,49],
                    mean_lai_m2m2 = output[,,49],
                    mean_annual_lai_m2m2 = output_annual[,,49],
                    gsi = output[,,50],
                    mean_gsi = output[,,50],
                    mean_annual_gsi = output_annual[,,50],
                    gsi_itemp = output[,,51],
                    mean_gsi_itemp = output[,,51],
                    mean_annual_gsi_itemp = output_annual[,,51],
                    gsi_iphoto = output[,,52],
                    mean_gsi_iphoto = output[,,52],
                    mean_annual_gsi_iphoto = output_annual[,,52],
                    gsi_ivpd = output[,,53],
                    mean_gsi_ivpd = output[,,53],
                    mean_annual_gsi_ivpd = output_annual[,,53],
                    # Photosynthesis / C~water coupling related
                    gs_demand_supply_ratio = output[,,54],
                    mean_gs_demand_supply_ratio = output[,,54],
                    mean_annual_gs_demand_supply_ratio = output_annual[,,54],
                    gs_mmolH2Om2s = output[,,55],
                    mean_gs_mmolH2Om2s = output[,,55],
                    mean_annual_gs_mmolH2Om2s = output_annual[,,55],
                    APAR_MJm2day = output[,,56],
                    mean_APAR_MJm2day = output[,,56],
                    mean_annual_APAR_MJm2day = output_annual[,,56],
                    gb_mmolH2Om2s = output[,,57],
                    mean_gb_mmolH2Om2s = output[,,57],
                    mean_annual_gb_mmolH2Om2s = output_annual[,,57],
                    CiCa = output[,,58],
                    mean_CiCa = output[,,58],
                    mean_annual_CiCa = output_annual[,,58],
                    # Misc
                    RootDepth_m = output[,,59],
                    mean_RootDepth_m = output[,,59],
                    mean_annual_RootDepth_m = output_annual[,,59],
                    ## Aggregated variables
                    # Mean Transit times
                    MTT_labile_years = MTT_years[,1],
                    MTT_foliage_years = MTT_years[,2],
                    MTT_roots_years = MTT_years[,3],
                    MTT_wood_years = MTT_years[,4],
                    MTT_litter_years = MTT_years[,5],
                    MTT_woodlitter_years = MTT_years[,6],
                    MTT_som_years = MTT_years[,7],
                    # Steady state estimates
                    SS_labile_gCm2 = SS_gCm2[,1],
                    SS_foliage_gCm2 = SS_gCm2[,2],
                    SS_roots_gCm2 = SS_gCm2[,3],
                    SS_wood_gCm2 = SS_gCm2[,4],
                    SS_litter_gCm2 = SS_gCm2[,5],
                    SS_woodlitter_gCm2 = SS_gCm2[,6],
                    SS_som_gCm2 = SS_gCm2[,7])
    # Determine the NPP fraction of expressed NPP
    # i.e. actual growth not GPP-Ra
    NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                         states_all$alloc_roots_gCm2day +
                         states_all$alloc_wood_gCm2day,1,mean)
    NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day,1,mean),
                         apply(states_all$alloc_roots_gCm2day,1,mean),
                         apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
    states_all$NPP_foliage_fraction = NPP_fraction[,1]
    states_all$NPP_roots_fraction = NPP_fraction[,2]
    states_all$NPP_wood_fraction = NPP_fraction[,3]
    # Tidy up variables
    rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.M2.#") { 
    output_dim = 49 ; MTT_dim = 5 ; SS_dim = 5
    # Load the required dalec shared object
    dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
    tmp=.Fortran( "rdalec16",output_dim=as.integer(output_dim)
                            ,MTT_dim=as.integer(MTT_dim)
                            ,SS_dim = as.integer(SS_dim)
                            ,met=as.double(t(met))
                            ,pars=as.double(pars_in)
                            ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                            ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                            ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                            ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                            ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                  
                            ,lat=as.double(lat)
                            ,nopars=as.integer(PROJECT$model$nopars[site])
                            ,nomet=as.integer(dim(met)[2])
                            ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                            ,nopools=as.integer(PROJECT$model$nopools[site])
                            ,nodays=as.integer(dim(met)[1])
                            ,nos_years=as.integer(noyears)
                            ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                            ,nos_iter=as.integer(nos_iter))
    # Extract the different output variables
    output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
    output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))              
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    mean_gpp_gCm2day = output_mean[,1],
                    mean_annual_gpp_gCm2day = output_annual[,,1],
                    rauto_gCm2day = output[,,2],
                    mean_rauto_gCm2day = output_mean[,2],
                    mean_annual_rauto_gCm2day = output_annual[,,2],
                    rhet_litter_gCm2day = output[,,3],
                    mean_rhet_litter_gCm2day = output_mean[,3],
                    mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                    rhet_som_gCm2day = output[,,4],
                    mean_rhet_som_gCm2day = output_mean[,4],
                    mean_annual_rhet_som_gCm2day = output_annual[,,4],
                    fire_gCm2day = output[,,5],
                    mean_fire_gCm2day = output_mean[,5],
                    mean_annual_fire_gCm2day = output_annual[,,5],
                    harvest_gCm2day = output[,,6],
                    mean_harvest_gCm2day = output_mean[,6],
                    mean_annual_harvest_gCm2day = output_annual[,,6],
                    grazing_gCm2day = output[,,7],
                    mean_grazing_gCm2day = output_mean[,7],
                    mean_annual_grazing_gCm2day = output_annual[,,7],
                    # Internal fluxes
                    alloc_foliage_gCm2day = output[,,8],
                    mean_alloc_foliage_gCm2day = output_mean[,8],
                    mean_annual_alloc_foliage_gCm2day = output_annual[,,8],
                    alloc_labile_gCm2day = output[,,9],
                    mean_alloc_labile_gCm2day = output_mean[,9],
                    mean_annual_alloc_labile_gCm2day = output_annual[,,9],
                    alloc_roots_gCm2day = output[,,10],
                    mean_alloc_roots_gCm2day = output_mean[,10],
                    mean_annual_alloc_roots_gCm2day = output_annual[,,10],
                    labile_to_foliage_gCm2day = output[,,11],
                    mean_labile_to_foliage_gCm2day = output_mean[,11],
                    mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                    foliage_to_litter_gCm2day = output[,,12],
                    mean_foliage_to_litter_gCm2day = output_mean[,12],
                    mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                    roots_to_litter_gCm2day = output[,,13],
                    mean_roots_to_litter_gCm2day = output_mean[,13],
                    mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                    litter_to_som_gCm2day = output[,,14],
                    mean_litter_to_som_gCm2day = output_mean[,14],
                    mean_annual_litter_to_som_gCm2day = output_annual[,,14],
                    # Disturbance fluxes
                    FIREemiss_labile_gCm2day = output[,,15],
                    mean_FIREemiss_labile_gCm2day = output_mean[,15],
                    mean_annual_FIREemiss_labile_gCm2day = output_annual[,,15],
                    FIRElitter_labile_gCm2day = output[,,16],
                    mean_FIRElitter_labile_gCm2day = output_mean[,16],
                    mean_annual_FIRElitter_labile_gCm2day = output_annual[,,16],
                    FIREemiss_foliage_gCm2day = output[,,17],
                    mean_FIREemiss_foliage_gCm2day = output_mean[,17],
                    mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,17],
                    FIRElitter_foliage_gCm2day = output[,,18],
                    mean_FIRElitter_foliage_gCm2day = output_mean[,18],
                    mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,18],
                    FIREemiss_roots_gCm2day = output[,,19],
                    mean_FIREemiss_roots_gCm2day = output_mean[,19],
                    mean_annual_FIREemiss_roots_gCm2day = output_annual[,,19],
                    FIRElitter_roots_gCm2day = output[,,20],
                    mean_FIRElitter_roots_gCm2day = output_mean[,20],
                    mean_annual_FIRElitter_roots_gCm2day = output_annual[,,20],
                    FIREemiss_litter_gCm2day = output[,,21],
                    mean_FIREemiss_litter_gCm2day = output_mean[,21],
                    mean_annual_FIREemiss_litter_gCm2day = output_annual[,,21],
                    FIRElitter_litter_gCm2day = output[,,22],
                    mean_FIRElitter_litter_gCm2day = output_mean[,22],
                    mean_annual_FIRElitter_litter_gCm2day = output_annual[,,22],
                    FIREemiss_som_gCm2day = output[,,23],
                    mean_FIREemiss_som_gCm2day = output_mean[,23],
                    mean_annual_FIREemiss_som_gCm2day = output_annual[,,23],
                    HARVESTextracted_labile_gCm2day = output[,,24],
                    mean_HARVESTextracted_labile_gCm2day = output_mean[,24],
                    mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,24],
                    HARVESTextracted_foliage_gCm2day = output[,,25],
                    mean_HARVESTextracted_foliage_gCm2day = output_mean[,25],
                    mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,25],
                    HARVESTextracted_roots_gCm2day = output[,,26],
                    mean_HARVESTextracted_roots_gCm2day = output_mean[,26],
                    mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,26],
                    HARVESTlitter_labile_gCm2day = output[,,27],
                    mean_HARVESTlitter_labile_gCm2day = output_mean[,27],
                    mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,27],
                    HARVESTlitter_foliage_gCm2day = output[,,28],
                    mean_HARVESTlitter_foliage_gCm2day = output_mean[,28],
                    mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,28],
                    HARVESTlitter_roots_gCm2day = output[,,29],
                    mean_HARVESTlitter_roots_gCm2day = output_mean[,29],
                    mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,29],
                    GRAZINGextracted_labile_gCm2day = output[,,30],
                    mean_GRAZINGextracted_labile_gCm2day = output_mean[,30],
                    mean_annual_GRAZINGextracted_labile_gCm2day = output_annual[,,30],
                    GRAZINGextracted_foliage_gCm2day = output[,,31],
                    mean_GRAZINGextracted_foliage_gCm2day = output_mean[,31],
                    mean_annual_GRAZINGextracted_foliage_gCm2day = output_annual[,,31],
                    GRAZINGextracted_roots_gCm2day = output[,,32],
                    mean_GRAZINGextracted_roots_gCm2day = output_mean[,32],
                    mean_annual_GRAZINGextracted_roots_gCm2day = output_annual[,,32],
                    GRAZINGlitter_labile_gCm2day = output[,,33],
                    mean_GRAZINGlitter_labile_gCm2day = output_mean[,33],
                    mean_annual_GRAZINGlitter_labile_gCm2day = output_annual[,,33],
                    GRAZINGlitter_foliage_gCm2day = output[,,34],
                    mean_GRAZINGlitter_foliage_gCm2day = output_mean[,34],
                    mean_annual_GRAZINGlitter_foliage_gCm2day = output_annual[,,34],
                    GRAZINGlitter_roots_gCm2day = output[,,35],
                    mean_GRAZINGlitter_roots_gCm2day = output_mean[,35],
                    mean_annual_GRAZINGlitter_roots_gCm2day = output_annual[,,35],
                    # Animal output (gC/m2/day)
                    animal_manure_to_soil_gCm2day = output[,,36],
                    mean_animal_manure_to_soil_gCm2day = output_mean[,36],
                    mean_annual_animal_manure_to_soil_gCm2day = output_annual[,,36],
                    animal_respiration_gCm2day = output[,,37],
                    mean_animal_respiration_gCm2day = output_mean[,37],
                    mean_annual_animal_respiration_gCm2day = output_annual[,,37],
                    animal_methane_gCm2day = output[,,38],
                    mean_animal_methane_gCm2day = output_mean[,38],
                    mean_annual_animal_methane_gCm2day = output_annual[,,38],
                    # C pools (gC/m2)
                    labile_gCm2 = output[,,39],
                    mean_labile_gCm2 = output_mean[,39],
                    mean_annual_labile_gCm2 = output_annual[,,39],
                    foliage_gCm2 = output[,,40],
                    mean_foliage_gCm2 = output_mean[,40],
                    mean_annual_foliage_gCm2 = output_annual[,,40],
                    roots_gCm2 = output[,,41],
                    mean_roots_gCm2 = output[,,41],
                    mean_annual_roots_gCm2 = output_annual[,,41],
                    litter_gCm2 = output[,,42],
                    mean_litter_gCm2 = output_mean[,42],
                    mean_annual_litter_gCm2 = output_annual[,,42],
                    som_gCm2 = output[,,43],
                    mean_som_gCm2 = output_mean[,43],
                    mean_annual_som_gCm2 = output_annual[,,43],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,44],
                    mean_lai_m2m2 = output_mean[,44],
                    mean_annual_lai_m2m2 = output_annual[,,44],
                    # Photosynthesis diagnostic
                    CiCa = output[,,45], 
                    mean_CiCa = output_mean[,45], 
                    mean_annual_CiCa = output_annual[,,45], 
                    # Canopy Phenology 
                    gsi = output[,,46],
                    mean_gsi = output_mean[,46],
                    mean_annual_gsi = output_annual[,,46],
                    gsi_itemp = output[,,47],
                    mean_gsi_itemp = output_mean[,47],
                    mean_annual_gsi_itemp = output_annual[,,47],
                    gsi_iphoto = output[,,48],
                    mean_gsi_iphoto = output_mean[,48],
                    mean_annual_gsi_iphoto = output_annual[,,48],
                    gsi_ivpd = output[,,49],
                    mean_gsi_ivpd = output_mean[,49],
                    mean_annual_gsi_ivpd = output_annual[,,49],
                    ## Aggregated variables
                    # Mean Transit times
                    MTT_labile_years = MTT_years[,1],
                    MTT_foliage_years = MTT_years[,2],
                    MTT_roots_years = MTT_years[,3],
                    MTT_litter_years = MTT_years[,4],
                    MTT_som_years = MTT_years[,5],
                    # Steady state estimates
                    SS_labile_gCm2 = SS_gCm2[,1],
                    SS_foliage_gCm2 = SS_gCm2[,2],
                    SS_roots_gCm2 = SS_gCm2[,3],
                    SS_litter_gCm2 = SS_gCm2[,4],
                    SS_som_gCm2 = SS_gCm2[,5])
    # Determine the NPP fraction of expressed NPP
    # i.e. actual growth not GPP-Ra                    
    NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                         states_all$alloc_roots_gCm2day +
                         states_all$alloc_foliage_gCm2day,1,mean)
    NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                         apply(states_all$alloc_roots_gCm2day,1,mean)) / NPP_fraction
    states_all$NPP_foliage_fraction = NPP_fraction[,1]
    states_all$NPP_roots_fraction = NPP_fraction[,2]
    # Tidy up variables
    rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A3.H2.M2.#") { 
    output_dim = 63 ; MTT_dim = 5 ; SS_dim = 5
    # Load the required dalec shared object
    dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
    tmp=.Fortran( "rdalec17",output_dim=as.integer(output_dim)
                            ,MTT_dim=as.integer(MTT_dim)
                            ,SS_dim = as.integer(SS_dim)
                            ,met=as.double(t(met))
                            ,pars=as.double(pars_in)
                            ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                            ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                            ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                            ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                            ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                  
                            ,lat=as.double(lat)
                            ,nopars=as.integer(PROJECT$model$nopars[site])
                            ,nomet=as.integer(dim(met)[2])
                            ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                            ,nopools=as.integer(PROJECT$model$nopools[site])
                            ,nodays=as.integer(dim(met)[1])
                            ,nos_years=as.integer(noyears)
                            ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                            ,nos_iter=as.integer(nos_iter))
    # Extract the different output variables
    output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
    output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))              
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    mean_gpp_gCm2day = output_mean[,1],
                    mean_annual_gpp_gCm2day = output_annual[,,1],
                    rauto_gCm2day = output[,,2],
                    mean_rauto_gCm2day = output_mean[,2],
                    mean_annual_rauto_gCm2day = output_annual[,,2],
                    rhet_litter_gCm2day = output[,,3],
                    mean_rhet_litter_gCm2day = output_mean[,3],
                    mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                    rhet_som_gCm2day = output[,,4],
                    mean_rhet_som_gCm2day = output_mean[,4],
                    mean_annual_rhet_som_gCm2day = output_annual[,,4],
                    fire_gCm2day = output[,,5],
                    mean_fire_gCm2day = output_mean[,5],
                    mean_annual_fire_gCm2day = output_annual[,,5],
                    harvest_gCm2day = output[,,6],
                    mean_harvest_gCm2day = output_mean[,6],
                    mean_annual_harvest_gCm2day = output_annual[,,6],
                    grazing_gCm2day = output[,,7],
                    mean_grazing_gCm2day = output_mean[,7],
                    mean_annual_grazing_gCm2day = output_annual[,,7],
                    # Internal fluxes
                    alloc_foliage_gCm2day = output[,,8],
                    mean_alloc_foliage_gCm2day = output_mean[,8],
                    mean_annual_alloc_foliage_gCm2day = output_annual[,,8],
                    alloc_labile_gCm2day = output[,,9],
                    mean_alloc_labile_gCm2day = output_mean[,9],
                    mean_annual_alloc_labile_gCm2day = output_annual[,,9],
                    alloc_roots_gCm2day = output[,,10],
                    mean_alloc_roots_gCm2day = output_mean[,10],
                    mean_annual_alloc_roots_gCm2day = output_annual[,,10],
                    labile_to_foliage_gCm2day = output[,,11],
                    mean_labile_to_foliage_gCm2day = output_mean[,11],
                    mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                    foliage_to_litter_gCm2day = output[,,12],
                    mean_foliage_to_litter_gCm2day = output_mean[,12],
                    mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                    roots_to_litter_gCm2day = output[,,13],
                    mean_roots_to_litter_gCm2day = output_mean[,13],
                    mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                    litter_to_som_gCm2day = output[,,14],
                    mean_litter_to_som_gCm2day = output_mean[,14],
                    mean_annual_litter_to_som_gCm2day = output_annual[,,14],
                    # Disturbance fluxes
                    FIREemiss_labile_gCm2day = output[,,15],
                    mean_FIREemiss_labile_gCm2day = output_mean[,15],
                    mean_annual_FIREemiss_labile_gCm2day = output_annual[,,15],
                    FIRElitter_labile_gCm2day = output[,,16],
                    mean_FIRElitter_labile_gCm2day = output_mean[,16],
                    mean_annual_FIRElitter_labile_gCm2day = output_annual[,,16],
                    FIREemiss_foliage_gCm2day = output[,,17],
                    mean_FIREemiss_foliage_gCm2day = output_mean[,17],
                    mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,17],
                    FIRElitter_foliage_gCm2day = output[,,18],
                    mean_FIRElitter_foliage_gCm2day = output_mean[,18],
                    mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,18],
                    FIREemiss_roots_gCm2day = output[,,19],
                    mean_FIREemiss_roots_gCm2day = output_mean[,19],
                    mean_annual_FIREemiss_roots_gCm2day = output_annual[,,19],
                    FIRElitter_roots_gCm2day = output[,,20],
                    mean_FIRElitter_roots_gCm2day = output_mean[,20],
                    mean_annual_FIRElitter_roots_gCm2day = output_annual[,,20],
                    FIREemiss_litter_gCm2day = output[,,21],
                    mean_FIREemiss_litter_gCm2day = output_mean[,21],
                    mean_annual_FIREemiss_litter_gCm2day = output_annual[,,21],
                    FIRElitter_litter_gCm2day = output[,,22],
                    mean_FIRElitter_litter_gCm2day = output_mean[,22],
                    mean_annual_FIRElitter_litter_gCm2day = output_annual[,,22],
                    FIREemiss_som_gCm2day = output[,,23],
                    mean_FIREemiss_som_gCm2day = output_mean[,23],
                    mean_annual_FIREemiss_som_gCm2day = output_annual[,,23],
                    HARVESTextracted_labile_gCm2day = output[,,24],
                    mean_HARVESTextracted_labile_gCm2day = output_mean[,24],
                    mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,24],
                    HARVESTextracted_foliage_gCm2day = output[,,25],
                    mean_HARVESTextracted_foliage_gCm2day = output_mean[,25],
                    mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,25],
                    HARVESTextracted_roots_gCm2day = output[,,26],
                    mean_HARVESTextracted_roots_gCm2day = output_mean[,26],
                    mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,26],
                    HARVESTlitter_labile_gCm2day = output[,,27],
                    mean_HARVESTlitter_labile_gCm2day = output_mean[,27],
                    mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,27],
                    HARVESTlitter_foliage_gCm2day = output[,,28],
                    mean_HARVESTlitter_foliage_gCm2day = output_mean[,28],
                    mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,28],
                    HARVESTlitter_roots_gCm2day = output[,,29],
                    mean_HARVESTlitter_roots_gCm2day = output_mean[,29],
                    mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,29],
                    GRAZINGextracted_labile_gCm2day = output[,,30],
                    mean_GRAZINGextracted_labile_gCm2day = output_mean[,30],
                    mean_annual_GRAZINGextracted_labile_gCm2day = output_annual[,,30],
                    GRAZINGextracted_foliage_gCm2day = output[,,31],
                    mean_GRAZINGextracted_foliage_gCm2day = output_mean[,31],
                    mean_annual_GRAZINGextracted_foliage_gCm2day = output_annual[,,31],
                    GRAZINGextracted_roots_gCm2day = output[,,32],
                    mean_GRAZINGextracted_roots_gCm2day = output_mean[,32],
                    mean_annual_GRAZINGextracted_roots_gCm2day = output_annual[,,32],
                    GRAZINGlitter_labile_gCm2day = output[,,33],
                    mean_GRAZINGlitter_labile_gCm2day = output_mean[,33],
                    mean_annual_GRAZINGlitter_labile_gCm2day = output_annual[,,33],
                    GRAZINGlitter_foliage_gCm2day = output[,,34],
                    mean_GRAZINGlitter_foliage_gCm2day = output_mean[,34],
                    mean_annual_GRAZINGlitter_foliage_gCm2day = output_annual[,,34],
                    GRAZINGlitter_roots_gCm2day = output[,,35],
                    mean_GRAZINGlitter_roots_gCm2day = output_mean[,35],
                    mean_annual_GRAZINGlitter_roots_gCm2day = output_annual[,,35],
                    # Animal output (gC/m2/day)
                    animal_manure_to_soil_gCm2day = output[,,36],
                    mean_animal_manure_to_soil_gCm2day = output_mean[,36],
                    mean_annual_animal_manure_to_soil_gCm2day = output_annual[,,36],
                    animal_respiration_gCm2day = output[,,37],
                    mean_animal_respiration_gCm2day = output_mean[,37],
                    mean_annual_animal_respiration_gCm2day = output_annual[,,37],
                    animal_methane_gCm2day = output[,,38],
                    mean_animal_methane_gCm2day = output_mean[,38],
                    mean_annual_animal_methane_gCm2day = output_annual[,,38],
                    # C pools (gC/m2)
                    labile_gCm2 = output[,,39],
                    mean_labile_gCm2 = output_mean[,39],
                    mean_annual_labile_gCm2 = output_annual[,,39],
                    foliage_gCm2 = output[,,40],
                    mean_foliage_gCm2 = output_mean[,40],
                    mean_annual_foliage_gCm2 = output_annual[,,40],
                    roots_gCm2 = output[,,41],
                    mean_roots_gCm2 = output[,,41],
                    mean_annual_roots_gCm2 = output_annual[,,41],
                    litter_gCm2 = output[,,42],
                    mean_litter_gCm2 = output_mean[,42],
                    mean_annual_litter_gCm2 = output_annual[,,42],
                    som_gCm2 = output[,,43],
                    mean_som_gCm2 = output_mean[,43],
                    mean_annual_som_gCm2 = output_annual[,,43],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,44],
                    mean_lai_m2m2 = output_mean[,44],
                    mean_annual_lai_m2m2 = output_annual[,,44],
                    # Photosynthesis diagnostic
                    CiCa = output[,,45], 
                    mean_CiCa = output_mean[,45], 
                    mean_annual_CiCa = output_annual[,,45], 
                    gs_demand_supply_ratio = output[,,46],
                    mean_gs_demand_supply_ratio = output_mean[,46],
                    mean_annual_gs_demand_supply_ratio = output_annual[,,46],
                    gs_mmolH2Om2s = output[,,47],
                    mean_gs_mmolH2Om2s = output_mean[,47],
                    mean_annual_gs_mmolH2Om2s = output_annual[,,47],
                    APAR_MJm2day = output[,,48],
                    mean_APAR_MJm2day = output_mean[,48],
                    mean_annual_APAR_MJm2day = output_annual[,,48],
                    gb_mmolH2Om2s = output[,,49],
                    mean_gb_mmolH2Om2s = output_mean[,49],
                    mean_annual_gb_mmolH2Om2s = output_annual[,,49],
                    # Water cycle related
                    ET_kgH2Om2day = output[,,50],
                    mean_ET_kgH2Om2day = output_mean[,50],
                    mean_annual_ET_kgH2Om2day = output_annual[,,50],
                    Etrans_kgH2Om2day = output[,,51],
                    mean_Etrans_kgH2Om2day = output_mean[,51],
                    mean_annual_Etrans_kgH2Om2day = output_annual[,,51],
                    Esoil_kgH2Om2day = output[,,52],
                    mean_Esoil_kgH2Om2day = output_mean[,52],
                    mean_annual_Esoil_kgH2Om2day = output_annual[,,52],
                    Ewetcanopy_kgH2Om2day = output[,,53],
                    mean_Ewetcanopy_kgH2Om2day = output_mean[,53],
                    mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,53],
                    runoff_kgH2Om2day = output[,,54],
                    mean_runoff_kgH2Om2day = output_mean[,54],
                    mean_annual_runoff_kgH2Om2day = output_annual[,,54],
                    underflow_kgH2Om2day = output[,,55],
                    mean_underflow_kgH2Om2day = output_mean[,55],
                    mean_annual_underflow_kgH2Om2day = output_annual[,,55],
                    SurfWater_kgH2Om2 = output[,,56],
                    mean_SurfWater_kgH2Om2 = output_mean[,56],
                    mean_annual_SurfWater_kgH2Om2 = output_annual[,,56],
                    wSWP_MPa = output[,,57],
                    mean_wSWP_MPa = output_mean[,57],
                    mean_annual_wSWP_MPa = output_annual[,,57],
                    snow_kgH2Om2 = output[,,58],
                    mean_snow_kgH2Om2 = output_mean[,58],
                    mean_annual_snow_kgH2Om2 = output_annual[,,58],
                    # Misc
                    RootDepth_m = output[,,59],
                    mean_RootDepth_m = output_mean[,59],
                    mean_annual_RootDepth_m = output_annual[,,59],
                    # Canopy Phenology 
                    gsi = output[,,60],
                    mean_gsi = output_mean[,60],
                    mean_annual_gsi = output_annual[,,60],
                    gsi_itemp = output[,,61],
                    mean_gsi_itemp = output_mean[,61],
                    mean_annual_gsi_itemp = output_annual[,,61],
                    gsi_iphoto = output[,,62],
                    mean_gsi_iphoto = output_mean[,62],
                    mean_annual_gsi_iphoto = output_annual[,,62],
                    gsi_ivpd = output[,,63],
                    mean_gsi_ivpd = output_mean[,63],
                    mean_annual_gsi_ivpd = output_annual[,,63],
                    ## Aggregated variables
                    # Mean Transit times
                    MTT_labile_years = MTT_years[,1],
                    MTT_foliage_years = MTT_years[,2],
                    MTT_roots_years = MTT_years[,3],
                    MTT_litter_years = MTT_years[,4],
                    MTT_som_years = MTT_years[,5],
                    # Steady state estimates
                    SS_labile_gCm2 = SS_gCm2[,1],
                    SS_foliage_gCm2 = SS_gCm2[,2],
                    SS_roots_gCm2 = SS_gCm2[,3],
                    SS_litter_gCm2 = SS_gCm2[,4],
                    SS_som_gCm2 = SS_gCm2[,5])
    # Determine the NPP fraction of expressed NPP
    # i.e. actual growth not GPP-Ra                    
    NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                         states_all$alloc_roots_gCm2day +
                         states_all$alloc_foliage_gCm2day,1,mean)
    NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                         apply(states_all$alloc_roots_gCm2day,1,mean)) / NPP_fraction
    states_all$NPP_foliage_fraction = NPP_fraction[,1]
    states_all$NPP_roots_fraction = NPP_fraction[,2]
    # Tidy up variables
    rm(output,MTT_years,SS_gCm2)
  } else {
      stop(paste("Model choice (",model_name,") does not have corresponding R interface",sep=""))
  }

  # return state variable means
  return(states_all) ; gc(verbose=FALSE)

} # end of function
## Use byte compile
simulate_all<-cmpfun(simulate_all)
