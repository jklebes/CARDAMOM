
###
## Function to extract obs needed for CARDAMOM
###

extract_obs<-function(latlon_wanted,lai_all,Csom_all,forest_all,Cwood_all,sand_clay_all,crop_man_all
                     ,burnt_all,soilwater_all,ctessel_pft,site_name,start_year,end_year,timestep_days
                     ,spatial_type,resolution,grid_type,modelname) {

    ###
    ## Get some LAI information (m2/m2)
    ###

    if (lai_source == "MODIS") {

      # Extract lai and uncertainty information
      # NOTE: assume default uncertainty (+/- scale)
      output = extract_modis_lai(timestep_days,spatial_type,resolution,grid_type,latlon_wanted,lai_all,as.numeric(start_year):as.numeric(end_year))
      lai = output$lai ; lai_unc = output$lai_unc

    } else if (lai_source == "COPERNICUS") {

      # Extract lai and uncertainty information
      # NOTE: assume default uncertainty (+/- scale)
      output = extract_copernicus_lai(timestep_days,spatial_type,resolution,grid_type,latlon_wanted,lai_all,as.numeric(start_year):as.numeric(end_year))
      lai = output$lai ; lai_unc = output$lai_unc

    } else if (lai_source == "site_specific") {

      # read from .csv or netcdf
      infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
      lai = read_site_specific_obs("LAI",infile) ; lai_unc = read_site_specific_obs("LAI_unc",infile)
      if (max(lai_unc) == -9999) {
        lai_unc = rep(-9999,times = length(lai))
        # borrowed linear approximtion of uncertainty form Copernicus
        lai_unc[which(lai != -9999)] = 0.25 #lai[which(lai != -9999)] * 0.14028508 + 1.0
      }

    } else {

      lai = -9999
      lai_unc = -9999

    }
    # Add model structural uncertainty to the uncertainty estimate if present
    # 0.5 estimates from comparison of LAI uncertainties trials at 0.5 and 0.25, resultant CI in both instances is +/- ~0.5
    lai_unc[lai_unc > 0] = lai_unc[lai_unc > 0] + 0.5

    ###
    ## Get some Cfoliage information (stock; gC/m2)
    ###

    if (Cfol_stock_source == "site_specific") {
      infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
      Cfol_stock=read_site_specific_obs("Cfol_stock",infile)
      Cfol_stock_unc=read_site_specific_obs("Cfol_stock_unc",infile)
      if (length(Cfol_stock_unc) == 1) {
        # on the other hand if not then we have no uncertainty info, so use default
        Cfol_stock_unc = rep(-9999,times = length(Cfol_stock))
        # See Smallman et al., (2017) for uncertainty estimate
        Cfol_stock_unc[which(Cfol_stock > 0)] = 0.38 * Cfol_stock[which(Cfol_stock > 0)]
      }
    } else {
      # assume no data available
      Cfol_stock=-9999 ; Cfol_stock_unc=-9999
    }

    ###
    ## Get some initial Csom (gC/m2) information
    ###

    if (Csom_source == "HWSD") {
        # could add other variables such as SOM (gC.m-2)
        Csom_initial = extract_hwsd_Csom(spatial_type,resolution,grid_type,latlon_wanted,Csom_all)
        Csom_initial_unc = Csom_initial * 0.47 # see papers assessing uncertainty of HWSD, ~47 %
    } else if (Csom_source == "SoilGrids") {
        Csom_info = extract_soilgrid_Csom(spatial_type,resolution,grid_type,latlon_wanted,Csom_all)
        Csom_initial = Csom_info$Csom_initial ; Csom_initial_unc = Csom_info$Csom_initial_unc
    } else if (Csom_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        Csom_initial = read_site_specific_obs("Csom_initial",infile)
        Csom_initial_unc = read_site_specific_obs("Csom_initial_unc",infile)
        if (Csom_initial_unc == -9999 & Csom_initial > 0) {
            # on the other hand if not then we have no uncertainty info, so use default
            Csom_initial_unc = 0.24 * Csom_initial
        }
    } else {
        # assume no data available
        Csom_initial = -9999 ; Csom_initial_unc = -9999
    }
    # Now assuming we have actual information we need to add the model structural uncertainty.
    # A structural of uncertainty has been estimates at ~ 1000 gC/m2 based on Smallman et al., (2017)
    if (Csom_initial != 9999) {Csom_initial_unc = Csom_initial_unc + 1000}

    ###
    ## Get some sand / clay information (%)
    ###

    if (sand_clay_source == "HWSD") {
        sand_clay=extract_hwsd_sand_clay(spatial_type,resolution,grid_type,latlon_wanted,sand_clay_all)
        top_sand = sand_clay$top_sand ; bot_sand = sand_clay$bot_sand
        top_clay = sand_clay$top_clay ; bot_clay = sand_clay$bot_clay
    } else if (sand_clay_source == "SoilGrids") {
        sand_clay=extract_soilgrid_sand_clay(spatial_type,resolution,grid_type,latlon_wanted,sand_clay_all)
        top_sand = sand_clay$top_sand ; bot_sand = sand_clay$bot_sand
        top_clay = sand_clay$top_clay ; bot_clay = sand_clay$bot_clay 
    } else if (sand_clay_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        top_sand = read_site_specific_obs("top_sand_initial",infile)
        bot_sand = read_site_specific_obs("bot_sand_initial",infile)
        top_clay = read_site_specific_obs("top_clay_initial",infile)
        bot_clay = read_site_specific_obs("bot_clay_initial",infile)
    } else {
        # assume no data available
        top_sand = 40 ; bot_sand = 40
        top_clay = 15 ; bot_clay = 15
    }

    ###
    ## Get some crop management information (day)
    ###

    if (ctessel_pft == 1 & crop_management_source == "sacks_crop_calendar") {
      # could add other variables such as SOM (gC.m-2)
      crop_dates=extract_sacks_crop_info(spatial_type,resolution,grid_type,latlon_wanted,crop_man_all)
      plant = crop_dates$plant ; plant_range = crop_dates$plant_range
      harvest = crop_dates$harvest ; harvest_range = crop_dates$harvest_range
    } else if (crop_management_source == "site_specific") {
      infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
      plant = read_site_specific_obs("plant_initial",infile)
      plant_range = read_site_specific_obs("plant_range_initial",infile)
      harvest = read_site_specific_obs("harvest_initial",infile)
      harvest_range = read_site_specific_obs("harvest_range_initial",infile)
    } else {
      # assume no data available
      plant = 304 ; plant_range = 15 # days
      harvest = 208 ; harvest_range = 15 # days
    }

    ###
    ## Get some Wood increment information (gC/m2/yr; time series)
    ###

    if (woodinc_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        woodinc=read_site_specific_obs("woodinc",infile)
        woodinc_unc=read_site_specific_obs("woodinc_unc",infile)
        if (length(woodinc_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            woodinc_unc = rep(-9999,times = length(woodinc))
            woodinc_unc[which(woodinc > 0)] = 0.25 * woodinc[which(woodinc > 0)]
        }
    } else {
      # assume no data available
      woodinc=-9999 ; woodinc_unc=-9999
    }

    ###
    ## Get some GPP information (time series; gC/m2/day)
    ###

    if (GPP_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")}
        #	if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater.csv",sep="")}
        if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater_copy.csv",sep="")}
        GPP = read_site_specific_obs("GPP",infile)
        GPP_unc = read_site_specific_obs("GPP_unc",infile)
        if (length(GPP_unc) == 1) {
            GPP_unc = rep(-9999,times = length(GPP))
            GPP_unc[which(GPP > 0)] = 0.5 * GPP[which(GPP > 0)] + 0.5
            if (modelname == "ACM") {GPP_unc = rep(mean(GPP*0.5),times=length(GPP))}
        }
    } else {
        # assume no data available
        GPP = -9999 ; GPP_unc = -9999
    }

    ###
    ## Get some Evapotranspiration information (time series; kgH2O/m2/day)
    ###

    if (Evap_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")}
        #if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater.csv",sep="")}
        if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater_copy.csv",sep="")}
        Evap = read_site_specific_obs("Evap",infile)
        Evap_unc = read_site_specific_obs("Evap_unc",infile)
        if (length(Evap_unc) == 1) {
            Evap_unc = rep(-9999,times = length(Evap))
            Evap_unc[which(Evap > 0)] = 0.5 * Evap[which(Evap > 0)]
        }
        if (modelname == "ACM") {
            # borrow woody increment for soil evaporation in ACM_ET recalibration
            woodinc = read_site_specific_obs("soilevap",infile)
            # borrow Cfol_stock for wet canopy evaporation in ACM_ET recalibration
            Cfol_stock = read_site_specific_obs("wetevap",infile)
            # actually lets make uncertainty half mean of total ET
            Evap_unc = rep(mean(Evap+woodinc+Cfol_stock)*0.5,times=length(Evap))
            woodinc_unc = Evap_unc
            Cfol_stock_unc = Evap_unc
        }
    } else {
        # assume no data available
        Evap = -9999 ; Evap_unc = -9999
    }

    ###
    ## Get some Reco information (time series; gC/m2/day)
    ###

    if (Reco_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Reco=read_site_specific_obs("Reco",infile)
        Reco_unc=read_site_specific_obs("Reco_unc",infile)
        if (length(Reco_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Reco_unc = rep(-9999,times = length(Reco))
            Reco_unc[which(Reco > 0)] = 0.50 * Reco[which(Reco > 0)] + 0.5
        }
    } else {
        # assume no data available
        Reco = -9999 ; Reco_unc = -9999
    }

    ###
    ## Get some NEE information (time series; gC/m2/day)
    ###

    if (NEE_source == "site_specific") {
      infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
      NEE=read_site_specific_obs("NEE",infile)
      NEE_unc=read_site_specific_obs("NEE_unc",infile)
      if (length(NEE_unc) == 1) {
        # on the other hand if not then we have no uncertainty info, so use default
        NEE_unc = rep(-9999,times = length(NEE))
        NEE_unc[which(NEE != -9999)] = 1.5
      }
    } else {
      # assume no data available
      NEE = -9999 ; NEE_unc = -9999
    }

    ###
    ## Get some Cfoliage information (initial conditions)
    ###

    if (Cfol_initial_source == "site_specific") {
      infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
      Cfol_initial=read_site_specific_obs("Cfol_initial",infile)
      Cfol_initial_unc=read_site_specific_obs("Cfol_initial_unc",infile)
      if (Cfol_initial_unc == -9999 & Cfol_initial > 0) {
        # on the other hand if not then we have no uncertainty info, so use default
        Cfol_initial_unc = 0.25 * Cfol_initial
      }
    } else {
      # assume no data available
      Cfol_initial = -9999 ; Cfol_initial_unc = -9999
    }

    ###
    ## Get some Cwood information (initial conditions)
    ###

    if (Cwood_initial_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        Cwood_initial=read_site_specific_obs("Cwood_initial",infile)
        Cwood_initial_unc=read_site_specific_obs("Cwood_initial_unc",infile)
        if (Cwood_initial_unc == -9999 & Cwood_initial > 0) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cwood_initial_unc = 0.25 * Cwood_initial
        }
    } else if (Cwood_initial_source == "GlobBIOMASS") {
        # get Cwood
        output = extract_globbiomass_biomass(timestep_days,spatial_type,resolution,grid_type,latlon_wanted,Cwood_all)
        Cwood_initial = output$Cwood_stock
        Cwood_initial_unc = output$Cwood_stock_unc
    } else if (Cwood_initial_source == "Avitabile") {
        # get Cwood
        output = extract_avitabile_biomass(timestep_days,spatial_type,resolution,grid_type,latlon_wanted,Cwood_all)
        Cwood_initial = output$Cwood_stock
        Cwood_initial_unc = output$Cwood_stock_unc
    } else if (Cwood_initial_source == "mpi_biomass") {
        # get Cwood
        output = extract_mpi_biomass(timestep_days,spatial_type,resolution,grid_type,latlon_wanted,Cwood_all)
        # locate time time data
        tmp = which(output$Cwood_stock == max(output$Cwood_stock,na.rm=TRUE))[1]
        if (length(tmp) > 0) {
            Cwood_initial = output$Cwood_stock[tmp]
            Cwood_initial_unc = output$Cwood_stock_unc[tmp]
        } else {
            Cwood_initial = -9999 ; Cwood_initial_unc = -9999
        }
    } else if (Cwood_initial_source == "UoL") {
        # this is a very bespoke modification so leave it here to avoid getting lost
        agb = raster(paste(path_to_biomass,"Kenya_0.25deg_AGB_stable_forest_2015_2017.tif", sep=""))
        unc = raster(paste(path_to_biomass,"Kenya_0.25deg_AGB_std_stable_forest_2015_2017.tif", sep=""))
        # extract dimension information for the grid, note the axis switching between raster and actual array
        xdim = dim(agb)[2] ; ydim = dim(agb)[1]
        # extract the lat / long information needed
        longitude = coordinates(agb)[,1] ; latitude = coordinates(agb)[,2]
        # restructure into correct orientation
        longitude = array(longitude, dim=c(xdim,ydim))
        latitude = array(latitude, dim=c(xdim,ydim))
        # break out from the rasters into arrays which we can manipulate
        agb = array(as.vector(unlist(agb)), dim=c(xdim,ydim))
        unc = array(as.vector(unlist(unc)), dim=c(xdim,ydim))
        output = closest2d(1,latitude,longitude,latlon_wanted[1],latlon_wanted[2],2)
        i1 = unlist(output)[1] ; j1 = unlist(output)[2]
        if (is.na(agb[i1,j1]) | is.na(unc[i1,j1])) {
           Cwood_initial = -9999 ; Cwood_initial_unc = -9999
        } else  {
           Cwood_initial = agb[i1,j1] ; Cwood_initial_unc = unc[i1,j1]
        }
    } else {
        # assume no data available
        Cwood_initial=-9999 ; Cwood_initial_unc=-9999
    }

    ###
    ## Get some Croots information (initial conditions)
    ###

    if (Croots_initial_source == "site_specific") {
      infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
      Croots_initial=read_site_specific_obs("Croots_initial",infile)
      Croots_initial_unc=read_site_specific_obs("Croots_initial_unc",infile)
      if (Croots_initial_unc == -9999 & Croots_initial > 0) {
        # on the other hand if not then we have no uncertainty info, so use default
        Croots_initial_unc = 0.44 * Croots_initial
      }
    } else {
      # assume no data available
      Croots_initial = -9999 ; Croots_initial_unc = -9999
    }

    ###
    ## Get some Clitter information (initial conditions)
    ###

    if (Clit_initial_source == "site_specific") {
      infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
      Clit_initial=read_site_specific_obs("Clit_initial",infile)
      Clit_initial_unc=read_site_specific_obs("Clit_initial_unc",infile)
      if (Clit_initial_unc == -9999 & Clit_initial > 0) {
        # on the other hand if not then we have no uncertainty info, so use default
        Clit_initial_unc = 0.25 * Clit_initial
      }
    } else {
      # assume no data available
      Clit_initial = -9999 ; Clit_initial_unc = -9999
    }

    ###
    ## Get some Cwood information (stock)
    ###

    if (Cwood_stock_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Cwood_stock=read_site_specific_obs("Cwood_stock",infile)
        Cwood_stock_unc=read_site_specific_obs("Cwood_stock_unc",infile)
        if (length(Cwood_stock_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cwood_stock_unc = rep(-9999,times = length(Cwood_stock))
            Cwood_stock_unc[which(Cwood_stock != -9999)] = abs(0.25 * Cwood_stock[which(Cwood_stock != -9999)])
        }
    } else if (Cwood_stock_source == "GlobBIOMASS") {
        # declare output variable
        Cwood_stock=array(-9999, dim=Cwood_all$step_of)
        Cwood_stock_unc=array(-9999, dim=Cwood_all$step_of)
        # only bother with this if 2010 is within time period
        if (Cwood_all$obs_step > 0) {
            # get Cwood
            output = extract_globbiomass_biomass(timestep_days,spatial_type,resolution,grid_type,latlon_wanted,Cwood_all)
            # and insert the extracted value into the correct location
            Cwood_stock[Cwood_all$obs_step]=output$Cwood_stock
            Cwood_stock_unc[Cwood_all$obs_step]=output$Cwood_stock_unc
        }

    } else if (Cwood_stock_source == "mpi_biomass") {
        # get Cwood
        output=extract_mpi_biomass(timestep_days,spatial_type,resolution,grid_type,latlon_wanted,Cwood_all)
        Cwood_stock=output$Cwood_stock
        Cwood_stock_unc=output$Cwood_stock_unc
    } else {
        # assume no data available
        Cwood_stock = -9999 ; Cwood_stock_unc = -9999
    }

    ###
    ## Get some Cagb information (stock)
    ###

    if (Cagb_stock_source == "site_specific") {
      infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
      Cagb_stock=read_site_specific_obs("Cagb_stock",infile)
      Cagb_stock_unc=read_site_specific_obs("Cagb_stock_unc",infile)
      if (length(Cagb_stock_unc) == 1) {
        # on the other hand if not then we have no uncertainty info, so use default
        Cagb_stock_unc = rep(-9999,times = length(Cagb_stock))
        Cagb_stock_unc[which(Cagb_stock != -9999)] = abs(0.25 * Cagb_stock[which(Cagb_stock != -9999)])
      }
    } else {
      # assume no data available
      Cagb_stock = -9999 ; Cagb_stock_unc = -9999
    }

    ###
    ## Get some Croots information (stock)
    ###

    if (Croots_stock_source == "site_specific") {
      infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
      Croots_stock=read_site_specific_obs("Croots_stock",infile)
      Croots_stock_unc=read_site_specific_obs("Croots_stock_unc",infile)
      if (length(Croots_stock_unc) == 1) {
        # on the other hand if not then we have no uncertainty info, so use default
        Croots_stock_unc = rep(-9999,times = length(Croots_stock))
        Croots_stock_unc[which(Croots_stock != -9999)] = abs(0.44 * Croots_stock[which(Croots_stock != -9999)])
      }
    } else {
      # assume no data available
      Croots_stock = -9999 ; Croots_stock_unc = -9999
    }

    ###
    ## Get some Clitter information (stock)
    ###
    #    print("checking Clitter")
    if (Clit_stock_source == "site_specific") {
      infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
      Clit_stock=read_site_specific_obs("Clit_stock",infile)
      Clit_stock_unc=read_site_specific_obs("Clit_stock_unc",infile)
      if (length(Clit_stock_unc) == 1) {
        # on the other hand if not then we have no uncertainty info, so use default
        Clit_stock_unc = rep(-9999,times = length(Clit_stock))
        Clit_stock_unc[which(Clit_stock != -9999)] = abs(0.38 * Clit_stock[which(Clit_stock != -9999)])
      }
    } else {
      # assume no data available
      Clit_stock = -9999 ; Clit_stock_unc = -9999
    }

    ###
    ## Get some Csom information (stock)
    ###
    #    print("checking Csom")
    if (Csom_stock_source == "site_specific") {
      infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
      Csom_stock=read_site_specific_obs("Csom_stock",infile)
      Csom_stock_unc=read_site_specific_obs("Csom_stock_unc",infile)
      if (length(Csom_stock_unc) == 1) {
        # on the other hand if not then we have no uncertainty info, so use default
        Csom_stock_unc = rep(-9999,times = length(Csom_stock))
        Csom_stock_unc[which(Csom_stock != -9999)] = abs(0.24 * Csom_stock[which(Csom_stock != -9999)])
      }
    } else {
      # assume no data available
      Csom_stock = -9999 ; Csom_stock_unc = -9999
    }
    # Now assuming we have actual information we need to add the model structural uncertainty.
    # A structural of uncertainty has been estimates at ~ 1000 gC/m2 based on Smallman et al., (2017)
    if (Csom_stock != 9999) {Csom_stock_unc = Csom_stock_unc + 1000}

    ###
    ## Get some Cstem information (stock)
    ###

    if (Cstem_stock_source == "site_specific") {
      infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
      Cstem_stock=read_site_specific_obs("Cstem_stock",infile)
      Cstem_stock_unc=read_site_specific_obs("Cstem_stock_unc",infile)
      if (length(Cstem_stock_unc) == 1) {
        # on the other hand if not then we have no uncertainty info, so use default
        Cstem_stock_unc = rep(-9999,times = length(Cstem_stock))
        Cstem_stock_unc[which(Cstem_stock != -9999)] = abs(0.24 * Cstem_stock[which(Cstem_stock != -9999)])
      }
    } else {
      # assume no data available
      Cstem_stock = -9999 ; Cstem_stock_unc = -9999
    }

    ###
    ## Get some Cbranch information (stock)
    ###

    if (Cbranch_stock_source == "site_specific") {
      infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
      Cbranch_stock=read_site_specific_obs("Cbranch_stock",infile)
      Cbranch_stock_unc=read_site_specific_obs("Cbranch_stock_unc",infile)
      if (length(Cbranch_stock_unc) == 1) {
        # on the other hand if not then we have no uncertainty info, so use default
        Cbranch_stock_unc = rep(-9999,times = length(Cbranch_stock))
        Cbranch_stock_unc[which(Cbranch_stock != -9999)] = abs(0.24 * Cbranch_stock[which(Cbranch_stock != -9999)])
      }
    } else {
      # assume no data available
      Cbranch_stock=-9999 ; Cbranch_stock_unc=-9999
    }

    ###
    ## Get some Ccoarseroot information (stock)
    ###

    if (Ccoarseroot_stock_source == "site_specific") {
      infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
      Ccoarseroot_stock=read_site_specific_obs("Ccoarseroot_stock",infile)
      Ccoarseroot_stock_unc=read_site_specific_obs("Ccoarseroot_stock_unc",infile)
      if (length(Ccoarseroot_stock_unc) == 1) {
        # on the other hand if not then we have no uncertainty info, so use default
        Ccoarseroot_stock_unc = rep(-9999,times = length(Ccoarseroot_stock))
        Ccoarseroot_stock_unc[which(Ccoarseroot_stock != -9999)] = abs(0.24 * Ccoarseroot_stock[which(Ccoarseroot_stock != -9999)])
      }
    } else {
      # assume no data available
      Ccoarseroot_stock = -9999 ; Ccoarseroot_stock_unc = -9999
    }

    ###
    ## Get some Cfolmax information (stock)
    ###

    if (Cfolmax_stock_source == "site_specific") {
      infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
      Cfolmax_stock=read_site_specific_obs("Cfolmax_stock",infile)
      Cfolmax_stock_unc=read_site_specific_obs("Cfolmax_stock_unc",infile)
      if (length(Cfolmax_stock_unc) == 1) {
        # on the other hand if not then we have no uncertainty info, so use default
        Cfolmax_stock_unc = rep(-9999,times = length(Cfolmax_stock))
        Cfolmax_stock_unc[which(Cfolmax_stock != -9999)] = abs(0.24 * Cfolmax_stock[which(Cfolmax_stock != -9999)])
      }
    } else {
      # assume no data available
      Cfolmax_stock = -9999 ; Cfolmax_stock_unc = -9999
    }

    ###
    ## Get some deforestation information (fraction time series)
    ###
    #    print("checking deforestation")
    if (deforestation_source == "site_specific") {
      infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
      deforestation = read_site_specific_obs("deforestation",infile)
      forest_management = read_site_specific_obs("management",infile)
      yield_class = 0 #read_site_specific_obs("yield_class",infile)
      age = read_site_specific_obs("age",infile)
      age = age[1] # we only want the age at the beginning of the simulation
    } else if (deforestation_source == "combined_dataset" | deforestation_source == "GFW") {
      output = extract_forestry_information(timestep_days,spatial_type,resolution,grid_type,latlon_wanted,forest_all,start_year,end_year,ctessel_pft)
      ctessel_pft = output$ctessel_pft
      deforestation = output$deforestation
      yield_class = output$yield_class
      age = output$age
      forest_management = 2
    } else {
      # assume no data available
      deforestation = 0
      forest_management = 1
      yield_class = 0
      age = -9999
      #	if (modelname == "DALEC_GSI_FR" | modelname == "DALECN_GSI_FR") {age=30}
    }

    ###
    ## Get some burnt area information (fraction time series)
    ###

    if (burnt_area_source == "site_specific") {
      infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
      burnt_area=read_site_specific_obs("burnt_area",infile)
    } else if (burnt_area_source == " "){
      # assume no data available
      burnt_area = 0
    } else {
      burnt_area = extract_burnt_area_information(latlon_wanted,timestep_days,spatial_type,grid_type,resolution,start_year,end_year,burnt_all)
    }

    ###
    ## Get some snow water equivalent (mm/day)
    ###

    if (snow_source == "site_specific") {
      infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
      SWE=read_site_specific_obs("snow_water",infile)
      SWE_unc=read_site_specific_obs("snow_water_unc",infile)
      if (length(which(SWE_unc != -9999)) == 0) {
        # on the other hand if not then we have no uncertainty info, so use default
        SWE_unc=rep(sd(SWE[which(SWE != -9999)],na.rm=TRUE),length.out=length(SWE))
      }
    } else {
      # assume no data available
      SWE = -9999 ; SWE_unc = -9999
    }

    ###
    ## Get initial soil water fraction prior (initial conditions)
    ###

    if (soilwater_initial_source == "site_specific") {
      infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
      soilwater=read_site_specific_obs("soil_water",infile)
      soilwater_unc=read_site_specific_obs("soil_water_unc",infile)
      if (soilwater_unc == -9999 & soilwater > 0) {
        # on the other hand if not then we have no uncertainty info, so use default
        soilwater_unc = 0.10 * soilwater
      }
    } else if (soilwater_initial_source == "GLEAM") {
      output = extract_soilwater_initial(spatial_type,resolution,grid_type,latlon_wanted,soilwater_all)
      soilwater = output$soil_water ; soilwater_unc = output$soil_water_unc
    } else {
      # assume no data available
      soilwater = -9999 ; soilwater_unc = -9999
    }

    # return output now
    return(list(LAT = latlon_wanted[1], LAI = lai, LAI_unc = lai_unc, GPP = GPP, GPP_unc = GPP_unc
      ,Evap = Evap, Evap_unc = Evap_unc, NEE = NEE, NEE_unc = NEE_unc, Reco = Reco, Reco_unc = Reco_unc
      ,woodinc = woodinc, woodinc_unc = woodinc_unc, Cfol_stock = Cfol_stock, Cfol_stock_unc = Cfol_stock_unc
      ,Cwood_stock = Cwood_stock,Cwood_stock_unc = Cwood_stock_unc, Cagb_stock=Cagb_stock, Cagb_stock_unc = Cagb_stock_unc
      ,Croots_stock = Croots_stock, Croots_stock_unc = Croots_stock_unc, Clit_stock = Clit_stock, Clit_stock_unc = Clit_stock_unc
      ,Csom_stock = Csom_stock, Csom_stock_unc = Csom_stock_unc, Cstem_stock = Cstem_stock, Cstem_stock_unc = Cstem_stock_unc
      ,Cbranch_stock = Cbranch_stock, Cbranch_stock_unc = Cbranch_stock_unc, Ccoarseroot_stock = Ccoarseroot_stock
      ,Ccoarseroot_stock_unc = Ccoarseroot_stock_unc, Cfolmax_stock = Cfolmax_stock, Cfolmax_stock_unc = Cfolmax_stock_unc
      ,Csom_initial = Csom_initial, Csom_initial_unc = Csom_initial_unc, Cfol_initial = Cfol_initial, Cfol_initial_unc = Cfol_initial_unc
      ,Cwood_initial = Cwood_initial, Cwood_initial_unc = Cwood_initial_unc, Croots_initial = Croots_initial
      ,Croots_initial_unc = Croots_initial_unc, Clit_initial = Clit_initial, Clit_initial_unc = Clit_initial_unc
      ,deforestation = deforestation, burnt_area = burnt_area, ctessel_pft = ctessel_pft, yield_class = yield_class
      ,age = age, forest_management = forest_management, top_sand = top_sand, bot_sand = bot_sand, top_clay = top_clay
      ,bot_clay = bot_clay, plant = plant, plant_range = plant_range, harvest = harvest, harvest_range = harvest_range
      ,SWE = SWE, SWE_unc = SWE_unc, soilwater = soilwater, soilwater_unc = soilwater_unc))


    }
