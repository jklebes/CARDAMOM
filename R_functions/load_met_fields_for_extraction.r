
###
## Function to load met data from global field ECMWF data
## subsequently extracted in extract_met_drivers.txt
###

load_met_fields_for_extraction<-function(latlon_in,met_source,modelname,startyear,endyear) {

    # let the user know this might take some time
    print("Loading global met fields for subsequent sub-setting ...")
    # declare timing variables
    t_grid = 0

    if (met_source == "site_specific") {

        # contruct output
        met_all=list(site_specific=TRUE)

    } else {

        if (met_source == "CRUJRA") {

            # declare variable ids needed to select files / infile variables
            varid = c("dswrf","tmax","pre","spfh","pres","tmin","wsp")
            infile_varid = c("dswrf","tmax","pre","spfh","pres","tmin","wsp")

            # open first ecmwf file to extract needed information
            input_file_1 = paste(path_to_met_source,"/",varid[1],"/crujra.V1.1.5d.",varid[1],".",startyear,".365d.noc.nc",sep="")
            data1 = nc_open(input_file_1)

            # Get timing variable
            # The native time step is 6 hours, but we directly average in the load_met_function.r
            steps_in_day = 1

            # extract location variables
            lat = ncvar_get(data1, "lat") ; long = ncvar_get(data1, "lon")
            # expand the one directional values here into 2 directional
            lat_dim = length(lat) ; long_dim = length(long)
            long = array(long,dim=c(long_dim,lat_dim))
            lat = array(lat,dim=c(lat_dim,long_dim)) ; lat=t(lat)

            # If we are using the GSI model we need the 21 days (or month) before the start date of the simulation, so we need to check if we have this information
            # NOTE that this section of code is duplicated for each of the available datasets because of differences in storage and file name
            extra_year = FALSE
            present = 0
            for (lag in seq(1,length(varid))) {
                 input_file_1 = paste(path_to_met_source,"/",varid[1],"/crujra.V1.1.5d.",varid[1],".",as.character(as.numeric(startyear)-1),".365d.noc.nc",sep="")
                 # check whether the files exist or not
                 if (file.exists(input_file_1)) {present = present+1}
            }
            # if all the files exist then we will use them
            if (present == length(varid)) {extra_year=TRUE}

        } else if (met_source == "ERA") {

            # declare variable ids needed to select files / infile variables
            varid=c("sw_radiation_daily_mean","airt_daily_max","precipitation_daily_mean","vpd_daily_mean","sfc_pressure_daily_mean","airt_daily_min","wind_spd_daily_mean")
            infile_varid=c("daily_swrad","airt_max","daily_precip","vpd_mean","sfc_pressure_mean","airt_min","wind_spd")

            # open first ecmwf file to extract needed information
            input_file_1=paste(path_to_met_source,varid[1],"_",startyear,"01.nc",sep="")
            data1=nc_open(input_file_1)

            # get timing variable
            steps_in_day=1

            # extract location variables
            lat = ncvar_get(data1, "Latitude") ; long = ncvar_get(data1, "Longitude")
            # expand the one directional values here into 2 directional
            lat_dim = length(lat) ; long_dim = length(long)
            long = array(long,dim=c(long_dim,lat_dim))
            lat = array(lat,dim=c(lat_dim,long_dim)) ; lat=t(lat)

            # If we are using the GSI model we need the 21 days (or month) before the start date of the simulation, so we need to check if we have this information
            # NOTE that this section of code is duplicated for each of the available datasets because of differences in storage and file name
            extra_year = FALSE ; present = 0
            for (lag in seq(1,length(varid))) {
                 input_file_1 = paste(path_to_met_source,varid[1],"_",as.character(as.numeric(startyear)-1),"01.nc",sep="")
                 # check whether the files exist or not
                 if (file.exists(input_file_1)) {present = present+1}
            }
            # if all the files exist then we will use them
            if (present == length(varid)) {extra_year=TRUE}

        } else if (met_source == "CHESS") {

            # declare variable ids needed to select files / infile variables
            varid=c("rsds","tas","precip","huss","psurf","dtr","sfcWind") ; infile_varid=c("rsds","tas","precip","huss","psurf","dtr","sfcWind")
            # open first ecmwf file to extract needed information
            input_file_1 = paste(path_to_met_source,"chess_",varid[1],"_",startyear,"01.nc",sep="")
            data1 = nc_open(input_file_1)

            # get timing variable
            steps_in_day = 1

            # extract location variables
            lat = ncvar_get(data1, "lat") ; long = ncvar_get(data1, "lon")

            # If we are using the GSI model we need the 21 days (or month) before the start date of the simulation, so we need to check if we have this information
            # NOTE that this section of code is duplicated for each of the available datasets because of differences in storage and file name
            extra_year = FALSE ; present = 0
            for (lag in seq(1,length(varid))) {
                 input_file_1 = paste(path_to_met_source,"chess_",varid[lag],"_",as.character(as.numeric(startyear)-1),"12.nc",sep="")
                 # check whether the files exist or not
                 if (file.exists(input_file_1)) {present=present+1}
            }
            # if all the files exist then we will use them
            if (present == length(varid)) {extra_year=TRUE}

        } # met_source

        # now find out which pixels we will filter through
        # if lat / long is 2 dimensional rather than vector we need to adjust for this
        lat_dim = dim(lat)[2] ; long_dim = dim(long)[1]
        # convert input data long to conform to what we need
        check1 = which(long > 180) ; if (length(check1) > 0) { long[check1] = long[check1]-360 }

        # which locations are within the desired zone
        remove_lat = intersect(which(lat < (max(latlon_in[,1])+1.0)),which(lat > (min(latlon_in[,1])-1.0)))
        remove_long = intersect(which(long < (max(latlon_in[,2])+1.0)),which(long > (min(latlon_in[,2])-1.0)))
        if (length(check1) > 0) { long[check1] = long[check1]+360 }
        # now find common where out in both contexts
        remove_lat = intersect(remove_lat,remove_long)
        # update both variables because of common matrix
        remove_long = remove_lat
        # adjust for matrix rather than vector arrangement
        remove_lat = remove_lat/dim(lat)[1]
        remove_long = (remove_long-(floor(remove_lat)*dim(lat)[1]))+1
        remove_lat = ceiling(remove_lat)

        # update new dimensions
        lat_dim = length(min(remove_lat):max(remove_lat)) ; long_dim = length(min(remove_long):max(remove_long))
        lat = lat[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat)] ; long = long[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat)]
        # next deal with extracting the wheat from the chaff
        # this should always be the swrad variable as we can easily put a lower limit
        tmp1 = ncvar_get(data1,infile_varid[1])
        # Also this should be applied to the first time step only as we want spatial pattern not temporal
        # Also restrict by the target area
        tmp1 = tmp1[,,1] ; tmp1 = tmp1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat)]
        # Include any pixels where the shortwave radiation is un-realistic.
        # We could expand this to other variables to ensure good realism but sw alone is generally sufficient.
        tmp1[tmp1 < 0] = NA
        # this section is key as near land sea borders it is possible for the nearest lat/long location to actually be a sea pixel
        wheat_from_chaff = which(is.na(tmp1) == FALSE)

        # convert input data long to conform to what we need
        check1 = which(long > 180) ; if (length(check1) > 0) { long[check1]=long[check1]-360 }
        # now filter through the reduced dataset for the specific locations
        # NOTE: this selection by met_in$wheat vectorises the lat and long variables therefore we need to use closest2 option 1 (I think...)
        output = lapply(1:dim(latlon_in)[1],FUN=closest2d,lat=lat[wheat_from_chaff],long=long[wheat_from_chaff],lat_in=latlon_in[,1],long_in=latlon_in[,2],nos_dim=1)
        var1_out = unlist(output, use.names=FALSE) ; rm(output)
        # select the correct location values for the original vector form the wheat_from_chaff
        wheat_from_chaff = wheat_from_chaff[var1_out]
        # return long to 0-360
        if (length(check1) > 0) { long[check1]=long[check1]+360 }

        # tidy
        nc_close(data1)

        # define timing variables
        years_to_load = as.numeric(startyear):as.numeric(endyear)
        if (extra_year) {
            load_years = c(as.character(as.numeric(startyear)-1),as.numeric(startyear):as.numeric(endyear))
        } else {
            load_years = as.numeric(startyear):as.numeric(endyear)
        }
        print("have finished determining met extraction parameters, now for the main load")
        # extract the data into lists which we will contruct.
        if (use_parallel) {

            # NOTE: that the use of mclapply() is due to reported improved efficiency over creating a virtual cluster.
            # However, mclapply does not (at the time of typing) work on Windows, i.e. Linux and Mac only
            cl <- min(length(load_years),numWorkers)
            var1_out_list=mclapply(load_years,FUN=load_met_function,varid=varid[1],infile_varid=infile_varid[1],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff,mc.cores = cl)
            print("...met load 14 %")
            var2_out_list=mclapply(load_years,FUN=load_met_function,varid=varid[2],infile_varid=infile_varid[2],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff,mc.cores = cl)
            print("...met load 29 %")
            var3_out_list=mclapply(load_years,FUN=load_met_function,varid=varid[3],infile_varid=infile_varid[3],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff,mc.cores = cl)
            print("...met load 43 %")
            var4_out_list=mclapply(load_years,FUN=load_met_function,varid=varid[4],infile_varid=infile_varid[4],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff,mc.cores = cl)
            print("...met load 57 %")
            var5_out_list=mclapply(load_years,FUN=load_met_function,varid=varid[5],infile_varid=infile_varid[5],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff,mc.cores = cl)
            print("...met load 71 %")
            var6_out_list=mclapply(load_years,FUN=load_met_function,varid=varid[6],infile_varid=infile_varid[6],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff,mc.cores = cl)
            print("...met load 86 %")
            wind_out_list=mclapply(load_years,FUN=load_met_function,varid=varid[7],infile_varid=infile_varid[7],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff,mc.cores = cl)
            print("...met load 100 %")

        } else {

            # or use serial
            var1_out_list=lapply(load_years,FUN=load_met_function,varid=varid[1],infile_varid=infile_varid[1],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
            print("...met load 14 %")
            var2_out_list=lapply(load_years,FUN=load_met_function,varid=varid[2],infile_varid=infile_varid[2],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
            print("...met load 29 %")
            var3_out_list=lapply(load_years,FUN=load_met_function,varid=varid[3],infile_varid=infile_varid[3],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
            print("...met load 43 %")
            var4_out_list=lapply(load_years,FUN=load_met_function,varid=varid[4],infile_varid=infile_varid[4],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
            print("...met load 57 %")
            var5_out_list=lapply(load_years,FUN=load_met_function,varid=varid[5],infile_varid=infile_varid[5],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
            print("...met load 71 %")
            var6_out_list=lapply(load_years,FUN=load_met_function,varid=varid[6],infile_varid=infile_varid[6],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
            print("...met load 86 %")
            wind_out_list=lapply(load_years,FUN=load_met_function,varid=varid[7],infile_varid=infile_varid[7],remove_lat=remove_lat,remove_long=remove_long,path_to_met_source=path_to_met_source,met_source=met_source,wheat=wheat_from_chaff)
            print("...met load 100 %")

        } # parallel option

        # user update
        print("...beginning restructuring of meteorological datasets")
        var1_out = 0 ; var2_out = 0 ; var3_out = 0 ; var4_out = 0 ; var5_out = 0 ; var6_out = 0 ; wind_out = 0 ; tmp_out = 0; t_grid = 0
        for (i in seq(1, length(var1_out_list))) {var1_out=append(var1_out,var1_out_list[[i]]$var_out) ; t_grid=append(t_grid,var1_out_list[[i]]$t_grid)}
        rm(var1_out_list)
        for (i in seq(1, length(var2_out_list))) {var2_out=append(var2_out,var2_out_list[[i]]$var_out)}
        rm(var2_out_list)
        for (i in seq(1, length(var3_out_list))) {var3_out=append(var3_out,var3_out_list[[i]]$var_out)}
        rm(var3_out_list)
        for (i in seq(1, length(var4_out_list))) {var4_out=append(var4_out,var4_out_list[[i]]$var_out)}
        rm(var4_out_list)
        for (i in seq(1, length(var5_out_list))) {var5_out=append(var5_out,var5_out_list[[i]]$var_out)}
        rm(var5_out_list)
        for (i in seq(1, length(var6_out_list))) {var6_out=append(var6_out,var6_out_list[[i]]$var_out)}
        rm(var6_out_list)
        for (i in seq(1, length(wind_out_list))) {wind_out=append(wind_out,wind_out_list[[i]]$var_out)}
        rm(wind_out_list)

        # remove initial value
        var1_out = var1_out[-1] ; var2_out = var2_out[-1]
        var3_out = var3_out[-1] ; var4_out = var4_out[-1]
        var5_out = var5_out[-1] ; var6_out = var6_out[-1]
        wind_out = wind_out[-1] ; gc()
        # we want total t_grid only
        t_grid = sum(t_grid)

        #
        # Unit conversions based on specific datasets
        #

        # convert ERA from W/m2 to MJ/m2/day
        if (met_source == "CHESS" | met_source == "ERA" | met_source == "CRUJRA") {
            var1_out = var1_out * 86400 * 1e-6
        }

        # Mean temperature and temperature range are provided in CHESS,
        # We want to make these into max / min temperature
        if (met_source == "CHESS") {
            tmp1 = var2_out + var6_out # max
            tmp2 = var2_out - var6_out # min
            var2_out = tmp1 ; var6_out = tmp2
        }

        # convert precipitation from kgH2O/m2/day -> per second
        if (met_source == "ERA") {
            var3_out = var3_out / 86400
        }

        # CHESS and CRUJRA provide specific humidity not vapour pressure deficit.
        # Best change that!
        if (met_source == "CHESS" | met_source == "CRUJRA") {
            var4_out = sp_humidity_to_vpd(var4_out,var5_out,(((var2_out + var6_out)*0.5)-273.15))
        }

        # generate some additional timing and data
        # NOTE that if we are using an extra of data for the purposes of a GSI model run then the number of days should be the number intended for simulation
        # not to match with the met files just at the moment.
        # This will be corrected later on when the location specific data are extracted
        print("...generating day of year variables")
        for (yr in seq(1,length(years_to_load))) {
             # if this includes the extra year add it on to the beginning
             if (extra_year) {extra_nos_days = nos_days_in_year(load_years[1])}
             # is current year a leap or not
             nos_days = nos_days_in_year(years_to_load[yr])

             # reconstruct arrays before doing site level
            if (yr == 1) {
                doy = seq(1,nos_days)
#                co2=rep(c(450,445,440,440,440,430,420,400,390,380,360,350,340,330,330,340,355,370,390,400,420,440,445,450),length.out=(nos_days*steps_in_day), each=floor((nos_days*steps_in_day)/24))
                # mass of dry air = 28.97(g/mol) ; mass of co2 = 44.01 (g/mol); *1e-6 scale from umol -> mol
                co2 = rep(380,length.out=(nos_days*steps_in_day)) # default ppm 570 / 370 face
                if (extra_year) {co2 = append(rep(380,length.out=(extra_nos_days*steps_in_day)),co2)}
            } else {
                doy = append(doy,seq(1,nos_days))
#                co2=append(co2,rep(c(450,445,440,440,440,430,420,400,390,380,360,350,340,330,330,340,355,370,390,400,420,440,445,450),length.out=(nos_days*steps_in_day), each=floor((nos_days*steps_in_day)/24))+(2*yr))
                co2 = append(co2,rep(380,length.out=(nos_days*steps_in_day)))
            }
        } # end of years loop

        # create day of run variable
        run_day = seq(1:length(doy))

        # user update
        print("...preparing final met output variables")
        # restructure to output variables
        var1_out = array(var1_out, dim=c(length(wheat_from_chaff),t_grid)) # SW Rad (MJ/m2/day)
        var2_out = array(var2_out, dim=c(length(wheat_from_chaff),t_grid)) # Max T (K)
        var3_out = array(var3_out, dim=c(length(wheat_from_chaff),t_grid)) # Precipitation (kgH2O/m2/s)
        var4_out = array(var4_out, dim=c(length(wheat_from_chaff),t_grid)) # VPD (Pa)
        var5_out = array(var5_out, dim=c(length(wheat_from_chaff),t_grid)) # Pressure (Pa)
        var6_out = array(var6_out, dim=c(length(wheat_from_chaff),t_grid)) # Min T (K)
        wind_out = array(wind_out, dim=c(length(wheat_from_chaff),t_grid)) # Wind (m/s)

        # output variables
        met_all = list(lat=lat,long=long,run_day=run_day,wheat=wheat_from_chaff,t_grid=t_grid
                      ,maxt=var2_out,swrad=var1_out,co2=co2,doy=doy,precip=var3_out,vpd=var4_out
                      ,pressure=var5_out,mint=var6_out,wind_spd=wind_out,extra_year=extra_year)

        # quick sanity check
        if (min(as.vector(met_all$swrad)) < -1) {stop(paste("SW_RAD summary: ",summary(as.vector(met_all$swrad)),sep="")) }
        met_all$swrad[which(met_all$swrad < 0)] = 0

        # clean up loose memory
        rm(lat_dim,long_dim,remove_lat,remove_long,years_to_load,load_years,tmp1,check1,nos_days,
           lat,long,run_day,wheat_from_chaff,t_grid,var2_out,var1_out,
           co2,doy,var3_out,var4_out,var5_out,extra_year,var6_out,wind_out)

    } # site_specific or not

    # final tidy up and return
    gc(reset=TRUE,verbose=FALSE)
    return(met_all)

} # function end
## Use byte compile
load_met_fields_for_extraction<-cmpfun(load_met_fields_for_extraction)
