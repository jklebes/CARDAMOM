
###
## Load met function
###

# This function assumed ERA-Interim data files are being used

load_met_function<- function (year_to_do,varid,infile_varid,remove_lat,remove_long,path_to_met_source,met_source,wheat) {

    if (met_source == "CHESS") {

        # open first file in the sequence
        m = 1
        input_file_1 = paste(path_to_met_source,"chess_",varid[1],"_",year_to_do,"0",m,".nc",sep="")
        # open netcdf files
        data1 = nc_open(input_file_1)

        # read the met drivers
        var1 = ncvar_get(data1, infile_varid[1]) ; var1=var1[,,1:(dim(var1)[3])]

        # close files after use
        nc_close(data1)

        # keep count of time steps
        t_grid = dim(var1)[3]

        # filter spatial extent
        var1 = var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),]
        # assign correct error value
        var1[is.na(var1)] = -9999

        # move through time removing the "chaff"
        var1_out = as.vector(var1[,,1])[wheat]
        for (i in seq(2, t_grid)) {var1_out = append(var1_out,as.vector(var1[,,i])[wheat])}

        # loop through months in the year
        for (m in seq(2, 12)) {
             # select correct file names
             if (m > 9) {
                 input_file_1 = paste(path_to_met_source,"chess_",varid[1],"_",year_to_do,m,".nc",sep="")
             } else {
                 input_file_1 = paste(path_to_met_source,"chess_",varid[1],"_",year_to_do,"0",m,".nc",sep="")
             }
             # open netcdf files
             data1 = nc_open(input_file_1)
             # read the met drivers
             var1 = ncvar_get(data1, infile_varid) ; var1=var1[,,1:(dim(var1)[3])]
             # In February 2013 some of the CHESS drivers incorectly have 29 days present, remove the final one!
             if (as.numeric(year_to_do) == 2013 & dim(var1)[3] == 29) {var1=var1[,,1:(dim(var1)[3]-1)]}
             tmp_t = dim(var1)[3] ; t_grid = t_grid + tmp_t
             # close files after use
             nc_close(data1)

             # filter spatial extent
             var1 = var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),]
             # assign correct error value
             var1[is.na(var1)] = -9999

             # move through time removing the "chaff"
             tmp = as.vector(var1[,,1])[wheat]
             for (i in seq(2, tmp_t)) {tmp = append(tmp,as.vector(var1[,,i])[wheat])}
             # now append to the output variable
             var1_out = append(var1_out,tmp)

        } # looping through dataset

        # clean up
        rm(var1,tmp,i,m,tmp_t) ; gc(reset=TRUE,verbose=FALSE)

    } else if (met_source == "ERA") {

        # open first file in the sequence
        m = 1
        input_file_1 = paste(path_to_met_source,varid[1],"_",year_to_do,"0",m,".nc",sep="")
        # open netcdf files
        data1 = nc_open(input_file_1)

        # read the met drivers
        var1 = ncvar_get(data1, infile_varid[1]) ; var1 = var1[,,1:(dim(var1)[3])]

        # close files after use
        nc_close(data1)

        # keep count of time steps
        t_grid = dim(var1)[3]

        # filter spatial extent
        var1 = var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),]
        # assign correct error value
        var1[is.na(var1)] = -9999

        # move through time removing the "chaff"
        tmp = as.vector(var1[,,1])[wheat]
        for (i in seq(2, t_grid)) {
             tmp=append(tmp,as.vector(var1[,,i])[wheat])
        }
        # now append to the output variable
        var1_out = tmp

        # loop through months in the year
        for (m in seq(2, 12)) {
             # select correct file names
             if (m > 9) {
                 input_file_1=paste(path_to_met_source,varid[1],"_",year_to_do,m,".nc",sep="")
             } else {
                 input_file_1=paste(path_to_met_source,varid[1],"_",year_to_do,"0",m,".nc",sep="")
             }
             # open netcdf files
             data1=nc_open(input_file_1)
             # read the met drivers
             var1=ncvar_get(data1, infile_varid) ; var1=var1[,,1:(dim(var1)[3])]
             tmp_t = dim(var1)[3] ; t_grid = t_grid + tmp_t
             # close files after use
             nc_close(data1)

             # update to the correct arrays
             # filter spatial extent
             var1=var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),]
             # assign correct error value
             var1[is.na(var1)] = -9999

             # move through time removing the "chaff"
             tmp=as.vector(var1[,,1])[wheat]
             for (i in seq(2, tmp_t)) {
                  tmp=append(tmp,as.vector(var1[,,i])[wheat])
             }
             # now append to the output variable
             var1_out=append(var1_out,tmp)

        } # looping through dataset

        # clean up
        rm(var1,tmp,i,m,tmp_t) ; gc(reset=TRUE,verbose=FALSE)

    } else if (met_source == "CRUJRA") {

        # open first file in the sequence
        input_file_1 = paste(path_to_met_source,"/",varid[1],"/crujra.V1.1.5d.",varid[1],".",year_to_do,".365d.noc.nc",sep="")
        # open netcdf files
        data1 = nc_open(input_file_1)
        # read the met drivers
        var1 = ncvar_get(data1, infile_varid[1]) ; var1 = var1[,,1:(dim(var1)[3])]
        # close files after use
        nc_close(data1)

        # keep count of time steps
        t_grid = dim(var1)[3]

        # filter spatial extent
        var1 = var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),]
        # assign correct error value
        var1[is.na(var1)] = -9999

        # Check whether this is a leap year or not
        nos_days = nos_days_in_year(year_to_do)

        # Based on the number of days aggregate to a daily time step.
        # Unlike in other data sources we will need to do something different
        # depending on the data type (i.e. temperature vs rainfall).
        aggregate_by = t_grid / 365 #  divide by 365 as this dataset neglects leap years, which are dealt with below

        if (varid[1] == "dswrf") {
            func <- function(var_in, aggregate_by) {
                       answer = rollapply(var_in, FUN = sum, by = aggregate_by, width = aggregate_by, align = "left")
                       answer = answer * 1.157407e-05 # (1/86400)
                       return(answer)
                    }
        } else if (varid[1] == "tmax") {
            # daily maximum temperature K
            func <- function(var_in, aggregate_by) {
                       answer = rollapply(var_in, FUN = max, by = aggregate_by, width = aggregate_by, align = "left")
                       return(answer)
                    }
        } else if (varid[1] == "pre") {
            # mm/6h -> kgH2O/m2/s
            func <- function(var_in, aggregate_by) {
                       answer = rollapply(var_in, FUN = sum, by = aggregate_by, width = aggregate_by, align = "left")
                       answer = answer * 1.157407e-05 # (1/86400)
                       return(answer)
                    }
        } else if (varid[1] == "spfh") {
            # calculate daily mean specific humidity (kg/kg)
            func <- function(var_in, aggregate_by) {
                       answer = rollapply(var_in, FUN = mean, by = aggregate_by, width = aggregate_by, align = "left")
                       return(answer)
                    }
        } else if (varid[1] == "pres") {
            # calculate daily mean pressure (Pa)
            func <- function(var_in, aggregate_by) {
                       answer = rollapply(var_in, FUN = mean, by = aggregate_by, width = aggregate_by, align = "left")
                       return(answer)
                    }
        } else if (varid[1] == "tmin") {
            # daily maximum temperature K
            func <- function(var_in, aggregate_by) {
                       answer = rollapply(var_in, FUN = min, by = aggregate_by, width = aggregate_by, align = "left")
                       return(answer)
                    }
        } else if (varid[1] == "wsp") {
            # daily mean wind speed m/s
            func <- function(var_in, aggregate_by) {
                       answer = rollapply(var_in, FUN = mean, by = aggregate_by, width = aggregate_by, align = "left")
                       return(answer)
                    }
        }
        # Use the correct function, in the averaging...
        var1 = apply(var1, c(1,2), func, aggregate_by = aggregate_by)
        t_grid = dim(var1)[1]

        # ...unfortunately this outputs an array which the wrong order for the dimension
        # which we need to fix now...
        # ...fortunately we can do so while we move through time removing the "chaff"
        tmp = as.vector(var1[i,,])[wheat]
        for (i in seq(2, t_grid)) {
             tmp = append(tmp,as.vector(var1[i,,])[wheat])
        }
        # if 366 days in the year we need to add a final day as the dataset only contains 365 days...sad-face
        if (nos_days == 366) {tmp = append(tmp,as.vector(var1[t_grid,,])[wheat]) ; t_grid = t_grid + 1}
        # now pass to the output variable
        var1_out = tmp

        # clean up
        rm(var1,tmp,i,m,tmp_t) ; gc(reset=TRUE,verbose=FALSE)

    } # end data source selection

    # return back to the user
    return(list(var_out=var1_out,t_grid=t_grid))

} # end function
## Use byte compile
load_met_function<-cmpfun(load_met_function)
