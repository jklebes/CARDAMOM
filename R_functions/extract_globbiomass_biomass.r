

extract_globbiomass_biomass<- function(timestep_days,spatial_type,resolution,grid_type,latlon_in,Cwood_all) {

            # find the nearest location
            output = closest2d(1,Cwood_all$lat_2010,Cwood_all$long_2010,latlon_in[1],latlon_in[2],2)
            i1_2010 = unlist(output)[1] ; j1_2010 = unlist(output)[2]
            output = closest2d(1,Cwood_all$lat_2017,Cwood_all$long_2017,latlon_in[1],latlon_in[2],2)
            i1_2017 = unlist(output)[1] ; j1_2017 = unlist(output)[2]

            print(paste("GlobBIOMASS data extracted for current location ",Sys.time(),sep=""))

            # work out number of pixels to average over
            if (spatial_type == "grid") {
                if (grid_type == "wgs84") {
                    # resolution of the product
                    product_res = abs(Cwood_all$lat_2010[1,2]-Cwood_all$lat_2010[1,1])+abs(Cwood_all$long_2010[2,1]-Cwood_all$long_2010[1,1])
                    product_res = product_res * 0.5
                    # radius is ceiling of the ratio of the product vs analysis ratio
                    radius = ceiling(resolution / product_res)
                } else if (grid_type == "UK") {
                    radius = max(0,floor(1*resolution*1e-3*0.5))
                } else {
                    stop("have not specified the grid used in this analysis")
                }
            } else {
                radius = 0
            }

            # Work out average areas for 2010
            average_i = max(1,(i1_2010-radius)):min(dim(Cwood_all$biomass_2010)[1],(i1_2010+radius)) 
            average_j=max(1,(j1_2010-radius)):min(dim(Cwood_all$biomass_2010)[2],(j1_2010+radius))
            # Carry out averaging
            # 2010 first...
            Cwood = mean(Cwood_all$biomass_2010[average_i,average_j], na.rm=TRUE)
            Cwood_unc = mean(Cwood_all$biomass_unc_2010[average_i,average_j], na.rm=TRUE)
      
            # Work out average areas for 2017
            average_i = max(1,(i1_2017-radius)):min(dim(Cwood_all$biomass_2010)[1],(i1_2017+radius)) 
            average_j=max(1,(j1_2017-radius)):min(dim(Cwood_all$biomass_2010)[2],(j1_2017+radius))
            # ... append 2017 next
            Cwood = append(Cwood,mean(Cwood_all$biomass_2017[average_i,average_j], na.rm=TRUE))
            Cwood_unc = append(Cwood_unc,mean(Cwood_all$biomass_unc_2017[average_i,average_j], na.rm=TRUE))

            # convert missing data back to -9999
            Cwood[which(is.na(Cwood))] = -9999.0 ; Cwood_unc[which(is.na(Cwood_unc))] = -9999.0

            # pass the information back
            return(list(Cwood_stock = Cwood, Cwood_stock_unc = Cwood_unc))


} # end of function
