
###
## Extracts location specific information on the leaf carbon per leaf area (gCm2) from a loaded gridded dataset
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_lca_prior<- function(spatial_type,resolution,grid_type,latlon_in,lca_all) {

   # Update the user
   print(paste("LCA prior extracted for current location ",Sys.time(),sep=""))

   # find the nearest location
   output = closest2d(1,lca_all$lat,lca_all$long,latlon_in[1],latlon_in[2],2)
   i1 = unlist(output)[1] ; j1 = unlist(output)[2]

   # If resolution has been provides as single value then adjust this here
   if (length(resolution) == 1 & spatial_type == "grid") {tmp_res = resolution * c(1,1)} else {tmp_res = resolution}

   # work out number of pixels to average over
   if (spatial_type == "grid") {
       # resolution of the product
       product_res = c(abs(lca_all$long[2,1]-lca_all$long[1,1]),abs(lca_all$lat[1,2]-lca_all$lat[1,1]))
       if (grid_type == "wgs84") {
           # radius is ceiling of the ratio of the product vs analysis ratio
           radius = floor(0.5*(resolution / product_res))
        } else if (grid_type == "UK") {
           # Estimate radius for UK grid assuming radius is determine by the longitude size
           # 6371e3 = mean earth radius (m)
           radius = round(rad2deg(sqrt((resolution / 6371e3**2))) / product_res, digits=0)
          #radius = max(0,floor(1*resolution*1e-3*0.5))
        } else {
           stop("have not specified the grid used in this analysis")
        }
   } else {
        radius = c(0,0)
   }

   # Work out average areas
   average_i = (i1-radius[1]):(i1+radius[1]) ; average_j = (j1-radius[2]):(j1+radius[2])
   average_i = max(1,(i1-radius[1])):min(dim(lca_all$lca_gCm2)[1],(i1+radius[1]))
   average_j = max(1,(j1-radius[2])):min(dim(lca_all$lca_gCm2)[2],(j1+radius[2]))
   # Carry out averaging
   lca_gCm2 = mean(lca_all$lca_gCm2[average_i,average_j], na.rm=TRUE)
   lca_unc_gCm2 = mean(lca_all$lca_uncertainty_gCm2[average_i,average_j], na.rm=TRUE)

   # Convert any NaN to missing data flag -9999
   lca_gCm2[which(is.na(lca_gCm2))] = -9999 ; lca_unc_gCm2[which(is.na(lca_unc_gCm2))] = -9999

   # pass the information back
   return(list(lca_gCm2 = lca_gCm2, lca_unc_gCm2 = lca_unc_gCm2))

} # end function extract_lca_prior

## Use byte compile
extract_lca_prior<-cmpfun(extract_lca_prior)
