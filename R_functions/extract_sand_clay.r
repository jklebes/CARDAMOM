
###
## Function to extract location specific information on soil texture from the gridded SoilGrids database
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_sand_clay<- function(spatial_type,resolution,grid_type,latlon_in,sand_clay_all) {

  # Update the user
  print(paste("Sand/clay data extracted for current location ",Sys.time(),sep=""))

  # find the nearest location
  output=closest2d_2(1,sand_clay_all$lat,sand_clay_all$long,latlon_in[1],latlon_in[2])
  j1=unlist(output, use.names=FALSE)[2];i1=unlist(output, use.names=FALSE)[1]

  # Extract the correct value
  top_sand = sand_clay_all$top_sand[i1,j1]
  bot_sand = sand_clay_all$bot_sand[i1,j1]
  top_clay = sand_clay_all$top_clay[i1,j1]
  bot_clay = sand_clay_all$bot_clay[i1,j1]

  # Guard against NaN values
  if (is.na(top_sand) | is.infinite(top_sand)) {top_sand = 40}
  if (is.na(bot_sand) | is.infinite(bot_sand)) {bot_sand = 40}
  if (is.na(top_clay) | is.infinite(top_clay)) {top_clay = 15}
  if (is.na(bot_clay) | is.infinite(bot_clay)) {bot_clay = 15}

  # just to check because when averaging sometimes the sand / clay combinations can be > 100 %
  # 94 % chosesn as this is the highest total % found in the HWSD dataset
  if ((top_sand+top_clay) > 94) {
      tmp1 = top_sand / (top_sand + top_clay + 6) # 6 % is implicit in the 94 % max value for silt / gravel
      tmp2 = top_clay / (top_sand + top_clay + 6) # 6 % is implicit in the 94 % max value for silt / gravel
      top_sand = tmp1*100 ; top_clay = tmp2*100
  }
  if ((bot_sand+bot_clay) > 94) {
       tmp1 = bot_sand / (bot_sand + bot_clay + 6) # 6 % is implicit in the 94 % max value for silt / gravel
       tmp2 = bot_clay / (bot_sand + bot_clay + 6) # 6 % is implicit in the 94 % max value for silt / gravel
       top_sand = tmp1*100 ; top_clay = tmp2*100
  }

  # pass the information back
  return(list(top_sand=top_sand,bot_sand=bot_sand,top_clay=top_clay,bot_clay=bot_clay))

} # end function extract_sand_clay

## Use byte compile
extract_sand_clay<-cmpfun(extract_sand_clay)
