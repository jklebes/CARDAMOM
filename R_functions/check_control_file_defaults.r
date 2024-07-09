
###
## Function to check the for the presence of options from control files,
## where appropriate set default values and provide warnings for obvious errors
###

# Author: T. Luke Smallman (02/05/2024)

check_control_file_defaults<-function() {

  ## Set defaults incase missing, NOTE: <<- to assign global
  # Analysis options
  if (exists("select_country") == FALSE)                {select_country <<- FALSE}
  if (exists("met_interp") == FALSE)                    {met_interp <<- FALSE}
  if (exists("pft_specific_parameters") == FALSE)       {pft_specific_parameters <<- FALSE}
  # Site combined forcings and assimilated data path
  if (exists("path_to_site_obs") == FALSE)              {path_to_site_obs <<- " "}
  # Forcings data paths 
  if (exists("path_to_landsea") == FALSE)               {path_to_landsea <<- "default"}
  if (exists("path_to_met_source") == FALSE)            {path_to_met_source <<- " "}
  if (exists("path_to_sand_clay") == FALSE)             {path_to_sand_clay <<- " "}
  if (exists("path_to_forestry") == FALSE)              {path_to_forestry <<- " "}
  if (exists("path_to_burnt_area") == FALSE)            {path_to_burnt_area <<- " "}
  if (exists("path_to_co2") == FALSE)                   {path_to_co2 <<- "./R_functions/"}
  # Assimilated data paths 
  if (exists("path_to_lai") == FALSE)                   {path_to_lai <<- " "}
  if (exists("path_to_fapar") == FALSE)                 {path_to_fapar <<- " "}
  if (exists("path_to_crop_management") == FALSE)       {path_to_crop_management <<- " "}
  if (exists("path_to_Csom") == FALSE)                  {path_to_Csom <<- " "}
  if (exists("path_to_Cwood_inc") == FALSE)             {path_to_Cwood_inc <<- " "}
  if (exists("path_to_Cwood_mortality") == FALSE)       {path_to_Cwood_mortality <<- "v"}
  if (exists("path_to_Cwood") == FALSE)                 {path_to_Cwood <<- " "}
  if (exists("path_to_Cwood_initial") == FALSE)         {path_to_Cwood_initial <<- " "}
  if (exists("path_to_Cwood_potential") == FALSE)       {path_to_Cwood_potential <<- " "}
  if (exists("path_to_gleam") == FALSE)                 {path_to_gleam <<- " "}
  if (exists("path_to_nbe") == FALSE)                   {path_to_nbe <<- " "}
  if (exists("path_to_gpp") == FALSE)                   {path_to_gpp <<- " "}
  if (exists("path_to_fire") == FALSE)                  {path_to_fire <<- " "}
  if (exists("path_to_lca") == FALSE)                   {path_to_lca <<- " "}
  # Forcings data options
  if (exists("met_source") == FALSE)                    {met_source <<- " "}
  if (exists("burnt_area_source") == FALSE)             {burnt_area_source <<- " "}
  if (exists("deforestation_source") == FALSE)          {deforestation_source <<- " "}
  if (exists("sand_clay_source") == FALSE)              {sand_clay_source <<- " "}
  # Assimilated data options 
  if (exists("lai_source") == FALSE)                    {lai_source <<- " "}
  if (exists("fapar_source") == FALSE)                  {fapar_source <<- " "}
  if (exists("Csom_source") == FALSE)                   {Csom_source <<- " "}
  if (exists("soilwater_initial_source") == FALSE)      {soilwater_initial_source <<- " "}
  if (exists("Evap_source") == FALSE)                   {Evap_source <<- " "}
  if (exists("Cwood_inc_source") == FALSE)              {Cwood_inc_source <<- " "}
  if (exists("Cwood_mortality_source") == FALSE)        {Cwood_mortality_source <<- " "}
  if (exists("GPP_source") == FALSE)                    {GPP_source <<- " "}
  if (exists("fire_source") == FALSE)                   {fire_source <<- " "}
  if (exists("Reco_source") == FALSE)                   {Reco_source <<- " "}
  if (exists("NEE_source") == FALSE)                    {NEE_source <<- " "}
  if (exists("nbe_source") == FALSE)                    {nbe_source <<- " "}
  if (exists("harvest_source") == FALSE)                {harvest_source <<- " "}
  if (exists("foliage_to_litter_source") == FALSE)      {foliage_to_litter_source <<- " "}
  if (exists("Cfol_initial_source") == FALSE)           {Cfol_initial_source <<- " "}
  if (exists("Cwood_initial_source") == FALSE)          {Cwood_initial_source <<- " "}
  if (exists("Croots_initial_source") == FALSE)         {Croots_initial_source <<- " "}
  if (exists("Clit_initial_source") == FALSE)           {Clit_initial_source <<- " "}
  if (exists("Cfol_stock_source") == FALSE)             {Cfol_stock_source <<- " "}
  if (exists("Cfolmax_stock_source") == FALSE)          {Cfolmax_stock_source <<- " "}
  if (exists("Cwood_stock_source") == FALSE)            {Cwood_stock_source <<- " "}
  if (exists("Cstem_stock_source") == FALSE)            {Cstem_stock_source <<- " "}
  if (exists("Cbranch_stock_source") == FALSE)          {Cbranch_stock_source <<- " "}
  if (exists("Cagb_stock_source") == FALSE)             {Cagb_stock_source <<- " "}
  if (exists("Ccoarseroot_stock_source") == FALSE)      {Ccoarseroot_stock_source <<- " "}
  if (exists("Croots_stock_source") == FALSE)           {Croots_stock_source <<- " "}
  if (exists("Clit_stock_source") == FALSE)             {Clit_stock_source <<- " "}
  if (exists("lca_source") == FALSE)                    {lca_source <<- " "}
  if (exists("frac_Cwood_coarse_root_source") == FALSE) {frac_Cwood_coarse_root_source <<- " "}
  if (exists("minLWP_source") == FALSE)                 {minLWP_source <<- " "}
  if (exists("Cwood_potential_source") == FALSE)        {Cwood_potential_source <<- " "}
  if (exists("crop_management_source") == FALSE)        {crop_management_source <<- " "}
  if (exists("snow_source") == FALSE)                   {snow_source <<- " "}
  # Algorithm options
  if (exists("request_nos_chains") == FALSE)            {request_nos_chains <<- 3}
  if (exists("request_nos_samples") == FALSE)           {request_nos_samples <<- 10e6}
  if (exists("request_nos_subsamples") == FALSE)        {request_nos_subsamples <<- 1e3}
  if (exists("request_use_EDCs") == FALSE)              {request_use_EDCs <<- TRUE}
  if (exists("request_extended_mcmc") == FALSE)         {request_extended_mcmc <<- 10e6}
  if (exists("request_cost_function_scaling") == FALSE) {request_cost_function_scaling <<- 0}
  # Computer options
  if (exists("request_use_server") == FALSE)            {request_use_server <<- FALSE}
  if (exists("request_use_local_slurm") == FALSE)       {request_use_local_slurm <<- FALSE} 
  if (exists("request_runtime") == FALSE)               {request_runtime <<- 48}
  if (exists("request_compile_server") == FALSE)        {request_compile_server <<- FALSE}
  if (exists("request_compile_local") == FALSE)         {request_compile_local <<- TRUE}

  ## Check for obvious combination errors
  # Forcings datasets
  if (met_source != "site_specific" & path_to_met_source == " ")                                       {stop(paste("specified 'met_source' and 'path_to_met_source' incompatible"))}
  if (burnt_area_source != "site_specific" & burnt_area_source != " " & path_to_burnt_area == " ")     {stop(paste("specified 'burnt_area_source' and 'path_to_burnt_area' incompatible"))}
  if (deforestation_source != "site_specific" & deforestation_source != " " & path_to_forestry == " ") {stop(paste("specified 'deforestation_source' and 'path_to_forestry' incompatible"))}
  if (sand_clay_source != "site_specific" & sand_clay_source != " " & path_to_sand_clay == " ")        {stop(paste("specified 'sand_clay_source' and 'path_to_sand_clay' incompatible"))}
  # Assimilated datasets, note these only consider those attached to gridded datasets
  if (lai_source != "site_specific" & lai_source != " " & path_to_lai == " ")                                    {stop(paste("specified 'lai_source' and 'path_to_lai' incompatible"))}
  if (fapar_source != "site_specific" & fapar_source != " " & path_to_fapar == " ")                              {stop(paste("specified 'fapar_source' and 'path_to_fapar' incompatible"))}
  if (crop_management_source != "site_specific" & crop_management_source != " " & path_to_crop_management == " "){stop(paste("specified 'crop_management_source' and 'path_to_crop_management' incompatible"))}
  if (Csom_source != "site_specific" & Csom_source != " " & path_to_Csom == " ")                                 {stop(paste("specified 'Csom_source' and 'path_to_Csom' incompatible"))}
  if (Cwood_inc_source != "site_specific" & Cwood_inc_source != " " & path_to_Cwood_inc == " ")                  {stop(paste("specified 'Cwood_inc_source' and 'path_to_Cwood_inc' incompatible"))}
  if (Cwood_mortality_source != "site_specific" & Cwood_mortality_source != " " & path_to_Cwood_mortality == " "){stop(paste("specified 'Cwood_mortality_source' and 'path_to_Cwood_mortality' incompatible"))}
  if (Cwood_stock_source != "site_specific" & Cwood_stock_source != " " & path_to_Cwood == " ")                  {stop(paste("specified 'Cwood_stock_source' and 'path_to_Cwood' incompatible"))}
  if (Cwood_initial_source != "site_specific" & Cwood_initial_source != " " & path_to_Cwood_initial == " ")      {stop(paste("specified 'Cwood_initial_source' and 'path_to_Cwood_initial' incompatible"))}
  if (Cwood_potential_source != "site_specific" & Cwood_potential_source != " " & path_to_Cwood_potential == " "){stop(paste("specified 'Cwood_potential_source' and 'path_to_Cwood_potential' incompatible"))}
  if (soilwater_initial_source != "site_specific" & soilwater_initial_source != " " & path_to_gleam == " ")      {stop(paste("specified 'soilwater_initial_source' and 'path_to_gleam' incompatible"))}
  if (nbe_source != "site_specific" & nbe_source != " " & path_to_nbe == " ")                                    {stop(paste("specified 'nbe_source' and 'path_to_nbe' incompatible"))}
  if (GPP_source != "site_specific" & GPP_source != " " & path_to_gpp == " ")                                    {stop(paste("specified 'GPP_source' and 'path_to_gpp' incompatible"))}
  if (fire_source != "site_specific" & fire_source != " " & path_to_fire == " ")                                 {stop(paste("specified 'fire_source' and 'path_to_fire' incompatible"))}
  if (lca_source != "site_specific" & lca_source != " " & path_to_lca == " ")                                    {stop(paste("specified 'fire_source' and 'path_to_lca' incompatible"))}

  # NOTE: current assimilated variables which have currently only been done on site scale, without gridded datasets
#  if (exists("Evap_source") == FALSE)                   {Evap_source <<- " "}
#  if (exists("fire_source") == FALSE)                   {fire_source <<- " "}
#  if (exists("Reco_source") == FALSE)                   {Reco_source <<- " "}
#  if (exists("NEE_source") == FALSE)                    {NEE_source <<- " "}
#  if (exists("harvest_source") == FALSE)                {harvest_source <<- " "}
#  if (exists("foliage_to_litter_source") == FALSE)      {foliage_to_litter_source <<- " "}
#  if (exists("Cfol_initial_source") == FALSE)           {Cfol_initial_source <<- " "}
#  if (exists("Croots_initial_source") == FALSE)         {Croots_initial_source <<- " "}
#  if (exists("Clit_initial_source") == FALSE)           {Clit_initial_source <<- " "}
#  if (exists("Cfol_stock_source") == FALSE)             {Cfol_stock_source <<- " "}
#  if (exists("Cfolmax_stock_source") == FALSE)          {Cfolmax_stock_source <<- " "}
#  if (exists("Cstem_stock_source") == FALSE)            {Cstem_stock_source <<- " "}
#  if (exists("Cbranch_stock_source") == FALSE)          {Cbranch_stock_source <<- " "}
#  if (exists("Cagb_stock_source") == FALSE)             {Cagb_stock_source <<- " "}
#  if (exists("Ccoarseroot_stock_source") == FALSE)      {Ccoarseroot_stock_source <<- " "}
#  if (exists("Croots_stock_source") == FALSE)           {Croots_stock_source <<- " "}
#  if (exists("Clit_stock_source") == FALSE)             {Clit_stock_source <<- " "}
#  if (exists("frac_Cwood_coarse_root_source") == FALSE) {frac_Cwood_coarse_root_source <<- " "}
#  if (exists("minLWP_source") == FALSE)                 {minLWP_source <<- " "}
#  if (exists("crop_management_source") == FALSE)        {crop_management_source <<- " "}
#  if (exists("snow_source") == FALSE)                   {snow_source <<- " "}

  return(print("Standard checks on control file completed"))

} # end function check_control_file_defaults

## Use byte compile
check_control_file_defaults<-cmpfun(check_control_file_defaults)


