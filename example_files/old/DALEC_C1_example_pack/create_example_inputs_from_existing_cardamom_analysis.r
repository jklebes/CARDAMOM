load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_LU_FIRES_MHMCMC/CSSP_BR-Sa1_NEElocal/RESULTS_PROCESSED/BR-Sa1.RData")
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_LU_FIRES_MHMCMC/CSSP_BR-Sa1_NEElocal/infofile.RData")

# Write drivers to file
drivers_out = data.frame(run_day = drivers$met[,1], mint_C = drivers$met[,2], maxt_C = drivers$met[,3],
                         swrad_MJm2day = drivers$met[,4], co2_ppm = drivers$met[,5], doy = drivers$met[,6],
                         precip_kgm2s = drivers$met[,7], deforestation_fraction = drivers$met[,8],
                         burnt_area_fraction = drivers$met[,9], mint_21day = drivers$met[,10],
                         dayl_21day = drivers$met[,11], vpd_21day = drivers$met[,12],
                         management_type = drivers$met[,13], avgTemp_C = drivers$met[,14],
                         wind_spd_ms = drivers$met[,15], vpd_Pa = drivers$met[,16])
write.table(drivers_out, file = "./inputs/drivers.csv", sep=",", row.names=FALSE, append=FALSE)

# Write LAI observations to file
observations_out = data.frame(run_day = drivers$met[,1], LAI_m2m2 = drivers$obs[,3], LAI_unc_m2m2 = drivers$obs[,4],
                                                         NEE_gCm2day = drivers$obs[,5], NEE_unc_gCm2day = drivers$obs[,6])

write.table(observations_out, file = "./inputs/observations.csv", sep=",", row.names=FALSE, append=FALSE)
# Write parameters to file
parameters_out = data.frame(lit_decomp = as.vector(parameters[1,,]), ragpp = as.vector(parameters[2,,]),
                            fgpp = as.vector(parameters[3,,]), rgpp = as.vector(parameters[4,,]),
                            LL_yrs = as.vector(parameters[5,,]), wood_turn = as.vector(parameters[6,,]),
                            root_turn = as.vector(parameters[7,,]), lit_Rhet = as.vector(parameters[8,,]),
                            som_Rhet = as.vector(parameters[9,,]), Tair_coef = as.vector(parameters[10,,]),
                            Ceff = as.vector(parameters[11,,]), MaxDayBudBurst = as.vector(parameters[12,,]),
                            lgpp = as.vector(parameters[13,,]), BudBurstDays = as.vector(parameters[14,,]),
                            MaxDayLeafFall = as.vector(parameters[15,,]), LeafFallDays = as.vector(parameters[16,,]),
                            LCA_gCm2 = as.vector(parameters[17,,]), Init_lab = as.vector(parameters[18,,]),
                            Init_fol = as.vector(parameters[19,,]), Init_root = as.vector(parameters[20,,]),
                            Init_wood = as.vector(parameters[21,,]), Init_lit = as.vector(parameters[22,,]),
                            Init_som = as.vector(parameters[23,,]))
write.table(parameters_out, file = "./inputs/parameters.csv", sep=",", row.names=FALSE, append=FALSE)
