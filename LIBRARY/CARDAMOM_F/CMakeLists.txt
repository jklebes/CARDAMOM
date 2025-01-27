set(FILES misc/math_functions.f90
                  misc/oksofar.f90
                  misc/brent_zero/brent_zero.f90
                  general/cardamom_structures.f90
                  method/MHMCMC/MCMC_FUN/MHMCMC_STRUCTURES.f90
                  #model/DALEC_1005a/src/DALEC_1005a_PARS.f90
                  model/DALEC.A1.C1.D2.F2.H2.P1/src/DALEC.A1.C1.D2.F2.H2.P1_PARS.f90
                  method/MHMCMC/MCMC_FUN/MHMCMC_StressTests.f90
                  general/cardamom_io.f90
                  general/model_likelihood_wrapper.f90
                  #general/model_likelihood_wrapper_for_R.f90
                  method/MHMCMC/MCMC_FUN/MHMCMC.f90
                  model/DALEC.A1.C1.D2.F2.H2.P1/likelihood/MODEL_LIKELIHOOD.f90
                  #model/DALEC_1005a/src/DALEC_1005a.f90
                  model/DALEC.A1.C1.D2.F2.H2.P1/src/DALEC.A1.C1.D2.F2.H2.P1.f90
                )
#try shared for R
add_library(cardamom-lib SHARED ${FILES})
#add_executable(cardamom general/cardamom_main.f90)
#target_link_libraries(cardamom cardamom-lib)
