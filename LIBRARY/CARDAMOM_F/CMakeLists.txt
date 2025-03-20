set(FILES1 misc/math_functions.f90
                  misc/oksofar.f90
                  misc/brent_zero/brent_zero.f90
                )

set(FILES2  general/cardamom_structures.f90
                  #model/DALEC_1005a/src/DALEC_1005a_PARS.f90
                  model/DALEC.A1.C1.D2.F2.H2.P1/src/DALEC.A1.C1.D2.F2.H2.P1_PARS.f90
                  general/cardamom_io.f90
                  general/model_likelihood_wrapper.f90
                  #general/model_likelihood_wrapper_for_R.f90
                  model/DALEC.A1.C1.D2.F2.H2.P1/likelihood/MODEL_LIKELIHOOD.f90
                  #model/DALEC_1005a/src/DALEC_1005a.f90
                  model/DALEC.A1.C1.D2.F2.H2.P1/src/DALEC.A1.C1.D2.F2.H2.P1.f90
                )
add_executable(cardamom general/cardamom_main.f90)
add_library(cardamom-lib SHARED ${FILES1} ${FILES2})
add_dependencies(cardamom-lib cardamom-samplers-lib) 
target_link_libraries(cardamom PUBLIC cardamom-samplers-lib cardamom-lib)
#target_link_libraries(cardamom-lib PUBLIC cardamom-samplers-lib)
target_include_directories(cardamom-lib PUBLIC "method")

add_subdirectory("method")

target_sources(cardamom-lib PUBLIC ${FILES2})

