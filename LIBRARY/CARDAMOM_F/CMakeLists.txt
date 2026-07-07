set(FILES1 misc/math_functions.f90
                  misc/brent_zero/brent_zero.f90
                  general/cardamom_structures.f90
                )
set(MODELPARS
                  model/DALEC.A1.C1.D2.F2.H2.P1.004/src/DALEC.A1.C1.D2.F2.H2.P1.004_PARS.f90
                )
set(FILES2
                  model/model_shared.f90
                  general/cardamom_io.f90
                  general/model_likelihood_wrapper.f90
                  general/cardamom_main_utils.f90
                  general/MHMCMC_StressTests.f90
                )
if(${CMAKE_BUILD_TYPE} STREQUAL "RELEASE")
  list(APPEND FILES2 general/wrappers_for_R.f90)
endif()

set(MODELLIKELIHOODFILE
                  model/DALEC.A1.C1.D2.F2.H2.P1.004/likelihood/MODEL_LIKELIHOOD.f90
                )
set(MODELFILE
                  model/DALEC.A1.C1.D2.F2.H2.P1.004/src/DALEC.A1.C1.D2.F2.H2.P1.004.f90
                  model/DALEC.A1.C1.D2.F2.H2.P1.004/src/DALEC.A1.C1.D2.F2.H2.P1.004_R_interface.f90
                )

set(ALLMODELNAMES
  #DALEC_1005 
  #DALEC_1005a
DALEC.A1.C1.D2.F2.H1.P1.003
DALEC.A1.C1.D2.F2.H2.P1.004
DALEC.A1.C1.D2.F2.H2.P1.R1.005
DALEC.A1.C1.D2.F2.H2.P2.018
DALEC.A1.C1.D2.F2.H2.P5.021
DALEC.A1.C1.D2.F2.H2.P6.022
#DALEC.A1.C1.D2.F2.H3.P1.029
DALEC.A1.C1.D2.F2.H5.P1.037
DALEC.A1.C1.D2.F2.H6.P1.R5.032
DALEC.A1.C1.D2.F2.H4.P1.024
DALEC.A1.C2.D2.F2.H1.P3.R1.008
DALEC.A1.C2.D2.F2.H1.P4.R2.010
DALEC.A1.C2.D2.F2.H2.P1.R1.006
DALEC.A1.C2.D2.F2.H2.P2.R1.007
DALEC.A1.C2.D2.F2.H2.P2.R3.019
DALEC.A1.C2.D2.F2.H2.P3.R1.009
DALEC.A1.C2.D2.F2.H2.P4.R2.011
DALEC.A1.C2.D2.F2.H2.P7.R2.023
#DALEC.A2.C1.D2.F2.H2.P1.020
DALEC.A3.C1.D2.F2.H2.P1.030
DALEC.A3.C3.H2.M1.015 
DALEC.A3.H2.M2.017 
DALEC.A4.C6.D2.F2.H3.P10.026
DALEC.A4.C6.D2.F2.H2.P11.031
DALEC.A4.C6.D2.F2.H3.P12.033
DALEC.C1.D1.F2.P1.002
#DALEC.C3.M1.014 # unfinished edit - template on 15
DALEC.C4.D1.F2.012
DALEC.C5.D1.F2.P1.013
DALEC.D1.F2.001
DALEC.M2.016
)

add_library(cardamom-samplers-lib
              method/math.f90
              method/random_uniform.f90
              method/samplers_io/samplers_io.f90
              method/samplers_shared.f90
              method/MHMCMC/cardamom-MHMCMC.f90
              method/MHMCMC/MHMCMC.f90
              method/DEMCz.f90
)

set_target_properties(
  "cardamom-samplers-lib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
)
target_include_directories(
  "cardamom-samplers-lib"
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)
if(NOT EXISTS "${PROJECT_BINARY_DIR}/include")
  make_directory("${PROJECT_BINARY_DIR}/include")
endif()

  find_package(OpenMP COMPONENTS Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(cardamom-samplers-lib PUBLIC OpenMP::OpenMP_Fortran)
endif()

if(TEST_COMPILE_ALL)
  #  loop through compiling with each model in list above, to check that they compile
  foreach(n IN LISTS ALLMODELNAMES)
    set(FILES ${FILES1})
    list(APPEND FILES "model/${n}/src/${n}_PARS.f90")
    list(APPEND FILES ${FILES2})
    list(APPEND FILES "model/${n}/src/${n}.f90")
    list(APPEND FILES "model/${n}/src/${n}_R_interface.f90")
    list(APPEND FILES "model/${n}/likelihood/MODEL_LIKELIHOOD.f90")
    add_library("cardamom-lib-${n}" STATIC ${FILES})
    add_dependencies("cardamom-lib-${n}" cardamom-samplers-lib) 
    target_link_libraries("cardamom-lib-${n}" PUBLIC cardamom-samplers-lib)
    add_executable("cardamom-${n}" general/cardamom_main.f90)
    target_link_libraries("cardamom-${n}" PUBLIC "cardamom-lib-${n}")
  endforeach()
endif()

# now with the model we chose to run with
add_library(cardamom-lib ${FILES1} ${MODELPARS} ${FILES2} ${MODELLIKELIHOODFILE} ${MODELFILE})
add_dependencies(cardamom-lib cardamom-samplers-lib) 
target_link_libraries(cardamom-lib PUBLIC cardamom-samplers-lib)
add_executable(cardamom general/cardamom_main.f90)
target_link_libraries(cardamom PUBLIC cardamom-lib)
add_executable(cardamom-diffev general/cardamom_main_diffev.f90)
target_link_libraries(cardamom-diffev PUBLIC cardamom-lib)

