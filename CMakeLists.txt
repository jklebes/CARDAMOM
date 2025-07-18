cmake_minimum_required(VERSION 3.28)

project(CARDAMOM VERSION 1.0 #TODO
	DESCRIPTION "University of Edinburgh / NCEO (UK) CARbon DAta MOdel FraMework" 
  LANGUAGES Fortran)
enable_language(Fortran)

include(CTest)

set (
  CMAKE_Fortran_FLAGS
  "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -O2 -g -pg -Wall -no-pie" 
  #"${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -O3" 
)
# I put a lot of extra flags here to try to make gprof work, but it turned out to be
# related to shared libraries - gprof can't do shared libraries
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg -g")

## Code coverage analysis
#add_compile_options(--coverage -fprofile-arcs -ftest-coverage) 
#add_link_options(-lgcov --coverage) 

option(BUILD_SHARED_LIBS "Build using shared libraries" OFF) #different needs for grpof and R

add_subdirectory(LIBRARY/CARDAMOM_F)

set_target_properties(
  "cardamom-lib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
)
target_include_directories(
  "cardamom-lib"
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)
if(NOT EXISTS "${PROJECT_BINARY_DIR}/include")
  make_directory("${PROJECT_BINARY_DIR}/include")
endif()

#install(TARGETS cardamom DESTINATION "bin")

if(ENABLE_TESTING)
  # add the testsuite
  enable_testing()
  add_subdirectory("test")
endif()
