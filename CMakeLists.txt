cmake_minimum_required(VERSION 3.28)

project(CARDAMOM VERSION 1.0 #TODO
	DESCRIPTION "University of Edinburgh / NCEO (UK) CARbon DAta MOdel FraMework" 
  LANGUAGES Fortran)
enable_language(Fortran)

add_compile_options(-ffree-line-length-0 -g -pg -Wuninitialized -Wmaybe-uninitialized)

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

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
