cmake_minimum_required(VERSION 3.28)

project(CARDAMOM VERSION 1.0 #TODO
	DESCRIPTION "University of Edinburgh / NCEO (UK) CARbon DAta MOdel FraMework" 
  LANGUAGES Fortran)
enable_language(Fortran)

add_subdirectory(LIBRARY/CARDAMOM_F)

install(TARGETS cardamom DESTINATION "bin")
