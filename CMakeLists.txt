# vim: set sts=0 ts=4 sw=4 noet tw=0 :
cmake_minimum_required(VERSION 3.0)
cmake_policy(SET CMP0011 NEW)
cmake_policy(SET CMP0054 NEW)

project(MedialAxisButchery)

option(VWC_USE_CGAL_MEDIAL_AXIS "Compile with CGAL and use it to compute the medial axis." ON)

option(VWC_USE_BOOST_MEDIAL_AXIS "Compile with BOOST and use it to compute the medial axis." ON)
option(VWC_USE_EXTERNAL_BOOST "Compile with some system-installed BOOST, not the local one." OFF)

option(VWC_USE_CAIRO_PDF "Compile with libcairo to enable PDF output." ON)

option(VWC_MAKE_ICESL_PLUGIN "Add a second build target for building a plugin for IceSL" OFF)

#SET(ICESL_PROJECT OFF)
IF(ICESL_PROJECT)
	SET(VWC_USE_CGAL_MEDIAL_AXIS OFF)
	SET(VWC_USE_BOOST_MEDIAL_AXIS ON)
	SET(VWC_USE_EXTERNAL_BOOST OFF)
	SET(VWC_USE_CAIRO_PDF OFF)
	SET(VWC_MAKE_ICESL_PLUGIN ON)
ENDIF()

if( NOT WIN32 )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif()

# = = = = = = = = = = = = = = = = MUST HAVE AT LEAST ONE WAY TO COMPUTE MEDIAL AXIS

set(ICESL_PLUGIN_VALUE 0)
if( VWC_MAKE_ICESL_PLUGIN )
	set(VWC_USE_BOOST_MEDIAL_AXIS ON)
	if( NOT ICESL_PROJECT )
		ADD_SUBDIRECTORY(icesl-interface)
	endif()
    set(ICESL_PLUGIN_VALUE 1)
endif()

if( (NOT VWC_USE_CGAL_MEDIAL_AXIS) AND (NOT VWC_USE_BOOST_MEDIAL_AXIS) )
    message(FATAL_ERROR "Please choose at least one way to compute the medial-axis.")
endif()

# = = = = = = = = = = = = = = = = SOURCES

set(BASE_SOURCES
src/MATGraph.cpp
src/utils.cpp
src/VariableWidthContouringSampling.cpp
src/Sampling.cpp
src/VariableWidthContouring.cpp
src/VariableWidthContouringLabeling.cpp
src/vec.cpp)

include_directories("${PROJECT_SOURCE_DIR}")

set(PLUGIN_SOURCES ${BASE_SOURCES}
    src/VWCInfillerPlugin.cpp)

set(FILL_SOURCES ${BASE_SOURCES}
    src/main.cpp)

# = = = = = = = = = = = = = = = = CLIPPER

find_library(CLIPPER_LIB polyclipping)

# = = = = = = = = = = = = = = = = CGAL

set(CGAL_C_VALUE 0)
if(VWC_USE_CGAL_MEDIAL_AXIS)
    set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE true)
    if(WIN32)
      set(CGAL_DIR "C:/dev/CGAL-4.10/")
      LIST(APPEND CMAKE_MODULE_PATH "C:/dev/CGAL-4.10/cmake/modules/")
    endif()
    find_package(CGAL COMPONENTS Core)
    if(NOT CGAL_FOUND)
            error("CGAL is missing")
    endif(NOT CGAL_FOUND)
    include (${CGAL_USE_FILE})
    include_directories(${CGAL_INCLUDE_DIRS})
    message(STATUS "CGAL include: ${CGAL_INCLUDE_DIRS}")
    set(FILL_LIBS ${CGAL_LIBRARIES})
    set(CGAL_C_VALUE 1)
    set(FILL_SOURCES ${FILL_SOURCES} src/MATGraphBuilder_with_CGAL.cpp)
endif()

# = = = = = = = = = = = = = = = = BOOST

set(BOOST_C_VALUE 0)
if(VWC_USE_BOOST_MEDIAL_AXIS)
	if( NOT VWC_USE_CGAL_MEDIAL_AXIS ) # Because CGAL already instructs CMake to use Boost
		if( VWC_USE_EXTERNAL_BOOST )
			set(Boost_USE_STATIC_LIBS OFF)
			set(Boost_USE_MULTITHREADED ON)
			set(Boost_USE_STATIC_RUNTIME OFF)
			find_package(Boost 1.50)
			include_directories(${Boost_INCLUDE_DIRS})
		else()
			set(Boost_USE_STATIC_LIBS ON)
			set(Boost_USE_MULTITHREADED OFF)
			set(Boost_USE_STATIC_RUNTIME ON)
		endif()
	endif( NOT VWC_USE_CGAL_MEDIAL_AXIS )
	set(BOOST_C_VALUE 1)
	set(FILL_SOURCES ${FILL_SOURCES} src/MATGraphBuilder_with_BOOST.cpp)
	set(PLUGIN_SOURCES ${PLUGIN_SOURCES} src/MATGraphBuilder_with_BOOST.cpp)
endif(VWC_USE_BOOST_MEDIAL_AXIS)

# = = = = = = = = = = = = = = = = CAIRO

set(NO_CAIRO_C_VALUE 1)
if(VWC_USE_CAIRO_PDF)
	set(NO_CAIRO_C_VALUE 0)
	find_library(CAIRO_LIBRARY NAMES cairo)
	set(FILL_LIBS ${FILL_LIBS} ${CAIRO_LIBRARY})
	set(FILL_SOURCES ${FILL_SOURCES} src/CairoSurface.cpp)
endif(VWC_USE_CAIRO_PDF)

# configure a header file to pass some of the CMake settings to the source code
configure_file(
  "${PROJECT_SOURCE_DIR}/fill_config.h.in"
  "${PROJECT_BINARY_DIR}/fill_config.h"
)
include_directories(${PROJECT_BINARY_DIR})

# = = = = = = = = = = = = = = = = FILL EXECUTABLE

if(NOT VWC_MAKE_ICESL_PLUGIN)
	add_executable(fill ${FILL_SOURCES})
	target_link_libraries(fill ${FILL_LIBS})
	if(VWC_USE_BOOST_MEDIAL_AXIS)
		target_link_libraries(fill ${CLIPPER_LIB})
	endif()
endif()

# = = = = = = = = = = = = = = = = ICESL PLUGIN

if(VWC_MAKE_ICESL_PLUGIN)
	set(PNAME variableWidthContouring)
	add_library(${PNAME} MODULE ${PLUGIN_SOURCES})
	if( NOT VWC_USE_EXTERNAL_BOOST )
		target_include_directories(${PNAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
	endif()
	if( ICESL_PROJECT )
		SET_TARGET_PROPERTIES(${PNAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../../icesl-plugins")
		SET_TARGET_PROPERTIES(${PNAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY_DEBUG "${CMAKE_CURRENT_SOURCE_DIR}/../../icesl-plugins-d")
		SET_TARGET_PROPERTIES(${PNAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY_RELEASE "${CMAKE_CURRENT_SOURCE_DIR}/../../icesl-plugins")
		SET_TARGET_PROPERTIES(${PNAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO "${CMAKE_CURRENT_SOURCE_DIR}/../../icesl-plugins")
		SET_TARGET_PROPERTIES(${PNAME} PROPERTIES DEBUG_POSTFIX "-d")
	endif()
	target_link_libraries(${PNAME} icesl_interface)
endif()
