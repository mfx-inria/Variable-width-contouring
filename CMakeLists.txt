cmake_minimum_required(VERSION 3.0)
cmake_policy(SET CMP0011 NEW)
cmake_policy(SET CMP0054 NEW)

project(MedialAxisButchery)

# = = = = = = = = = = = = = = = = CGAL

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

set(LIBS ${LIBS} ${CGAL_LIBRARIES})

# = = = = = = = = = = = = = = = = CAIRO

option(NO_CAIRO "Disable PDF output. No link with libcairo" OFF)

set(NO_CAIRO_C_VALUE 1)
if(NOT NO_CAIRO)
set(NO_CAIRO_C_VALUE 0)
find_library(CAIRO_LIBRARY NAMES cairo)
set(LIBS ${LIBS} ${CAIRO_LIBRARY})
endif(NOT NO_CAIRO)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file(
  "${PROJECT_SOURCE_DIR}/fill_config.h.in"
  "${PROJECT_BINARY_DIR}/fill_config.h"
)

if( NOT WIN32 )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif()

include_directories("${PROJECT_SOURCE_DIR}")
file(GLOB_RECURSE FILL_SOURCES      src/*.cpp)
if(NO_CAIRO) # remove CairoSurface.*
	list(REMOVE_ITEM FILL_SOURCES ${PROJECT_SOURCE_DIR}/src/CairoSurface.cpp)
endif(NO_CAIRO)

# ================================

set(LIBS ${LIBS} "/opt/local/lib/libpolyclipping.dylib")

#set_source_files_properties(${FILL_SOURCES} PROPERTIES COMPILE_FLAGS "${CPP11_FLAG}")

set(SOURCES ${FILL_SOURCES} ${lib_sources})

set(APP_NAME fill)

add_executable(${APP_NAME} ${SOURCES})
target_link_libraries(${APP_NAME} ${LIBS})

message( STATUS "FILL_SOURCES: " ${FILL_SOURCES})
message( STATUS "LIBS: " ${LIBS} )

