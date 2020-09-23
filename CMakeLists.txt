cmake_minimum_required(VERSION 3.0)
cmake_policy(SET CMP0011 NEW)
cmake_policy(SET CMP0054 NEW)

project(MedialAxisButchery)

option(USE_CGAL_MEDIAL_AXIS "Compile with CGAL and use it to compute the medial axis." ON)

option(USE_BOOST_MEDIAL_AXIS "Compile with BOOST and use it to compute the medial axis." ON)

option(USE_CAIRO_PDF "Compile with libcairo to enable PDF output." ON)

if( NOT WIN32 )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif()

# = = = = = = = = = = = = = = = = MUST HAVE AT LEAST ONE WAY TO COMPUTE MEDIAL AXIS

if( (NOT USE_CGAL_MEDIAL_AXIS) AND (NOT USE_BOOST_MEDIAL_AXIS) )
    message(FATAL_ERROR "Please choose at least one way to compute the medial-axis.")
endif()

# = = = = = = = = = = = = = = = = SOURCES

include_directories("${PROJECT_SOURCE_DIR}")
file(GLOB_RECURSE FILL_SOURCES      src/*.cpp)

# = = = = = = = = = = = = = = = = CGAL

if(USE_CGAL_MEDIAL_AXIS)

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

    set(LIBS ${LIBS} ${CGAL_LIBRARIES})

    set(CGAL_C_VALUE 1)
else(USE_CGAL_MEDIAL_AXIS)
    list(REMOVE_ITEM FILL_SOURCES ${PROJECT_SOURCE_DIR}/src/MATGraphBuilder_with_CGAL.cpp)
    set(CGAL_C_VALUE 0)
endif()

# = = = = = = = = = = = = = = = = BOOST

if(USE_BOOST_MEDIAL_AXIS)
    if( NOT USE_CGAL_MEDIAL_AXIS )
        set(Boost_USE_STATIC_LIBS OFF)
        set(Boost_USE_MULTITHREADED ON)
        set(Boost_USE_STATIC_RUNTIME OFF)
        find_package(Boost 1.50)# REQUIRED COMPONENTS filesystem)
        include_directories(${Boost_INCLUDE_DIRS})
    endif( NOT USE_CGAL_MEDIAL_AXIS )
    set(BOOST_C_VALUE 1)
else(USE_BOOST_MEDIAL_AXIS)
    list(REMOVE_ITEM FILL_SOURCES ${PROJECT_SOURCE_DIR}/src/MATGraphBuilder_with_BOOST.cpp)
    list(REMOVE_ITEM FILL_SOURCES ${PROJECT_SOURCE_DIR}/src/VariableWidthContouringForIceSL.cpp)
    set(BOOST_C_VALUE 0)
endif(USE_BOOST_MEDIAL_AXIS)

# = = = = = = = = = = = = = = = = CAIRO

set(NO_CAIRO_C_VALUE 1)
if(USE_CAIRO_PDF)
    set(NO_CAIRO_C_VALUE 0)
    find_library(CAIRO_LIBRARY NAMES cairo)
    set(LIBS ${LIBS} ${CAIRO_LIBRARY})
else(USE_CAIRO_PDF)
	list(REMOVE_ITEM FILL_SOURCES ${PROJECT_SOURCE_DIR}/src/CairoSurface.cpp)
endif(USE_CAIRO_PDF)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file(
  "${PROJECT_SOURCE_DIR}/fill_config.h.in"
  "${PROJECT_BINARY_DIR}/fill_config.h"
)
include_directories(${PROJECT_BINARY_DIR})

# ================================

if(APPLE)
	set(LIBS ${LIBS} "/opt/local/lib/libpolyclipping.dylib")
else(APPLE)
	set(LIBS ${LIBS} "/usr/local/lib/libpolyclipping.so")
endif(APPLE)

#set_source_files_properties(${FILL_SOURCES} PROPERTIES COMPILE_FLAGS "${CPP11_FLAG}")

set(SOURCES ${FILL_SOURCES} ${lib_sources})

set(APP_NAME fill)

add_executable(${APP_NAME} ${SOURCES})
target_link_libraries(${APP_NAME} ${LIBS})

message( STATUS "FILL_SOURCES: " ${FILL_SOURCES})
message( STATUS "LIBS: " ${LIBS} )

