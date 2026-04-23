cmake_minimum_required(VERSION 3.20)

project(parallel_plate_scan LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(PARALLEL_PLATE_SHARED_GARFIELD_HINT "${CMAKE_SOURCE_DIR}/../../local/garfield")

if(DEFINED ENV{ROOTSYS})
  list(PREPEND CMAKE_PREFIX_PATH "$ENV{ROOTSYS}")
endif()
if(DEFINED ENV{GARFIELD_INSTALL} AND NOT "$ENV{GARFIELD_INSTALL}" STREQUAL "")
  list(PREPEND CMAKE_PREFIX_PATH "$ENV{GARFIELD_INSTALL}")
elseif(EXISTS "${PARALLEL_PLATE_SHARED_GARFIELD_HINT}/lib/cmake/Garfield/GarfieldConfig.cmake")
  list(PREPEND CMAKE_PREFIX_PATH "${PARALLEL_PLATE_SHARED_GARFIELD_HINT}")
endif()

find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(Garfield CONFIG QUIET)

if(NOT Garfield_FOUND)
  message(FATAL_ERROR
    "Garfield++ was not found. Set GARFIELD_INSTALL or CMAKE_PREFIX_PATH to a Garfield installation, "
    "or source ../../local/garfield/share/Garfield/setupGarfield.sh when using the shared workspace install."
  )
endif()

add_executable(parallel_plate_scan
  src/parallel_plate_scan.cc
)

target_include_directories(parallel_plate_scan PRIVATE
  ${CMAKE_SOURCE_DIR}/third_party
)

target_link_libraries(parallel_plate_scan PRIVATE
  Garfield::Garfield
  ROOT::Core
  ROOT::RIO
  ROOT::Hist
  ROOT::Graf
  ROOT::Gpad
)

target_compile_options(parallel_plate_scan PRIVATE
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
    -Wall
    -Wextra
    -Wpedantic
    -Wshadow
  >
)

set(_runtime_paths "")
if(DEFINED ENV{GARFIELD_INSTALL} AND NOT "$ENV{GARFIELD_INSTALL}" STREQUAL "")
  list(APPEND _runtime_paths "$ENV{GARFIELD_INSTALL}/lib")
elseif(EXISTS "${PARALLEL_PLATE_SHARED_GARFIELD_HINT}/lib")
  list(APPEND _runtime_paths "${PARALLEL_PLATE_SHARED_GARFIELD_HINT}/lib")
endif()
if(DEFINED ENV{ROOTSYS})
  list(APPEND _runtime_paths "$ENV{ROOTSYS}/lib")
endif()

set_target_properties(parallel_plate_scan PROPERTIES
  BUILD_RPATH "${_runtime_paths}"
  INSTALL_RPATH "${_runtime_paths}"
)

include(CTest)

if(BUILD_TESTING)
  add_test(
    NAME smoke_parallel_plate_scan
    COMMAND parallel_plate_scan
      --config ${CMAKE_SOURCE_DIR}/config/smoke_parallel_plate.json
      --scan d
      --out ${CMAKE_BINARY_DIR}/smoke-output
  )
  set_tests_properties(smoke_parallel_plate_scan PROPERTIES
    ENVIRONMENT
      "GARFIELD_INSTALL=${PARALLEL_PLATE_SHARED_GARFIELD_HINT};HEED_DATABASE=${PARALLEL_PLATE_SHARED_GARFIELD_HINT}/share/Heed/database"
  )
endif()
