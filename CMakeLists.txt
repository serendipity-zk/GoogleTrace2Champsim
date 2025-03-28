cmake_minimum_required(VERSION 3.16)
project(newParser)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_PREFIX_PATH "/home/kanzhu/ally-artifact/DynamoRIO-Linux-11.3.0-1/cmake")

find_package(DynamoRIO)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -O3")
set(CMAKE_BUILDE_TYPE "Release")

# enable fPIC for DynamoRIO 

# set(DynamoRIO_RPATH ON)
add_definitions(-DLINUX)
add_definitions(-DX86_64)

set(dynamorio_dir "/home/kanzhu/ally-artifact/dynamorio/clients/drcachesim")
add_executable(test main.cpp ${dynamorio_dir}/tools/common/decode_cache.cpp  ${dynamorio_dir}/tracer/raw2trace_shared.cpp)
set_target_properties(test PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/"
)

find_package(Threads REQUIRED)


configure_DynamoRIO_standalone(test)
target_include_directories(test PRIVATE ${PROJECT_SOURCE_DIR} "/home/kanzhu/ally-artifact/DynamoRIO-Linux-11.3.0-1/tools/include" "/home/kanzhu/ally-artifact/DynamoRIO-Linux-11.3.0-1/include")

target_link_directories(test PRIVATE "/home/kanzhu/ally-artifact/DynamoRIO-Linux-11.3.0-1/tools/lib64/release/")


function(add_all_subdirs_include_dirs BASE_DIR)
  # List all entries in the base directory.
  file(GLOB SUB_DIRS RELATIVE "${BASE_DIR}" "${BASE_DIR}/*")
  foreach(subdir ${SUB_DIRS})
    set(FULL_PATH "${BASE_DIR}/${subdir}")
    # Check if it's a directory.
    if(IS_DIRECTORY "${FULL_PATH}")
      message(STATUS "Adding include directory: ${FULL_PATH}")
      include_directories("${FULL_PATH}")
    endif()
  endforeach()
endfunction()


include_directories(${dynamorio_dir}/simulator)
include_directories(${dynamorio_dir}/common)
include_directories(${dynamorio_dir}/reader)
include_directories(${dynamorio_dir}/tracer)
include_directories(${dynamorio_dir}/scheduler)
include_directories(${dynamorio_dir}/tools/common)
include_directories(${dynamorio_dir}/tools)
include_directories("/home/kanzhu/ally-artifact/dynamorio/core/ir/")

include_directories(${dynamorio_dir})
include_directories("/home/kanzhu/ally-artifact/dynamorio/include")
include_directories("/home/kanzhu/ally-artifact/dynamorio/tools/include")
add_all_subdirs_include_dirs("/home/kanzhu/ally-artifact/dynamorio/ext/")
target_link_directories(test PRIVATE "/home/kanzhu/ally-artifact/DynamoRIO-Linux-11.3.0-1/ext/lib64/release/")
target_link_libraries(test Threads::Threads drmemtrace drmemtrace_analyzer drcovlib_static)