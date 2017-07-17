cmake_minimum_required(VERSION 2.6)

project(epic)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")    #enable c++ 11

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
link_directories(${OpenCV_LIBRARY_DIRS})

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})

include_directories(${CMAKE_SOURCE_DIR})
link_directories(${CMAKE_INSTALL_PREFIX}/lib)

file(GLOB HEADERS "cpp/*.hpp" "cpp/*.h")
file(GLOB FILES "cpp/*.cpp" "cpp/*.c")

add_library(epic SHARED ${HEADERS} ${FILES})
target_link_libraries(epic png jpeg lapack ${PNG_LIBRARY})

file(GLOB exes "epicflow.cpp")
foreach(exe IN ITEMS ${exes})
  get_filename_component(script_name ${exe} NAME_WE)
  
  # Add executable
  add_executable(${script_name} ${exe})
  target_link_libraries(${script_name}
    ${PROJECT_NAME}
    ${OpenCV_LIBS})
endforeach()
