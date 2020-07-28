cmake_minimum_required(VERSION 3.13)

project(Integral_DT VERSION 2.0 LANGUAGES CXX)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/)

set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -O0")

include_directories(./include/)

file(GLOB CPPS "./src/*.cpp")

add_executable(${PROJECT_NAME} ${CPPS})

install(TARGETS ${PROJECT_NAME} DESTINATION bin)