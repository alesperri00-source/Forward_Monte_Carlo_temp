cmake_minimum_required(VERSION 3.6)
project(MC_with_T)

set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

include_directories(/home/alessandro/alessandro/PhD_Alessandro/first_project/MaxEnt-Chromosome-Caulobacter-0.1/Forward_Monte_Carlo_with_T/eigen-3.4.0)
set(SOURCE_FILES main.cpp Initialize.h Moves.h Energy_changes.h global.h Functions.h)
add_executable(MC_with_T ${SOURCE_FILES})