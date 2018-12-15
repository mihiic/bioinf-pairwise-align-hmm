cmake_minimum_required(VERSION 3.12)
project(bionf)

set(CMAKE_CXX_STANDARD 17)

add_executable(bionf main.cpp alignment_model.cpp alignment_model.h)