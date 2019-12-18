cmake_minimum_required(VERSION 3.15)
project(MSE760 C)

set(CMAKE_C_STANDARD 99)

add_executable(MSE760 main.c)
add_executable(Homework3 Homework3/hw3.c)
target_link_libraries(Homework3 m)

add_executable(Homework2 Homework2/hw2.c Homework2/hw2Extented.h)
target_link_libraries(Homework2 m)