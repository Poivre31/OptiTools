cmake_minimum_required(VERSION 3.15)
project(example LANGUAGES CXX)

find_package(Python COMPONENTS Interpreter Development REQUIRED)

add_subdirectory(pybind11)
pybind11_add_module(example module.cpp)

add_custom_command(TARGET example POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    $<TARGET_FILE:example>
    ${CMAKE_SOURCE_DIR}/python/example/_example.pyd
)