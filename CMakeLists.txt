cmake_minimum_required(VERSION 4.0)
set(DEV_MODULE Development.Module)
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/install" CACHE PATH "" FORCE)

project(opti_tools LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(optiToolsLib STATIC)
target_include_directories(optiToolsLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_sources(
        optiToolsLib
        PRIVATE src/solver/euler.cpp
        PUBLIC FILE_SET
        HEADERS
        BASE_DIRS
        src
        FILES
        src/solver/euler.h
        src/math/vec3d.h)

include(FetchContent)

FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt
        GIT_TAG e69e5f977d458f2650bb346dadf2ad30c5320281) # 10.2.1
FetchContent_MakeAvailable(fmt)

target_link_libraries(optiToolsLib fmt::fmt)


if (${BUILD_FOR} MATCHES "PYTHON")
    find_package(
            Python
            COMPONENTS Interpreter Development
            REQUIRED)

    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/nanobind)
    nanobind_add_module(opti_tools src/module/module.cpp)
    set_property(TARGET opti_tools PROPERTY EXPORT_NAME opti_tools)

    target_link_libraries(opti_tools PRIVATE optiToolsLib)
    target_link_libraries(opti_tools PRIVATE -static -lpthread -static-libgcc
            -static-libstdc++)
    set_target_properties(opti_tools PROPERTIES COMPILE_FLAGS -fvisibility=hidden)
    set_target_properties(opti_tools PROPERTIES OUTPUT_NAME opti_tools)
    set_target_properties(opti_tools PROPERTIES SUFFIX ".pyd")

    nanobind_add_stub(
            opti_tools_stub
            MODULE
            opti_tools
            OUTPUT
            ${CMAKE_CURRENT_BINARY_DIR}/opti_tools.pyi
            PYTHON_PATH
            $<TARGET_FILE_DIR:opti_tools>
            DEPENDS
            opti_tools)

    install(TARGETS opti_tools DESTINATION .)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/opti_tools.pyi DESTINATION .)
else ()
    add_executable(cpp_tests src/test.cpp)
    target_link_libraries(cpp_tests PRIVATE optiToolsLib)
endif ()
