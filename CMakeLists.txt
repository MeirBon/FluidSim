cmake_minimum_required(VERSION 3.10)
if (DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
            CACHE STRING "")
endif ()

project(FluidSim)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/CMakeModules")

find_package(OpenGL REQUIRED)

include_directories("${PROJECT_SOURCE_DIR}/include")

add_executable(${PROJECT_NAME}
        main.cpp
        Shader.cpp
        Shader.h
        Particle.cpp
        Particle.h
        Simulator.cpp
        Simulator.h
        Plane.cpp
        Plane.h
        Camera.cpp 
        Camera.h 
        Timer.h
        )

target_link_libraries(${PROJECT_NAME} OpenGL::GL)

if (NOT WIN32)
    find_package(glfw3 REQUIRED)
    find_package(GLM REQUIRED)
    find_package(GLEW REQUIRED)

    target_link_libraries(${PROJECT_NAME}
            glfw
            GLEW::GLEW)
else ()
    if (CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(LIB_DIR "x64")
    elseif (CMAKE_SIZEOF_VOID_P EQUAL 4)
        set(LIB_DIR "Win32")
    endif ()

    target_link_libraries(${PROJECT_NAME}
            "${PROJECT_SOURCE_DIR}/lib/${LIB_DIR}/glew32.lib"
            "${PROJECT_SOURCE_DIR}/lib/${LIB_DIR}/glfw3.lib"
            )

    file(GLOB DLLS "${PROJECT_SOURCE_DIR}/dll/${LIB_DIR}/*.dll")
    file(COPY ${DLLS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif ()