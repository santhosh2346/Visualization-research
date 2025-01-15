cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(VectorFieldGenerator)

find_package(VTK COMPONENTS 
  CommonCore
  CommonDataModel
  IOLegacy
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "Unable to find the VTK build folder.")
endif()

add_executable(VectorFieldGenerator VectorFieldGenerator.cxx)
target_link_libraries(VectorFieldGenerator PRIVATE ${VTK_LIBRARIES})

# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS VectorFieldGenerator
  MODULES ${VTK_LIBRARIES}
)
