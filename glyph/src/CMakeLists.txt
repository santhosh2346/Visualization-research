cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(SquidGlyph LANGUAGES CXX)

# Find VTK-m installation
find_package(VTKm REQUIRED COMPONENTS cont io)

if(TARGET vtkm::cont)
  add_executable(squid_glyph main.cxx)

  target_link_libraries(squid_glyph
    PRIVATE 
    vtkm::cont   # Core data structures
    vtkm::io     # For writing VTK files
  )

  # Configure for GPU/CPU portability
  vtkm_add_target_information(squid_glyph
    DROP_UNUSED_SYMBOLS ON
    MODIFY_CUDA_FLAGS ON
    DEVICE_SOURCES main.cxx
  )

  # Enforce C++17 standard
  set_target_properties(squid_glyph PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
  )

else()
  message(FATAL_ERROR "VTK-m not found. Set CMAKE_PREFIX_PATH to VTK-m install dir.")
endif()