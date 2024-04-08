cmake_minimum_required(VERSION 3.20)
project(SKFlatAnalyzer VERSION 1.0.0)

# Find ROOT packages
find_package(ROOT REQUIRED COMPONENTS)
include(${ROOT_USE_FILE})

# Add subdirectories
add_subdirectory(DataFormats)
add_subdirectory(AnalyzerTools)
add_subdirectory(external/GEScaleSyst)
add_subdirectory(Analyzers)
