# find LHAPDF
set(LHAPDF_INCLUDE_DIR $ENV{LHAPDFDIR}/include)
set(LHAPDF_LIBRARY $ENV{LHAPDFDIR}/lib/libLHAPDF.so)
set(GEScaleSyst_DIR $ENV{SKFlat_WD}/external/GEScaleSyst)

# include ROOT
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${LHAPDF_INCLUDE_DIR})
include_directories(${GEScaleSyst_DIR})

# Generate the ROOT dictionary before defining the library
file(GLOB HEADERS include/*.h)
ROOT_GENERATE_DICTIONARY(G__Analyzers 
    ${HEADERS}
    LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/Analyzers_LinkDef.hpp
)

# Define the library with source files and the generated dictionary
file(GLOB SOURCES src/*.C)
add_library(Analyzers SHARED
    ${SOURCES}
    G__Analyzers.cxx
)

# Specify where the library should look for header files
target_include_directories(Analyzers PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Link ROOT libraries
target_link_libraries(Analyzers ${ROOT_LIBRARIES})
target_link_libraries(Analyzers ${LHAPDF_LIBRARY})
target_link_libraries(Analyzers DataFormats)
target_link_libraries(Analyzers AnalyzerTools)
target_link_libraries(Analyzers GEScaleSyst)

# Specify where to install the library
install(TARGETS Analyzers DESTINATION lib)

# Iinstall .rootmap and .pcm files to the lib
set(ROOTMAP "${CMAKE_CURRENT_BINARY_DIR}/libAnalyzers.rootmap")
set(PCM "${CMAKE_CURRENT_BINARY_DIR}/libAnalyzers_rdict.pcm")
message(STATUS "ROOTMAP: ${ROOTMAP}")
message(STATUS "PCM: ${PCM}")
install(FILES ${ROOTMAP} ${PCM} DESTINATION lib)
