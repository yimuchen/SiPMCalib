file(GLOB SiPMCalc_src src/*.cc)
add_library(SiPMCalc SHARED ${SiPMCalc_src})
target_link_libraries(SiPMCalc
  Common
  PlotUtils
  ${ROOT_LIBRARIES}
)

## Function for compiling unit tests
function(make_sipmcalc_bin binfile)
  get_filename_component( binname ${binfile} NAME_WE )
  set( exebin "SiPM_${binname}" )
  add_executable( ${exebin} ${binfile} )
  set_target_properties( ${exebin} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/bin/SiPMCalc )
  target_link_libraries( ${exebin}
    SiPMCalc
    Common
    PlotUtils
    MathUtils
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES}
     )
endfunction()

## Listing all unit test functions to perform
file(GLOB sipmcalc_bin_files "bin/*.cc")
foreach(sipmcalc_bin_file ${sipmcalc_bin_files})
  make_sipmcalc_bin( ${sipmcalc_bin_file} )
endforeach()
