
#export VTK_DIR=/home/qian/Software/vtk/build
find_package(VTK PATHS /home/qian/Software/vtk/build 
COMPONENTS
IOLegacy
IOXML
IOImage
FiltersCore
)

add_executable(ErrorControlledResample ErrorControlledResample.cpp AdpG.cpp)
target_link_libraries(ErrorControlledResample ${VTK_LIBRARIES} mgard::mgard)

