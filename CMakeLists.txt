cmake_minimum_required (VERSION 3.2)
project(UMM LANGUAGES CXX)

#export VTK_DIR=/home/qian/Software/vtk/build

find_package(VTK PATHS /home/qian/Software/vtk/build 
COMPONENTS
IOLegacy
IOXML
IOImage
FiltersCore
)

add_executable(addVar addVar.cpp)
target_link_libraries(addVar ${VTK_LIBRARIES})

find_package(mgard)
add_executable(resampleCompress resampleCompress.cpp)
target_link_libraries(resampleCompress ${VTK_LIBRARIES} mgard::mgard)
