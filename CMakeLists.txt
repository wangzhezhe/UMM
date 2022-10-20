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

#set (CMAKE_CXX_FLAGS "-DDEBUG_BUILD")
add_executable(interp interp.cpp)
target_link_libraries(interp ${VTK_LIBRARIES})

find_package(MGARD REQUIRED)

#find_package(mgard PATHS /home/qian/Software/MGARD/build_scripts/install-serial/lib/cmake/mgard/)
add_executable(resampleCompress resampleCompress.cpp)
target_link_libraries(resampleCompress ${VTK_LIBRARIES} mgard::mgard)
add_executable(resampleCompressRect resampleCompressRect.cpp)
target_link_libraries(resampleCompressRect ${VTK_LIBRARIES} mgard::mgard)
