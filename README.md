UMM represents unstructured medsh for mgard, this is the experimental project to explore
how to leverage the existing method to let the mgard proceeds the unstructured data.


This is an exmaple to build the source code
Assuming there is a build dir in the src dir

```
cmake .. -DVTK_DIR=/home/zw/cworkspace/build/vtk
```

Example for mgard compression example on local env

```
cmake .. -Dmgard_DIR=/home/zw/cworkspace/build/MGARD/install -Dzstd_DIR=~/cworkspace/build/zstd/install/lib/cmake/zstd/
```


### Examples

```
$./addVar ../datasets/sampleUmesh.vtk 
$./resampleCompress ../datasets/sampleUmeshWithVar.vtk
```