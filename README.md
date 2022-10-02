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


questions:

naive way to process the unstructructed data, how to put them into mgard
unstrucutred supported for mgard?

construct the baseline, derived properties examples?

how about the 2d or 3d case, coordinates how to do that?

data set, do not too symetric

try to validate if the new added points are zero? from the mgard's perspective