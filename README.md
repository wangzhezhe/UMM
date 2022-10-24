UMM represents unstructured medsh for mgard, this is the experimental project to explore
how to leverage the existing method to let the mgard proceeds the unstructured data.


install mgard and vtk, then build the project by executing

```
$ cd installscripts
$ /bin/bash summit_install.sh
```

the datasets is located in 

```
/gpfs/alpine/proj-shared/csc143/zhewang/datasets/UMM_datasets
```

### Examples

```
$ ./addVar ../../../datasets/large_case4.vtk 
add field: v_center_dist
add field: v_sin
create vtk file: ../../../datasets/large_case4WithVar.vtk

$ ./resampleCompress ../../../datasets/large_case4WithVar.vtk v_center_dist 5000
$ ./interp ../../../datasets/twoholes v_center_dist
or
$ ./resampleCompress ../../../datasets/large_case4WithVar.vtk v_sin 5000
$ ./interp ../../../datasets/large_case4 v_sin 
load files: ../../../datasets/large_case4.vtk,../../../datasets/large_case4WithVar.vtk,../../../datasets/large_case4WithVarResample.vtk
add field:v_sin_interp
create vtk file: ../../../datasets/large_case4Interp.vtk
accumulated error 8.87303
max error 0.00150419
avg error 1.5308e-05

```



