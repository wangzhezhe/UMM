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
$ ./addVar ../../../datasets/twoholes.vtk 
add field: v_center_dist
add field: v_sin
create vtk file: ../../../datasets/twoholesWithVar.vtk

$ ./resampleCompress ../../../datasets/twoholesWithVar.vtk v_center_dist 200
$ ./interp ../../../datasets/twoholes v_center_dist
or
$ ./resampleCompress ../../../datasets/twoholesWithVar.vtk v_sin 200
$ ./interp ../../../datasets/twoholes v_sin

```



