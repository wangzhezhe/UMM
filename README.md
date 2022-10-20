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
$ ./addVar ../UMM_datasets/twoholes.vtk 
add field: v_center_dist
add field: v_sin
generate vtk file: ../UMM_datasets/twoholesWithVar.vtk


$ ./resampleCompress ../UMM_datasets/twoholesWithVar.vtk v_center_dist 100

$ ./resampleCompress ../UMM_datasets/twoholesWithVar.vtk v_sin 100

$ ./interp ../datasets/sampleUmeshWithVarResample.vtk ../datasets/sampleUmesh.vtk 
```



