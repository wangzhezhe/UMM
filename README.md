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
$./addVar sampleUmesh.vtk 
$./resampleCompress sampleUmeshWithVar.vtk 100
```

