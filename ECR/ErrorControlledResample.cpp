

/*
input: unstructured data grid (Usg), fieldname(field), error bound(e)
return: uniform grid with the controlled error for a specific data field

// get the uniform grid to make sure the number of vertex of unstructured grid is less then k
AdaptG = InitUniGrid(Usg, k)

// refine the adaptive grid to make sure the interp error of original data is less then a bounud 
Refine(Usg, field, e, AdapG)

// create the uniform grid based on the AdapG
Unig = AdapG.getUniGrid()
*/

