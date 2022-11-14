

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

#include <iostream>
#include <vtkUnstructuredGridReader.h>
#include <vtkDataSetWriter.h>
#include "AdpG.hpp"

int main(int argc, char *argv[])
{
  // Parse command line arguments
  if (argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " <Filename> <FieldName>" << std::endl;
    return EXIT_FAILURE;
  }

  std::string vtkUStructuredFile = argv[1];

  std::string fieldName = argv[2];

  // load the vtk file
  vtkSmartPointer<vtkUnstructuredGridReader> reader =
      vtkSmartPointer<vtkUnstructuredGridReader>::New();
  reader->SetFileName(vtkUStructuredFile.c_str());
  reader->Update();

  // get the specific polydata and check the results
  vtkUnstructuredGrid *unsGridData = reader->GetOutput();

  ADPG AdpG(unsGridData, 20, 2);

  AdpG.BuildAdpGrid();

  std::cout << "---ok to build the adp grid" << std::endl;

  vtkSmartPointer<vtkUnstructuredGrid> adpUnsGridData = AdpG.ConvertToUnstructuredGrid();

  if (adpUnsGridData != nullptr)
  {
    adpUnsGridData->PrintSelf(std::cout, vtkIndent(1));
    std::cout << "---ok to convert the adp grid" << std::endl;
  }
  else
  {
    std::cout << "adpUnsGridData is null" << std::endl;
  }

  vtkSmartPointer<vtkDataSetWriter> writer =
      vtkSmartPointer<vtkDataSetWriter>::New();
  std::string fileSuffix = vtkUStructuredFile.substr(0, vtkUStructuredFile.size() - 4);
  std::string outputFileName = fileSuffix + std::string("Adp.vtk");

  writer->SetFileName(outputFileName.c_str());

  // get the specific polydata and check the results
  writer->SetInputData(adpUnsGridData);
  // Optional - set the mode. The default is binary.
  // writer->SetDataModeToBinary();
  // writer->SetDataModeToAscii();
  writer->Write();
}