#include <iostream>
#include <vtkCellLocator.h>
#include <vtkDataSetReader.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>

int main(int argc, char *argv[]) {
  // Parse command line arguments
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0]
              << " <uniformGrid file> <unstructuredGrid file>" << std::endl;
    return EXIT_FAILURE;
  }

  std::string structuredGridFile = argv[1];

  std::string unstructuredGridFile = argv[2];

  // load the uniform grid
  vtkSmartPointer<vtkStructuredPointsReader> imgReader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();
  imgReader->SetFileName(structuredGridFile.c_str());
  imgReader->Update();

  vtkStructuredPoints *struPoints = imgReader->GetOutput();

  // std::cout << "---check loaded struGrid---" << std::endl;
  // struGrid->Print(std::cout);
  std::string fieldName = "v_center_dist";

  vtkPointData *pointData = struPoints->GetPointData();
  auto pointDataArray = pointData->GetScalars(fieldName.c_str());
  // pointDataArray->Print(std::cout);

  // load the unstructured grid
  vtkSmartPointer<vtkUnstructuredGridReader> reader =
      vtkSmartPointer<vtkUnstructuredGridReader>::New();
  reader->SetFileName(unstructuredGridFile.c_str());
  reader->Update();

  // get the specific polydata and check the results
  vtkUnstructuredGrid *unsGridData = reader->GetOutput();

  // std::cout << "---check loaded unsGridData---" << std::endl;
  // unsGridData->Print(std::cout);

  // range each point in the unsGridData
  // find the assocaited three points in the struPoints
  // then using the interpolation to do associated operation

  // build the cell locator based on stuctured points
  vtkNew<vtkCellLocator> cellLocator;
  cellLocator->SetDataSet(struPoints);
  cellLocator->BuildLocator();

  double pointcoord[3];
  for (vtkIdType id = 0; id < unsGridData->GetNumberOfPoints(); id++) {

    unsGridData->GetPoint(id, pointcoord);

    // use this pointcoord to find the specific cell in the structured grid
    vtkIdType interpCellId = cellLocator->FindCell(pointcoord);
    vtkCell *interpCell = struPoints->GetCell(interpCellId);

    int interpPointsNum = interpCell->GetNumberOfPoints();

    if(interpCellId == -1){
      // do not find the cell id
      // use the closetPoint to find
      double closestPoint[3];
      vtkIdType closestCellId;
      vtkGenericCell*closestCell;
      int subid=0;
      double dist=0;
      cellLocator->FindClosestPoint(pointcoord,closestPoint,closestCellId,subid,dist);
      interpCellId = closestCellId;
    }

    std::cout << "id " << id << " x " << pointcoord[0] << " y " << pointcoord[1]
              << " z " << pointcoord[2] << " interpCellId is " << interpCellId
              << " interpPointsNum is " << interpPointsNum << std::endl;

    // TODO, how to process the case that the points at boundry returns -1 for
    // the FindCell
  }
}