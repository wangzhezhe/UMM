#include <iostream>
#include <vector>

#include <vtkCellLocator.h>
#include <vtkDataSetReader.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkDoubleArray.h>
#include <vtkUnstructuredGridWriter.h>

#ifdef DEBUG_BUILD
#define DEBUG(x) std::cout << x << std::endl;
#else
#define DEBUG(x) \
  do             \
  {              \
  } while (0)
#endif

// input coordinates of 4 points
// input fieldvalue for 4 points

// input coordinates of dedicated point
// get the interpolated calue of dedicated point

double Interpolate2D(std::vector<std::vector<double>> inputCoor, std::vector<double> fieldValues, double dedicatedPoint[3])
{
  if (inputCoor.size() != fieldValues.size())
  {
    throw std::runtime_error("fields value should equal to the input coordinates");
  }

  // the point sequence is bottom left corner 0 (0,0), bottom right corner 1 (1,0)
  // upper left corner 2 (0,1) , upper right corner 3 (1,1)
  // https://www.omnicalculator.com/math/bilinear-interpolation

  double x1 = inputCoor[0][0];
  double y1 = inputCoor[0][1];

  double x2 = inputCoor[3][0];
  double y2 = inputCoor[3][1];

  // look at the bound, make sure the targeted point is in the bounds
  // there are tolerance issues
  double tolerance = 0.000001;

  // dedicatedPoint[0] should within x1 and x2
  // dedicatedPoint[1] should within y1 and y2
  double tx = dedicatedPoint[0];
  double ty = dedicatedPoint[1];
  if ((tx < x1 - tolerance) && (tx > x2 + tolerance))
  {
    DEBUG(tx << " , " << x1 << " , " << x2);
    throw std::runtime_error("tx is out of bounds");
  }

  if ((ty < y1 - tolerance) && (ty > y2 + tolerance))
  {
    DEBUG(ty << " , " << y1 << " , " << y2);
    throw std::runtime_error("ty is out of bounds");
  }

  double q11 = fieldValues[0];
  double q21 = fieldValues[1];
  double q12 = fieldValues[2];
  double q22 = fieldValues[3];

  DEBUG("x1 y1 x2 y2 " << x1 << "," << y1 << "," << x2 << "," << y2);
  DEBUG("q11 q21 q12 q22 " << q11 << "," << q21 << "," << q12 << "," << q22);

  double r1 = q11 * ((x2 - tx) / (x2 - x1)) + q21 * ((tx - x1) / (x2 - x1));
  double r2 = q12 * ((x2 - tx) / (x2 - x1)) + q22 * ((tx - x1) / (x2 - x1));

  double p = r1 * ((y2 - ty) / (y2 - y1)) + r2 * ((ty - y1) / (y2 - y1));
  DEBUG("interpolated p " << p);

  return p;
}

// compute the diff between the interpolated value and the original value
void computeDiff(vtkDoubleArray *interpFieldArray, vtkDataArray *originalFieldArray)
{
  if (interpFieldArray->GetNumberOfComponents() != originalFieldArray->GetNumberOfComponents())
  {
    throw std::runtime_error("field array have different components");
  }

    if (interpFieldArray->GetNumberOfTuples() != originalFieldArray->GetNumberOfTuples())
  {
    throw std::runtime_error("field array have different tuples");
  }
  double accuError = 0;
  double maxError = 0;
  double avgError = 0;
  for (int i = 0; i < interpFieldArray->GetNumberOfTuples(); i++)
  {
    double *v1 = interpFieldArray->GetTuple(i);
    double *v2 = originalFieldArray->GetTuple(i);

    double absError = abs(*v1 - *v2);
    if(absError>0.5){
    std::cout << "id " << i 
              << " interp " << *v1 << " original " << *v2 << std::endl;
    }
    accuError = accuError + absError;
    
    maxError = std::max (maxError, abs(*v1 - *v2));
  }

  std::cout << "accumulated error " << accuError << std::endl;
  std::cout << "max error " << maxError << std::endl;
  std::cout << "avg error " << accuError/interpFieldArray->GetNumberOfTuples() << std::endl;
}

int main(int argc, char *argv[])
{
  // Parse command line arguments
  if (argc != 3)
  {
    std::cerr << "Usage: " << argv[0]
              << " <dirname+suffix> <fieldName>" << std::endl;
    return EXIT_FAILURE;
  }

  std::string datasetDirSuffix = argv[1];

  std::string fieldName = argv[2];

  std::string unstructuredMeshRaw = datasetDirSuffix + ".vtk";
  std::string unstructuredWithVarFile = datasetDirSuffix + "WithVar.vtk";
  std::string structuredReampleFile = datasetDirSuffix + "WithVarResample.vtk";

  std::cout << "load files: " << unstructuredMeshRaw << "," << unstructuredWithVarFile << "," << structuredReampleFile << std::endl;

  // load the uniform grid
  vtkSmartPointer<vtkStructuredPointsReader> imgReader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();
  imgReader->SetFileName(structuredReampleFile.c_str());
  imgReader->Update();

  vtkStructuredPoints *struPoints = imgReader->GetOutput();

  // std::cout << "---check loaded struPoints---" << std::endl;
  // struPoints->Print(std::cout);

  vtkPointData *pointData = struPoints->GetPointData();
  auto pointDataArray = pointData->GetScalars(fieldName.c_str());
  // pointDataArray->Print(std::cout);

  // load the unstructured grid
  vtkSmartPointer<vtkUnstructuredGridReader> reader =
      vtkSmartPointer<vtkUnstructuredGridReader>::New();
  reader->SetFileName(unstructuredMeshRaw.c_str());
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
  auto fieldArray = struPoints->GetPointData()->GetArray(fieldName.c_str());
  DEBUG("check field");
  DEBUG(fieldArray->GetNumberOfTuples());
  DEBUG(fieldArray->GetNumberOfComponents());

  vtkNew<vtkDoubleArray> interpFieldArray;

  std::string interpFieldName = fieldName + "_interp";
  interpFieldArray->SetName(interpFieldName.c_str());
  interpFieldArray->SetNumberOfComponents(1);
  vtkIdType pointsNum = unsGridData->GetNumberOfPoints();
  interpFieldArray->SetNumberOfTuples(pointsNum);

  for (vtkIdType id = 0; id < pointsNum; id++)
  {

    unsGridData->GetPoint(id, pointcoord);

    // use this pointcoord to find the specific cell in the structured grid
    vtkIdType interpCellId = cellLocator->FindCell(pointcoord);
    vtkCell *interpCell = struPoints->GetCell(interpCellId);

    int interpPointsNum = interpCell->GetNumberOfPoints();

    if (interpCellId == -1)
    {
      // do not find the cell id
      // use the closetPoint to find
      double closestPoint[3];
      vtkIdType closestCellId;
      vtkGenericCell *closestCell;
      int subid = 0;
      double dist = 0;
      cellLocator->FindClosestPoint(pointcoord, closestPoint, closestCellId, subid, dist);
      interpCellId = closestCellId;
      interpCell = struPoints->GetCell(interpCellId);
      DEBUG("original is -1 new id is " << closestCellId);
    }

    DEBUG("id " << id << " x " << pointcoord[0] << " y " << pointcoord[1]
                << " z " << pointcoord[2] << " interpCellId is " << interpCellId
                << " interpPointsNum is " << interpPointsNum);

    // get the cell from the cellid from struPoints and do the interpolation
    vtkPoints *pointsArray = interpCell->GetPoints();
    int pointNum = pointsArray->GetNumberOfPoints();
    DEBUG(pointNum);
    vtkIdList *pointsIds = interpCell->GetPointIds();

    // std::cout << *(pointsArray->GetPoint(i)) << ",";
    std::vector<std::vector<double>> inputCoor;
    inputCoor.clear();
    std::vector<double> fieldValues;
    fieldValues.clear();
    for (int i = 0; i < pointNum; i++)
    {
      double strucCoors[3];
      pointsArray->GetPoint(i, strucCoors);

      DEBUG("point id " << pointsIds->GetId(i) << " coord " << strucCoors[0] << " " << strucCoors[1] << " " << strucCoors[2]);
      //  what are the value of assocaited points?
      //  get the field value for specific points
      double *fieldValue = pointDataArray->GetTuple(pointsIds->GetId(i));
      DEBUG("check data " << *fieldValue);
      //  use the strucCoors to interpolate pointcoord
      inputCoor.push_back({strucCoors[0], strucCoors[1], strucCoors[2]});
      fieldValues.push_back(*fieldValue);
    }

    double interPolatedValue = Interpolate2D(inputCoor, fieldValues, pointcoord);
    interpFieldArray->SetValue(id, interPolatedValue);
  }

  auto unsDataset = unsGridData->GetPointData();
  unsDataset->AddArray(interpFieldArray);
  std::cout << "add field:" << interpFieldName << std::endl;

  // output the data into the vtk file
  vtkSmartPointer<vtkUnstructuredGridWriter> writer =
      vtkSmartPointer<vtkUnstructuredGridWriter>::New();
  std::string fileSuffix = unstructuredMeshRaw.substr(0, unstructuredMeshRaw.size() - 4);
  std::string outputFileName = fileSuffix + std::string("Interp.vtk");
  std::cout << "create vtk file: " << outputFileName << std::endl;
  writer->SetFileName(outputFileName.c_str());
  // get the specific polydata and check the results
  writer->SetInputData(unsGridData);
  writer->Write();

  // Load the data with var and then extract associated file (value with var)
  // and then compare it with the interpolated array
  // maybe udpate the code here, just provide the name of the original mesh surfix

  vtkSmartPointer<vtkUnstructuredGridReader> readerWithVar =
      vtkSmartPointer<vtkUnstructuredGridReader>::New();
  readerWithVar->SetFileName(unstructuredWithVarFile.c_str());
  readerWithVar->Update();

  // get the specific polydata and check the results
  vtkUnstructuredGrid *unsGridDataWithVar = readerWithVar->GetOutput();

  vtkPointData *pointDataWithVar = unsGridDataWithVar->GetPointData();
  auto origianlDataArray = pointDataWithVar->GetScalars(fieldName.c_str());
  //TODO, why there are large error aound the edge of the hole? still need to investigate
  computeDiff(interpFieldArray, origianlDataArray);
}