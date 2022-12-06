
#include <cmath>
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkImageReader.h>
#include <vtkDataSetWriter.h>
#include <vtkImageData.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsWriter.h>

struct Block
{
  Block() = default;
  Block(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
  {
    bounds[0] = xmin;
    bounds[1] = xmax;
    bounds[2] = ymin;
    bounds[3] = ymax;
    bounds[4] = zmin;
    bounds[5] = zmax;
  };
  double bounds[6] = {0, 0, 0, 0, 0, 0};
  void printBound()
  {
    for (int i = 0; i < 6; i++)
    {
      std::cout << bounds[i] << ",";
    }
    std::cout << std::endl;
  };

  bool contain(double x, double y, double z)
  {
    if (x >= bounds[0] && x <= bounds[1] && y >= bounds[2] && y <= bounds[3] && z >= bounds[4] && z <= bounds[5])
    {
      return true;
    }
    return false;
  }
  ~Block() = default;
};

std::vector<Block> buildBlockList(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
  int blockNumPerDim = 4;
  std::vector<Block> blockList;

  double dx = (xmax - xmin) / blockNumPerDim;
  double dy = (ymax - ymin) / blockNumPerDim;
  double dz = (zmax - zmin) / blockNumPerDim;

  for (int i = 0; i < blockNumPerDim; i++)
  {
    for (int j = 0; j < blockNumPerDim; j++)
    {
      for (int k = 0; k < blockNumPerDim; k++)
      {
        // i j k is the index of x y z
        double bxmin = xmin + i * dx;
        double bxmax = xmin + (i + 1) * dx;
        double bymin = ymin + j * dy;
        double bymax = ymin + (j + 1) * dy;
        double bzmin = zmin + k * dz;
        double bzmax = zmin + (k + 1) * dz;
        Block b(bxmin, bxmax, bymin, bymax, bzmin, bzmax);
        //b.printBound();
        blockList.push_back(b);
      }
    }
  }

  return blockList;
}

Block findBlock(double x, double y, double z, std::vector<Block> &blockList)
{
  // if the point is on the boundry
  // return the block that first contains this point
  for (int i = 0; i < blockList.size(); i++)
  {
    Block b = blockList[i];
    if (b.contain(x, y, z))
    {
      return b;
    }
  }
  // empty block
  return Block();
};

// refer to
// https://stackoverflow.com/questions/38937139/how-to-store-a-vector-field-with-vtk-c-vtkwriter
void addVariableVector(vtkSmartPointer<vtkUnstructuredGrid> inputDataset)
{
  // the variable is the distance between the center
  double testpoint[3];
  const double *bounds = inputDataset->GetBounds();
  double xmin = bounds[0];
  double xmax = bounds[1];
  double ymin = bounds[2];
  double ymax = bounds[3];
  double zmin = bounds[4];
  double zmax = bounds[5];

  double center[3] = {(xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2};

  double xydiameter = sqrt(pow((xmin - center[0]) / 2, 2) +
                           pow((ymin - center[1]) / 2, 2));

  vtkNew<vtkDoubleArray> VarArray;
  VarArray->SetName("velocity");
  VarArray->SetNumberOfComponents(3);
  for (vtkIdType id = 0; id < inputDataset->GetNumberOfPoints(); id++)
  {

    inputDataset->GetPoint(id, testpoint);
    double xydist = sqrt(pow((testpoint[0] - center[0]), 2) +
                         pow((testpoint[1] - center[1]), 2));

    if (fabs(xydist - xydiameter) > 0.2 * xydiameter)
    {
      // insert vector
      // testpoint to center, then another vector
      double dx = testpoint[0] - center[0];
      double dy = testpoint[1] - center[1];
      VarArray->InsertNextTuple3(-dy, dx, 0);
    }
    else
    {
      VarArray->InsertNextTuple3(0, 0, 0);
    }

    // std::cout << "Point: " << testpoint[0] << ", " << testpoint[1] << ", "
    //           << testpoint[2] << " center dist " << dist << endl;
    //  Add variable into the data set
  }

  vtkPointData *pointDataset = inputDataset->GetPointData();
  pointDataset->AddArray(VarArray);
}

void addVariableVector(vtkSmartPointer<vtkStructuredPoints> inputDataset)
{
  // the variable is the distance between the center
  double testpoint[3];
  const double *bounds = inputDataset->GetBounds();
  double xmin = bounds[0];
  double xmax = bounds[1];
  double ymin = bounds[2];
  double ymax = bounds[3];
  double zmin = bounds[4];
  double zmax = bounds[5];

  for (int i = 0; i < 6; i++)
  {
    std::cout << bounds[i] << " ";
  }
  std::cout << std::endl;

  double center[3] = {(xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2};

  double xydiameter = sqrt(pow((xmin - center[0]), 2) +
                           pow((ymin - center[1]), 2));

  vtkNew<vtkDoubleArray> VarArray;
  VarArray->SetName("velocity");
  VarArray->SetNumberOfComponents(3);
  for (vtkIdType id = 0; id < inputDataset->GetNumberOfPoints(); id++)
  {

    inputDataset->GetPoint(id, testpoint);
    double xydist = sqrt(pow((testpoint[0] - center[0]), 2) +
                         pow((testpoint[1] - center[1]), 2));

    if (xydist > 0.1 * xydiameter)
    {
      // insert vector
      // testpoint to center, then another vector
      double dx = testpoint[0] - center[0];
      double dy = testpoint[1] - center[1];
      VarArray->InsertNextTuple3(dy, -dx, 0);
    }
    else
    {
      VarArray->InsertNextTuple3(0, 0, 0);
    }

    // std::cout << "Point: " << testpoint[0] << ", " << testpoint[1] << ", "
    //           << testpoint[2] << " center dist " << dist << endl;
    //  Add variable into the data set
  }

  vtkPointData *pointDataset = inputDataset->GetPointData();
  pointDataset->AddArray(VarArray);
}

void addVariableVectorUnEven(vtkSmartPointer<vtkStructuredPoints> inputDataset)
{
  // the variable is the distance between the center
  double testpoint[3];
  const double *bounds = inputDataset->GetBounds();
  double xmin = bounds[0];
  double xmax = bounds[1];
  double ymin = bounds[2];
  double ymax = bounds[3];
  double zmin = bounds[4];
  double zmax = bounds[5];

  for (int i = 0; i < 6; i++)
  {
    std::cout << bounds[i] << " ";
  }
  std::cout << std::endl;

  double center[3] = {(xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2};

  double diameter = sqrt(pow((xmin - center[0]) / 2, 2) +
                         pow((ymin - center[1]) / 2, 2));

  vtkNew<vtkDoubleArray> VarArray;
  VarArray->SetName("velocity");
  VarArray->SetNumberOfComponents(3);
  for (vtkIdType id = 0; id < inputDataset->GetNumberOfPoints(); id++)
  {

    inputDataset->GetPoint(id, testpoint);
    double dist = sqrt(pow((testpoint[0] - center[0]), 2) +
                       pow((testpoint[1] - center[1]), 2) +
                       pow((testpoint[2] - center[2]), 2));

    if ((testpoint[2] - center[2]) > 0.5 * diameter)
    {
      // insert vector
      // testpoint to center, then another vector
      double dx = testpoint[0] - center[0];
      double dy = testpoint[1] - center[1];
      VarArray->InsertNextTuple3(dy, dx, 0);
    }
    else
    {
      VarArray->InsertNextTuple3(0, 0, 1);
    }

    // std::cout << "Point: " << testpoint[0] << ", " << testpoint[1] << ", "
    //           << testpoint[2] << " center dist " << dist << endl;
    //  Add variable into the data set
  }

  vtkPointData *pointDataset = inputDataset->GetPointData();
  pointDataset->AddArray(VarArray);
}

void addVectorUnstructured(std::string vtkFile)
{

  // create the unstructured instance
  // Read all the data from the file
  vtkSmartPointer<vtkUnstructuredGridReader> reader =
      vtkSmartPointer<vtkUnstructuredGridReader>::New();
  reader->SetFileName(vtkFile.c_str());
  reader->Update();

  // get the specific unstructureGridData and check the results
  vtkSmartPointer<vtkUnstructuredGrid> unstructureGridData =
      reader->GetOutput();

  // add variables
  // addVariableDistCenter(unstructureGridData);
  // addVariableSin(unstructureGridData);
  addVariableVector(unstructureGridData);
  // output

  // unstructureGridData->Print(std::cout);

  vtkSmartPointer<vtkUnstructuredGridWriter> writer =
      vtkSmartPointer<vtkUnstructuredGridWriter>::New();
  std::string fileSuffix = vtkFile.substr(0, vtkFile.size() - 4);
  std::string outputFileName = fileSuffix + std::string("WithVarVector.vtk");
  std::cout << "create vtk file: " << outputFileName << std::endl;
  writer->SetFileName(outputFileName.c_str());
  // get the specific polydata and check the results
  // visit only support the 42 version
  // otherwise, the data can not be loaded properly
  writer->SetFileVersion(42);
  writer->SetInputData(unstructureGridData);
  writer->Write();
}

void addVectorImage(std::string vtkFile)
{
  // create the unstructured instance
  // Read all the data from the file
  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(vtkFile.c_str());
  reader->Update();

  // get the specific unstructureGridData and check the results
  vtkSmartPointer<vtkStructuredPoints> imageData = reader->GetOutput();

  std::cout << "Number of points: " << imageData->GetNumberOfPoints()
            << std::endl;
  std::cout << "Number of cells: " << imageData->GetNumberOfCells()
            << std::endl;

  addVariableVector(imageData);
  // addVariableVectorUnEven(imageData);

  vtkSmartPointer<vtkStructuredPointsWriter>
      writer =
          vtkSmartPointer<vtkStructuredPointsWriter>::New();
  std::string fileSuffix = vtkFile.substr(0, vtkFile.size() - 4);
  std::string outputFileName = fileSuffix + std::string("WithVarVector.vtk");
  std::cout << "create vtk file: " << outputFileName << std::endl;
  writer->SetFileName(outputFileName.c_str());
  // get the specific polydata and check the results
  // visit only support the 42 version
  // otherwise, the data can not be loaded properly
  writer->SetFileVersion(42);
  writer->SetInputData(imageData);
  writer->Write();
}

void addVectorImageComplex(std::string vtkFile)
{
  // create the unstructured instance
  // Read all the data from the file
  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();
  reader->SetFileName(vtkFile.c_str());
  reader->Update();

  // get the specific unstructureGridData and check the results
  vtkSmartPointer<vtkStructuredPoints> imageData = reader->GetOutput();

  std::cout << "Number of points: " << imageData->GetNumberOfPoints()
            << std::endl;
  std::cout << "Number of cells: " << imageData->GetNumberOfCells()
            << std::endl;

  // the variable is the distance between the center
  double testpoint[3];
  const double *bounds = imageData->GetBounds();
  double xmin = bounds[0];
  double xmax = bounds[1];
  double ymin = bounds[2];
  double ymax = bounds[3];
  double zmin = bounds[4];
  double zmax = bounds[5];

  for (int i = 0; i < 6; i++)
  {
    std::cout << bounds[i] << " ";
  }
  std::cout << std::endl;

  // create the list of the centers/blocks
  std::vector<Block> blocklist =
      buildBlockList(xmin, xmax, ymin, ymax, zmin, zmax);

  // go through each point, add the vector according to their positions

  vtkNew<vtkDoubleArray> VarArray;
  VarArray->SetName("velocity");
  VarArray->SetNumberOfComponents(3);
  for (vtkIdType id = 0; id < imageData->GetNumberOfPoints(); id++)
  {

    imageData->GetPoint(id, testpoint);
    Block b = findBlock(testpoint[0], testpoint[1], testpoint[2], blocklist);
    if (fabs(b.bounds[0] - 0) <= 0.00001 && fabs(b.bounds[1] - 0) <= 0.00001)
    {
      throw std::runtime_error("did not find bounds");
    }

    // find bounds
    double bcenterx = (b.bounds[0] + b.bounds[1]) / 2;
    double bcentery = (b.bounds[2] + b.bounds[3]) / 2;
    double diameter = (b.bounds[3] - b.bounds[2]) / 2;
    double xydistToCenter = sqrt(pow((testpoint[0] - bcenterx), 2) +
                                 pow((testpoint[1] - bcentery), 2));
    if (xydistToCenter > 0.1 * diameter)
    {
      // insert vector
      // testpoint to center, then another vector
      double dx = testpoint[0] - bcenterx;
      double dy = testpoint[1] - bcentery;
      VarArray->InsertNextTuple3(dy, -dx, 0);
    }
    else
    {
      VarArray->InsertNextTuple3(0, 0, 0);
    }

    // std::cout << "Point: " << testpoint[0] << ", " << testpoint[1] << ", "
    //           << testpoint[2] << " center dist " << dist << endl;
    //  Add variable into the data set
  }

  vtkPointData *pointDataset = imageData->GetPointData();
  pointDataset->AddArray(VarArray);

  // ok to add variables

  vtkSmartPointer<vtkStructuredPointsWriter>
      writer =
          vtkSmartPointer<vtkStructuredPointsWriter>::New();
  std::string fileSuffix = vtkFile.substr(0, vtkFile.size() - 4);
  std::string outputFileName = fileSuffix + std::string("WithComplexVector.vtk");
  std::cout << "create vtk file: " << outputFileName << std::endl;
  writer->SetFileName(outputFileName.c_str());
  // get the specific polydata and check the results
  // visit only support the 42 version
  // otherwise, the data can not be loaded properly
  writer->SetFileVersion(42);
  writer->SetInputData(imageData);
  writer->Write();
}

int main(int argc, char *argv[])
{
  // Parse command line arguments
  if (argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " Filename "
              << " type(unstructured/image) " << std::endl;
    return EXIT_FAILURE;
  }

  std::string dataType(argv[2]);
  if (dataType == "unstructured")
  {
    addVectorUnstructured(argv[1]);
  }
  else if (dataType == "image")
  {
    addVectorImage(argv[1]);
    addVectorImageComplex(argv[1]);
  }
  else
  {
    std::cout << "unsupported type" << std::endl;
  }

  return 0;
}
