
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGridWriter.h>

struct Point {
  std::vector<double> m_coor;
  double m_v1;
};

double randNum(double rmin, double rmax) {
  double f = (double)rand() / RAND_MAX;
  return rmin + f * (rmax - rmin);
};

// generates a random point list in a region
// refer to this
// https://kitware.github.io/vtk-examples/site/Cxx/UnstructuredGrid/UGrid/ refer
// to this
// https://discourse.vtk.org/t/unstructured-grid-from-randomly-disributed-points/7222

void generatePointLists() {

  // create bounding box and put points here
  int pointNum = 100;
  // for 2d cases
  std::vector<double> lb = {0, 0};
  std::vector<double> ub = {5, 5};

  std::vector<Point> pointLists(pointNum);

  for (int i = 0; i < pointNum; i++) {
    // create the coordinates randomly
    double x = randNum(lb[0], ub[0]);
    double y = randNum(lb[1], ub[1]);
    pointLists[i].m_coor.push_back(x);
    pointLists[i].m_coor.push_back(y);
  }

  // output to xyz format
  std::ofstream pointsFile;
  pointsFile.open("pointsLists.xyz");
  pointsFile << pointNum << std::endl;
  pointsFile << "x y z" << std::endl;

  for (int i = 0; i < pointNum; i++) {
    pointsFile << i << " " << pointLists[i].m_coor[0] << " "
               << pointLists[i].m_coor[1] << " " << 0 << std::endl;
  }

  pointsFile.close();
}
//refer to https://stackoverflow.com/questions/10847007/using-the-gaussian-probability-density-function-in-c
double normalPdf(double x, double m, double s)
{
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

// input an plan unstructured mesh
// add varaible to the mesh
// and output the data
void addVariables(std::string vtkFile) {
  // create the unstructured instance
  //  Read all the data from the file
  vtkSmartPointer<vtkUnstructuredGridReader> reader =
      vtkSmartPointer<vtkUnstructuredGridReader>::New();
  reader->SetFileName(vtkFile.c_str());
  reader->Update();

  // get the specific unstructureGridData and check the results
  vtkSmartPointer<vtkUnstructuredGrid> unstructureGridData =
      reader->GetOutput();

  // TODO, go through each points and adding variable here
  // the variable is the distance between the center
  double testpoint[3];
  const double *bounds = unstructureGridData->GetBounds();
  double xmin = bounds[0];
  double xmax = bounds[1];
  double ymin = bounds[2];
  double ymax = bounds[3];
  double zmin = bounds[4];
  double zmax = bounds[5];

  double center[3] = {(xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2};

  std::cout << "center " << center[0] << "," << center[1] << "," << center[2]
            << std::endl;

  vtkNew<vtkDoubleArray> VarArray;
  VarArray->SetName("v_center_dist");
  VarArray->SetNumberOfComponents(1);

  for (vtkIdType id = 0; id < unstructureGridData->GetNumberOfPoints(); id++) {

    unstructureGridData->GetPoint(id, testpoint);
    double dist = sqrt(pow((testpoint[0] - center[0]), 2) +
                  pow((testpoint[1] - center[1]), 2) +
                  pow((testpoint[2] - center[2]), 2));
    double normalPdfValue = normalPdf(dist, 0.0, 0.15);
    //the array bonds with the point based on the id
    if(dist<0.5){
      normalPdfValue=normalPdfValue*5;
    }
    VarArray->InsertNextValue(normalPdfValue);

    // std::cout << "Point: " << testpoint[0] << ", " << testpoint[1] << ", "
    //           << testpoint[2] << " center dist " << dist << endl;
    //  Add variable into the data set
  }

  auto dataset = unstructureGridData->GetPointData();
  dataset->AddArray(VarArray);

  unstructureGridData->Print(std::cout);

  //output
    vtkSmartPointer<vtkUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    std::string fileSuffix = vtkFile.substr(0, vtkFile.size() - 4);
    std::string outputFileName = fileSuffix+std::string("WithVar.vtk");
    writer->SetFileName(outputFileName.c_str());
    // get the specific polydata and check the results
    writer->SetInputData(unstructureGridData);
    writer->Write();
}

int main(int argc, char *argv[]) {
  // Parse command line arguments
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " Filename" << std::endl;
    return EXIT_FAILURE;
  }

  addVariables(argv[1]);
  return 0;
}