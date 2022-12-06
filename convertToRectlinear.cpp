#include <vector>
#include <stdio.h>
#include <math.h>
#include <vtkUnstructuredGridReader.h>


void addCoord(std::vector<double> &coords, double p)
{
    double tol = 0.00001;
    for (int i = 0; i < coords.size(); i++)
    {
        if (fabs(coords[i] - p) < tol)
        {
            // exist
            return;
        }
    }

    // did not find
    // push p into the list
    coords.push_back(p);
}

int main(int argc, char *argv[])
{


  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " Filename" << std::endl;
    return EXIT_FAILURE;
  }

  std::string vtkFile = std::string(argv[1]);

    // input a unstructed data
      vtkSmartPointer<vtkUnstructuredGridReader> reader =
      vtkSmartPointer<vtkUnstructuredGridReader>::New();
  reader->SetFileName(vtkFile.c_str());
  reader->Update();

  // get the specific unstructureGridData and check the results
  vtkSmartPointer<vtkUnstructuredGrid> unstructureGridData =
      reader->GetOutput();

    // output the rectlinear grid based on the unstructured data
    // the original point can be kept, so there is zero error
    // we also added some new points accidentally

    // refer to https://kitware.github.io/vtk-examples/site/Cxx/RectilinearGrid/RGrid/

    std::vector<double> xcoords;
    std::vector<double> ycoords;
    std::vector<double> zcoords;

    // load the unstructred data

    // go throught the unstructured data

    // for each point
    // call the xcood ycood and zcood separately to test the existance
    // for x y and z dim


    //when we get the rectangle grid

    //then call the resample

    //instead of calling the resample direactly
}