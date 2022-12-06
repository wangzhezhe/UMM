
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
#include <vtkNew.h>
#include <vtkCell.h>
#include <vtkImageData.h>
#include <vtkDataSetWriter.h>
#include <vtkImageWriter.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>

void createSampleImage()
{
    std::string fileNameSuffix = "sample_image";
    int dimx = 100;
    int dimy = 100;
    int dimz = 100;
    // Specify the size of the image data
    vtkNew<vtkStructuredPoints> imageData;
    imageData->SetDimensions(dimx, dimy, dimz);
    imageData->SetSpacing(1.0, 1.0, 1.0);
    imageData->SetOrigin(0.0, 0.0, 0.0);

    std::cout << "Number of points: " << imageData->GetNumberOfPoints()
              << std::endl;
    std::cout << "Number of cells: " << imageData->GetNumberOfCells()
              << std::endl;

    char temp[128];
    sprintf(temp, "_%d_%d_%d", dimx, dimy, dimz);
    std::string fileName = fileNameSuffix + std::string(temp) + ".vtk";

    vtkSmartPointer<vtkStructuredPointsWriter> writer =
        vtkSmartPointer<vtkStructuredPointsWriter>::New();
    std::cout << "create vtk file: " << fileName << std::endl;
    writer->SetFileName(fileName.c_str());
    // get the specific polydata and check the results
    // visit only support the 42 version
    // otherwise, the data can not be loaded properly
    writer->SetFileVersion(42);
    writer->SetInputData(imageData);
    writer->Write();
}

int main(int argc, char *argv[])
{

    createSampleImage();
    return 0;
}