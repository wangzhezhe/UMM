

#include "AdpG.hpp"

ADPG::ADPG(vtkDataSet *dataset, int maxDepth, int vertexNum)
{
  this->m_Top = nullptr;
  this->m_LeafNodeList = nullptr;
  this->m_dataset = dataset;
  this->m_maxDepth = maxDepth;
  this->m_vertexNumThreshold = vertexNum;
}

ADPG::ADPG(vtkDataSet *dataset)
{
  this->m_Top = nullptr;
  this->m_LeafNodeList = nullptr;
  this->m_dataset = dataset;
  this->m_maxDepth = 10;
  this->m_vertexNumThreshold = 10;
}

void ADPG::DivideRegion(vtkOctreePointLocatorNode *node, int *regionPointIds, int depth)
{

  if (node == nullptr)
  {
    return;
  }

  // get to the maximal depth
  // std::cout << "debug " << depth << std::endl;
  if (depth > this->m_maxDepth)
  {
    return;
  }

  // point number in region is smaller then threshold
  if (node->GetNumberOfPoints() < this->m_vertexNumThreshold)
  {
    return;
  }

  // did not hit the return condition, continue to refine things

  node->CreateChildNodes();

  int numberOfPoints = node->GetNumberOfPoints();

  std::cout << "debug DivideRegion " << numberOfPoints << std::endl;

  double bounds[6];
  // TODO there are some issues for getting bounds here
  node->GetBounds(bounds);
  std::cout << "debug DivideRegion bounds ";

  for (int i = 0; i < 6; i++)
  {
    std::cout << bounds[i] << " ";
  }
  std::cout << std::endl;

  vtkDataSet *ds = this->m_dataset;

  // ds->Print(std::cout);

  // there are 7 sub domains
  // each subdomain contains one list
  std::vector<int> points[7];
  int i;
  // number of points in each subdivision
  int subOctantNumberOfPoints[8] = {0, 0, 0, 0, 0, 0, 0, 0};

  // go through all points in current region
  for (i = 0; i < numberOfPoints; i++)
  {
    // check the current point, and to see it is located in which quadrant
    double *temp = ds->GetPoint(regionPointIds[i]);
    // std::cout << "debug point:" << regionPointIds[i] << " " << temp[0] << " " << temp[1] << " " << temp[2] << std::endl;
    int index = node->GetSubOctantIndex(ds->GetPoint(regionPointIds[i]), 0);
    // std::cout << "index " << index << std::endl;
    if (index)
    {
      points[index - 1].push_back(regionPointIds[i]);
    }
    else if (index == 0)
    {
      // the first one?
      regionPointIds[subOctantNumberOfPoints[0]] = regionPointIds[i];
    }
    else
    {
      // index is a negative value
      throw std::runtime_error("index is negative");
    }
    subOctantNumberOfPoints[index]++;
  }

  int counter = 0;
  int sizeOfInt = sizeof(int);
  for (i = 0; i < 8; i++)
  {
    // how many points for each subregion
    counter += subOctantNumberOfPoints[i];
    if (!points[i].empty())
    {
      // put the list into the subregion
      memcpy(regionPointIds + counter, points[i].data(), subOctantNumberOfPoints[i + 1] * sizeOfInt);
    }
  }

  counter = 0;
  for (i = 0; i < 8; i++)
  {
    // go through each subregion
    if (subOctantNumberOfPoints[i] == 0)
    {
      continue;
    }
    node->GetChild(i)->SetNumberOfPoints(subOctantNumberOfPoints[i]);
    std::cout << "depth " << depth << " debug point number " << i << " " << subOctantNumberOfPoints[i] << std::endl;
    double childBounds[6];
    node->GetChild(i)->GetBounds(childBounds);
    std::cout << "depth " << depth << " bounds ";
    for (int j = 0; j < 6; j++)
    {
      std::cout << childBounds[j] << " ";
    }
    std::cout << std::endl;
    this->DivideRegion(node->GetChild(i), regionPointIds + counter, depth + 1);
    counter += subOctantNumberOfPoints[i];
  }
}

void ADPG::BuildAdpGrid()
{

  // In VTK, there are some boundry processing operations, some tolerance thing?

  vtkOctreePointLocatorNode *node = this->m_Top = vtkOctreePointLocatorNode::New();

  double bounds[6];
  this->m_dataset->GetBounds(bounds);
  int numPoints = this->m_dataset->GetNumberOfPoints();

  node->SetNumberOfPoints(numPoints);

  for (int i = 0; i < 6; i++)
  {
    std::cout << "debug bound " << i << " " << bounds[i] << std::endl;
  }

  // SetBounds will set the location for the data
  // SetDataBounds will set the data contained by the point
  node->SetBounds(bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

  this->m_subRegionIds = new int[numPoints];
  this->m_subRegionPoints = new float[3 * numPoints];

  if (this->m_subRegionIds == nullptr || this->m_subRegionPoints == nullptr)
  {
    throw std::runtime_error("failed to allocated locatorids");
  }

  for (int i = 0; i < numPoints; i++)
  {
    this->m_subRegionIds[i] = i;
  }

  this->DivideRegion(node, this->m_subRegionIds, 0);
}

void ADPG::traverse(vtkOctreePointLocatorNode *node, vtkUnstructuredGrid *ugrid, std::vector<std::array<double, 3>> &pointsCoord, int depth)
{
  if (node == nullptr)
  {
    return;
  }

  // get the bounudry of the cell
  double bounds[6];
  node->GetBounds(bounds);

  // add new points into list
  int sizeOffset = pointsCoord.size();

  // only test 2d case currently

  vtkIdType pid0 = sizeOffset + 0;
  vtkIdType pid1 = sizeOffset + 1;
  vtkIdType pid2 = sizeOffset + 2;
  vtkIdType pid3 = sizeOffset + 3;

  // make sure there is no crossing line
  vtkIdType cellIds[4] = {pid0, pid1, pid2, pid3};

  // std::cout << "insert cell id:" << pid0 << " " << pid1 << " " << pid2 << " " << pid3 << std::endl;
  ugrid->InsertNextCell(VTK_QUAD, 4, cellIds);

  // create the points coordinates
  // create new points id and cell

  // insert into the grid
  // refer to https://kitware.github.io/vtk-examples/site/Cxx/UnstructuredGrid/UGrid/
  // only test 2d case here

  // std::cout << "debug traverse " << bounds[0] << " " << bounds[1] << std::endl;
  double xmin = bounds[0];
  double xmax = bounds[1];
  double ymin = bounds[2];
  double ymax = bounds[3];
  double zmin = bounds[4];
  double zmax = bounds[5];

  std::cout << "traverse depth " << depth << " " << xmin << " " << xmax << " " << ymin << " " <<ymax << std::endl; 

  // TODO change this for the 3d case
  pointsCoord.push_back({xmin, ymin, 0});
  pointsCoord.push_back({xmax, ymin, 0});
  pointsCoord.push_back({xmax, ymax, 0});
  pointsCoord.push_back({xmin, ymax, 0});

  // get children
  for (int i = 0; i < 8; i++)
  {
    // go through each subregion
    if(node->GetChild(i)==nullptr){
      continue;
    }
    if(node->GetChild(i)->GetNumberOfPoints()==0){
      continue;
    }
    traverse(node->GetChild(i), ugrid, pointsCoord, depth+1);
  }

  // call the reverse for childrena
}


//refer to https://kitware.github.io/vtk-examples/site/Cxx/RectilinearGrid/RGrid/
vtkSmartPointer<vtkRectilinearGrid> ADPG::ConvertToRectlinearGrid(){
// go through the adp and conver to rectlinear grid

}

vtkSmartPointer<vtkUnstructuredGrid> ADPG::ConvertToUnstructuredGrid()
{
  if (this->m_Top == nullptr)
  {
    throw std::runtime_error("the m_Top is not supposed to be nullptr");
  }

  vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  std::vector<std::array<double, 3>> pointsCoord;

  // traverse the whole tree
  traverse(this->m_Top, ugrid, pointsCoord, 0);

  // add one cell when we access one node
  vtkNew<vtkPoints> points;
  for (int i = 0; i < pointsCoord.size(); i++)
  {
    double temp[3] = {pointsCoord[i][0], pointsCoord[i][1], pointsCoord[i][2]};
    std::cout << "insert point id " << i << " " << temp[0] << " " << temp[1] << " " << temp[2] << std::endl;
    points->InsertPoint(i, temp);
  }
  ugrid->SetPoints(points);

  return ugrid;
}

ADPG::~ADPG()
{
  delete[] this->m_LeafNodeList;
  this->m_LeafNodeList = nullptr;
}
