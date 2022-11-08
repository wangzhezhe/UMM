

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

  vtkDataSet *ds = this->m_dataset;

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
    int index = node->GetSubOctantIndex(ds->GetPoint(regionPointIds[i]), 0);
    if (index)
    {
      points[index - 1].push_back(regionPointIds[i]);
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
  for (i = 0; i < 7; i++)
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
    node->GetChild(i)->SetNumberOfPoints(subOctantNumberOfPoints[i]);
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
  node->SetDataBounds(bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

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

void ADPG::traverse(vtkOctreePointLocatorNode *node, vtkUnstructuredGrid *ugrid, std::vector<double[3]> &pointsCoord)
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

  vtkIdType cellIds[4] = {pid0, pid1, pid2, pid3};

  ugrid->InsertNextCell(VTK_QUAD, 4, cellIds);

  // create the points coordinates
  // create new points id and cell

  // insert into the grid
  // refer to https://kitware.github.io/vtk-examples/site/Cxx/UnstructuredGrid/UGrid/
  // only test 2d case here
  pointsCoord.push_back({bounds[0], bounds[0], 0});
  pointsCoord.push_back({bounds[0], bounds[1], 0});
  pointsCoord.push_back({bounds[1], bounds[0], 0});
  pointsCoord.push_back({bounds[1], bounds[1], 0});

  // get children
  for (int i = 0; i < 8; i++)
  {
    // go through each subregion
    traverse(node->GetChild(i), ugrid, pointsCoord);
  }

  // call the reverse for childrena
}

vtkUnstructuredGrid *ADPG::ConvertToUnstructuredGrid()
{
  if (this->m_Top == nullptr)
  {
    throw std::runtime_error("the m_Top is not supposed to be nullptr");
  }

  vtkNew<vtkUnstructuredGrid> ugrid;

  std::vector<double[3]> pointsCoord;

  // traverse the whole tree
  traverse(this->m_Top, ugrid, pointsCoord);

  // add one cell when we access one node
  vtkNew<vtkPoints> points;
  for (int i = 0; i < pointsCoord.size(); i++)
  {
    points->InsertPoint(i, pointsCoord[i]);
  }

  ugrid->SetPoints(points);

  return ugrid;
}

ADPG::~ADPG()
{
  delete[] this->m_LeafNodeList;
  this->m_LeafNodeList = nullptr;
}
