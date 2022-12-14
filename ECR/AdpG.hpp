#ifndef ADPG_h
#define ADPG_h

#include <vtkOctreePointLocatorNode.h>
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vector>
#include <array>
#include <vtkRectilinearGrid.h>

class ADPG
{
public:
    ADPG(vtkDataSet* dataset, int maxDepth, int vertexNum);
    ADPG(vtkDataSet* dataset);
    ~ADPG();
    void BuildAdpGrid();
    void DivideRegion(vtkOctreePointLocatorNode *node, int* subRegionIds, int depth);
    //TODO
    //ConvertToUniformGrid();

    vtkSmartPointer<vtkUnstructuredGrid> ConvertToUnstructuredGrid();
    //ECBuildAdpGrid() build the grid with error controlled
    vtkSmartPointer<vtkRectilinearGrid> ConvertToRectlinearGrid();


private:

    void traverse(vtkOctreePointLocatorNode *node, vtkUnstructuredGrid*ugrid, std::vector<std::array<double,3>> &pointsCoord, int depth);

    // the condition to control when to stop the division
    int m_maxDepth;
    int m_vertexNumThreshold;

    float* m_subRegionPoints;
    int* m_subRegionIds;


    vtkDataSet* m_dataset = nullptr;
    
    vtkOctreePointLocatorNode *m_Top=nullptr;
    vtkOctreePointLocatorNode **m_LeafNodeList=nullptr; // indexed by region/node ID
};

#endif