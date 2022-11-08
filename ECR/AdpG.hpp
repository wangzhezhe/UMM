#ifndef ADPG_h
#define ADPG_h

#include <vtkOctreePointLocatorNode.h>
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>

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

    vtkUnstructuredGrid* ConvertToUnstructuredGrid();
    //ECBuildAdpGrid() build the grid with error controlled

private:

    void traverse(vtkOctreePointLocatorNode *node, vtkUnstructuredGrid*ugrid, std::vector<double[3]> &pointsCoord);

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