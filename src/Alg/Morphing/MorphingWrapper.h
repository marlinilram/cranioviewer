#ifndef MorphingWrapper_H
#define MorphingWrapper_H

#include "Morphing.h"
#include "Mesh.h"

class MorphingWrapper
{
public:
    MorphingWrapper();
    ~MorphingWrapper();

    void loadMesh(std::string fName, vtkSmartPointer<vtkRenderer> renderer = nullptr);

private:
    std::vector<Mesh> meshes;
    Morphing *morphing_hander;
};

#endif