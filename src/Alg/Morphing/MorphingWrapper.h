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
    Morphing *getMorphingHandler() { return morphing_handler; };
    void doMorphing(double *paras, int n_paras);
    void doLinearMorphing(double *paras, int n_paras);
    Mesh *getMeshPtr(size_t i_mesh) { return meshes[i_mesh]; };
    void setCenterMesh(size_t i_mesh);
    inline size_t getMeshSize() { return meshes.size(); };

private:
    std::vector<Mesh *> meshes;
    Morphing *morphing_handler;
};

#endif