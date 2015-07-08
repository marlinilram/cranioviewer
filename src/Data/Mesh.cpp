#include "Mesh.h"

#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkLookupTable.h>
#include <vtkPointData.h>
#include <vtkPointPicker.h>
#include <vtkWorldPointPicker.h>

Mesh::Mesh()
{
    mapper = nullptr;
    actor = nullptr;
    renderer = nullptr;
    widget_transform = nullptr;
    widget_transform_mat = nullptr;
    kdtree_picker = nullptr;
}

Mesh::~Mesh()
{
    if (actor)
    {
        renderer->RemoveActor(actor);
        std::cout<<"remove mesh actor\n";
    }
}

Mesh::Mesh(std::string fName)
{
    setMeshData(fName);
}

Mesh::Mesh(std::string fName, vtkSmartPointer<vtkRenderer> mainWin_renderer)
{
    setMeshData(fName);
    setRenderer(mainWin_renderer);
}

void Mesh::setMeshData(std::string fName)
{
    vtkSmartPointer<vtkOBJReader> obj_reader = vtkSmartPointer< vtkOBJReader>::New();

    obj_reader->SetFileName(fName.c_str());
    obj_reader->Update();

    vtkSmartPointer<vtkLookupTable> colorLut = vtkSmartPointer<vtkLookupTable>::New();
    colorLut->SetNumberOfTableValues(10);
    colorLut->Build();

    // Fill in a few known colors, the rest will be generated if needed
    colorLut->SetTableValue(0     , 0     , 0     , 0, 1);  //Black
    colorLut->SetTableValue(1, 0.8900, 0.8100, 0.3400, 1); // Banana
    colorLut->SetTableValue(2, 1.0000, 0.3882, 0.2784, 1); // Tomato
    colorLut->SetTableValue(3, 0.9608, 0.8706, 0.7020, 1); // Wheat
    colorLut->SetTableValue(4, 0.9020, 0.9020, 0.9804, 1); // Lavender
    colorLut->SetTableValue(5, 1.0000, 0.4900, 0.2500, 1); // Flesh
    colorLut->SetTableValue(6, 0.5300, 0.1500, 0.3400, 1); // Raspberry
    colorLut->SetTableValue(7, 0.9804, 0.5020, 0.4471, 1); // Salmon
    colorLut->SetTableValue(8, 0.7400, 0.9900, 0.7900, 1); // Mint
    colorLut->SetTableValue(9, 0.2000, 0.6300, 0.7900, 1); // Peacock

    colorLut->SetRange(0,10);

    vtkSmartPointer<vtkUnsignedCharArray> cellColorId = vtkSmartPointer<vtkUnsignedCharArray>::New();
    cellColorId->SetNumberOfComponents(1);
    cellColorId->SetNumberOfTuples(obj_reader->GetOutput()->GetNumberOfCells());
    for (int i = 0; i < cellColorId->GetNumberOfTuples(); ++i)
    {
      cellColorId->InsertTuple1(i, 1);
    }

    cellColorId->SetLookupTable(colorLut);
    obj_reader->GetOutput()->GetCellData()->SetScalars(cellColorId);

    vtkSmartPointer<vtkUnsignedCharArray> ptColorId = vtkSmartPointer<vtkUnsignedCharArray>::New();
    ptColorId->SetNumberOfComponents(1);
    ptColorId->SetNumberOfTuples(obj_reader->GetOutput()->GetNumberOfPoints());
    for (int i = 0; i < ptColorId->GetNumberOfTuples(); ++i)
    {
      ptColorId->InsertTuple1(i, 9);
    }

    ptColorId->SetLookupTable(colorLut);
    obj_reader->GetOutput()->GetPointData()->SetScalars(ptColorId);

    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(obj_reader->GetOutput());

    mapper->SetScalarModeToUseCellData();
    mapper->SetColorModeToMapScalars();
    mapper->UseLookupTableScalarRangeOn();
    mapper->ScalarVisibilityOn();

    actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    double bounds[6];
    actor->GetBounds(bounds);
    computeCenter(bounds);

    widget_transform = vtkSmartPointer<vtkTransform>::New();
    widget_transform->PostMultiply();

    widget_transform_mat = vtkSmartPointer<vtkMatrix4x4>::New();
    widget_transform_mat->Identity();
}

void Mesh::computeCenter(double bounds[6])
{
    mesh_center[0] = (bounds[0] + bounds[1]) / 2;
    mesh_center[1] = (bounds[2] + bounds[3]) / 2;
    mesh_center[2] = (bounds[4] + bounds[5]) / 2;
    std::cout<<mesh_center[0]<<"\t"<<mesh_center[1]<<"\t"<<mesh_center[2]<<"\n";
}

void Mesh::setRenderer(vtkSmartPointer<vtkRenderer> mainWin_renderer)
{
    if (mainWin_renderer)
    {
        renderer = mainWin_renderer;
        renderer->AddActor(actor);

        renderer->ResetCamera();
    }
}

void Mesh::updateTransform(double *tr_vals)
{
    widget_transform->Identity();
    widget_transform->Translate(-mesh_center[0], -mesh_center[1], -mesh_center[2]);
    widget_transform->RotateX(tr_vals[3]);
    widget_transform->RotateY(tr_vals[4]);
    widget_transform->RotateZ(tr_vals[5]);
    widget_transform->Scale(tr_vals[6], tr_vals[6], tr_vals[6]);
    widget_transform->Translate(mesh_center[0], mesh_center[1], mesh_center[2]);
    widget_transform->Translate(tr_vals[0], tr_vals[1], tr_vals[2]);
    widget_transform->GetMatrix(widget_transform_mat);

    actor->SetUserTransform(widget_transform);
}

void Mesh::reloadTransform()
{
    actor->SetUserTransform(widget_transform);

    //vtkSmartPointer<vtkMatrix4x4> mat = vtkSmartPointer<vtkMatrix4x4>::New();
    //widget_transform->GetMatrix(mat);
    //mat->Print(std::cout);
}

void Mesh::setVisible(int state)
{
  if (state)
  {
    actor->VisibilityOn();
  }
  else
  {
    actor->VisibilityOff();
  }
}

void Mesh::resetMesh(vtkSmartPointer<vtkPolyData> new_mesh_data)
{
    renderer->RemoveActor(actor);
    mapper->SetInputData(new_mesh_data);
    actor->SetMapper(mapper);

    renderer->AddActor(actor);
}

void Mesh::setMeshCenter(double new_center[3])
{
    vtkSmartPointer<vtkCenterOfMass> c_mass = vtkSmartPointer<vtkCenterOfMass>::New();

    c_mass->SetInputData(mapper->GetInput());
    c_mass->SetUseScalarsAsWeights(false);
    c_mass->Update();
    c_mass->GetCenter(mesh_center);
    //mesh_center[0] = new_center[0];
    //mesh_center[1] = new_center[1];
    //mesh_center[2] = new_center[2];
}

void Mesh::setKDTreeLocator()
{
    kdtree_picker = vtkSmartPointer<vtkKdTreePointLocator>::New();
    kdtree_picker->SetDataSet(mapper->GetInput());
    kdtree_picker->BuildLocator();
}

int Mesh::getClosestPtID(double world_coord[3])
{
    if (kdtree_picker)
    {
        int vertex_id = -1;
        vertex_id = kdtree_picker->FindClosestPoint(world_coord);
        double *vertex_coord;
        vertex_coord = mapper->GetInput()->GetPoint(vertex_id);

        return vertex_id;
    }

    return -1;
}

void Mesh::getWorldCoord(double x, double y, double coord[3])
{
  vtkSmartPointer<vtkWorldPointPicker> picker = vtkSmartPointer<vtkWorldPointPicker>::New();
  picker->AddPickList(actor);
  picker->PickFromListOn();
  picker->Pick(x, y, 0.0, renderer);


  double pos[3];
  picker->GetPickPosition(pos);
  //std::cout<<"picked position: "<<pos[0]<<"\t"<<pos[1]<<"\t"<<pos[2]<<"\n";

  if (coord)
  {
    coord[0] = pos[0];
    coord[1] = pos[1];
    coord[2] = pos[2];
  }
}

int Mesh::pickVertex(double x, double y)
{
  double coord[3];
  this->getWorldCoord(x, y, coord);
  int pickedId = this->getClosestPtID(coord);

  mapper->GetInput()->GetPointData()->GetScalars()->InsertTuple1(pickedId, 1);
  mapper->GetInput()->Modified();

  return pickedId;
}

void Mesh::markRegion(vtkIdType centerId, int nRing)
{
  vtkPolyData* mesh = mapper->GetInput();
  mesh->GetCellData()->GetScalars()->FillComponent(0, 1);
  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
  vtkIdType nPts;
  vtkIdType* pts;
  std::vector<vtkIdType> centers;
  mapMarkReg.clear();
  markRegBound.clear();
  std::map<vtkIdType, int>::iterator iter_mapPtIds;
  mapMarkReg[centerId] = 0;
  centers.push_back(centerId);
  for (int i = 0; i < nRing; ++i)
  {
    std::vector<vtkIdType> newCenters;
    for (int j = 0; j < centers.size(); ++j)
    {
      mesh->GetPointCells(centers[j], cellIds);
      for (vtkIdType iCell = 0; iCell < cellIds->GetNumberOfIds(); ++iCell)
      {
        // draw all cell connecting to center vertex
        if (i != nRing)
        {
          mesh->GetCellData()->GetScalars()->InsertTuple1(cellIds->GetId(iCell), 5);
        }

        // search i-Ring vertex set
        // put them in to mapMarkReg
        mesh->GetCellPoints(cellIds->GetId(iCell), nPts, pts);

        for (vtkIdType iPt = 0; iPt < nPts; ++iPt)
        {
          iter_mapPtIds = mapMarkReg.find(pts[iPt]);
          if (iter_mapPtIds == mapMarkReg.end())
          {
            if (i != nRing)
            {
              mapMarkReg[pts[iPt]] = i + 1;
            }
            else
            {
              markRegBound.push_back(pts[iPt]);
            }
            newCenters.push_back(pts[iPt]);
            //mesh->GetPointData()->GetScalars()->InsertTuple1(pts[iPt], 5);
          }
        }
      }
    }
    centers = newCenters;
  }
  mapper->GetInput()->Modified();
}

std::map<vtkIdType, int>& Mesh::getMarkRegion()
{
  return mapMarkReg;
}

std::vector<int>& Mesh::getMarkRegionBound()
{
  return markRegBound;
}