#include "VTKEventHandler.h"

#include <vtkInteractorStyleDrawPolygon.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkInteractorStyleRubberBand3D.h>
#include <vtkTransform.h>
#include "SelectModeStyle.h"
#include "RegionMoveStyle.h"
#include "FullRegionMoveStyle.h"

VTKEventHandler::VTKEventHandler()
{
    event_connector = vtkSmartPointer<vtkEventQtSlotConnect>::New();

    for (size_t i = 0; i < 4; ++i)
    {
        qvtk_widgets[i] = nullptr;
    }

    for (size_t i = 0; i < 3; ++i)
    {
        imge_slice_widgets[i] = nullptr;
    }

    status_bar = nullptr;
    main_viewer = nullptr;
    cur_crsp_count = 0;

    picked_ids.clear();
    crsp_pos.clear();
    picked_vertices_actors.clear();

    select_mode = false;
    regionMove_mode = false;

    centerTrackball = vtkSmartPointer<vtkAxesActor>::New();
    centerTrackball->SetTotalLength(10,10,10);
}

VTKEventHandler::~VTKEventHandler()
{

}

void VTKEventHandler::setQVtkWidgets(QVTKWidget *mainWin_widgets[4])
{
    for (size_t i = 0; i < 4; ++i)
    {
        qvtk_widgets[i] = mainWin_widgets[i];
        event_connector->Connect(qvtk_widgets[i]->GetRenderWindow()->GetInteractor(), vtkCommand::MouseMoveEvent, this, SLOT(displayCoord(vtkObject *)));
        event_connector->Connect(qvtk_widgets[i]->GetRenderWindow()->GetInteractor(), vtkCommand::LeftButtonPressEvent, this, SLOT(getWorldCoord(vtkObject *)));
        event_connector->Connect(qvtk_widgets[i]->GetRenderWindow()->GetInteractor(), vtkCommand::KeyReleaseEvent, this, SLOT(setUserCrsp(vtkObject *)));
    }
}

void VTKEventHandler::setImageSliceWidgets(ImageSliceWidget *widgets[3])
{
    for (size_t i = 0; i < 3; ++i)
    {
        imge_slice_widgets[i] = widgets[i];
    }
}

void VTKEventHandler::displayCoord(vtkObject *obj)
{
    vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::SafeDownCast(obj);

    int event_pos[2];
    iren->GetEventPosition(event_pos);

    vtkSmartPointer<vtkRenderer> cur_renderer = iren->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
    double val = 0.0;
    for (size_t i = 0; i < 3; ++i)
    {
        //std::cout<<cur_renderer.GetPointer()<<"\t"<<imge_slice_widgets[i]->getImageSlice()->getRenderer().GetPointer()<<"\n";
        if (imge_slice_widgets[i]->getImageSlice() && cur_renderer.GetPointer() == imge_slice_widgets[i]->getImageSlice()->getRenderer().GetPointer())
        {
            val = imge_slice_widgets[i]->getImageSlice()->getPixelVal(event_pos[0], event_pos[1]);
        }
        if (imge_slice_widgets[i]->getIntersector() && cur_renderer.GetPointer() == imge_slice_widgets[i]->getIntersector()->getRenderer().GetPointer())
        {
            imge_slice_widgets[i]->getIntersector()->getWorldCoord(event_pos[0], event_pos[1]);
        }
    }


    QString str;
    str.sprintf("x=%d : y=%d : val=%f", event_pos[0], event_pos[1], val);
    
    status_bar->showMessage(str);
}

void VTKEventHandler::setStatusBar(QStatusBar *mainWin_statusBar)
{
    status_bar = mainWin_statusBar;
}

void VTKEventHandler::getWorldCoord(vtkObject *obj)
{
    vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::SafeDownCast(obj);

    int event_pos[2];
    iren->GetEventPosition(event_pos);

    double coord[3] = {0,0,0};

    vtkSmartPointer<vtkRenderer> cur_renderer = iren->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
    for (size_t i = 0; i < 3; ++i)
    {
        if (imge_slice_widgets[i]->getIntersector() && cur_renderer.GetPointer() == imge_slice_widgets[i]->getIntersector()->getRenderer().GetPointer())
        {
            imge_slice_widgets[i]->getIntersector()->getWorldCoord(event_pos[0], event_pos[1], coord);
            //main_viewer->getMeshData()->getClosestPtID(coord);
        }
    }
}

void VTKEventHandler::setMainViewer(MainViewer *viewer)
{
    main_viewer = viewer;
}

void VTKEventHandler::initUserCrsp()
{
    main_viewer->getMeshData()->setKDTreeLocator();

    cur_crsp_count = 0;

    for (size_t i = 0; i < picked_vertices_actors.size(); ++i)
    {
      main_viewer->getRenderer()->RemoveActor(picked_vertices_actors[i]);
    }

    picked_ids.clear();
    crsp_pos.clear();
    picked_vertices_actors.clear();

    emit(main_viewer->updateRenderers());
}

void VTKEventHandler::setUserCrsp(vtkObject *obj)
{
    vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::SafeDownCast(obj);
    std::string key = iren->GetKeySym();

    int event_pos[2];
    iren->GetEventPosition(event_pos);

    vtkSmartPointer<vtkRenderer> cur_renderer = iren->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
    for (size_t i = 0; i < 3; ++i)
    {
        //std::cout<<cur_renderer.GetPointer()<<"\t"<<imge_slice_widgets[i]->getImageSlice()->getRenderer().GetPointer()<<"\n";
        if (imge_slice_widgets[i]->getIntersector() && cur_renderer.GetPointer() == imge_slice_widgets[i]->getIntersector()->getRenderer().GetPointer())
        {
            double coord[3] = {0,0,0};
            imge_slice_widgets[i]->getIntersector()->getWorldCoord(event_pos[0], event_pos[1], coord);
            if (key == "1")
            {
                int picked_id = main_viewer->getMeshData()->getClosestPtID(coord);
                if (picked_id >= 0) picked_ids.push_back(picked_id);
                std::cout<<"picked id: "<<picked_id<<"\tCurrent Crsps: "<<picked_ids.size()<<"\n";
            }
            else if (key == "2")
            {
                crsp_pos.push_back(coord[0]);
                crsp_pos.push_back(coord[1]);
                crsp_pos.push_back(coord[2]);
                std::cout<<"crsp pos: "<<coord[0]<<"\t"<<coord[1]<<"\t"<<coord[2]<<"\tCurrent Crsps: "<<crsp_pos.size()/3<<"\n";
            }
        }
    }

    //Create a line
    if ((picked_ids.size() == crsp_pos.size()/3) && (picked_ids.size() == (cur_crsp_count + 1)))
    {
      ++cur_crsp_count;

      double *vertex_coord = main_viewer->getMeshData()->getMeshData()->GetPoint(picked_ids.back());

      vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
      lineSource->SetPoint1(vertex_coord);
      lineSource->SetPoint2(&crsp_pos[3*(cur_crsp_count-1)]);
      lineSource->Update();

      //  vtkSmartPointer<vtkSphereSource> sphereSource =
      //  vtkSmartPointer<vtkSphereSource>::New();
      //sphereSource->SetCenter(vertex_coord[0], vertex_coord[1], vertex_coord[2]);
      //sphereSource->SetRadius(1.0);

      //Create a mapper and actor
      vtkSmartPointer<vtkPolyDataMapper> pick_mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
      pick_mapper->SetInputConnection(lineSource->GetOutputPort());

      vtkSmartPointer<vtkActor> pick_actor =
        vtkSmartPointer<vtkActor>::New();
      pick_actor->SetMapper(pick_mapper);
      pick_actor->GetProperty()->SetDiffuseColor(1.0, 0.0, 0.0);
      pick_actor->GetProperty()->SetLineWidth(2);

      picked_vertices_actors.push_back(pick_actor);
      main_viewer->getRenderer()->AddActor(picked_vertices_actors.back());
      emit(main_viewer->updateRenderers());
    }
}

void VTKEventHandler::selectMode()
{
  if (select_mode)
  {
    qvtk_widgets[0]->GetRenderWindow()->GetInteractor()->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New());

    main_viewer->getRenderer()->RemoveActor(centerTrackball);

    event_connector->Disconnect(qvtk_widgets[0]->GetRenderWindow()->GetInteractor(), vtkCommand::StartInteractionEvent, this, SLOT(getSelectCenter(vtkObject *)));
    event_connector->Disconnect(qvtk_widgets[0]->GetRenderWindow()->GetInteractor(), vtkCommand::MouseMoveEvent, this, SLOT(getSelectRad(vtkObject *)));

    select_mode = false;
  }
  else
  {
    qvtk_widgets[0]->GetRenderWindow()->GetInteractor()->SetInteractorStyle(vtkSmartPointer<SelectModeStyle>::New());

    main_viewer->getMeshData()->getActor()->GetMapper()->SetScalarModeToUseCellData();

    //event_connector->Disconnect(qvtk_widgets[0]->GetRenderWindow()->GetInteractor(), vtkCommand::StartInteractionEvent, this, SLOT(initRegionMoveMode(vtkObject *)));
    event_connector->Disconnect(qvtk_widgets[0]->GetRenderWindow()->GetInteractor(), vtkCommand::MouseMoveEvent, this, SLOT(regionMoveCall(vtkObject *)));

    event_connector->Connect(qvtk_widgets[0]->GetRenderWindow()->GetInteractor(), vtkCommand::StartInteractionEvent, this, SLOT(getSelectCenter(vtkObject *)));
    event_connector->Connect(qvtk_widgets[0]->GetRenderWindow()->GetInteractor(), vtkCommand::MouseMoveEvent, this, SLOT(getSelectRad(vtkObject *)));

    select_mode = true;
    regionMove_mode = false;
  }
}

void VTKEventHandler::getSelectCenter(vtkObject *obj)
{
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::SafeDownCast(obj);
  SelectModeStyle* selectModeStyle = SelectModeStyle::SafeDownCast(iren->GetInteractorStyle());

  selectModeStyle->SetDrawPolygonPixels(false);

  int startPos[2];
  selectModeStyle->GetStartPos(startPos);
  
  if (main_viewer->getMeshData())
  {
    main_viewer->getMeshData()->setKDTreeLocator();

    selectModeStyle->SetStartPtId(main_viewer->getMeshData()->pickVertex(startPos[0], startPos[1]));

    centerTransformPtId = selectModeStyle->GetStartPtId();
    double selectCenter[3];
    main_viewer->getMeshData()->getMeshData()->GetPoint(selectModeStyle->GetStartPtId(), selectCenter);
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->Translate(selectCenter);
    centerTrackball->SetUserTransform(transform);
    main_viewer->getRenderer()->RemoveActor(centerTrackball);
    main_viewer->getRenderer()->AddActor(centerTrackball);
    std::cout<<"Selected Center Pos: "<<selectCenter[0]<<" "<<selectCenter[1]<<" "<<selectCenter[2]<<"\n";

    emit(main_viewer->updateRenderers());
  }

  std::cout<<"Start Pos: "<<startPos[0]<<", "<<startPos[1]<<"\n";
}

void VTKEventHandler::getSelectRad(vtkObject *obj)
{
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::SafeDownCast(obj);
  SelectModeStyle* selectModeStyle = SelectModeStyle::SafeDownCast(iren->GetInteractorStyle());

  if (!selectModeStyle->GetMovingTag())
  {
    return;
  }
  
  int endPos[2];
  selectModeStyle->GetEndPos(endPos);

  int *size = iren->GetRenderWindow()->GetSize();

  double selectRad = selectModeStyle->GetRadius();
  double selectRadMax = sqrt(size[0]*size[0]+size[1]*size[1]);

  int nRing = 100 * selectRad / selectRadMax;

  std::cout<<"nRing: "<<nRing<<"\n";

  int centerId = selectModeStyle->GetStartPtId();

  if (main_viewer->getMeshData())
  {
    main_viewer->getMeshData()->markRegion(centerId, nRing);

    emit(main_viewer->updateRenderers());
  }

  std::cout<<"End Pos: "<<endPos[0] <<", "<<endPos[1]<<"\n";

  std::cout<<"Select Radius: "<<selectModeStyle->GetRadius()<<"\n";
}

void VTKEventHandler::regionMoveMode()
{
  if (!regionMove_mode)
  {
    qvtk_widgets[0]->GetRenderWindow()->GetInteractor()->SetInteractorStyle(vtkSmartPointer<FullRegionMoveStyle>::New());

    FullRegionMoveStyle* fullRegionMove = FullRegionMoveStyle::SafeDownCast(
      qvtk_widgets[0]->GetRenderWindow()->GetInteractor()->GetInteractorStyle());

    event_connector->Disconnect(qvtk_widgets[0]->GetRenderWindow()->GetInteractor(), vtkCommand::StartInteractionEvent, this, SLOT(getSelectCenter(vtkObject *)));
    event_connector->Disconnect(qvtk_widgets[0]->GetRenderWindow()->GetInteractor(), vtkCommand::MouseMoveEvent, this, SLOT(getSelectRad(vtkObject *)));

    nonRigid->getNonRigid()->setSubMesh(
      main_viewer->getMeshData()->getMarkRegion(), 
      main_viewer->getMeshData()->getMarkRegionBound());

    // initialize group of center vertices, cache them for transform
    centerOldPos.clear();
    centerPtIds.clear();
    centerPtIds.push_back(centerTransformPtId);
    Mesh* mesh = main_viewer->getMeshData();
    std::map<vtkIdType, int>& markRegion = mesh->getMarkRegion();
    std::map<vtkIdType, int>::iterator it_markRegion;
    double tempPt[3];
    for (it_markRegion = markRegion.begin(); it_markRegion != markRegion.end(); ++it_markRegion)
    {
      if (it_markRegion->second != 0 && it_markRegion->second < 1)
      {
        centerPtIds.push_back(it_markRegion->first);
      }
    }
    // transform center stores at the first place
    for (int i = 0; i < centerPtIds.size(); ++i)
    {
      mesh->getMeshData()->GetPoint(centerPtIds[i], tempPt);
      centerOldPos.push_back(tempPt[0]);
      centerOldPos.push_back(tempPt[1]);
      centerOldPos.push_back(tempPt[2]);
    }

    main_viewer->getRenderer()->RemoveActor(centerTrackball);
    vtkTransform* centerTrackballTransform = vtkTransform::SafeDownCast(centerTrackball->GetUserTransform());
    centerTrackballTransform->Identity();
    centerTrackballTransform->Translate(
      centerOldPos[0],
      centerOldPos[1],
      centerOldPos[2]);
    main_viewer->getRenderer()->AddActor(centerTrackball);
    fullRegionMove->SetInteractionProp(centerTrackball);

    //event_connector->Connect(qvtk_widgets[0]->GetRenderWindow()->GetInteractor(), vtkCommand::StartInteractionEvent, this, SLOT(initRegionMoveMode(vtkObject *)));
    event_connector->Connect(qvtk_widgets[0]->GetRenderWindow()->GetInteractor(), vtkCommand::MouseMoveEvent, this, SLOT(regionMoveCall(vtkObject *)));

    select_mode = false;
    regionMove_mode = true;
  }
  else
  {
    qvtk_widgets[0]->GetRenderWindow()->GetInteractor()->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New());

    main_viewer->getRenderer()->RemoveActor(centerTrackball);

    //event_connector->Disconnect(qvtk_widgets[0]->GetRenderWindow()->GetInteractor(), vtkCommand::StartInteractionEvent, this, SLOT(initRegionMoveMode(vtkObject *)));
    event_connector->Disconnect(qvtk_widgets[0]->GetRenderWindow()->GetInteractor(), vtkCommand::MouseMoveEvent, this, SLOT(regionMoveCall(vtkObject *)));


    regionMove_mode = false;
  }
}

void VTKEventHandler::initRegionMoveMode(vtkObject *obj)
{
  //vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::SafeDownCast(obj);
  //RegionMoveStyle* regionMoveStyle = RegionMoveStyle::SafeDownCast(iren->GetInteractorStyle());

  //std::map<vtkIdType, int>& markRegion = main_viewer->getMeshData()->getMarkRegion();
  //std::map<vtkIdType, int>::iterator it_markRegion;
  //vtkIdType ptId = 0;
  //for (it_markRegion = markRegion.begin(); it_markRegion != markRegion.end(); ++it_markRegion)
  //{
  //  if (it_markRegion->second == 0)
  //  {
  //    ptId = it_markRegion->first;
  //    break;
  //  }
  //}
  //double StartWPosition[3];
  //main_viewer->getMeshData()->getMeshData()->GetPoint(ptId, StartWPosition);

  //std::cout<<"center pt: "<<ptId<<"\n";

  //regionMoveStyle->SetStartWPosition(StartWPosition);
}

void VTKEventHandler::regionMoveCall(vtkObject *obj)
{
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::SafeDownCast(obj);
  FullRegionMoveStyle* fullRegionMoveStyle = FullRegionMoveStyle::SafeDownCast(iren->GetInteractorStyle());

  if (!fullRegionMoveStyle->GetEventState())
  {
    return;
  }

  std::cout<<"test\n";

  Mesh* mesh = main_viewer->getMeshData();
  double tempPt[3];
  centerPos.resize(centerOldPos.size());
  for (int i = 0; i < centerPtIds.size(); ++i)
  {
    tempPt[0] = centerOldPos[3*i + 0] - centerOldPos[0];
    tempPt[1] = centerOldPos[3*i + 1] - centerOldPos[1];
    tempPt[2] = centerOldPos[3*i + 2] - centerOldPos[2];
    centerTrackball->GetUserTransform()->TransformPoint(tempPt, &centerPos[3*i + 0]);
    //mesh->getMeshData()->GetPoints()->SetPoint(
    //  centerPtIds[i],
    //  &centerPos[3*i + 0]
    //  );
  }
  //mesh->getMeshData()->GetPoints()->Modified();

  nonRigid->getNonRigid()->refineSubMesh(centerPos);
  if (nonRigid->getNonRigid()->getVisCrspLines())
  {
    nonRigid->drawCrspLines();
  }

  emit(main_viewer->updateRenderers());
}

void VTKEventHandler::setNonRigid(NonRigidWrapper* nonRigid)
{
  this->nonRigid = nonRigid;
}

void VTKEventHandler::setLocalSmooth(int state)
{
  this->nonRigid->getNonRigid()->setLocalSmoothMode(state);
}