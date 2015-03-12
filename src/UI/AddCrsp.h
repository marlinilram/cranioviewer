#ifndef ADDCRSP_H
#define ADDCRSP_H

#include <QMainWindow>
#include "ui_AddCrsp.h"

#include "vtkEventQtSlotConnect.h"
#include "vtkSmartPointer.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkInteractorStyleImage.h"
#include "vtkImageViewer2.h"
#include "vtkImageData.h"
#include "vtkImageActor.h"
#include "vtkImageSlice.h"
#include "vtkImageProperty.h"
#include "vtkMapper.h"
#include "vtkImageSliceMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkActorCollection.h"
#include "vtkOBJReader.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkTransform.h"
#include "vtkBoxWidget.h"
#include "vtkCutter.h"
#include "vtkPlane.h"
#include "vtkImplicitFunction.h"
#include "vtkStripper.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkCamera.h"
#include "vtkMatrix4x4.h"
#include "vtkLinearTransform.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkPoints.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkPointData.h"
#include "vtkOBJExporter.h"
#include "vtkIdList.h"
#include "vtkIdentityTransform.h"
#include "vtkStructuredPoints.h"
#include "vtkShortArray.h"
#include "vtkMarchingCubes.h"
#include "vtkUnsignedCharArray.h"
#include "vtkPointPicker.h"
#include "vtkCellPicker.h"
#include "vtkObjectFactory.h"
#include "vtkRendererCollection.h"
#include "vtkSphereSource.h"
#include "vtkCell.h"

#include <vector>
using std::vector;

// Define interaction style
class MouseInteractorStylePP : public vtkInteractorStyleTrackballCamera
{
  public:
    static MouseInteractorStylePP* New();
    vtkTypeMacro(MouseInteractorStylePP, vtkInteractorStyleTrackballCamera);

	virtual void OnKeyPress()
	{
		switch (this->Interactor->GetKeyCode()) {
			case '1':
				setCrsp(this->Interactor->GetEventPosition()[0], this->Interactor->GetEventPosition()[1]);
				break;
			default:
				break;
		}

		vtkInteractorStyleTrackballCamera::OnKeyPress();
	}

	vtkIdType get_picked_ID()
	{
		return picked_ID;
	}

private:
	void setCrsp(int event_position_x, int event_position_y)
	{
		std::cout << "Picking pixel: " << this->Interactor->GetEventPosition()[0] << " " << this->Interactor->GetEventPosition()[1] << std::endl;
		this->Interactor->GetPicker()->Pick(this->Interactor->GetEventPosition()[0],
			this->Interactor->GetEventPosition()[1],
			0,  // always zero.
			this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
		double picked[3];
		this->Interactor->GetPicker()->GetPickPosition(picked);
		std::cout << "Picked value: " << picked[0] << " " << picked[1] << " " << picked[2] << std::endl;

		//Cell picker
		vtkSmartPointer< vtkCellPicker > cell_picker = vtkSmartPointer< vtkCellPicker >::New();
		cell_picker->Pick(this->Interactor->GetEventPosition()[0],
			this->Interactor->GetEventPosition()[1],
			0,
			this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
		picked_ID = cell_picker->GetCellId();
		double *picked_from_actor;
		double dist = std::numeric_limits < double >::max();
		vtkSmartPointer< vtkPolyData > mesh = vtkPolyDataMapper::SafeDownCast(cell_picker->GetActor()->GetMapper())->GetInput();
		vtkSmartPointer< vtkIdList > picked_pt_list = vtkSmartPointer< vtkIdList >::New();
		mesh->GetCellPoints(picked_ID, picked_pt_list);
		for (int i = 0; i != 3; ++i) {
			picked_from_actor = mesh->GetPoint(picked_pt_list->GetId(i));
			double temp = (picked_from_actor[0] - picked[0])*(picked_from_actor[0] - picked[0]) + (picked_from_actor[1] - picked[1])*(picked_from_actor[1] - picked[1]) + (picked_from_actor[2] - picked[2])*(picked_from_actor[2] - picked[2]);
			if (temp < dist) dist = temp, picked_ID = picked_pt_list->GetId(i);
		}
		picked_from_actor = mesh->GetPoint(picked_ID);

		//Create a sphere
		vtkSmartPointer<vtkSphereSource> sphereSource =
			vtkSmartPointer<vtkSphereSource>::New();
		sphereSource->SetCenter(picked_from_actor[0], picked_from_actor[1], picked_from_actor[2]);
		sphereSource->SetRadius(1.0);

		//Create a mapper and actor
		vtkSmartPointer<vtkPolyDataMapper> mapper =
			vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(sphereSource->GetOutputPort());

		vtkSmartPointer<vtkActor> actor =
			vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetDiffuseColor(0.0, 1.0, 1.0);

		this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(actor);
		this->Interactor->GetRenderWindow()->Render();
	};

	vtkIdType picked_ID;
};
 


class AddCrsp : public QWidget, public Ui::AddCrsp
{
	Q_OBJECT
	
public:
	AddCrsp();
	~AddCrsp();

	void passPara(vtkSmartPointer< vtkActor > &source_actor, vtkSmartPointer< vtkActor > &target_actor);
	void passCrsp(std::vector<vtkIdType> &idx_U, std::vector<vtkIdType> &idx_T_U);

private slots:
	void getCrsp();

private:
	vtkSmartPointer< vtkRenderer > source_renderer;
	vtkSmartPointer< vtkRenderer > target_renderer;
	vtkSmartPointer< MouseInteractorStylePP > source_pointpicker;
	vtkSmartPointer< MouseInteractorStylePP > target_pointpicker;
	vector<vtkIdType> source_user_crsp;
	vector<vtkIdType> target_user_crsp;

};

#endif