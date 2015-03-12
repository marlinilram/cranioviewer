#include "AddCrsp.h"

vtkStandardNewMacro(MouseInteractorStylePP);

AddCrsp::AddCrsp()
{
	setupUi(this);
	this->show();

	source_renderer = vtkSmartPointer< vtkRenderer >::New();
	m_SourceViewer->GetRenderWindow()->AddRenderer(source_renderer);
	source_renderer->SetBackground(0.3, 0.3, 0.3);
	target_renderer = vtkSmartPointer< vtkRenderer >::New();
	m_TargetViewer->GetRenderWindow()->AddRenderer(target_renderer);
	target_renderer->SetBackground(0.3, 0.3, 0.3);

	source_pointpicker = vtkSmartPointer< MouseInteractorStylePP >::New();
	target_pointpicker = vtkSmartPointer< MouseInteractorStylePP >::New();

	connect( m_PushButtonGetCrspPt, SIGNAL( clicked() ), this, SLOT( getCrsp()) );

	source_user_crsp.clear();
	target_user_crsp.clear();
}

AddCrsp::~AddCrsp()
{
}

void AddCrsp::passPara(vtkSmartPointer< vtkActor > &source_actor, vtkSmartPointer< vtkActor > &target_arcor)
{
	source_renderer->AddActor(source_actor);
	m_SourceViewer->GetRenderWindow()->GetInteractor()->SetInteractorStyle(source_pointpicker);
	target_renderer->AddActor(target_arcor);
	m_TargetViewer->GetRenderWindow()->GetInteractor()->SetInteractorStyle(target_pointpicker);
}

void AddCrsp::getCrsp()
{
	cout << "Point ID in source: " << source_pointpicker->get_picked_ID() << "\tPoint ID intarget: " << target_pointpicker->get_picked_ID() << endl;
	source_user_crsp.push_back(source_pointpicker->get_picked_ID());
	target_user_crsp.push_back(target_pointpicker->get_picked_ID());
}

void AddCrsp::passCrsp(vector<vtkIdType> &idx_U, vector<vtkIdType> &idx_T_U)
{
	if (source_user_crsp.size() == target_user_crsp.size()) {
		idx_U = source_user_crsp;
		idx_T_U = target_user_crsp;
	}
	else cout << "Number of user define correspondences in Source and Target unmatch" << endl;
}