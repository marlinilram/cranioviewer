#include "MorphingViewer.h"

MorphingViewer::MorphingViewer()
{
    setupUi(this);
    morphing_renderer = nullptr;
    morphing_wrapper = new MorphingWrapper;

    morphing_renderer = vtkSmartPointer<vtkRenderer>::New();
    morphing_renderer->SetBackground(0.0,0.0,0.0);
    qvtkWidgetMorphing->GetRenderWindow()->AddRenderer(morphing_renderer);

    connect( pushButtonAddMesh, SIGNAL( clicked() ), this, SLOT( addMesh() ) );
    connect( pushButtonSaveMesh, SIGNAL( clicked() ), this, SLOT( saveMesh() ) );
    connect( sliderMorphing, SIGNAL( sliderMoved( int ) ), this, SLOT( updateCurMorph( int ) ) );
    connect( sliderMorphing, SIGNAL( valueChanged( int ) ), this, SLOT( updateCurMorph( int ) ) );
}

MorphingViewer::~MorphingViewer()
{
    delete morphing_wrapper;
}

void MorphingViewer::show()
{
    QWidget::show();
}

void MorphingViewer::addMesh()
{
    QDir dir;
    //QString fileName = QFileDialog::getOpenFileName( this, QString(tr("Open Image")), dir.absolutePath() , filter );
    QString fileDir = QFileDialog::getExistingDirectory(this, QString(tr("Choose Dir")),dir.absolutePath(), QFileDialog::ShowDirsOnly
        | QFileDialog::DontResolveSymlinks);
    if ( fileDir.isEmpty() == true ) return;

    // 支持带中文路径的读取
    std::string temp = fileDir.toStdString();

    QStringList nameFilter("*.obj");
    QDir directory(fileDir);
    QStringList txtFilesAndDirectories = directory.entryList(nameFilter);
    QStringList::iterator qstring_iter = txtFilesAndDirectories.begin();
    for (; qstring_iter != txtFilesAndDirectories.end(); ++qstring_iter)
    {
        morphing_wrapper->loadMesh(temp+"/"+(*qstring_iter).toStdString(), morphing_renderer);
    }
    morphing_renderer->ResetCamera();
    morphing_wrapper->setCenterMesh(0);
    morphing_wrapper->getMorphingHandler()->computeTransform();
    updateRenderer();
}

void MorphingViewer::updateRenderer()
{
    qvtkWidgetMorphing->GetRenderWindow()->Render();
}

void MorphingViewer::updateCurMorph(int val)
{
    double morph_paras[2];
    morph_paras[1] = (double)val / sliderMorphing->maximum();
    morph_paras[0] = 1-morph_paras[1];

    morphing_wrapper->doMorphing(morph_paras, 2);
    updateRenderer();
}

void MorphingViewer::saveMesh()
{
    vtkSmartPointer<vtkOBJWriter> writer = vtkSmartPointer<vtkOBJWriter>::New();
    writer->SetInputData(morphing_wrapper->getMeshPtr(0)->getMeshData());
    writer->SetFileName("morph_result.obj");
    writer->Update();
}