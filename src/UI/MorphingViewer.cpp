#include "MorphingViewer.h"

MorphingViewer::MorphingViewer()
{
    setupUi(this);
    morphing_renderer = nullptr;
    morphing_wrapper = new MorphingWrapper;

    connect( pushButtonAddMesh, SIGNAL( clicked() ), this, SLOT( addMesh() ) );
}

MorphingViewer::~MorphingViewer()
{
    delete morphing_wrapper;
}

void MorphingViewer::show()
{
    QWidget::show();

    if (morphing_renderer)
        qvtkWidgetMorphing->GetRenderWindow()->RemoveRenderer(morphing_renderer);

    morphing_renderer = vtkSmartPointer<vtkRenderer>::New();
    morphing_renderer->SetBackground(0.0,0.0,0.0);
    qvtkWidgetMorphing->GetRenderWindow()->AddRenderer(morphing_renderer);
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
        morphing_wrapper->loadMesh(temp+"/"+(*qstring_iter).toStdString());
    }
    updateRenderer();
}

void MorphingViewer::updateRenderer()
{
    qvtkWidgetMorphing->GetRenderWindow()->Render();
}