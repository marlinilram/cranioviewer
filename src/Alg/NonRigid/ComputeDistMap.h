#ifndef ComputeDistMap_H
#define ComputeDistMap_H

#include "itkFastMarchingImageToNodePairContainerAdaptor.h"
#include "itkFastMarchingImageFilterBase.h"
#include "itkFastMarchingThresholdStoppingCriterion.h"
#include "itkImageFileWriter.h"
#include "itkImageToVTKImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"

#include <vtkImageActor.h>
#include <vtkImageViewer.h>
#include <vtkRenderWindow.h>
#include <vtkRendererCollection.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLImageDataReader.h>

#include <ctime>

//#define DEBUG_PRINT

class ComputeDistMap
{
    typedef float PixelType;
    static const unsigned int Dimension = 3;
    typedef itk::Image< PixelType, Dimension > FloatImageType;
    typedef itk::Image< short, Dimension > ShortImageType;
    typedef itk::FastMarchingThresholdStoppingCriterion< FloatImageType, FloatImageType > CriterionType;
    typedef itk::FastMarchingImageFilterBase< FloatImageType, FloatImageType > FastMarchingType;

    typedef FastMarchingType::NodeType NodeType;
    typedef FastMarchingType::NodePairType NodePairType;
    typedef FastMarchingType::NodePairContainerType NodePairContainerType;

    typedef itk::ImageToVTKImageFilter<FloatImageType> ConnectorType;
    typedef itk::ImageToVTKImageFilter<ShortImageType> ShortConnectorType;
    typedef itk::FastMarchingImageToNodePairContainerAdaptor< FloatImageType, FloatImageType, FloatImageType > AdaptorType;
public:
    ComputeDistMap();
    ~ComputeDistMap();

    void runFastMarching();
    vtkSmartPointer<vtkImageData> computeOuterDist();
    vtkSmartPointer<vtkImageData> computeInnerDist();
    void computeFinalDistMap();
    void setVTKImg(vtkSmartPointer<vtkImageData> img_data, double mid_interval = 550, double wid_interval = 450);
    void setActiveImg(FloatImageType::Pointer active_ptr, bool out_tag = true);
    void setTrialImg(FloatImageType::Pointer trial_ptr, FloatImageType::Pointer active_ptr);
    vtkSmartPointer<vtkImageData> getDistMap() { return dist_map; };
    void gaussianSmooth(vtkSmartPointer<vtkImageData> img);

private:
    CriterionType::Pointer criterion;
    FloatImageType::SizeType size;
    FloatImageType::RegionType region;
    FloatImageType::SpacingType output_spacing;

    vtkSmartPointer<vtkImageData> bone_img;
    vtkSmartPointer<vtkImageData> dist_map;
    double max_interval;
    double min_interval;
};
#endif