#include "NiiLoader.h"
#include <iostream>
using std::cout;

NiiLoader::NiiLoader()
{
}

NiiLoader::NiiLoader(const char* niiFile)
{
  loadNii(niiFile);
}

void NiiLoader::loadNii(const char* niiFile)
{
  //binary nii format part is copied from my ..\NiiLoader\ project which is adapted from http://sourceforge.net/projects/niftilib/

  cout << "Voxel grid initializing (to " << niiFile << ")..\n";

  nifti_1_header hdr;
  FILE *fp;
  int ret,i;
  double totalIntensity;
  MY_DATATYPE *data = NULL;

  /********** open and read header */
  fp = fopen(niiFile,"rb"); //seeing this took my 2 hours :(
  if (fp == NULL) {
    fprintf(stderr, "\nError opening header file %s\n",niiFile);
    exit(1);
  }
  ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
  if (ret != 1) {
    fprintf(stderr, "\nError reading header file %s\n",niiFile);
    exit(1);
  }

  fclose(fp); //close here 'cos i'll reopen it below with a jump to the start of the image data

  /********** print a little header information */
  fprintf(stderr, "\nnii file header information:");
  fprintf(stderr, "\nXYZT dimensions: %d %d %d %d",hdr.dim[1],hdr.dim[2],hdr.dim[3],hdr.dim[4]);
  fprintf(stderr, "\nDatatype code and bits/pixel: %d %d",hdr.datatype,hdr.bitpix); //for my CTSample\ folder, datatype = 4 which is NIFTI_TYPE_INT16
  fprintf(stderr, "\nScaling slope and intercept: %.6f %.6f",hdr.scl_slope,hdr.scl_inter);
  fprintf(stderr, "\nByte offset to data in datafile: %ld",(long)(hdr.vox_offset));
  fprintf(stderr, "\nVoxel spacing: %f, %f, %f\n", hdr.pixdim[1], hdr.pixdim[2], hdr.pixdim[3]);

  /********** open the datafile, jump to data offset */
  fp = fopen(niiFile,"rb");
  if (fp == NULL) {
    fprintf(stderr, "\nError opening data file %s\n",niiFile);
    exit(1);
  }

  ret = fseek(fp, (long)(hdr.vox_offset), SEEK_SET); //jump
  if (ret != 0) {
    fprintf(stderr, "\nError doing fseek() to %ld in data file %s\n",(long)(hdr.vox_offset), niiFile);
    exit(1);
  }

  /********** allocate buffer and read first 3D volume from data file */
  data = (MY_DATATYPE *) malloc(sizeof(MY_DATATYPE) * hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
  if (data == NULL) {
    fprintf(stderr, "\nError allocating data buffer for %s\n",niiFile);
    exit(1);
  }

  //	printf("\n\nmemo allocated for data: %f x %d\n", (float) sizeof(MY_DATATYPE), hdr.dim[1]*hdr.dim[2]*hdr.dim[3]); //(float) and (long) crucial o/w 0 printed
  //	printf("# bytes jumped to skip header: %ld\n", (long) hdr.vox_offset);
  ret = fread(data, sizeof(MY_DATATYPE), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], fp);
  //	printf("\n# elements successfully read: %d == %d\n", ret, hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
  if (ret != hdr.dim[1]*hdr.dim[2]*hdr.dim[3]) {
    fprintf(stderr, "\nError reading volume 1 from %s (%d)\n",niiFile,ret);
    exit(1);
  }
  fclose(fp);

  /********** scale the data buffer  */
  if (hdr.scl_slope != 0) {
    for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; ++i)
      data[i] = (data[i] * hdr.scl_slope) + hdr.scl_inter; //slope is 1.00f and inter = -1024.00f so even if when MY_DATATYPE is unsigned short, i don't lose data
  }


  //find mean of data and more importantly initialize each voxel
  totalIntensity = 0;
  width = hdr.dim[1];
  height = hdr.dim[2];
  depth = hdr.dim[3]; //of each cell to be added to voxels[] below
  voxelWidth = hdr.pixdim[1];
  voxelHeight = hdr.pixdim[2];
  voxelDepth = hdr.pixdim[3];
  int toRight = 0, toUp = 0, toFront = 0;
  double x, y, z;
  //all voxels that includes bones, brains, eyes, .. everything
  //Voxel** voxelsAll = new Voxel*[512*512*365]; and [i] = new Voxel() below memory crashes on 6GB desktop;
  //int* voxelsAll = new int[512*512*365]; no memory crash 'cos int is way cheaper than Voxel*; already keeping them in data[] anyway

  niiImg = vtkSmartPointer<vtkImageData>::New();
  niiImg->SetDimensions(width, height, depth);
  niiImg->SetSpacing(voxelWidth, voxelHeight, voxelDepth);
  //niiImg->SetScalarTypeToShort();
  //niiImg->SetScalarTypeToFloat();
  //niiImg->SetNumberOfScalarComponents(1);
  niiImg->AllocateScalars(VTK_SHORT, 1);
  signed short *ptr = (signed short *)niiImg->GetScalarPointer();
  //fArray = new float[width*height*depth];
  int nSkullPt = 0;


  for (i = 0; i < hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; ++i)
  {


    //id, intensity, and coord of this current voxel
    x = toRight * hdr.pixdim[1]; //voxel spacing, space b/w x coords, i.e. pixdim[1]

    y = toUp * hdr.pixdim[2]; //voxel spacing, space b/w y coords, i.e. pixdim[2]
    z = toFront * hdr.pixdim[3]; //voxel spacing, space b/w z coords, i.e. pixdim[3]

    //load data
    //fArray[i] = data[i];
    *ptr ++ = data[i];
    if ( isSkullVoxel(data[i]) )
    {
      //if (toRight % 5 == 0 && toUp % 5 == 0 && toFront % 5 == 0) //extra thresholding to downsample the bone voxels by drawing every fifth in x-direction and in y-direction and in z-direction
      //if (nSkullPt % 50 == 0)
      {
        //skullVoxel.push_back( x ); //downsampled active/bone voxels
        //skullVoxel.push_back( y );
        //skullVoxel.push_back( z );	
        //skullVoxelVal.push_back( data[i] );
      }
      nSkullPt += 1;
    }
    //else fArray[i] = 0;//data[i];//fArray[i] = 0;

    ++toRight;
    //direction/pointer adjustments
    if (toRight >= hdr.dim[1])
    {
      toRight = 0; //reset to the left end
      ++toUp; //1 level up
      if (toUp >= hdr.dim[2])
      {
        ++toFront; //1 slice towards me, that is towards front
        toUp = 0; //start from the bottom level
      }
    }
  } //end of i
  free(data);
  cout << "total " << nSkullPt << " skull voxels\n";
}

vtkSmartPointer< vtkImageData > NiiLoader::getData()
{
  // return fArray
  return niiImg;
}

void NiiLoader::clearData()
{
  niiImg = vtkSmartPointer<vtkImageData>::New();
  vector<double>().swap(skullVoxel);
  vector<double>().swap(skullVoxelForIso);
  vector<signed short>().swap(skullVoxelVal);
}

std::vector<double>& NiiLoader::extractSkullVertex(int v, int t, int& num)
{
  signed short* data = (signed short*)niiImg->GetScalarPointer();
  vector<signed short> fArray;
  int nSkullPt = 0;
  int i;
  int toRight = 0, toUp = 0, toFront = 0;
  skullVoxelForIso.clear();

  for (i = 0; i < width*height*depth; i++)
  {
    if (data[i] - v <= t && data[i] - v >= -t) {
      skullVoxelForIso.push_back(toRight * voxelWidth);
      skullVoxelForIso.push_back(toUp * voxelHeight);
      skullVoxelForIso.push_back(toFront * voxelDepth);
    }

    //direction/pointer adjustments
    toRight++;
    if (toRight >= width)
    {
      toRight = 0; //reset to the left end
      toUp++; //1 level up
      if (toUp >= height)
      {
        toFront++; //1 slice towards me, that is towards front
        toUp = 0; //start from the bottom level
      }
    }
  } //end of i
  num = skullVoxelForIso.size()/3;
  return skullVoxelForIso;
};

void NiiLoader::getUnique(vector<signed short> &u)
{
  u = skullVoxelVal;
  std::sort(u.begin(), u.end());
  vector<signed short>::iterator iter = std::unique(u.begin(), u.end());
  u.erase(iter, u.end());
}

vtkSmartPointer< vtkPolyData > NiiLoader::mcSkullVertex(vector<double> &T_Pts, vector<double> &T_Norms, double isovalue, int sample_ratio)
{
  int downsample_ritio = sample_ratio;
  int d_width = width%downsample_ritio ? 1+width/downsample_ritio : width/downsample_ritio;
  int d_height = height%downsample_ritio ? 1+height/downsample_ritio : height/downsample_ritio;
  int d_depth = depth%downsample_ritio ? 1+depth/downsample_ritio : depth/downsample_ritio; //of each cell to be added to voxels[] below
  float d_voxelWidth = voxelWidth*downsample_ritio;
  float d_voxelHeight = voxelHeight*downsample_ritio;
  float d_voxelDepth = voxelDepth*downsample_ritio;
  int toRight = 0, toUp = 0, toFront = 0;

  signed short* data = (signed short*)niiImg->GetScalarPointer();
  vector<signed short> fArray;
  int nSkullPt = 0;
  int i;

  for (i = 0; i < width*height*depth; i++)
  {
    if (toRight % downsample_ritio == 0 && toUp % downsample_ritio == 0 && toFront % downsample_ritio == 0) fArray.push_back(data[i]);

    //direction/pointer adjustments
    toRight++;
    if (toRight >= width)
    {
      toRight = 0; //reset to the left end
      toUp++; //1 level up
      if (toUp >= height)
      {
        toFront++; //1 slice towards me, that is towards front
        toUp = 0; //start from the bottom level
      }
    }
  } //end of i

  vtkSmartPointer<vtkStructuredPoints> sPoints = vtkSmartPointer<vtkStructuredPoints>::New();
  sPoints->SetDimensions(d_width,d_height,d_depth);
  sPoints->SetSpacing(d_voxelWidth,d_voxelHeight,d_voxelDepth);
  //sPoints->SetScalarTypeToShort();
  //sPoints->AllocateScalars(VTK_SHORT, 1);

  vtkSmartPointer<vtkShortArray> vtkFArray = vtkSmartPointer<vtkShortArray>::New();
  //vtkFArray->SetNumberOfValues(depth*height*width);
  vtkFArray->SetArray(&fArray[0], d_depth*d_height*d_width, 1);

  sPoints->GetPointData()->SetScalars(vtkFArray);
  //sPoints->Update();

  /* create iso surface with marching cubes algorithm */
  fprintf(stderr, "\nrunning marching cubes");
  vtkSmartPointer<vtkMarchingCubes> mcSource = vtkSmartPointer<vtkMarchingCubes>::New();
  mcSource->SetInputData(sPoints);
  mcSource->ComputeNormalsOn();
  mcSource->SetNumberOfContours(1);
  mcSource->SetValue(0,isovalue);
  mcSource->Update();

  fprintf(stderr, "\nrunning connectivity filter\n");
  vtkSmartPointer<vtkPolyDataConnectivityFilter> cfilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
  cfilter->SetInputConnection(mcSource->GetOutputPort());
  cfilter->SetExtractionModeToLargestRegion();
  //cfilter->SetExtractionModeToPointSeededRegions();
  //cfilter->SetClosestPoint(x,y,z);
  cfilter->Update();

  vtkSmartPointer< vtkPolyData > temp = vtkSmartPointer< vtkPolyData >::New();
  temp = cfilter->GetOutput();
  //cout << "Total points: " << temp->GetNumberOfPoints() << endl;
  //cout << temp->GetPoints()->GetPoint(0)[0] << "\t" << temp->GetPoints()->GetPoint(0)[1] << "\t" << temp->GetPoints()->GetPoint(0)[2] << endl;
  //cout << temp->GetPointData()->GetNormals()->GetTuple(0)[0] << "\t" << temp->GetPointData()->GetNormals()->GetTuple(0)[1] << "\t" << temp->GetPointData()->GetNormals()->GetTuple(0)[2] << endl;
  //system("pause");

  T_Pts.clear();
  T_Norms.clear();
  for (decltype(temp->GetNumberOfPoints()) i = 0; i != temp->GetNumberOfPoints(); ++i) {
    T_Pts.push_back(temp->GetPoints()->GetPoint(i)[0]);
    T_Pts.push_back(temp->GetPoints()->GetPoint(i)[1]);
    T_Pts.push_back(temp->GetPoints()->GetPoint(i)[2]);
    T_Norms.push_back(temp->GetPointData()->GetNormals()->GetTuple(i)[0]);
    T_Norms.push_back(temp->GetPointData()->GetNormals()->GetTuple(i)[1]);
    T_Norms.push_back(temp->GetPointData()->GetNormals()->GetTuple(i)[2]);
    //cout << i << endl;
  }
  return temp;
}