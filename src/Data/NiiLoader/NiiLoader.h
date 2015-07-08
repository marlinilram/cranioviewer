#ifndef _NIILOADER_
#define _NIILOADER_

#include "nifti1.h"
#include <vector>
#include <iterator>
#include <algorithm>
using std::vector;

#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkStructuredPoints.h"
#include "vtkShortArray.h"
#include "vtkMarchingCubes.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"

typedef signed short MY_DATATYPE;
#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352

class NiiLoader
{
public:
	int width, height, depth;
	float voxelWidth, voxelHeight, voxelDepth;

	NiiLoader();
	NiiLoader(const char* niiFile);
	void loadNii(const char* niiFile);

	vtkSmartPointer<vtkImageData> getData();

	void clearData();

	std::vector<double>& extractSkullVertex(int v, int t, int& num);
	void getUnique(vector<signed short> &u);
	vtkSmartPointer< vtkPolyData > mcSkullVertex(vector<double> &T_Pts, vector<double> &T_Norms, double isovalue, int sample_ratio = 3);

private:
	vtkSmartPointer<vtkImageData> niiImg;
	vector<double> skullVoxel;
	vector<double> skullVoxelForIso;
	vector<signed short> skullVoxelVal;
	//float *fArray;
	inline bool isSkullVoxel(signed short intensity)
	{
		if (intensity >= 0) return true;
		else return false;
	};
};

#endif