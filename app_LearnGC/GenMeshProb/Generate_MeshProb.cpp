#include <string.h>
#include <iostream> 
#include <time.h> 
#include <math.h>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>


#include "STLreader.h"
#include "zxhImageGipl.h" 
#include "zxhImageModelingLinear.h"

using namespace std;
vtkSmartPointer<vtkIdList> GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int id);

int main(int argc, char* argv[])
{
	string mainfold, Datafold, TargetImageName, MeshWallName, mesh_name, MeshWallLabelName, \
		WallDistName, ScarDistName, ScarProbImageName;

	mainfold = argv[1];
	string  casename = argv[2];
	string PathName = mainfold + "\\" + casename;
	TargetImageName = PathName + "\\enhanced.nii.gz";
	mesh_name = PathName + "\\LA_Mesh_M.stl";
	MeshWallLabelName = PathName + "\\LA_MeshWallLabel_M.nii.gz";
	WallDistName = PathName + "\\LA_MeshWall_421_dismap_M.nii.gz";
	ScarDistName = PathName + "\\LA_MeshWall_422_dismap_M.nii.gz";
	
	ScarProbImageName = PathName + "\\LA_MeshWall_422_prob_M.nii.gz";




	zxhImageData SourceImage, MeshWallLabelImage;
	zxhImageDataT<float>WallDist, ScarDist, ScarProbImage;
	zxh::OpenImageSafe(&SourceImage, TargetImageName);
	zxh::OpenImageSafe(&MeshWallLabelImage, MeshWallLabelName);	
	zxh::OpenImageSafe(&WallDist, WallDistName);
	zxh::OpenImageSafe(&ScarDist, ScarDistName);

	ScarProbImage.NewImage(MeshWallLabelImage.GetImageInfo());

	const int * Size = SourceImage.GetImageSize();
	for (int scz = 0; scz < Size[2]; scz++)
	{
		for (int scy = 0; scy < Size[1]; scy++)
		{
			for (int scx = 0; scx < Size[0]; scx++)
			{
				ScarProbImage.SetPixelByGreyscale(scx, scy, scz, 0,-1);
			}
		}
	}

	STLreader *stlreader = new STLreader(mesh_name);
	vtkSmartPointer<vtkPolyData> LAMesh = stlreader->decimated;	//load in the mesh
	const int iNumOfMeshPoints = LAMesh->GetNumberOfPoints(); //the number of mesh point


	float ProbValue;
	for (int ptId = 0; ptId < iNumOfMeshPoints; ptId++)
	{
		float MeshNode_P2I_Coor[] = { LAMesh->GetPoint(ptId)[0], LAMesh->GetPoint(ptId)[1], LAMesh->GetPoint(ptId)[2], 0 };
		SourceImage.GetImageInfo()->PhysicalToImage(MeshNode_P2I_Coor);//物理坐标转成图像坐标
		int scx = zxh::round(MeshNode_P2I_Coor[0]);
		int scy = zxh::round(MeshNode_P2I_Coor[1]);
		int scz = zxh::round(Size[2] - MeshNode_P2I_Coor[2]);		
		bool bIsInsideImage = SourceImage.InsideImage(scx, scy, scz, 0);
		if (!bIsInsideImage)
		{
			std::cout<<"error: node "<<ptId<<" not inside image\n";
			continue;
		} 

		//------------------------------Get New input world coordinate-------------------------------- 
		float WallDistValue = WallDist.GetPixelGreyscaleClosest(scx, scy, scz, 0);
		float ScarDistValue = ScarDist.GetPixelGreyscaleClosest(scx, scy, scz, 0);
		if ((WallDistValue > 0) && (ScarDistValue < 0))
		{
			ProbValue = 1 - 1 / (1 + exp(ScarDistValue));
		}

		else if ((ScarDistValue > 0) && (WallDistValue < 0))
		{
			ProbValue = 1 / (1 + exp(WallDistValue));
		}
		else
		{
			ProbValue = -1;//this node is mitral valve
		}

		ScarProbImage.SetPixelByGreyscale(scx, scy, scz, 0, ProbValue);



	}

	const char *SegResultName = ScarProbImageName.data();
	zxh::SaveImage(&ScarProbImage, SegResultName);

	return 1;
}

