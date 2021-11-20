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


int main(int argc, char* argv[])
{
	string mainfold, Datafold, TargetImageName, TargetLabelName, MeshWallName, mesh_name, LABloodPoolGaussianBlurName, \
		LALabelName, GoldLALabelName, MeshWallLabelName;

	mainfold = argv[1];
	string  casename = argv[2];

	//mainfold = "E:\\2018MedAI_ScarSeg\\Data\\TwoTimePointData";
	//string  casename = "P01025_6";

	string PathName = mainfold + "\\" + casename;
	TargetImageName = PathName + "\\enhanced.nii.gz";
	TargetLabelName = PathName + "\\en_seg_msp_M.nii.gz";
	mesh_name = PathName + "\\LA_Mesh_M.stl";
	MeshWallName = PathName + "\\LA_MeshWall_M.nii.gz";
	LABloodPoolGaussianBlurName = PathName + "\\LA_label_GauiisanBlur_M.nii.gz"; // never use LA_MeshWall_GauiisanBlur
	GoldLALabelName = PathName + "\\scarSegImgM.nii.gz";

	MeshWallLabelName = PathName + "\\LA_PatchLabel_M.nii.gz";



	zxhImageData SourceImage, SourceLabel, MeshWallImage, LABloodPoolGaussianBlur, GoldLALabel, MeshWallLabelImage;
	zxh::OpenImageSafe(&SourceImage, TargetImageName);
	zxh::OpenImageSafe(&SourceLabel, TargetLabelName);
	zxh::OpenImageSafe(&MeshWallImage, MeshWallName);
	zxh::OpenImageSafe(&LABloodPoolGaussianBlur, LABloodPoolGaussianBlurName);
	zxh::OpenImageSafe(&GoldLALabel, GoldLALabelName);


	MeshWallLabelImage.CloneFrom(&MeshWallImage);//copy from wall
	const int * Size = SourceImage.GetImageSize();

	STLreader *stlreader = new STLreader(mesh_name);
	vtkSmartPointer<vtkPolyData> LAMesh = stlreader->decimated;	//load in the mesh
	const int iNumOfMeshPoints = LAMesh->GetNumberOfPoints(); //the number of mesh point

	zxhImageModelingLinear GradientMod;
	GradientMod.SetImage(&LABloodPoolGaussianBlur); ///////-- 错误，不能使用wall做smooth后图像求法向量---------------
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
			std::cout << "error: node " << ptId << " not inside image\n";
			continue;
		}


		ZXHPixelTypeDefault MeshWallSurfaceValue = MeshWallImage.GetPixelGreyscaleClosest(scx, scy, scz, 0);
		if (MeshWallSurfaceValue == 0)//PV
			continue;

		//------------------------------Get New input world coordinate-------------------------------- 
		float InputWorldCoord[] = { MeshNode_P2I_Coor[0], MeshNode_P2I_Coor[1], Size[2] - MeshNode_P2I_Coor[2], 0 };
		SourceImage.ImageToWorld(InputWorldCoord);

		float  pwGrad[4] = { 0 };
		GradientMod.GetPixelGradientByWorld(pwGrad, InputWorldCoord[0], InputWorldCoord[1], InputWorldCoord[2], 0);
		float ix = pwGrad[0];
		float iy = pwGrad[1];
		float iz = pwGrad[2];
		float mag = sqrt(ix*ix + iy*iy + iz*iz);
		if (mag<ZXH_FloatInfinitesimal)
		{
			std::cout << "error: magnitude " << mag << " too small for node " << ptId << "\n";
			return -1;
		}
		float Ia = 0.5*ix / mag; // 每步只有0.5mm步长，防止错过scar厚度小于1mm的点
		float Ib = 0.5*iy / mag;
		float Ic = 0.5*iz / mag;

		int MVIndex = 0;
		int RAIndex = 0;
		int DOIndex = 0;
		for (int step = 0; step <= 15; step++) // 0.5mm for each step
		{
			// forward step
			float  NewInputWorldCoord[4] = { InputWorldCoord[0] + step*Ia, InputWorldCoord[1] + step*Ib, InputWorldCoord[2] + step*Ic, 0 };
			float NewInputImageCoord[] = { NewInputWorldCoord[0], NewInputWorldCoord[1], NewInputWorldCoord[2], 0 };
			SourceImage.WorldToImage(NewInputImageCoord);
			int cx = zxh::round(NewInputImageCoord[0]);
			int cy = zxh::round(NewInputImageCoord[1]);
			int cz = zxh::round(NewInputImageCoord[2]);
			bool bIsInsideImage = SourceImage.InsideImage(cx, cy, cz, 0);
			float SourceLabelValue;
			if (bIsInsideImage)
			{

				//remove PV and Mitral valve
				/*if ((LALabelValue == 1) || (LALabelValue == 500))*/ //remove PV and Mitral valve
				SourceLabelValue = SourceLabel.GetPixelGreyscaleClosest(cx, cy, cz, 0);
				if (SourceLabelValue == 500)
				{
					MVIndex++;//mitral valve
				}
				else if (SourceLabelValue == 550)
				{
					RAIndex++;//ra
				}
				if (SourceLabelValue == 820)
				{
					DOIndex++;//DO
				}
			}
			// backward step
			float  BackInputWorldCoord[4] = { InputWorldCoord[0] - step*Ia, InputWorldCoord[1] - step*Ib, InputWorldCoord[2] - step*Ic, 0 };
			float BackInputImageCoord[] = { BackInputWorldCoord[0], BackInputWorldCoord[1], BackInputWorldCoord[2], 0 };
			SourceImage.WorldToImage(BackInputImageCoord);
			cx = zxh::round(BackInputImageCoord[0]);
			cy = zxh::round(BackInputImageCoord[1]);
			cz = zxh::round(BackInputImageCoord[2]);
			bIsInsideImage = SourceImage.InsideImage(cx, cy, cz, 0);
			if (!bIsInsideImage)
				continue;

			
			SourceLabelValue = SourceLabel.GetPixelGreyscaleClosest(cx, cy, cz, 0);
			if (SourceLabelValue == 500)
			{
				MVIndex++;//mitral valve
			}
			else if (SourceLabelValue == 550)
			{
				RAIndex++;//RA
			}
			else if (SourceLabelValue == 820)
			{
				DOIndex++;//DO
			}


		}

		if (RAIndex > 0)
		{
			MeshWallLabelImage.SetPixelByGreyscale(scx, scy, scz, 0, 550);//RA
		}
		if (DOIndex > 0)
		{
			MeshWallLabelImage.SetPixelByGreyscale(scx, scy, scz, 0, 820);//DO
		} 
		if(MVIndex > 0)
		{
			MeshWallLabelImage.SetPixelByGreyscale(scx, scy, scz, 0, 500);//mitral valve
		}
		


	}

	zxh::SaveImage(&MeshWallLabelImage, MeshWallLabelName);

	return 0;
}

