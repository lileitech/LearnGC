#include <string.h>
#include <iostream> 

#include <time.h> 
#include <math.h>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "zxhImageGipl.h" 
#include "zxhImageModelingLinear.h"

using namespace std;


int main(int argc, char* argv[])
{
	string mainfold, Datafold, TargetImageName, LABloodPoolGaussianBlurName, \
		BloodLabelName,LALabelName, GoldLALabelName, MeshWallLabelName;

	//mainfold = argv[1];
	//string  casename = argv[2];
	//string wallname = argv[3];
	//string surfacename = argv[4];

	mainfold = "E:\\zxh_work\\MICCAI2018_ScarSeg_zxh\\LA_Segmentation";
	string  casename = "C45_Pre_V_2_4"; //C19_Post  C45_Pre  C36_Post
	string wallname = "method3comp_Sub3D_To_MultiImage0_scar.nii.gz";
	string surfacename = "method3comp_Sub3D_To_MultiImage0_scar_surface_blood.nii.gz";

	string PathName = mainfold + "\\" + casename;
	TargetImageName = PathName + "\\enhanced.nii.gz";
	LABloodPoolGaussianBlurName = PathName + "\\LA_label_GauiisanBlur.nii.gz"; // never use LA_MeshWall_GauiisanBlur
	LALabelName = PathName + "\\LA_label.nii.gz";
	BloodLabelName = PathName + "\\en_seg_msp_la.nii.gz";
	GoldLALabelName = PathName + "\\" + wallname;
	MeshWallLabelName = PathName + "\\" + surfacename;



	zxhImageData SourceImage,LABloodPoolGaussianBlur, GoldLALabel, LALabel, BloodLabel, MeshWallLabelImage;
	zxh::OpenImageSafe(&SourceImage, TargetImageName);
	zxh::OpenImageSafe(&LABloodPoolGaussianBlur, LABloodPoolGaussianBlurName);
	zxh::OpenImageSafe(&GoldLALabel, GoldLALabelName);
	zxh::OpenImageSafe(&LALabel, LALabelName);
	zxh::OpenImageSafe(&BloodLabel, BloodLabelName);



	MeshWallLabelImage.CloneFrom(&BloodLabel);//copy from wall
	const int * Size = SourceImage.GetImageSize();


	for (int scx = 0; scx < Size[0];scx++)
	{
		for (int scy = 0; scy < Size[0]; scy++)
		{
			for (int scz = 0; scz < Size[0]; scz++)
			{
				ZXHPixelTypeDefault BloodLabelValue = BloodLabel.GetPixelGreyscaleClosest(scx, scy, scz, 0);
				if (BloodLabelValue == 0)//PV
					continue;

				//------------------------------Get New input world coordinate-------------------------------- 
				float InputWorldCoord[] = { scx, scy, scz, 0 };
				SourceImage.ImageToWorld(InputWorldCoord);

				zxhImageModelingLinear GradientMod;
				GradientMod.SetImage(&LABloodPoolGaussianBlur); ///////-- 错误，不能使用wall做smooth后图像求法向量---------------
				float  pwGrad[4] = { 0 };
				GradientMod.GetPixelGradientByWorld(pwGrad, InputWorldCoord[0], InputWorldCoord[1], InputWorldCoord[2], 0);
				float ix = pwGrad[0];
				float iy = pwGrad[1];
				float iz = pwGrad[2];
				float mag = sqrt(ix*ix + iy*iy + iz*iz);
				/*if (mag<ZXH_FloatInfinitesimal)
				{
					std::cout << "error: magnitude " << mag << " too small"<< "\n";
					return -1;
				}*/
				float Ia = 0.5*ix / mag; // 每步只有0.5mm步长，防止错过scar厚度小于1mm的点
				float Ib = 0.5*iy / mag;
				float Ic = 0.5*iz / mag;

				int ScarIndex = 0;
				int PVIndex = 0;
				for (int step = 0; step <= 20; step++) // 0.5mm for each step
				{
					// forward step
					float  NewInputWorldCoord[4] = { InputWorldCoord[0] + step*Ia, InputWorldCoord[1] + step*Ib, InputWorldCoord[2] + step*Ic, 0 };
					float NewInputImageCoord[] = { NewInputWorldCoord[0], NewInputWorldCoord[1], NewInputWorldCoord[2], 0 };
					SourceImage.WorldToImage(NewInputImageCoord);
					int cx = zxh::round(NewInputImageCoord[0]);
					int cy = zxh::round(NewInputImageCoord[1]);
					int cz = zxh::round(NewInputImageCoord[2]);
					bool bIsInsideImage = SourceImage.InsideImage(cx, cy, cz, 0);
					if (!bIsInsideImage)
						continue;//放弃这样的点



					float GoldLALabelValue = GoldLALabel.GetPixelGreyscaleClosest(cx, cy, cz, 0);
					if (GoldLALabelValue > 0)
					{
						ScarIndex++;
						break;
					}
					//remove PV and Mitral valve
					float LALabelValue = LALabel.GetPixelGreyscaleClosest(cx, cy, cz, 0);
					if (LALabelValue == 500)
					{
						PVIndex++;
						break;
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

					GoldLALabelValue = GoldLALabel.GetPixelGreyscaleClosest(cx, cy, cz, 0);
					if (GoldLALabelValue > 0)
					{
						ScarIndex++;
						break;
					}

					LALabelValue = LALabel.GetPixelGreyscaleClosest(cx, cy, cz, 0);
					/*if ((LALabelValue == 1) || (LALabelValue == 500))*/ //remove PV and Mitral valve
					if (LALabelValue == 500) //remove Mitral valve
					{
						PVIndex++;
						break;
					}

				}

				if (PVIndex > 0)//只要出现PV,就赋为0
				{
					MeshWallLabelImage.SetPixelByGreyscale(scx, scy, scz, 0, 0);
				}

				else if ((ScarIndex > 0) && (PVIndex == 0))//只要出现scar，且不在PV处，就赋为422
				{
					MeshWallLabelImage.SetPixelByGreyscale(scx, scy, scz, 0, 422);
				}


				else if ((ScarIndex == 0) && (PVIndex == 0))//只要出现scar，就赋为421
				{
					MeshWallLabelImage.SetPixelByGreyscale(scx, scy, scz, 0, 421);
				}

			}
		}
	}

	const char *SegResultName = MeshWallLabelName.data();
	zxh::SaveImage(&MeshWallLabelImage, SegResultName);

	return 1;
}

