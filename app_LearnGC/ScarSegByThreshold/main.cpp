
#include <iostream> 
#include <string>   
#include <iomanip> 
#include <sstream>  

#include "STLreader.h"

#include "zxhImageGipl.h" 
#include "zxhImageModelingLinear.h"
#include <time.h> 
#include <math.h>


#include <cstdlib>  
#include <limits>  
#include <cmath> 
#include <vtkExtractEdges.h>
#include <fstream>


using namespace std;


zxhImageData SourceImage, LAlabelImg, segMask, MeshWallLabel, LAGaussianBlurLabel;
zxhImageDataT<float> probimg;
vtkSmartPointer<vtkPolyData> LAMesh;


int main(int argc, char* argv[])
{
	string mainfold, casename, mesh_name, image_name, strSaveFilename, LABloodPoolGaussianBlurName, MeshWallLabelName;
	
	mainfold = argv[1];
	casename = argv[2];
	int HalfPatchLength = atoi(argv[3]);//5
	int N = atoi(argv[4]);//2

	///*mainfold = "E:\\LA_Segmentation\\TestData\\CASA_AF_May20Data_43_22\\";*/
	//mainfold = "J:\\lilei\\SJTU\\LA_Segmentation\\TestData\\TwoTimePointData_20_10\\";
	//casename = "P01025_6";
	//int HalfPatchLength = 5;
	//int N = 2;

	string PathName = mainfold + casename;
	mesh_name = PathName + "\\LA_Mesh.stl";
	image_name = PathName + "\\enhanced.nii.gz";
	LABloodPoolGaussianBlurName = PathName + "\\LA_label_GauiisanBlur.nii.gz";
	MeshWallLabelName = PathName + "\\LA_MeshWallLabel_fixed.nii.gz";

	//patch的fold
	string PSize = to_string(2 * N + 1);
	string PLength = to_string(2 * HalfPatchLength + 1);
	string sep = "_";
	string PInfo = "Patch_" + PSize + sep + PSize + sep + PLength;
	string PatchFold = PathName + "\\" + PInfo + "\\";

	strSaveFilename = PatchFold + "\\ScarSegThreshold.nii.gz";

	string T_SaveFoldScar = PatchFold + "Patch" + PSize + sep + PSize + sep + PLength + "_T_info_Predict.txt";
	string T_SaveFoldNorm = PatchFold + "Patch" + PSize + sep + PSize + sep + PLength + "_T_info_norm_Predict.txt";
	string N_SaveFold = PatchFold + "Patch" + PSize + sep + PSize + sep + PLength + "_N_info_Predict.txt";


	STLreader *stlreader = new STLreader(mesh_name);
	vtkSmartPointer<vtkPolyData> LAMesh = stlreader->decimated;	//load in the mesh

	zxh::OpenImageSafe(&SourceImage, image_name);
	zxh::OpenImageSafe(&MeshWallLabel, MeshWallLabelName);


	segMask.NewImage(SourceImage.GetImageInfo());// this is where we store the results
	const int * Size = SourceImage.GetImageSize();
	const int iNumOfMeshPoints = LAMesh->GetNumberOfPoints();

	//-----------------读取DL训练的nlink weight， fg t-linkweight----------------------------
	std::ifstream ifs_scar, ifs_norm, ifs_nlink;
	const int static_buffer_size_of_base = 1024;
	char buffer[static_buffer_size_of_base];
	std::string sLine, sContent, sComment;
	int nodeid, cellid; float weight;

	ifs_scar.open(T_SaveFoldScar, std::ios_base::in);
	ifs_scar.getline(buffer, static_buffer_size_of_base); // ignor first line
	

	//-----------------------------------------------------------------------------------------------

	for (int ptId = 0; ptId < iNumOfMeshPoints; ptId++)
	{
		float MeshNode_P2I_Coor[] = { LAMesh->GetPoint(ptId)[0], LAMesh->GetPoint(ptId)[1], LAMesh->GetPoint(ptId)[2], 0 };
		SourceImage.GetImageInfo()->PhysicalToImage(MeshNode_P2I_Coor);//物理坐标转成图像坐标
		int scx = zxh::round(MeshNode_P2I_Coor[0]);
		int scy = zxh::round(MeshNode_P2I_Coor[1]);
		int scz = zxh::round(Size[2] - MeshNode_P2I_Coor[2]);

		bool bNonInsterestNode = SourceImage.InsideImage(scx, scy, scz, 0) == false;
		if (bNonInsterestNode == false && MeshWallLabel.GetPixelGreyscaleClosest(scx, scy, scz)<420)
			bNonInsterestNode = true;
		if (bNonInsterestNode)
			continue; // currently non interested node will NOT be stored for test 

		if (MeshWallLabel.GetPixelGreyscaleClosest(scx, scy, scz) < 420)
			continue; //this node is Mitral valve, ignor (notice, PV boundaries are considered as normal myo)

		//------------------------------------------load tlink weight----------------------------------------------------
		//scar
		ifs_scar.getline(buffer, static_buffer_size_of_base);
		sLine = buffer;
		zxh::ParseStringLine(sContent, sComment, sLine);
		zxh::trim_both(sContent);
		if (sContent.empty() || ifs_scar.fail())
		{
			std::cout << "error: please find it code-1;\n";
			return 1;
		}
		std::istringstream istr(sContent);
		istr >> nodeid;
		istr >> weight;
		float SourceStrength = weight; //--------------output source weight from txt file

		if (MeshWallLabel.GetPixelGreyscaleClosest(scx, scy, scz)<420)
		{
			segMask.SetPixelByGreyscale(scx, scy, scz, 0, 0);
		}

		//if it is foreground 
		if (SourceStrength >= 0.5)
		{
			segMask.SetPixelByGreyscale(scx, scy, scz, 0, 422);
		}
		// if it is background 
		else
		{
			segMask.SetPixelByGreyscale(scx, scy, scz, 0, 421);
		}
	}


	zxh::SaveImage(&segMask, strSaveFilename); // save label image

	return 0;
}

