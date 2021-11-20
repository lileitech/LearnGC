/**************************************************************************
Program:   ZXH PatchExtraction Software
Author: lei li
Module:    ......   $
Language:  C++
Date:      $Date: From  2018-02-05 $
Version:   $Revision: v 1.0 $
**************************************************************************/

#include <string.h>
#include <iostream> 
#include <time.h> 
#include <math.h>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "zxhImageGipl.h" 
#include "PatchExtractByWorldCoordinate.h"
#include "STLreader.h"
#include <vtkExtractEdges.h>

#include <io.h>  
#include <direct.h> //创建文件夹

using namespace std;

vtkSmartPointer<vtkIdList> GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int id);//寻找node的领域点函数

int main(int argc, char* argv[])
{
		string mainfold, Datafold, TargetImageName, LABloodPoolGaussianBlurName, MeshWallLabelName, ScarProbImageName, mesh_name;

		//mainfold = "J:\\lilei\\SJTU\\2018MedAI_ScarSeg\\Data\\TestData";
		//string  casename = "P01040_6";
		//int HalfPatchLength = 9;
		//int N = 7;

		mainfold = argv[1];
		string  casename = argv[2];
		int HalfPatchLength = atoi(argv[3]);//5
		int N = atoi(argv[4]);//2

		srand((unsigned)time(NULL));
		int SiglePatchSize = (N * 2 + 1)*(N * 2 + 1)*(HalfPatchLength * 2 + 1);
		int PatchInfo[4] = { N, HalfPatchLength, 0, 0 };

		string PathName = mainfold + "\\" + casename;
		TargetImageName = PathName + "\\enhanced.nii.gz";

		LABloodPoolGaussianBlurName = PathName + "\\LA_label_GauiisanBlur_M.nii.gz";
		MeshWallLabelName = PathName + "\\LA_MeshWallLabel_M.nii.gz";
		ScarProbImageName = PathName + "\\LA_MeshWall_422_prob_M.nii.gz";
		mesh_name = PathName + "\\LA_Mesh_M.stl";


		//Patch fold
		string PSize = to_string(2 * N + 1);
		string PLength = to_string(2 * HalfPatchLength + 1);
		string sep = "_";
		string PInfo = "Patch_" + PSize + sep + PSize + sep + PLength;
		string PatchFold = PathName + "\\" + PInfo + "\\";


		int ftyp = _access(PatchFold.c_str(), 0);
		if (0 != ftyp)  //if this is not a directory!  
		{
			const char *PatchFoldName = PatchFold.data();
			_mkdir(PatchFoldName);
		}

		string T_SaveFoldScar = PatchFold + "Patch_T_info.txt";
		string T_SaveFoldNorm = PatchFold + "Patch_T_info_norm.txt";

		//--------------------load image-----------------------------------------------------------
		zxhImageData SourceImage, LABloodPoolGaussianBlur, MeshWallLabel;
		zxhImageDataT<float> ScarProbImage;
		zxh::OpenImageSafe(&SourceImage, TargetImageName);
		zxh::OpenImageSafe(&LABloodPoolGaussianBlur, LABloodPoolGaussianBlurName);
		// wall label (>=420) generated from mesh and WHS label which doesnot have mitral valve boundary, not only from WHS label @Li Lei 
		zxh::OpenImageSafe(&MeshWallLabel, MeshWallLabelName);
		zxh::OpenImageSafe(&ScarProbImage, ScarProbImageName);

		//--------------------load mesh-----------------------------------------------------------
		STLreader *stlreader = new STLreader(mesh_name);
		vtkSmartPointer<vtkPolyData> LAMesh = stlreader->decimated;	//load in the mesh

		//为了求取mesh相连的边，先把边萃取出来
		vtkSmartPointer<vtkExtractEdges> extractEdges = vtkSmartPointer<vtkExtractEdges>::New();
		extractEdges->SetInputData(LAMesh);
		extractEdges->Update();
		vtkSmartPointer<vtkPolyData> mesh = extractEdges->GetOutput();

		const int iNumOfMeshPoints = LAMesh->GetNumberOfPoints(); //the number of mesh point
		const int iNumOfMeshSurfaces = LAMesh->GetNumberOfPolys(); //the number of mesh surface
		int iNumOfMeshEdges = 0;
		for (int ptId = 0; ptId < iNumOfMeshPoints; ptId++)
		{
			vtkSmartPointer<vtkIdList> connectedVertices = GetConnectedVertices(mesh, ptId);
			int NumClist = connectedVertices->GetNumberOfIds();
			iNumOfMeshEdges += NumClist;
		}
		// save node weight info txt files
		ofstream outfile_scar(T_SaveFoldScar + "_____temp.txt", ios::beg);// first line number of nodes;  
		outfile_scar << -1 << "\n"; // the number of saved nodes is not known yet
		ofstream outfile_norm(T_SaveFoldNorm + "_____temp.txt", ios::beg);// first line number of nodes; 
		outfile_norm << -1 << "\n";

		//-------------------------------------------------------------------------------

		const int * Size = SourceImage.GetImageSize();

		//-----------Generate Patch Image with the real patch size ----------------------
		float spacing111[] = { 1, 1, 1, 1 }; // mm
		//定义patch library的size
		int T_newsize[] = { SiglePatchSize, 30000, iNumOfMeshPoints / 30000 + 1, 1 };//TO DO


		zxhImageData T_PatchImage1, T_PatchImage2, T_PatchImage3, T_PatchImage4;
		T_PatchImage1.NewImage(3, T_newsize, spacing111, SourceImage.GetImageInfo());
		T_PatchImage2.NewImage(3, T_newsize, spacing111, SourceImage.GetImageInfo());
		T_PatchImage3.NewImage(3, T_newsize, spacing111, SourceImage.GetImageInfo());
		T_PatchImage4.NewImage(3, T_newsize, spacing111, SourceImage.GetImageInfo());

		zxhImageData* ImageArray_T1[6] = { &SourceImage, &LABloodPoolGaussianBlur, &T_PatchImage1, &T_PatchImage2, &T_PatchImage3, &T_PatchImage4 };

		int nEdgeSaved = 0;
		int nNodeSaved = 0;
		for (int ptId = 0; ptId < iNumOfMeshPoints; ptId++)
		{
			float MeshNode_P2I_Coor[] = { LAMesh->GetPoint(ptId)[0], LAMesh->GetPoint(ptId)[1], LAMesh->GetPoint(ptId)[2], 0 };
			SourceImage.GetImageInfo()->PhysicalToImage(MeshNode_P2I_Coor);//物理坐标转成图像坐标
			int scx = zxh::round(MeshNode_P2I_Coor[0]);
			int scy = zxh::round(MeshNode_P2I_Coor[1]);
			int scz = zxh::round(Size[2] - MeshNode_P2I_Coor[2]);

			bool bNonInsterestNode = SourceImage.InsideImage(scx, scy, scz, 0) == false;
			if (bNonInsterestNode == false && MeshWallLabel.GetPixelGreyscaleClosest(scx, scy, scz) < 420)
				bNonInsterestNode = true;
			if (bNonInsterestNode)
				continue; // currently non interested node will NOT be stored for test 

			if (MeshWallLabel.GetPixelGreyscaleClosest(scx, scy, scz) < 420)
				continue; //this node is Mitral valve, ignor (notice, PV boundaries are considered as normal myo)

			float ScarProbValue = ScarProbImage.GetPixelGreyscaleClosest(scx, scy, scz, 0);

			//------------------------------------Generate mesh node patch-----------------------------------------------------	
			float InputWorldCoord[] = { MeshNode_P2I_Coor[0], MeshNode_P2I_Coor[1], Size[2] - MeshNode_P2I_Coor[2], 0 };
			SourceImage.GetImageInfo()->ImageToWorld(InputWorldCoord);
			if (PatchExtractByWorldCoordinate::GeneratePatchExtractByWorldWithRandOffset_test(InputWorldCoord, PatchInfo, 0, ImageArray_T1, nNodeSaved) == false)
				continue;
			nNodeSaved++;
			//Save t-link weight: nodeid T_weight for each line
			outfile_scar << ptId << " " << ScarProbValue << "\n";
			outfile_norm << ptId << " " << 1 - ScarProbValue << "\n";
			//------------------------------------------------------------------------------------------------------------------------------------------
		}
		outfile_scar.close();
		outfile_norm.close();


		//-----------------------correct the number of nodes/edges in the info text files --------------------@Li Lei: 请参考这些代码，看如何快速从txt中读取数字
		outfile_scar.open(T_SaveFoldScar, ios::beg);
		outfile_norm.open(T_SaveFoldNorm, ios::beg);
		outfile_scar << nNodeSaved << "\n";
		outfile_norm << nNodeSaved << "\n";

		std::ifstream ifs_scar, ifs_norm, ifs_nlink;
		const int static_buffer_size_of_base = 1024;
		char buffer[static_buffer_size_of_base];
		std::string sLine, sContent, sComment;
		int nodeid, cellid; float weight; float currDist;

		ifs_scar.open(T_SaveFoldScar + "_____temp.txt", std::ios_base::in);
		ifs_norm.open(T_SaveFoldNorm + "_____temp.txt", std::ios_base::in);

		ifs_scar.getline(buffer, static_buffer_size_of_base); // ignor first line
		ifs_norm.getline(buffer, static_buffer_size_of_base);


		// 1 T link weights
		for (int id = 0; id < nNodeSaved; ++id)
		{
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
			outfile_scar << nodeid << " " << weight << "\n";
			//normal
			ifs_norm.getline(buffer, static_buffer_size_of_base);
			sLine = buffer;
			zxh::ParseStringLine(sContent, sComment, sLine);
			zxh::trim_both(sContent);
			if (sContent.empty() || ifs_norm.fail())
			{
				std::cout << "error: please find it code-2;\n";
				return 2;
			}
			std::istringstream istr2(sContent);
			istr2 >> nodeid;
			istr2 >> weight;
			outfile_norm << nodeid << " " << weight << "\n";
		}
		outfile_scar.close();
		outfile_norm.close();
		ifs_scar.close();
		ifs_norm.close();

		//-----------------------save patch with PatchInfo--------------------

		zxh::SaveImage(&T_PatchImage1, PatchFold + "Patch_T_scale1.nii.gz");
		zxh::SaveImage(&T_PatchImage2, PatchFold + "Patch_T_scale2.nii.gz");
		zxh::SaveImage(&T_PatchImage3, PatchFold + "Patch_T_scale3.nii.gz");
		zxh::SaveImage(&T_PatchImage4, PatchFold + "Patch_T_scale4.nii.gz");


	return 0;
}


vtkSmartPointer<vtkIdList> GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int id)
{
	vtkSmartPointer<vtkIdList> connectedVertices =
		vtkSmartPointer<vtkIdList>::New();

	//get all cells that vertex 'id' is a part of
	vtkSmartPointer<vtkIdList> cellIdList =
		vtkSmartPointer<vtkIdList>::New();
	mesh->GetPointCells(id, cellIdList);

	for (vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
	{

		vtkSmartPointer<vtkIdList> pointIdList =
			vtkSmartPointer<vtkIdList>::New();
		mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);

		if (pointIdList->GetId(0) != id)
		{

			connectedVertices->InsertNextId(pointIdList->GetId(0));
		}
		else
		{
			connectedVertices->InsertNextId(pointIdList->GetId(1));
		}
	}

	return connectedVertices;
}

