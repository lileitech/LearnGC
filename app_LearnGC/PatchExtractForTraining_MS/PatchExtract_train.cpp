/**************************************************************************
Program:   ZXH PatchExtraction Software
Author: lei li
Module:    ......   $
Aim:    for traning   $
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

	//mainfold = "C:\\lilei\\2018MedAI_ScarSeg\\Data\\repeat_exp";
	//string  casename = "C02_Post";
	//int HalfPatchLength = 8;
	//int N = 6;
	//int Offset = 8;

	mainfold = argv[1];
	string  casename = argv[2];
	int HalfPatchLength = atoi(argv[3]);//5
	int N = atoi(argv[4]);//2
	int Offset = atoi(argv[5]);

	srand((unsigned)time(NULL));
	int SiglePatchSize = (N * 2 + 1)*(N * 2 + 1)*(HalfPatchLength * 2 + 1);

	string PathName = mainfold + "\\" + casename;
	TargetImageName = PathName + "\\enhanced.nii.gz";
	LABloodPoolGaussianBlurName = PathName + "\\LA_label_GauiisanBlur.nii.gz";
	MeshWallLabelName = PathName + "\\LA_MeshWallLabel_fixed.nii.gz";
	ScarProbImageName = PathName + "\\LA_MeshWall_422_prob_fixed.nii.gz";
	mesh_name = PathName + "\\LA_Mesh.stl";

	//string PathName = mainfold + "\\" + casename;
	//TargetImageName = PathName + "\\enhanced.nii.gz";
	//LABloodPoolGaussianBlurName = PathName + "\\LA_label_GauiisanBlur_M.nii.gz";
	//MeshWallLabelName = PathName + "\\LA_MeshWallLabel_M.nii.gz";
	//ScarProbImageName = PathName + "\\LA_MeshWall_422_prob_M.nii.gz";
	//mesh_name = PathName + "\\LA_Mesh_M.stl";

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
	string N_SaveFold = PatchFold + "Patch_N_info.txt";

	//--------------------load image-----------------------------------------------------------
	zxhImageData SourceImage, LABloodPoolGaussianBlur, MeshWallLabel;
	zxhImageDataT<float> ScarProbImage;
	zxh::OpenImageSafe(&SourceImage, TargetImageName);
	zxh::OpenImageSafe(&LABloodPoolGaussianBlur, LABloodPoolGaussianBlurName);
	zxh::OpenImageSafe(&MeshWallLabel, MeshWallLabelName);
	zxh::OpenImageSafe(&ScarProbImage, ScarProbImageName);

	//--------------------load mesh-----------------------------------------------------------
	STLreader *stlreader = new STLreader(mesh_name);
	vtkSmartPointer<vtkPolyData> LAMesh = stlreader->decimated;	//load in the mesh
	const int iNumOfMeshPoints = LAMesh->GetNumberOfPoints(); //the number of mesh point


	//为了求取mesh相连的边，先把边萃取出来
	vtkSmartPointer<vtkExtractEdges> extractEdges = vtkSmartPointer<vtkExtractEdges>::New();
	extractEdges->SetInputData(LAMesh);
	extractEdges->Update();
	vtkSmartPointer<vtkPolyData> mesh = extractEdges->GetOutput();
	//-------------------------------------------------------------------------------

	const int * Size = SourceImage.GetImageSize();

	int ScarAndBoundaryNodesNum = 0;
	for (int ptId = 0; ptId < iNumOfMeshPoints; ptId++)
	{
		float MeshNode_P2I_Coor[] = { LAMesh->GetPoint(ptId)[0], LAMesh->GetPoint(ptId)[1], LAMesh->GetPoint(ptId)[2], 0 };
		SourceImage.GetImageInfo()->PhysicalToImage(MeshNode_P2I_Coor);//物理坐标转成图像坐标
		int scx = zxh::round(MeshNode_P2I_Coor[0]);
		int scy = zxh::round(MeshNode_P2I_Coor[1]);
		int scz = zxh::round(Size[2] - MeshNode_P2I_Coor[2]);
		bool bIsInsideImage = SourceImage.InsideImage(scx, scy, scz, 0); // 超过图像边界的，不给予考虑，也就是说，默认为normal myo
		if (!bIsInsideImage)
		{
			std::cout << "warning: node " << ptId << " not inside image " << scx << "," << "\n"; //-----------------
			continue;
		}

		float ScarProbValue = ScarProbImage.GetPixelGreyscaleClosest(scx, scy, scz, 0);

		// ------------------------------------go over the mesh neighbors------------------------------------
		vtkSmartPointer<vtkIdList> connectedVertices = GetConnectedVertices(mesh, ptId);
		int NumClist = connectedVertices->GetNumberOfIds();
		int cellId = rand() % NumClist;//rand 
		int nNodeId = connectedVertices->GetId(cellId);
		float NMeshNode_P2I_Coor[] = { LAMesh->GetPoint(nNodeId)[0], LAMesh->GetPoint(nNodeId)[1], LAMesh->GetPoint(nNodeId)[2], 0 };
		SourceImage.GetImageInfo()->PhysicalToImage(NMeshNode_P2I_Coor);//物理坐标转成图像坐标
		int si = zxh::round(NMeshNode_P2I_Coor[0]);
		int sj = zxh::round(NMeshNode_P2I_Coor[1]);
		int sk = zxh::round(Size[2] - NMeshNode_P2I_Coor[2]);
		bool N_bIsInsideImage = SourceImage.InsideImage(si, sj, sk, 0);
		if (!N_bIsInsideImage)
		{
			std::cout << "warning: node not inside image " << si << "," << "\n"; //-----------------
			continue;
		}
		float N_ScarProbValue = ScarProbImage.GetPixelGreyscaleClosest(si, sj, sk, 0);
		if ((ScarProbValue < 0) || (N_ScarProbValue < 0))
			continue;//this node is PV or Maitral valve

		//---------------lei: using float data to calculate dist
		float wCurrNode[] = { MeshNode_P2I_Coor[0], MeshNode_P2I_Coor[1], Size[2] - MeshNode_P2I_Coor[2], 0 },
			wNeiNode[] = { NMeshNode_P2I_Coor[0], NMeshNode_P2I_Coor[1], Size[2] - NMeshNode_P2I_Coor[2], 0 };
		SourceImage.ImageToWorld(wCurrNode);
		SourceImage.ImageToWorld(wNeiNode);
		float currDist = zxh::VectorOP_Distance(wCurrNode, wNeiNode, 3);
		if (currDist > 5)
			continue;


		if (ScarProbValue > 0.1) // || N_ScarProbValue > 0.1 )  first node is scar 
			ScarAndBoundaryNodesNum++;
	}

	if (ScarAndBoundaryNodesNum > 20000)
	{
		ScarAndBoundaryNodesNum = 20000;//控制scar node的点的数目不要超过20000，因为nii img的最大存储就是30000
	}

	int TraningNodeNum = 1.5*ScarAndBoundaryNodesNum;//TO DO
	int Num_Scar = ScarAndBoundaryNodesNum; //prob>0.1
	int Num_Normal = 0.5*ScarAndBoundaryNodesNum; //prob <= 0.1 && prob>=0

	//-----------Generate Patch Image with the real patch size ----------------------
	float spacing111[] = { 1, 1, 1, 1 }; // mm
	//定义patch library的size
	int T_newsize[] = { SiglePatchSize, TraningNodeNum, 1, 1 };
	int N_newsize[] = { SiglePatchSize, TraningNodeNum, 1, 1 };
	//-------------------------------------------------------------------------------

	int ScarNumIdex = 0;
	int NormalNumIdex = 0;
	int iOutliers = 0;
	//--------------------------------rand mesh node-------------------------------------- 
	ofstream outfile_scar(T_SaveFoldScar, ios::beg);//output
	outfile_scar << TraningNodeNum << "\n";
	ofstream outfile_norm(T_SaveFoldNorm, ios::beg);//output
	outfile_norm << TraningNodeNum << "\n";

	zxhImageData T_PatchImage1, T_PatchImage2, T_PatchImage3;
	T_PatchImage1.NewImage(2, T_newsize, spacing111, SourceImage.GetImageInfo());
	T_PatchImage2.NewImage(2, T_newsize, spacing111, SourceImage.GetImageInfo());
	T_PatchImage3.NewImage(2, T_newsize, spacing111, SourceImage.GetImageInfo());



	zxhImageData* ImageArray_T1[5] = { &SourceImage, &LABloodPoolGaussianBlur, &T_PatchImage1, &T_PatchImage2, &T_PatchImage3 };


	for (int ptId = 0; (ScarNumIdex < Num_Scar) || (NormalNumIdex < Num_Normal); ptId++)
	{
		if ((ScarNumIdex > ScarAndBoundaryNodesNum) || (ptId >= iNumOfMeshPoints))
		{
			ptId = rand() % iNumOfMeshPoints;//rand for normal
		}
		float MeshNode_P2I_Coor[] = { LAMesh->GetPoint(ptId)[0], LAMesh->GetPoint(ptId)[1], LAMesh->GetPoint(ptId)[2], 0 };
		SourceImage.GetImageInfo()->PhysicalToImage(MeshNode_P2I_Coor);//物理坐标转成图像坐标
		int scx = zxh::round(MeshNode_P2I_Coor[0]);
		int scy = zxh::round(MeshNode_P2I_Coor[1]);
		int scz = zxh::round(Size[2] - MeshNode_P2I_Coor[2]);
		bool bIsInsideImage = SourceImage.InsideImage(scx, scy, scz, 0); // 超过图像边界的，不给予考虑，也就是说，默认为normal myo
		if (!bIsInsideImage)
		{
			std::cout << "warning: node " << ptId << " not inside image " << scx << "," << "\n"; //-----------------
			continue;
		}

		float ScarProbValue = ScarProbImage.GetPixelGreyscaleClosest(scx, scy, scz, 0);
		if (ScarProbValue < 0)//this node is Mitral valve, ignor (notice, PV boundaries are considered as normal myo)
			continue;

		// ------------------------------------go over the mesh neighbors------------------------------------
		vtkSmartPointer<vtkIdList> connectedVertices = GetConnectedVertices(mesh, ptId);
		int NumClist = connectedVertices->GetNumberOfIds();
		int cellId = rand() % NumClist;//rand 
		int nNodeId = connectedVertices->GetId(cellId);
		float NMeshNode_P2I_Coor[] = { LAMesh->GetPoint(nNodeId)[0], LAMesh->GetPoint(nNodeId)[1], LAMesh->GetPoint(nNodeId)[2], 0 };
		SourceImage.GetImageInfo()->PhysicalToImage(NMeshNode_P2I_Coor);//物理坐标转成图像坐标
		int si = zxh::round(NMeshNode_P2I_Coor[0]);
		int sj = zxh::round(NMeshNode_P2I_Coor[1]);
		int sk = zxh::round(Size[2] - NMeshNode_P2I_Coor[2]);
		bool N_bIsInsideImage = SourceImage.InsideImage(si, sj, sk, 0);
		if (!N_bIsInsideImage)
		{
			std::cout << "warning: node not inside image " << si << "," << "\n"; //-----------------
			continue;
		}
		float N_ScarProbValue = ScarProbImage.GetPixelGreyscaleClosest(si, sj, sk, 0);
		if ((ScarProbValue < 0) || (N_ScarProbValue < 0))
			continue;//this node is PV or Maitral valve

		//---------------lei: using float data to calculate dist
		float wCurrNode[] = { MeshNode_P2I_Coor[0], MeshNode_P2I_Coor[1], Size[2] - MeshNode_P2I_Coor[2], 0 },
			wNeiNode[] = { NMeshNode_P2I_Coor[0], NMeshNode_P2I_Coor[1], Size[2] - NMeshNode_P2I_Coor[2], 0 };
		SourceImage.ImageToWorld(wCurrNode);
		SourceImage.ImageToWorld(wNeiNode);
		float currDist = zxh::VectorOP_Distance(wCurrNode, wNeiNode, 3);
		if (currDist > 5)
		{
			iOutliers++;
			//std::cout << "error: dist:" << currDist << " is too large for node " << ptId << "\n";
			continue;
		}

		//Save n-link weight 
		bool bIsScarNode = false;
		if (ScarProbValue > 0.1) // || N_ScarProbValue > 0 )  first node is scar 
			bIsScarNode = true;

		if (ScarNumIdex < Num_Scar && bIsScarNode == false)//scar点没有取满，先就全部不跑normal点
			continue;

		if (ScarNumIdex == Num_Scar && bIsScarNode == true)
			continue;
		if (ScarNumIdex == Num_Normal && bIsScarNode == false)
			continue;


		//---------------lei: using float data to extract patch
		float InputWorldCoordRand[] = { MeshNode_P2I_Coor[0], MeshNode_P2I_Coor[1], Size[2] - MeshNode_P2I_Coor[2], 0 };
		SourceImage.GetImageInfo()->ImageToWorld(InputWorldCoordRand);

		//------------------------------------Generate neighboor node patch-----------------------------------------------------	
		int PatchInfo[4] = { N, HalfPatchLength, 0, 0 };//kk保证一样

		float randoffset = 0;
		if (Offset > 0)
		{
			randoffset = float((rand() % (200 * Offset)) - 100 * Offset) / 100.0;//rand shift	
		}
		if (PatchExtractByWorldCoordinate::GeneratePatchExtractByWorldWithRandOffset(InputWorldCoordRand, PatchInfo, randoffset, ImageArray_T1, ScarNumIdex + NormalNumIdex) == false)
			continue;//for t-link



		outfile_scar << ptId << " " << ScarProbValue << "\n";
		outfile_norm << ptId << " " << 1.0 - ScarProbValue << "\n";

		if (bIsScarNode == true)
			ScarNumIdex++;
		else //if ( bIsScarNode==false )
			NormalNumIdex++;

	}

	outfile_scar.close();
	outfile_norm.close();

	//-----------------------save patch with PatchInfo--------------------


	zxh::SaveImage(&T_PatchImage1, PatchFold + "Patch_T_scale1.nii.gz");
	zxh::SaveImage(&T_PatchImage2, PatchFold + "Patch_T_scale2.nii.gz");
	zxh::SaveImage(&T_PatchImage3, PatchFold + "Patch_T_scale3.nii.gz");



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

