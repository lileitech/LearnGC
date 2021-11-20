
#include <iostream> 
#include <string>   
#include <iomanip> 
#include <sstream>  

#include "graph.h"
#include "ComputeTime.h"
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

#define pi 3.1415926
#define max(a,b) (((a) > (b)) ? (a) : (b))

using namespace std;
//using namespace cv;
const float INT32_CONST = 1000;
const float HARD_CONSTRAINT_CONST = 1000;

//index
bool OneCutIdex = false;
bool GraphCutIdex = true;

//float WeightT = 1;
//float WeightN = 1;

//prob修正领域范围
float NN = 4;
float r = 0.6;

//for onecut
int numUsedBins = 0;
float bha_slope = 1;
int numBinsPerChannel = 100;


zxhImageData SourceImage, LAlabelImg, segMask, MeshWallLabel, LAGaussianBlurLabel;
zxhImageDataT<float> probimg;
vtkSmartPointer<vtkPolyData> LAMesh;

int  init(char * imgFileName, vtkSmartPointer<vtkPolyData>LAMesh, string T_SaveFoldScar, string T_SaveFoldNorm, string N_SaveFold, float WeightN);
vtkSmartPointer<vtkIdList> GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int id);
void getBinPerPixel(zxhImageDataT<float> & binPerPixelImg, string T_SaveFoldScar, vtkSmartPointer<vtkPolyData>LAMesh, int numBinsPerChannel, int & numUsedBins);

typedef Graph<int, int, int> GraphType;
GraphType *myGraph;

int main(int argc, char* argv[])
{
	string mainfold, casename, mesh_name, image_name, strSaveFilename, LABloodPoolGaussianBlurName, MeshWallLabelName;
	
	mainfold = argv[1];
	casename = argv[2];
	int HalfPatchLength = atoi(argv[3]);//5
	int N = atoi(argv[4]);//2
	float WeightN = atof(argv[5]);


	//mainfold = "E:\\LA_Segmentation\\TestData\\TwoTimePointData_20_10\\";
	///*mainfold = "J:\\lilei\\SJTU\\LA_Segmentation\\TestData\\TwoTimePointData_20_10\\";*/
	//casename = "PX0001_6";
	//int HalfPatchLength = 5;
	//int N = 2;

	//string PathName = mainfold + casename;
	//mesh_name = PathName + "\\LA_Mesh.stl";
	//image_name = PathName + "\\enhanced.nii.gz";
	//LABloodPoolGaussianBlurName = PathName + "\\LA_label_GauiisanBlur.nii.gz";
	//MeshWallLabelName = PathName + "\\LA_MeshWallLabel_fixed.nii.gz";

	string PathName = mainfold + casename;
	mesh_name = PathName + "\\LA_Mesh_M.stl";
	image_name = PathName + "\\enhanced.nii.gz";
	LABloodPoolGaussianBlurName = PathName + "\\LA_label_GauiisanBlur_M.nii.gz";
	MeshWallLabelName = PathName + "\\LA_MeshWallLabel_M.nii.gz";

	//patch的fold
	string PSize = to_string(2 * N + 1);
	string PLength = to_string(2 * HalfPatchLength + 1);
	string sep = "_";
	string PInfo = "Patch_" + PSize + sep + PSize + sep + PLength;
	string PatchFold = PathName + "\\" + PInfo + "\\";

	strSaveFilename = PatchFold + "\\ScarSegGraphCut_M.nii.gz";

	
	//string T_SaveFoldScar = PatchFold + "Patch" + PSize + sep + PSize + sep + PLength + "_T_info_Predict.txt";
	//string T_SaveFoldNorm = PatchFold + "Patch" + PSize + sep + PSize + sep + PLength + "_T_info_norm_Predict.txt";
	//string N_SaveFold = PatchFold + "Patch" + PSize + sep + PSize + sep + PLength + "_N_info_Predict.txt";

	string T_SaveFoldScar = PatchFold + "Patch" + PSize + sep + PSize + sep + PLength + "_T_info_Predict.txt";
	string T_SaveFoldNorm = PatchFold + "Patch" + PSize + sep + PSize + sep + PLength + "_T_info_norm_Predict.txt";
	string N_SaveFold = PatchFold + "Patch" + PSize + sep + PSize + sep + PLength + "_N_info_Predict.txt";



	STLreader *stlreader = new STLreader(mesh_name);
	vtkSmartPointer<vtkPolyData> LAMesh = stlreader->decimated;	//load in the mesh

	char * imgFileName = (char *)image_name.c_str();
	char * WallFileName = (char *)MeshWallLabelName.c_str();

	zxh::OpenImageSafe(&MeshWallLabel, MeshWallLabelName);


	//-------------------------------------------------start init---------------------------------------
	ComputeTime ct_init;
	ct_init.Begin();
	if (init(imgFileName, LAMesh, T_SaveFoldScar, T_SaveFoldNorm, N_SaveFold, WeightN) == -1)
	{
		cout << "Could not initialize" << endl;
		return -1;
	}
	cout << "初始化运行时间：  " << ct_init.End() << "ms" << endl;


	//------------------------------------------------start to segmentation---------------------------------------
	ComputeTime ct;
	ct.Begin();
	cout << "setting the hard constraints..." << endl;
	cout << "maxflow..." << endl;

	int flow = myGraph->maxflow();
	cout << "done maxflow..." << endl;

	segMask.NewImage(SourceImage.GetImageInfo());// this is where we store the results
	const int * Size = SourceImage.GetImageSize();
	const int iNumOfMeshPoints = LAMesh->GetNumberOfPoints();
	for (size_t ptId = 0; ptId < iNumOfMeshPoints; ptId++)
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


		if (MeshWallLabel.GetPixelGreyscaleClosest(scx, scy, scz)<420)
		{
			segMask.SetPixelByGreyscale(scx, scy, scz, 0, 0);
		}

		//if it is foreground 
		else if (myGraph->what_segment((GraphType::node_id)ptId) == GraphType::SOURCE)
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

	cout << "运行时间：  " << ct.End() << "ms" << endl;
	return 0;
}

int init(char * imgFileName, vtkSmartPointer<vtkPolyData>LAMesh, string T_SaveFoldScar, string T_SaveFoldNorm, string N_SaveFold, float WeightN)
{
	zxhImageData img, goldsurfacelabel;
	zxhImageDataT<float> SaltProbImg;
	zxh::OpenImageSafe(&img, std::string(imgFileName));

	const int iNumOfMeshPoints = LAMesh->GetNumberOfPoints(); //the number of mesh point
	const int * Size = img.GetImageSize();
	SourceImage.CloneFrom(&img);

	// Check for invalid input
	if (SourceImage.IsEmpty())
	{
		cout << "Could not open or find the image: " << imgFileName << std::endl;
		return -1;
	}

	zxhImageDataT<float> binPerPixelImg;
	binPerPixelImg.NewImage(img.GetImageInfo());

	if (OneCutIdex)
	{
		getBinPerPixel(binPerPixelImg, T_SaveFoldScar, LAMesh, numBinsPerChannel, numUsedBins);

	}


	myGraph = new GraphType(/*estimated # of mesh nodes*/  iNumOfMeshPoints + numUsedBins,
		/*estimated # of edges=11 spatial neighbors and one link to auxiliary*/ 8 * iNumOfMeshPoints);

	GraphType::node_id ptId = myGraph->add_node((int)iNumOfMeshPoints + numUsedBins);


	//-----------------为了求取mesh相连的边，先把边萃取出来----------------------------
	vtkSmartPointer<vtkExtractEdges> extractEdges =
		vtkSmartPointer<vtkExtractEdges>::New();
	extractEdges->SetInputData(LAMesh);
	extractEdges->Update();
	vtkSmartPointer<vtkPolyData> mesh = extractEdges->GetOutput();
	//-------------------------------------------------------------------------------



	//-----------------读取DL训练的nlink weight， fg t-linkweight， bg t-linkweight----------------------------
	std::ifstream ifs_scar, ifs_norm, ifs_nlink;
	const int static_buffer_size_of_base = 1024;
	char buffer[static_buffer_size_of_base];
	std::string sLine, sContent, sComment;
	int nodeid, cellid; float weight; float dist;

	ifs_scar.open(T_SaveFoldScar, std::ios_base::in);
	ifs_norm.open(T_SaveFoldNorm, std::ios_base::in);
	ifs_nlink.open(N_SaveFold, std::ios_base::in);
	ifs_scar.getline(buffer, static_buffer_size_of_base); // ignor first line
	ifs_norm.getline(buffer, static_buffer_size_of_base);
	ifs_nlink.getline(buffer, static_buffer_size_of_base);
	
	//-----------------------------------------------------------------------------------------------

	for (ptId = 0; ptId < iNumOfMeshPoints; ptId++)
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
		float SinkStrength = weight;//--------------output sink weight from txt file



		// go over the mesh neighbors
		vtkSmartPointer<vtkIdList> connectedVertices = GetConnectedVertices(mesh, ptId);
		int NumClist = connectedVertices->GetNumberOfIds();
		GraphType::node_id cellId;
		for (cellId = 0; cellId <NumClist; cellId++)
		{
			GraphType::node_id nNodeId = connectedVertices->GetId(cellId);
			float NMeshNode_P2I_Coor[] = { LAMesh->GetPoint(nNodeId)[0], LAMesh->GetPoint(nNodeId)[1], LAMesh->GetPoint(nNodeId)[2], 0 };
			SourceImage.GetImageInfo()->PhysicalToImage(NMeshNode_P2I_Coor);//物理坐标转成图像坐标
			int si = zxh::round(NMeshNode_P2I_Coor[0]);
			int sj = zxh::round(NMeshNode_P2I_Coor[1]);
			int sk = zxh::round(Size[2] - NMeshNode_P2I_Coor[2]);
			bool N_bNonInsterestNode = SourceImage.InsideImage(si, sj, sk, 0) == false;
			if (N_bNonInsterestNode == false && MeshWallLabel.GetPixelGreyscaleClosest(si, sj, sk)<420)
				bNonInsterestNode = true;
			if (N_bNonInsterestNode)
				continue; // currently non interested node will NOT be stored for test 

			if (MeshWallLabel.GetPixelGreyscaleClosest(si, sj, sk) < 420)
				continue; //this node is Mitral valve, ignor (notice, PV boundaries are considered as normal myo)

			float wCurrNode[] = { MeshNode_P2I_Coor[0], MeshNode_P2I_Coor[1], Size[2] - MeshNode_P2I_Coor[2], 0 },\
				wNeiNode[] = { NMeshNode_P2I_Coor[0], NMeshNode_P2I_Coor[1], Size[2] - NMeshNode_P2I_Coor[2], 0 };
			SourceImage.ImageToWorld(wCurrNode);
			SourceImage.ImageToWorld(wNeiNode);
			float currDist = zxh::VectorOP_Distance(wCurrNode, wNeiNode, 3);
			if (currDist>5)
				N_bNonInsterestNode = true;
			if (currDist < ZXH_FloatInfinitesimal)
			{
			std:cerr << "error: distance == 0 \n";
				currDist = 0.1;
			}
			if (N_bNonInsterestNode)
				continue; // currently non interested node will NOT be stored for test  


			// distance
			ifs_nlink.getline(buffer, static_buffer_size_of_base);
			sLine = buffer;
			zxh::ParseStringLine(sContent, sComment, sLine);
			zxh::trim_both(sContent);
			if (sContent.empty() || ifs_scar.fail())
			{
				std::cout << "error: please find it code-3;\n";
				return 3;
			}
			std::istringstream istr(sContent);
			istr >> nodeid;
			istr >> cellid;
			istr >> dist;
			

			// weight
			ifs_nlink.getline(buffer, static_buffer_size_of_base);
			sLine = buffer;
			zxh::ParseStringLine(sContent, sComment, sLine);
			zxh::trim_both(sContent);
			if (sContent.empty() || ifs_scar.fail())
			{
				std::cout << "error: please find it code-3;\n";
				return 3;
			}
			std::istringstream istr2(sContent);
			istr2 >> nodeid;
			istr2 >> cellid;
			istr2 >> weight;

			float Simialrty = weight;

			float currEdgeStrength = WeightN*0.5*((float)0.95 * Simialrty + (float)0.05);

			if (currEdgeStrength < 0) 
			{
				std::cout << "error: currEdgeStrength is negative;\n";
				continue;
			}

			if (ptId >= 0 && nNodeId >= 0 && (ptId < iNumOfMeshPoints) && (nNodeId < iNumOfMeshPoints) && (ptId != nNodeId))
			{
				myGraph->add_edge(ptId, nNodeId,  /* capacities */ (int)ceil(INT32_CONST*currEdgeStrength + 0.5), (int)ceil(INT32_CONST*currEdgeStrength + 0.5));
			}
		}

		//--------------------------------one cut--------------------------------------
		if (OneCutIdex)
		{
			//add the edge to the auxiliary node
			int currBin = (int)binPerPixelImg.GetPixelGreyscaleClosest(scx, scy, scz, 0);
			myGraph->add_edge(ptId, (GraphType::node_id)(currBin + iNumOfMeshPoints),
				/* capacities */ (int)ceil(INT32_CONST*bha_slope + 0.5), (int)ceil(INT32_CONST*bha_slope + 0.5));

		}

		//--------------------------------graph cut--------------------------------------
		/*float WeightT = 1 - WeightN; */
		float WeightT = 1;//还是用原来的reda的定义
		float  epsilon = pow(10, -8);
		if (GraphCutIdex)
		{
			float LogSourceStrength = -WeightT * log(SinkStrength + epsilon);
			float LogSinkStrength = -WeightT * log(SourceStrength + epsilon);
			if ((LogSourceStrength < 0) || (LogSinkStrength < 0))
			{
				std::cout << "error: LogSourceStrength or LogSinkStrength is negative;\n";
				continue;
			}

			myGraph->add_tweights(ptId, /* capacities */  (int)ceil(INT32_CONST * LogSourceStrength + 0.5), (int)ceil(INT32_CONST * LogSinkStrength + 0.5));

		}

	}

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

// get bin index for each image pixel, store it in binPerPixelImg
void getBinPerPixel(zxhImageDataT<float> & binPerPixelImg, string T_SaveFoldScar, vtkSmartPointer<vtkPolyData>LAMesh, int numBinsPerChannel, int & numUsedBins)
{
	// this vector is used to through away bins that were not used 计算x的y次幂。初值64*64*64空间中初值都是-1
	vector<int> occupiedBinNewIdx((int)pow((double)numBinsPerChannel, (double)3), -1);

	// go over the image
	int newBinIdx = 0;
	const int * Size = SourceImage.GetImageSize();
	const int iNumOfMeshPoints = LAMesh->GetNumberOfPoints();

	//#pragma omp parallel for
	vtkIdType ptId, cellId;
	for (ptId = 0; ptId < iNumOfMeshPoints; ptId++)
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

		float ProbValue; //-----------------TO DO----Load prob

		int bin = (int)(floor(ProbValue*(float)numBinsPerChannel)) + (float)numBinsPerChannel * floor((1 - ProbValue)*(float)numBinsPerChannel);


		if (occupiedBinNewIdx[bin] == -1)
		{
			occupiedBinNewIdx[bin] = newBinIdx;
			newBinIdx++;
		}
		float BinIdx = occupiedBinNewIdx[bin];
		binPerPixelImg.SetPixelByGreyscale(scx, scy, scz, 0, BinIdx);
	}


	double maxBin;
	maxBin = binPerPixelImg.GetPixelGreyscaleMax();//图像中的最大值
	numUsedBins = (int)maxBin + 1;

	occupiedBinNewIdx.clear();
	cout << "Num occupied bins:" << numUsedBins << endl;
}
