#include "PatchExtractByWorldCoordinate.h"

bool PatchExtractByWorldCoordinate::GeneratePatchExtractByWorldWithRandOffset(float InputOutputWorldCoord[], int PatchSize[], float randoffset, zxhImageData* ImageArray[], int PatchCorId)
{

	zxhImageData * pSourceImage = ImageArray[0];
	zxhImageData * pLAGaussianBlurLabel =  ImageArray[1];
	zxhImageData * pPatchImage_1 = ImageArray[2];
	zxhImageData * pPatchImage_2 = ImageArray[3];
	zxhImageData * pPatchImage_3 = ImageArray[4];
	zxhImageData * pPatchImage_4 = ImageArray[5];

	int N = PatchSize[0];
	int HalfPatchLength = PatchSize[1]; 

	//Get gradient for LA Label
	zxhImageModelingLinear GradientMod;
	GradientMod.SetImage(pLAGaussianBlurLabel);//using LA label image
	float  pwGrad[4] = { 0 };//assign the initial value
	GradientMod.GetPixelGradientByWorld(pwGrad, InputOutputWorldCoord[0], InputOutputWorldCoord[1], InputOutputWorldCoord[2], 0);
	 
	if (zxh::VectorOP_Magnitude(pwGrad,3) <ZXH_FloatInfinitesimal)
	{
		std::cout << "error: magnitude too small for node " << PatchCorId << "\n";
		return false ;
	}
	zxh::VectorOP_Normalise( pwGrad, 3 ) ;
	float Ia = 1.0*pwGrad[0]; // 每步只有1mm步长，没有防止错过scar厚度小于1mm的点
	float Ib = 1.0*pwGrad[1];
	float Ic = 1.0*pwGrad[2];
	 
	//Rand shift
	float InputWorldCoord_x = InputOutputWorldCoord[0] + randoffset*Ia;
	float InputWorldCoord_y = InputOutputWorldCoord[1] + randoffset*Ib;
	float InputWorldCoord_z = InputOutputWorldCoord[2] + randoffset*Ic;

	//get the tangent plane , @李雷，(1)请写出A B C三个单位向量计算的思路；（2）请double check这个代码，（3）同时，写个代码计算一下这三个向量是否相互垂直（（abs(cos)<ZXH_FloatInfinitesimal)，且mag都差不多等于1（abs(mag-1)<ZXH_FloatInfinitesimal))
	float ttt;
	ttt = -(InputWorldCoord_x*Ia + InputWorldCoord_y*Ib + InputWorldCoord_z*Ic);
	float ta = -InputWorldCoord_x - ttt*Ia;
	float tb = -InputWorldCoord_y - ttt*Ib;
	float tc = -InputWorldCoord_z - ttt*Ic;
	float tna = ta / (sqrt(ta*ta + tb*tb + tc*tc) + ZXH_FloatInfinitesimal);
	float tnb = tb / (sqrt(ta*ta + tb*tb + tc*tc) + ZXH_FloatInfinitesimal);
	float tnc = tc / (sqrt(ta*ta + tb*tb + tc*tc) + ZXH_FloatInfinitesimal);
	float tva = Ib*tnc - Ic*tnb;
	float tvb = Ic*tna - Ia*tnc;
	float tvc = Ia*tnb - Ib*tna;
	///  （3）同时，写个代码计算一下这三个向量是否相互垂直（（abs(cos)<ZXH_FloatInfinitesimal)，且mag都差不多等于1（abs(mag-1)<ZXH_FloatInfinitesimal))
	float vA[]={tna,tnb,tnc}, vB[]={tva,tvb,tvc}, vC[]={Ia,Ib,Ic};
	if( zxh::absf(zxh::VectorOP_Cosine(vA,vB,3))>1.0e-3 ||zxh::absf(zxh::VectorOP_Cosine(vA,vC,3))>1.0e-3 ||zxh::absf(zxh::VectorOP_Cosine(vB,vC,3))>1.0e-3 )
	{
		std::cerr<<"error: axises should be perpendicular\n";
		return false ;
	}
	if( zxh::absf(zxh::VectorOP_Magnitude(vA,3)-1.0 )>1.0e-3 || zxh::absf(zxh::VectorOP_Magnitude(vB,3)-1.0 )>1.0e-3 || zxh::absf(zxh::VectorOP_Magnitude(vC,3)-1.0 )>1.0e-3 )
	{
		std::cerr<<"error: axis vectors should be normalized\n";
		return false ;
	}


	float PatchPointWorldCoord[3] = { 0 };  //啊？为啥之前这个定义为int？？？？？
	//Computing the intensity of points in the tangent plane using interpolation	
	//Interpolation
	zxhImageModelingLinear InterpolationMod;
	InterpolationMod.SetImage(pSourceImage);
	for (int i = -N; i < N + 1; i++)
	{
		for (int j = -N; j < N + 1; j++)
		{
			for (int k = -HalfPatchLength; k < HalfPatchLength + 1; k++)//extend the patch along the gradient orientation
			{ 
				//scale 1
				PatchPointWorldCoord[0] = InputWorldCoord_x + (i*tna + j*tva + k*Ia);
				PatchPointWorldCoord[1] = InputWorldCoord_y + (i*tnb + j*tvb + k*Ib);
				PatchPointWorldCoord[2] = InputWorldCoord_z + (i*tnc + j*tvc + k*Ic);	
				float IntensityValue = InterpolationMod.GetPixelFloatValueWithCheckByWorld(PatchPointWorldCoord[0], PatchPointWorldCoord[1], PatchPointWorldCoord[2], 0);				 
				int W2I_Coor_SetPatchPixel = (i + N) + (j + N) * (2 * N + 1) + (k + HalfPatchLength)* (2 * N + 1)* (2 * N + 1);

				int Idex = PatchCorId / 30000;//取商
				if (pPatchImage_1->InsideImage(W2I_Coor_SetPatchPixel, PatchCorId - 30000 * Idex, Idex, 0) == false)
					std::cout << "error\n";
				pPatchImage_1->SetPixelByGreyscale(W2I_Coor_SetPatchPixel, PatchCorId - 30000 * Idex, Idex, 0, IntensityValue);

				//scale 2
				PatchPointWorldCoord[0] = InputWorldCoord_x + (i*tna + j*tva + k*Ia)*2;
				PatchPointWorldCoord[1] = InputWorldCoord_y + (i*tnb + j*tvb + k*Ib)*2;
				PatchPointWorldCoord[2] = InputWorldCoord_z + (i*tnc + j*tvc + k*Ic)*2;			
				IntensityValue = InterpolationMod.GetPixelFloatValueWithCheckByWorld(PatchPointWorldCoord[0], PatchPointWorldCoord[1], PatchPointWorldCoord[2], 0);
				W2I_Coor_SetPatchPixel = (i + N) + (j + N) * (2 * N + 1) + (k + HalfPatchLength)* (2 * N + 1)* (2 * N + 1);

				Idex = PatchCorId / 30000;//取商
				if (pPatchImage_2->InsideImage(W2I_Coor_SetPatchPixel, PatchCorId - 30000 * Idex, Idex, 0) == false)
					std::cout << "error\n";
				pPatchImage_2->SetPixelByGreyscale(W2I_Coor_SetPatchPixel, PatchCorId - 30000 * Idex, Idex, 0, IntensityValue);

				//scale 3
				PatchPointWorldCoord[0] = InputWorldCoord_x + (i*tna + j*tva + k*Ia)*3;
				PatchPointWorldCoord[1] = InputWorldCoord_y + (i*tnb + j*tvb + k*Ib)*3;
				PatchPointWorldCoord[2] = InputWorldCoord_z + (i*tnc + j*tvc + k*Ic)*3;			
				IntensityValue = InterpolationMod.GetPixelFloatValueWithCheckByWorld(PatchPointWorldCoord[0], PatchPointWorldCoord[1], PatchPointWorldCoord[2], 0);
				W2I_Coor_SetPatchPixel = (i + N) + (j + N) * (2 * N + 1) + (k + HalfPatchLength)* (2 * N + 1)* (2 * N + 1);

				Idex = PatchCorId / 30000;//取商
				if (pPatchImage_3->InsideImage(W2I_Coor_SetPatchPixel, PatchCorId - 30000 * Idex, Idex, 0) == false)
					std::cout << "error\n";
				pPatchImage_3->SetPixelByGreyscale(W2I_Coor_SetPatchPixel, PatchCorId - 30000 * Idex, Idex, 0, IntensityValue);

				//scale 4
				PatchPointWorldCoord[0] = InputWorldCoord_x + (i*tna + j*tva + k*Ia) * 4;
				PatchPointWorldCoord[1] = InputWorldCoord_y + (i*tnb + j*tvb + k*Ib) * 4;
				PatchPointWorldCoord[2] = InputWorldCoord_z + (i*tnc + j*tvc + k*Ic) * 4;
				IntensityValue = InterpolationMod.GetPixelFloatValueWithCheckByWorld(PatchPointWorldCoord[0], PatchPointWorldCoord[1], PatchPointWorldCoord[2], 0);
				W2I_Coor_SetPatchPixel = (i + N) + (j + N) * (2 * N + 1) + (k + HalfPatchLength)* (2 * N + 1)* (2 * N + 1);

				Idex = PatchCorId / 30000;//取商
				if (pPatchImage_4->InsideImage(W2I_Coor_SetPatchPixel, PatchCorId - 30000 * Idex, Idex, 0) == false)
					std::cout << "error\n";
				pPatchImage_4->SetPixelByGreyscale(W2I_Coor_SetPatchPixel, PatchCorId - 30000 * Idex, Idex, 0, IntensityValue);
			
			}
		}
	}
	 
	InputOutputWorldCoord[0] += randoffset*Ia; // return offset world
	InputOutputWorldCoord[1] += randoffset*Ib;
	InputOutputWorldCoord[2] += randoffset*Ic;
	
	return true ;

}


PatchExtractByWorldCoordinate::~PatchExtractByWorldCoordinate()
{

}
