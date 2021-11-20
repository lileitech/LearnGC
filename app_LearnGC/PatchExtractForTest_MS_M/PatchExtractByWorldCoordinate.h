#pragma once
#include <string.h>
#include <iostream> 
#include <time.h> 
#include <math.h>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <malloc.h>

#include "zxhImageGipl.h" 
#include "zxhImageModelingLinear.h"

#define pi 3.1415926 
using namespace std;

class PatchExtractByWorldCoordinate
{
public: // for methods
	PatchExtractByWorldCoordinate(){} ; 
	virtual ~PatchExtractByWorldCoordinate() ;
	
	//virtual void GeneraetPatchExtractByWorldCoordinate(float InputWorldCoord[], int PatchInfo[], zxhImageData ImageArray[], int PatchCorId);
	///
	static bool GeneratePatchExtractByWorldWithRandOffset_test(float InputWorldCoord[], int PatchSize[], float randoffset, zxhImageData *ImageArray[6], int PatchCorId);

public: // for variables
	//zxhImageData *pPatchImage;
};

