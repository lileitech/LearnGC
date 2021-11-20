//Copyright (c) 2014, Lena Gorelick
//All rights reserved.
//

//
//THIS SOFTWARE IMPLEMENTS THE OneCut ALGORITHM THAT USES SCRIBBLES AS HARD CONSTRAINTS.
//PLEASE USE THE FOLLOWING CITATION:
//
//@inproceedings{iccv2013onecut,
//	title	= {Grabcut in One Cut},
//	author	= {Tang, Meng and Gorelick, Lena and Veksler, Olga and Boykov, Yuri},
//	booktitle={International Conference on Computer Vision},
//	month	= {December},
//	year	= {2013}}
//
//THIS SOFTWARE USES maxflow/min-cut CODE THAT WAS IMPLEMENTED BY VLADIMIR KOLMOGOROV,
//THAT CAN BE DOWNLOADED FROM http://vision.csd.uwo.ca/code/.
//PLEASE USE THE FOLLOWING CITATION:
//
//@ARTICLE{Boykov01anexperimental,
//    author = {Yuri Boykov and Vladimir Kolmogorov},
//    title = {An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Vision},
//    journal = {IEEE TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE},
//    year = {2001},
//    volume = {26},
//    pages = {359--374}}
//
//
//
//THIS SOFTWARE USES OpenCV 2.4.3 THAT CAN BE DOWNLOADED FROM http://opencv.org
//
#pragma once
#include "graph.h"
#include <opencv2/imgproc/imgproc.hpp>  // Gaussian Blur
#include <opencv2/core/core.hpp>        // Basic OpenCV structures (cv::Mat, Scalar)
#include <opencv2/highgui/highgui.hpp>  // OpenCV window I/O

using namespace std;
using namespace cv;

class segmenter
{
public:
	segmenter(void);
	~segmenter(void);
	// init all images/vars
	int  init(char * imgFileName);
	// clear everything before closing
	void destroyAll();
	// mouse listener
	static void onMouse(int event, int x, int y, int, void*);
	// set bin index for each image pixel, store it in binPerPixelImg
	void getBinPerPixel(Mat & binPerPixelImg, Mat & inputImg, int numBinsPerChannel, int & numUsedBins);
	// compute the variance of image edges between neighbors
	void getEdgeVariance(Mat & inputImg, Mat & showEdgesImg, float & varianceSquared);

	typedef Graph<int, int, int> GraphType;
	GraphType *myGraph;
};

