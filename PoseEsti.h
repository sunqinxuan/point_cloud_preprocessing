#pragma once
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "Pointset.h"

#include <vector>
#include <string>
#include <time.h>
#include <iostream>
#include <math.h>

using namespace cv;
using namespace std;

namespace dgc
{
	namespace matching
	{
		typedef struct PointRansac {
			float x;
			float y;
			float z;
		} PointRansac;

		class PoseEsti
		{
		public:
			PoseEsti();
			~PoseEsti();
			void SetDebug(bool d) {debug_=d;}
			void SetDetectThresh(int th) {eps1=th;}
			void SetMatchDistThresh(float th) {eps2=th;}
			void SetRansacIteration(int i) {iter=i;}
			void SetInlierDistThresh(float th) {eps=th;}
			Mat GetRf(void) {return Rf;}
			Mat GetTf(void) {return Tf;}
			PointRansac GetPoint1(int i) {return points1_final[i];}
			PointRansac GetPoint2(int i) {return points2_final[i];}
			int GetSize(void) {return points1_final.size();}
			void LoadImg(Mat im1, Mat im2) {img1=Mat::Mat(im1); img2=Mat::Mat(im2);}
			void clear(void);
			void FeatureMatch(void);//20,20
			int RansacEsti(GBPPointset *ps1, GBPPointset *ps2);//400,50
			//Mat GetImg1(void) {return img1;}
			//Mat GetImg2(void) {return img2;}

		private:
			Mat Rf,Tf;
			Mat img1,img2;
			vector<KeyPoint> keyPoints1,keyPoints2;
			vector<DMatch> matches;
			vector<PointRansac> points1_final,points2_final;
			bool debug_;
			int eps1,iter;
			float eps2,eps;
		};
	}
}


//class KeyPointEsti
//{
//public:
//	KeyPointEsti();
//	~KeyPointEsti();
//	void Extract();
//	int RansacEsti();
//	Mat Rf;	
//	Mat Tf;
//private:
//	vector<KeyPoint> keyPoints1,keyPoints2;
//	vector<DMatch> matches;	
//	vector<PointRansac> points1,points2;
//	PointRansac Point;
//};