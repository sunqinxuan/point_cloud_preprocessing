#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "ANN.h"
#include "transform.h"

using namespace std;
using namespace cv;

namespace dgc
{
	namespace matching
	{
		typedef double mat_t[3][3];

		typedef struct GBPPoint
		{
			float x,y,z;
			unsigned char r,g,b;
			float nx,ny,nz;
			mat_t C;
			int index;
		}GBPPoint;

		class GBPPointset
		{
		public:
			GBPPointset();
			~GBPPointset();
			void AppendPoint(GBPPoint const &pt) {scan.push_back(pt);}//640*480
			int SizeOrg() {return scan.size();}
			int SizeDS() {return pointset.size();}
			GBPPoint GetPointOrg(int i,int j) {return scan[i*640+j];}//******************
			GBPPoint GetPointDS(int i) {return pointset[i];}//*********************
			void SetPointDS(int i, GBPPoint point) {pointset[i]=point;}
			void SetPlaneThresh(float th) {eps1=th;}
			void SetInlierThresh(float th) {eps2=th;}
			void SetColorThresh(float th) {eps3=th;}
			void SetGICPEpsilon(double eps) {gicp_epsilon=eps;}
			void SetDebug(bool d) {debug_=d;}

			GBPPoint & operator[](int i) { return pointset[i]; }
			GBPPoint const& operator[](int i) const { return pointset[i]; }
			GBPPoint & operator()(int i) { return scan[i]; }
			GBPPoint const& operator()(int i) const { return scan[i]; }

			void UniformSample(void);
			void Downsample(void);
			void ComputeMatrices();
			void Clear(void);
			//int AlignScan(GBPPointset *ds, dgc_transform_t base_t, dgc_transform_t t);

		private:
			void GridSplit(vector<GBPPoint> &gridpoint, int index, float xbar, float ybar, float zbar, float nx, float ny, float nz);
			bool toSplit(vector<GBPPoint> &gridpoint, int index, float xbar, float ybar, float zbar, float nx, float ny, float nz);
			void ComputeCov(vector<GBPPoint> gridpoint, GBPPoint &point_avg, Mat &Cov, Mat &Covrgb);

			vector<GBPPoint> scan;
			vector<GBPPoint> pointset;
			float eps1,eps2,eps3;
			int LengthPerGrid,Row,Col,PointsPerGrid,NumOfGrids;

			bool matrices_done;
			double gicp_epsilon;
			bool debug_;
		};
	}
}