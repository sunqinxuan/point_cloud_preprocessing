#include "Pointset.h"
#include "ScanMatching.h"

namespace dgc
{
	namespace matching
	{
		GBPPointset::GBPPointset()
		{
			eps1=660000;
			eps2=0.9;
			eps3=50;
			gicp_epsilon=0.0004;
			LengthPerGrid=16;
			Row=480/LengthPerGrid;
			Col=640/LengthPerGrid;
			PointsPerGrid=LengthPerGrid*LengthPerGrid;
			NumOfGrids=Row*Col;
			scan.clear();
			pointset.clear();
			debug_=false;
		}

		GBPPointset::~GBPPointset()
		{
		}

		void GBPPointset::UniformSample(void)
		{
			float nx,ny,nz;
			vector<GBPPoint> gridpoint;
			GBPPoint point_tmp;
			Mat Cov,Covrgb;
			GBPPoint point,point_avg;
			Mat eigenValues,eigenVectors;

			pointset.clear();
			//int num=0;
			for(int i=0;i<480;i+=9)
			{
				for(int j=0;j<640;j+=9)
				{
					if(scan[640*i+j].x==0&&scan[640*i+j].y==0&&scan[640*i+j].z==0)
						continue;
					else
					{
						point=scan[640*i+j];
						//pointset.push_back(scan[640*i+j]);
						//num++;
					}

					gridpoint.clear();
					for(int k=i-2;k<=i+2;k++)
					{
						for(int l=j-2;l<=j+2;l++)
						{
							point_tmp=scan[640*k+l];
							if(point_tmp.x==0&&point_tmp.y==0&&point_tmp.z==0)
								continue;
							else
							{
								gridpoint.push_back(point_tmp);
							}
						}
					}
					if(gridpoint.size()!=0)
					{
						ComputeCov(gridpoint,point_avg,Cov,Covrgb);
						eigen(Cov, eigenValues, eigenVectors);
						point.nx=eigenVectors.at<float>(2,0);
						point.ny=eigenVectors.at<float>(2,1);
						point.nz=eigenVectors.at<float>(2,2);
						for(int m=0;m<3;m++)
						{
							for(int n=0;n<3;n++)
							{
								point.C[m][n]=Cov.at<float>(m,n);
							}
						}
						pointset.push_back(point);
					}
					else
						continue;

				}
			}
			
			//ofstream fp;
			//fp.open("pointset.txt",ios::out);
			//for(int i=0;i<pointset.size();i++)
			//{
			//	fp<<pointset[i].x<<" "<<pointset[i].y<<" "<<pointset[i].z<<endl;
			//	for(int k=0;k<3;k++)
			//	{
			//		for(int l=0;l<3;l++)
			//		{
			//			fp<<pointset[i].C[k][l]<<" ";
			//		}
			//		fp<<endl;
			//	}
			//	fp<<endl;
			//}
			//fp.close();
		}

		void GBPPointset::Downsample(void)
		{
			float nx,ny,nz;
			vector<GBPPoint> gridpoint;
			GBPPoint point_tmp;
			Mat Cov,Covrgb;
			GBPPoint point_avg;
			Mat eigenValues,eigenVectors;
			Mat eigenValuesrgb,eigenVectorsrgb;
			vector<int> indexes;
			int index;
			vector<GBPPoint> gridpoint_tmp;

			for(int row=0;row<Row;row++)
			{
				for(int col=0;col<Col;col++)
				{
					gridpoint.clear();
					int first=row*Col*PointsPerGrid+col*LengthPerGrid;
					for(int i=0;i<LengthPerGrid;i++)
					{
						for(int j=0;j<LengthPerGrid;j++)
						{
							point_tmp=scan[first+640*i+j];
							if(point_tmp.x==0&&point_tmp.y==0&&point_tmp.z==0)
								continue;
							else
							{
								gridpoint.push_back(point_tmp);
							}
						}
					}
					if(gridpoint.size()!=0)
					{
						ComputeCov(gridpoint,point_avg,Cov,Covrgb);
					}
					else
						continue;

					eigen(Cov, eigenValues, eigenVectors);
					eigen(Covrgb, eigenValuesrgb, eigenVectorsrgb);

					indexes.clear();
					index=0;
					indexes.push_back(index);
					while(indexes.size()>0)
					{
						index=*(indexes.end()-1);
						indexes.pop_back();
						if(index>0)
						{
							gridpoint_tmp.clear();
							for(vector<GBPPoint>::iterator it=gridpoint.begin();it!=gridpoint.end();++it)
							{
								if(it->index==index)
								{
									gridpoint_tmp.push_back(*it);
								}
							}

							ComputeCov(gridpoint_tmp,point_avg,Cov,Covrgb);
							
							eigen(Cov, eigenValues, eigenVectors);
							eigen(Covrgb, eigenValuesrgb, eigenVectorsrgb);
						}
						//(nx,ny,nz): normal vector at point_avg
						nx=eigenVectors.at<float>(2,0);
						ny=eigenVectors.at<float>(2,1);
						nz=eigenVectors.at<float>(2,2);
						//**********************************************************eps3=50
						if(toSplit(gridpoint,index,point_avg.x,point_avg.y,point_avg.z,nx,ny,nz)
							&&index<=126&&eigenValuesrgb.at<float>(0)>eps3)//((eigenValuesn.at<double>(0)>70*zbar*zbar/700000)&&index<=14)//
						{
							GridSplit(gridpoint,index,point_avg.x,point_avg.y,point_avg.z,eigenVectors.at<float>(0,0),eigenVectors.at<float>(0,1),eigenVectors.at<float>(0,2));
							indexes.push_back(index*2+1);
							indexes.push_back(index*2+2);
						}
						else
						{
							point_avg.nx=nx;
							point_avg.ny=ny;
							point_avg.nz=nz;
							point_avg.index=index;
							for(int i=0;i<3;i++)
							{
								for(int j=0;j<3;j++)
								{
									point_avg.C[i][j]=Cov.at<float>(i,j);
								}
							}
							pointset.push_back(point_avg);
						}
					}
				}
			}
			if(debug_)
			{
				ofstream fp1;
				fp1.open("downsample.txt",ios::out);
				for(int k=0;k<pointset.size();k++)
				{
					fp1<<pointset[k].index<<endl;
					fp1<<pointset[k].x<<" "<<pointset[k].y<<" "<<pointset[k].y<<endl;
					fp1<<(int)pointset[k].r<<" "<<(int)pointset[k].g<<" "<<(int)pointset[k].b<<endl;
					fp1<<pointset[k].nx<<" "<<pointset[k].ny<<" "<<pointset[k].nz<<endl;
					for(int i=0;i<3;i++)
					{
						for(int j=0;j<3;j++)
						{
							fp1<<pointset[k].C[i][j]<<" ";
						}
						fp1<<endl;
					}
					fp1<<endl;
				}
				fp1.close();
			}
		}

		void GBPPointset::ComputeCov(vector<GBPPoint> gridpoint, GBPPoint &point_avg, Mat &Cov, Mat &Covrgb)
		{
			Cov = Mat::zeros(3, 3, CV_32FC1);
			Covrgb = Mat::zeros(3, 3, CV_32FC1);
			float xdelta,ydelta,zdelta;
			float rbar=0,gbar=0,bbar=0;
			point_avg.x=0;
			point_avg.y=0;
			point_avg.z=0;
			//point_avg.r=0;
			//point_avg.g=0;
			//point_avg.b=0;
			int num=gridpoint.size();
			for(vector<GBPPoint>::iterator it=gridpoint.begin();it!=gridpoint.end();++it)
			{
				point_avg.x+=it->x;
				point_avg.y+=it->y;
				point_avg.z+=it->z;
				rbar+=it->r;
				gbar+=it->g;
				bbar+=it->b;
			}
			point_avg.x/=num;
			point_avg.y/=num;
			point_avg.z/=num;
			rbar/=num;
			gbar/=num;
			bbar/=num;
			point_avg.r=(unsigned char)rbar;
			point_avg.g=(unsigned char)gbar;
			point_avg.b=(unsigned char)bbar;
			point_avg.nx=0;
			point_avg.ny=0;
			point_avg.nz=0;
			point_avg.index=0;
			for(vector<GBPPoint>::iterator it=gridpoint.begin();it!=gridpoint.end();++it)
			{
				xdelta=it->x-point_avg.x;
				ydelta=it->y-point_avg.y;
				zdelta=it->z-point_avg.z;
				Cov.at<float>(0, 0) += xdelta*xdelta;
				Cov.at<float>(0, 1) += xdelta*ydelta;
				Cov.at<float>(0, 2) += xdelta*zdelta;
				Cov.at<float>(1, 0) += ydelta*xdelta;
				Cov.at<float>(1, 1) += ydelta*ydelta;
				Cov.at<float>(1, 2) += ydelta*zdelta;
				Cov.at<float>(2, 0) += zdelta*xdelta;
				Cov.at<float>(2, 1) += zdelta*ydelta;
				Cov.at<float>(2, 2) += zdelta*zdelta;
				xdelta=(float)it->r-rbar;
				ydelta=(float)it->g-gbar;
				zdelta=(float)it->b-bbar;
				Covrgb.at<float>(0, 0) += xdelta*xdelta;
				Covrgb.at<float>(0, 1) += xdelta*ydelta;
				Covrgb.at<float>(0, 2) += xdelta*zdelta;
				Covrgb.at<float>(1, 0) += ydelta*xdelta;
				Covrgb.at<float>(1, 1) += ydelta*ydelta;
				Covrgb.at<float>(1, 2) += ydelta*zdelta;
				Covrgb.at<float>(2, 0) += zdelta*xdelta;
				Covrgb.at<float>(2, 1) += zdelta*ydelta;
				Covrgb.at<float>(2, 2) += zdelta*zdelta;
			}
			Cov.at<float>(0, 0) /=num;
			Cov.at<float>(0, 1) /=num;
			Cov.at<float>(0, 2) /=num;
			Cov.at<float>(1, 0) /=num;
			Cov.at<float>(1, 1) /=num;
			Cov.at<float>(1, 2) /=num;
			Cov.at<float>(2, 0) /=num;
			Cov.at<float>(2, 1) /=num;
			Cov.at<float>(2, 2) /=num;
			Covrgb.at<float>(0, 0) /=num;
			Covrgb.at<float>(0, 1) /=num;
			Covrgb.at<float>(0, 2) /=num;
			Covrgb.at<float>(1, 0) /=num;
			Covrgb.at<float>(1, 1) /=num;
			Covrgb.at<float>(1, 2) /=num;
			Covrgb.at<float>(2, 0) /=num;
			Covrgb.at<float>(2, 1) /=num;
			Covrgb.at<float>(2, 2) /=num;
		}

		void GBPPointset::GridSplit(vector<GBPPoint> &gridpoint, int index, float xbar, float ybar, float zbar, float nx, float ny, float nz)
		{
			//(nx,ny,nz) is the eigenvector correspond to the largest eigenvalue at (xbar,ybar,zbar)
			float d=nx*xbar+ny*ybar+nz*zbar;
			for(vector<GBPPoint>::iterator it=gridpoint.begin();it!=gridpoint.end();++it)
			{
				if(it->index==index)
				{
					if(nx*it->x+ny*it->y+nz*it->z<=d)
						it->index=it->index*2+1;
					else
						it->index=it->index*2+2;
				}
			}
		}

		bool GBPPointset::toSplit(vector<GBPPoint> &gridpoint, int index, float xbar, float ybar, float zbar, float nx, float ny, float nz)
		{
			//(nx,ny,nz) is the normal vector at (xbar,ybar,zbar)
			int inliers=0;
			int num=0;
			float d=nx*xbar+ny*ybar+nz*zbar;
			float level=1;
			if(index==0)
				level=1;
			else if(index==1||index==2)
				level=1.2;
			else if(index>=3&&index<=14)
				level=1.5;
			else if(index>=15&&index<=30)
				level=1.7;
			else
				level=2;
			for(vector<GBPPoint>::iterator it=gridpoint.begin();it!=gridpoint.end();++it)
			{
				if(it->index==index)
				{
					//**************************************************************eps1=700000
					if(fabs(nx*it->x+ny*it->y+nz*it->z-d)<zbar*zbar/eps1/level)
						inliers++;
					num++;
				}
			}
			//cout<<"num="<<num<<",inliers="<<inliers<<",ratio="<<(float)inliers/num<<endl;
			//****************************eps2=0.9
			if((float)inliers/num<eps2)
				return true;
			else
				return false;
		}

		void GBPPointset::ComputeMatrices()
		{
			if(matrices_done)
				return;
			matrices_done=true;

			gsl_vector *work = gsl_vector_alloc(3);
			if (work == NULL)
				return;
			gsl_vector *gsl_singulars = gsl_vector_alloc(3);
			if (gsl_singulars == NULL)
				return;
			gsl_matrix *gsl_v_mat = gsl_matrix_alloc(3, 3);
			if (gsl_v_mat == NULL)
				return;

			for(int i=0;i<pointset.size();i++)
			{
				mat_t &cov=pointset[i].C;
				gsl_matrix_view gsl_cov=gsl_matrix_view_array(&cov[0][0],3,3);

				//cout<<"i="<<i<<endl;
				//cout<<"gsl_cov:"<<endl;
				//for(int k=0;k<3;k++)
				//{
				//	for(int l=0;l<3;l++)
				//	{
				//		cout<<gsl_matrix_get(&gsl_cov.matrix,k,l)<<" ";
				//	}
				//	cout<<endl;
				//}
				//cout<<"A->size2="<<gsl_cov.matrix.size2<<endl;
				//cout<<"V->size1="<<gsl_v_mat->size1<<endl;
				//cout<<"gsl_v_mat before SVD:"<<endl;
				//for(int k=0;k<3;k++)
				//{
				//	for(int l=0;l<3;l++)
				//	{
				//		cout<<gsl_matrix_get(gsl_v_mat,k,l)<<" ";
				//	}
				//	cout<<endl;
				//}

				gsl_linalg_SV_decomp(&gsl_cov.matrix,gsl_v_mat,gsl_singulars,work);

				//cout<<"V->size1(after SVD)="<<gsl_v_mat->size1<<endl;
				//cout<<"gsl_v_mat:"<<endl;
				//for(int k=0;k<3;k++)
				//{
				//	for(int l=0;l<3;l++)
				//	{
				//		cout<<gsl_matrix_get(gsl_v_mat,k,l)<<" ";
				//	}
				//	cout<<endl;
				//}
				//cout<<"gsl_singulars:"<<endl;
				//for(int k=0;k<3;k++)
				//{
				//	cout<<gsl_vector_get(gsl_singulars,k)<<" ";
				//}
				//cout<<endl<<endl;

				for(int k=0;k<3;k++)
				{
					for(int l=0;l<3;l++)
					{
						cov[k][l]=0;
					}
				}
				for (int k = 0; k < 3; k++) {
					gsl_vector_view col = gsl_matrix_column(gsl_v_mat, k);
					double v = 1.; // biggest 2 singular values replaced by 1
					//****************************************************gicp_epsilon=0.0004
					if (k == 2) // smallest singular value replaced by gicp_epsilon
						v = gicp_epsilon;
					gsl_blas_dger(v, &col.vector, &col.vector, &gsl_cov.matrix);
				}
			}
			if (work != NULL)
				gsl_vector_free(work);
			if (gsl_v_mat != NULL)
				gsl_matrix_free(gsl_v_mat);
			if (gsl_singulars != NULL)
				gsl_vector_free(gsl_singulars);
		}

		void GBPPointset::Clear(void)
		{
			pointset.clear();
			scan.clear();
		}

		//int GBPPointset::AlignScan(GBPPointset *ds, dgc_transform_t base_t, dgc_transform_t t)
		//{

		//}
	}
}