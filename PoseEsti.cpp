#include "PoseEsti.h"

namespace dgc
{
	namespace matching
	{
		PoseEsti::PoseEsti()
		{
			eps1=20;
			eps2=15;
			eps=50;
			iter=500;
			debug_=false;
		}

		PoseEsti::~PoseEsti()
		{
		}

		void PoseEsti::FeatureMatch(void)
		{
			FastFeatureDetector fast(eps1);//参数1：越大，点数越少！15：1500；30：450
			fast.detect(img1,keyPoints1);
			fast.detect(img2,keyPoints2);
			
			BriefDescriptorExtractor extractor;
			Mat descriptors1, descriptors2;
			extractor.compute(img1, keyPoints1, descriptors1);
			extractor.compute(img2, keyPoints2, descriptors2);

			BruteForceMatcher<Hamming> matcher;
			matcher.match(descriptors1, descriptors2, matches);

			if(debug_)
			{
				cout<<"points number:"<<keyPoints1.size()<<","<<keyPoints2.size()<<endl;
				cout<<"matched points:"<<matches.size()<<endl;
			}
			for (vector<DMatch>::iterator it = matches.begin(); it != matches.end(); )
			{
				//********************************************************
				if (it->distance > eps2) //参数2：匹配特征点筛选阈值，10：300
				{
					matches.erase(it);
				}
				else
				{
					it++;
				}
			}

			if(debug_)
			{
				cout <<"matched points after filtering:" << matches.size() << endl;
				namedWindow("matches", 1);
				Mat img_matches;
				drawMatches(img1, keyPoints1, img2, keyPoints2, matches, img_matches);
				imshow("matches", img_matches);
			}

		}

		int PoseEsti::RansacEsti(GBPPointset *ps1, GBPPointset *ps2)
		{
			vector<PointRansac> points1,points2;
			PointRansac Point;
			points1.clear();
			points2.clear();
			int m1,n1,m2,n2;
			for(vector<DMatch>::iterator it=matches.begin();it!=matches.end();++it)
			{
				m1=keyPoints1.at((*it).queryIdx).pt.x;
				n1=keyPoints1.at((*it).queryIdx).pt.y;
				m2=keyPoints2.at((*it).trainIdx).pt.x;
				n2=keyPoints2.at((*it).trainIdx).pt.y;
				if(ps1->GetPointOrg(n1,m1).x+ps1->GetPointOrg(n1,m1).y+ps1->GetPointOrg(n1,m1).z==0||ps2->GetPointOrg(n2,m2).x+ps2->GetPointOrg(n2,m2).y+ps2->GetPointOrg(n2,m2).z==0)
					continue;
				else
				{
					Point.x=ps1->GetPointOrg(n1,m1).x;
					Point.y=ps1->GetPointOrg(n1,m1).y;
					Point.z=ps1->GetPointOrg(n1,m1).z;
					points1.push_back(Point);
					Point.x=ps2->GetPointOrg(n2,m2).x;
					Point.y=ps2->GetPointOrg(n2,m2).y;
					Point.z=ps2->GetPointOrg(n2,m2).z;
					points2.push_back(Point);//points1和points2是匹配好的特征点集在三维空间的投影
				}
			}

			Mat p11=Mat(3,1,CV_32FC1);
			Mat p12=Mat(3,1,CV_32FC1);
			Mat p13=Mat(3,1,CV_32FC1);
			Mat p21=Mat(3,1,CV_32FC1);
			Mat p22=Mat(3,1,CV_32FC1);
			Mat p23=Mat(3,1,CV_32FC1);
			Mat x1=Mat(3,1,CV_32FC1);
			Mat y1=Mat(3,1,CV_32FC1);
			Mat z1=Mat(3,1,CV_32FC1);
			Mat x2=Mat(3,1,CV_32FC1);
			Mat y2=Mat(3,1,CV_32FC1);
			Mat z2=Mat(3,1,CV_32FC1);
			Mat M1=Mat(3,3,CV_32FC1);
			Mat M2=Mat(3,3,CV_32FC1);
			Mat R=Mat(3,3,CV_32FC1);
			Mat T=Mat(3,1,CV_32FC1);
			Mat p1_bar=Mat(3,1,CV_32FC1);
			Mat p2_bar=Mat(3,1,CV_32FC1);
			Mat test1=Mat(3,1,CV_32FC1);
			Mat test2=Mat(3,1,CV_32FC1);
			Mat test3=Mat(3,1,CV_32FC1);
			int k1,k2,k3;
			vector<PointRansac> points1_in,points2_in;
			srand((int)time(0));
			float dist=0,pre_dist=99999999;
			int count=0;
			int count_pre=0;

			Rf=Mat(3,3,CV_32FC1);
			Tf=Mat(3,1,CV_32FC1);
			//cout<<"ransac start"<<endl;
			//*******************************************
			for(int j=0;j<iter;j++)//参数3：ransac迭代次数
			//while(count<100)
			{
				while(1)
				{
					k1=rand()%points1.size();
					k2=rand()%points1.size();
					k3=rand()%points1.size();
					if(k1!=k2&&k1!=k3&&k2!=k3)
						break;
				}
				//cout<<"k1:"<<k1<<endl;
				//cout<<"k2:"<<k2<<endl;
				//cout<<"k3:"<<k3<<endl;
				//fprintf(rt,"k1:%d\nk2:%d\nk3:%d\n",k1,k2,k3);
				p11.at<float>(0,0)=points1.at(k1).x;
				p11.at<float>(1,0)=points1.at(k1).y;
				p11.at<float>(2,0)=points1.at(k1).z;
				p12.at<float>(0,0)=points1.at(k2).x;
				p12.at<float>(1,0)=points1.at(k2).y;
				p12.at<float>(2,0)=points1.at(k2).z;
				p13.at<float>(0,0)=points1.at(k3).x;
				p13.at<float>(1,0)=points1.at(k3).y;
				p13.at<float>(2,0)=points1.at(k3).z;
				p21.at<float>(0,0)=points2.at(k1).x;
				p21.at<float>(1,0)=points2.at(k1).y;
				p21.at<float>(2,0)=points2.at(k1).z;
				p22.at<float>(0,0)=points2.at(k2).x;
				p22.at<float>(1,0)=points2.at(k2).y;
				p22.at<float>(2,0)=points2.at(k2).z;
				p23.at<float>(0,0)=points2.at(k3).x;
				p23.at<float>(1,0)=points2.at(k3).y;
				p23.at<float>(2,0)=points2.at(k3).z;
	
				x1=p12-p11;
				x1=x1/sqrt(x1.at<float>(0,0)*x1.at<float>(0,0)+x1.at<float>(1,0)*x1.at<float>(1,0)+x1.at<float>(2,0)*x1.at<float>(2,0));
				x2=p22-p21;
				x2=x2/sqrt(x2.at<float>(0,0)*x2.at<float>(0,0)+x2.at<float>(1,0)*x2.at<float>(1,0)+x2.at<float>(2,0)*x2.at<float>(2,0));
				y1=(p13-p11)-x1.dot(p13-p11)*x1;
				y1=y1/sqrt(y1.at<float>(0,0)*y1.at<float>(0,0)+y1.at<float>(1,0)*y1.at<float>(1,0)+y1.at<float>(2,0)*y1.at<float>(2,0));
				y2=(p23-p21)-x2.dot(p23-p21)*x2;
				y2=y2/sqrt(y2.at<float>(0,0)*y2.at<float>(0,0)+y2.at<float>(1,0)*y2.at<float>(1,0)+y2.at<float>(2,0)*y2.at<float>(2,0));
				z1=x1.cross(y1);
				z2=x2.cross(y2);

				M1.at<float>(0,0)=x1.at<float>(0,0);
				M1.at<float>(1,0)=x1.at<float>(1,0);
				M1.at<float>(2,0)=x1.at<float>(2,0);
				M1.at<float>(0,1)=y1.at<float>(0,0);
				M1.at<float>(1,1)=y1.at<float>(1,0);
				M1.at<float>(2,1)=y1.at<float>(2,0);
				M1.at<float>(0,2)=z1.at<float>(0,0);
				M1.at<float>(1,2)=z1.at<float>(1,0);
				M1.at<float>(2,2)=z1.at<float>(2,0);
				M2.at<float>(0,0)=x2.at<float>(0,0);
				M2.at<float>(1,0)=x2.at<float>(1,0);
				M2.at<float>(2,0)=x2.at<float>(2,0);
				M2.at<float>(0,1)=y2.at<float>(0,0);
				M2.at<float>(1,1)=y2.at<float>(1,0);
				M2.at<float>(2,1)=y2.at<float>(2,0);
				M2.at<float>(0,2)=z2.at<float>(0,0);
				M2.at<float>(1,2)=z2.at<float>(1,0);
				M2.at<float>(2,2)=z2.at<float>(2,0);
				R=M2*M1.t();

				p1_bar=(p11+p12+p13)/3;
				p2_bar=(p21+p22+p23)/3;
				T=p2_bar-R*p1_bar;

				count=0;
				//dist=0;
				points1_in.clear();
				points2_in.clear();
				for(int i=0;i<points1.size();i++)
				{
					test1.at<float>(0,0)=points1.at(i).x;
					test1.at<float>(1,0)=points1.at(i).y;
					test1.at<float>(2,0)=points1.at(i).z;
					test1=R*test1+T;
					test2.at<float>(0,0)=points2.at(i).x;
					test2.at<float>(1,0)=points2.at(i).y;
					test2.at<float>(2,0)=points2.at(i).z;
					dist=norm(test1,test2);
					//******************************
					if(dist<eps)//参数4：内点筛选阈值
					{
						count++;
						points1_in.push_back(points1.at(i));
						points2_in.push_back(points2.at(i));
					}
					//cout<<i<<":  "<<dist<<endl;
					//fprintf(rt,"%d: %f\n",i+1,dist);
				}
				//cout << endl;
				//cout << "count=" << count << endl;
				//cout << "points1_in.size()=" << points1_in.size() << endl;
				//waitKey(0);
		

				if(count>count_pre)
				{
					points1_final.clear();
					points2_final.clear();
					for(vector<PointRansac>::iterator it=points1_in.begin();it!=points1_in.end();++it)
					{
						points1_final.push_back(*it);
					}
					for(vector<PointRansac>::iterator it=points2_in.begin();it!=points2_in.end();++it)
					{
						points2_final.push_back(*it);
					}
					Rf=R.clone();
					Tf=T.clone();
					//cout<<"count_pre: "<<count_pre<<endl;
					//cout<<"count: "<<count<<endl<<endl;
					count_pre=count;
				}
			}


			if(debug_)
			{
				dist = 0;
				for (int i = 0; i<points1_final.size(); i++)
				{
					test1.at<float>(0, 0) = points1_final.at(i).x;
					test1.at<float>(1, 0) = points1_final.at(i).y;
					test1.at<float>(2, 0) = points1_final.at(i).z;
					test1 = Rf*test1 + Tf;
					test2.at<float>(0, 0) = points2_final.at(i).x;
					test2.at<float>(1, 0) = points2_final.at(i).y;
					test2.at<float>(2, 0) = points2_final.at(i).z;
					//cout << norm(test1, test2) << ",";
					dist += norm(test1, test2);
				}
				//cout << endl;
				dist /= count_pre;
				//dist = 0;
				//vector<PointRansac>::iterator it1 = points2_final.begin();
				//for (vector<PointRansac>::iterator it = points1_final.begin(); it != points1_final.end(); ++it)
				//{
				//	dist += sqrt((it->x - it1->x)*(it->x - it1->x) + (it->y - it1->y)*(it->y - it1->y) + (it->z - it1->z)*(it->z - it1->z));
				//	it1++;
				//}
				//dist /= count_pre;

				cout << "dist=" << dist << endl;
				cout << "inliers=" << count_pre << endl;
			}

			return count_pre;
		}

		void clear(void)
		{

		}
	}
}