#include "GL/freeglut.h"
#include "PoseEsti.h"
#include "transform.h"
#include "Pointset.h"
#include "ScanMatching.h"
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "gicp.h"

using namespace dgc::gicp;
using namespace dgc::matching;

#define Height 480
#define Width 640
//#define M_PI 3.14159265358979323846
#define uchar unsigned char
#define SIGN(x) ( (x)<0 ? -1:((x)>0?1:0 ) )

static string filename1="0_pointCloud.txt";
static string filename2="1_pointCloud.txt";
static string filename_t_base;
static bool debug = true;
static double gicp_epsilon = 1e-3;
static double max_distance = 5.;

int glWinWidth = 640, glWinHeight = 480;
double eyex, eyey, eyez, atx, aty, atz;
double sx = 1, sy = 1, sz = 1;
bool leftClickHold = false, rightClickHold = false;
bool wheelup = false, wheeldown = false;
int mx, my;
int ry = 90, rx = 90;
double mindepth = 200.0, maxdepth = 2000.0;
double radius = 3000.0;
extern bool flag_endloop = 0;
GBPPointset ds1, ds2;

void mouse(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON)
	{
		if (state == GLUT_DOWN)
		{
			leftClickHold = true;
		}
		else
		{
			leftClickHold = false;
		}
	}

	if (button == GLUT_RIGHT_BUTTON)
	{
		if (state == GLUT_DOWN)
		{
			rightClickHold = true;
		}
		else
		{
			rightClickHold = false;
		}
	}
}

void key(int key, int x, int y)
{
	switch (key)
	{
	case GLUT_KEY_UP:
		sx += 0.01;
		sy += 0.01;
		sz += 0.01;
		glutPostRedisplay();
		break;
	case GLUT_KEY_DOWN:
		sx -= 0.01;
		sy -= 0.01;
		sz -= 0.01;
		glutPostRedisplay();
		break;
	case GLUT_KEY_END:
		flag_endloop = true;
		break;
	default:
		sx = sx;
		sy = sy;
		sz = sz;
	}
}

void motion(int x, int y)
{
	int rstep = 5;
	if (leftClickHold == true)
	{
		if (abs(x - mx) > abs(y - my))
		{
			rx += SIGN(x - mx)*rstep;
		}
		else
		{
			ry -= SIGN(y - my)*rstep;
		}

		mx = x;
		my = y;
		glutPostRedisplay();
	}

	if (rightClickHold == true)
	{
		radius += SIGN(y - my)*100.0;
		radius = std::max(radius, 100.0);
		mx = x;
		my = y;
		glutPostRedisplay();
	}
}

void reshape(int w, int h)
{
	glWinWidth = w;
	glWinHeight = h;
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45, (GLfloat)w / (GLfloat)h, 1.0, 15000.0);
	glMatrixMode(GL_MODELVIEW);
}

void renderScene(void)
{
	// clear screen and depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// Reset the coordinate system before modifying 
	glLoadIdentity();
	// set the camera position
	atx = 0.0f;
	aty = 0.0f;
	atz = (mindepth - maxdepth) / 2.0f;
	eyex = atx + radius * sin(CV_PI * ry / 180.0f) * cos(CV_PI * rx / 180.0f);
	eyey = aty + radius * cos(CV_PI * ry / 180.0f);
	eyez = atz + radius * sin(CV_PI * ry / 180.0f) * sin(CV_PI * rx / 180.0f);
	gluLookAt(eyex, eyey, eyez, atx, aty, atz, 0.0, -1.0, 0.0);
	//gluLookAt(0.0f,0.0f,6000.0f,0.0f,0.0f,1000.0f,0.0,1.0,0.0);
	glRotatef(0, 0, 1, 0);
	glRotatef(-180, 1, 0, 0);
	glScalef(sx, sy, sz);
	
	glPointSize(3.0);
	glBegin(GL_POINTS);

	glColor3f(1, 0, 0);
	glVertex3f(0, 0, 0);
	//z:b
	glColor3f(0, 0, 1);
	glVertex3f(0, 0, 100);
	glColor3f(0, 0, 1);
	glVertex3f(0, 0, 500);
	glColor3f(0, 0, 1);
	glVertex3f(0, 0, 1000);
	//x:r
	glColor3f(1, 0, 0);
	glVertex3f(100, 0, 0);
	glColor3f(1, 0, 0);
	glVertex3f(500, 0, 0);
	glColor3f(1, 0, 0);
	glVertex3f(1000, 0, 0);
	//y:g
	glColor3f(0, 1, 0);
	glVertex3f(0, 100, 0);
	glColor3f(0, 1, 0);
	glVertex3f(0, 500, 0);
	glColor3f(0, 1, 0);
	glVertex3f(0, 1000, 0);

	//original points
	for (int i=0;i<480;i+=1)
	{
		for(int j=0;j<640;j+=1)
		{
			glColor3f(1,0,0);
			//glColor3f((float)ds1.GetPointOrg(i,j).r/255,(float)ds1.GetPointOrg(i,j).g/255,(float)ds1.GetPointOrg(i,j).b/255);
			glVertex3f(ds1.GetPointOrg(i,j).x,ds1.GetPointOrg(i,j).y,ds1.GetPointOrg(i,j).z);
			glColor3f(1,1,1);
			//glColor3f((float)ds2.GetPointOrg(i,j).r/255,(float)ds2.GetPointOrg(i,j).g/255,(float)ds2.GetPointOrg(i,j).b/255);
			glVertex3f(ds2.GetPointOrg(i,j).x,ds2.GetPointOrg(i,j).y,ds2.GetPointOrg(i,j).z);
		}
	}
	glEnd();

	//downsampled points
	//glPointSize(3.0);
	//glBegin(GL_POINTS);
	//glColor3f(0,1,0);
	//for(int i=0;i<ds1.SizeDS();i++)
	//{
	//	//glColor3f((float)it->r/255,(float)it->g/255,(float)it->b/255);
	//	glVertex3f(ds1[i].x,ds1[i].y,ds1[i].z);
	//}
	//glColor3f(1,0,0);
	//for(int i=0;i<ds2.SizeDS();i++)
	//{
	//	//glColor3f((float)it->r/255,(float)it->g/255,(float)it->b/255);
	//	glVertex3f(ds2[i].x,ds2[i].y,ds2[i].z);
	//}
	//glEnd();

	//normals
	//glBegin(GL_LINES);
	//glColor3f(1, 1, 1);
	//for (int i = 0; i < ds1.SizeDS(); i++)
	//{
	//	glVertex3f(ds1[i].x,ds1[i].y,ds1[i].z);
	//	glVertex3f(ds1[i].x+ds1[i].nx*200,ds1[i].y+ds1[i].ny*200,ds1[i].z+ds1[i].nz*200);
	//}	
	//glEnd();

	glFlush();
	glutSwapBuffers();
}

bool loadData(GBPPointset &ds, Mat &img, char *points, char *image)
{
	ifstream in(points);
	if(!in)
	{
		cout<<"Could not open "<<points<<endl;
		return true;
	}
	//Mat img1;
	img = imread(image);
	//imshow("img1",img1);
	//waitKey(0);
	uchar *pt1;
	GBPPoint point_tmp;
	point_tmp.nx=0;
	point_tmp.ny=0;
	point_tmp.nz=0;
	point_tmp.index=0;
	//cout<<point_tmp.nx<<" "<<point_tmp.ny<<endl;
	string line;
	//fp1 = fopen(points, "r");
	for (int i = 0; i<Height; i++)
	{
		pt1 = img.ptr<uchar>(i);
		for (int j = 0; j<Width; j++)
		{
			point_tmp.b = pt1[j * 3];
			point_tmp.g = pt1[j * 3 + 1];
			point_tmp.r = pt1[j * 3 + 2];
			getline(in,line);
			istringstream sin(line);
			sin>>point_tmp.x>>point_tmp.y>>point_tmp.z;
			ds.AppendPoint(point_tmp);
		}
	}
	in.close();
	return false;
	//fclose(fp1);
}

bool load_points(GICPPointSet *set, GBPPointset ds) 
{
	bool error = false;

	GICPPoint pt;
	pt.range = -1;
	for(int k = 0; k < 3; k++) 
	{
		for(int l = 0; l < 3; l++) 
		{
			pt.C[k][l] = (k == l)?1:0;
		}
	}
	for(int i=0;i<ds.SizeDS();i++)
	{
		pt.x=ds[i].x;
		pt.y=ds[i].y;
		pt.z=ds[i].z;
		set->AppendPoint(pt);
	}
	return false;  
}

int main(int argc, char *argv[])
{
	//Mat::Mat(int rows, int cols, int type, void* data, size_t step=AUTO_STEP)
	//data ¨C Pointer to the user data. Matrix constructors that take data and step parameters 
	//do not allocate matrix data. Instead, they just initialize the matrix header that points 
	//to the specified data, which means that no data is copied. This operation is very 
	//efficient and can be used to process external data using OpenCV functions. The 
	//external data is not automatically deallocated, so you should take care of it.
	//float aa[4][4];
	//Mat bb(4,4,CV_32FC1,aa);
	//for(int i=0;i<4;i++)
	//{
	//	for(int j=0;j<4;j++)
	//	{
	//		aa[i][j]=i+j*2;
	//	}
	//}
	//for(int i=0;i<4;i++)
	//{
	//	for(int j=0;j<4;j++)
	//	{
	//		cout<<bb.at<float>(i,j)<<" ";
	//	}
	//	cout<<endl;
	//}
	//return 0;

	clock_t start,end;
	bool error;
	char points[255], image[255];
	Mat img1, img2;

	int zhen = 43;
	sprintf(points, "data863_3/%d_pointCloud.txt", zhen);
	sprintf(image, "data863_3/%d_saveImage.jpg", zhen);
	cout<<"loading points for ds1..."<<endl;
	start=clock();
	error=loadData(ds1,img1,points,image);
	end=clock();
	if(error)
		return 1;
	cout<<"loaded "<<ds1.SizeOrg()<<" points for ds1 in "<<end-start<<"ms."<<endl;

	zhen = 44;
	sprintf(points, "data863_3/%d_pointCloud.txt", zhen);
	sprintf(image, "data863_3/%d_saveImage.jpg", zhen);
	cout<<"loading points for ds2..."<<endl;
	start=clock();
	error=loadData(ds2,img2,points,image);
	end=clock();
	if(error)
		return 1;
	cout<<"loaded "<<ds2.SizeOrg()<<" points for ds2 in "<<end-start<<"ms."<<endl;

	cout<<"downsampling the pointset for ds1..."<<endl;
	start=clock();
	ds1.SetDebug(debug);
	//ds1.Downsample();
	ds1.UniformSample();
	ds1.ComputeMatrices();
	end=clock();
	cout<<"time: "<<end-start<<"ms, there are "<<ds1.SizeDS()<<" points in ds1 after downsampling."<<endl;

	//return 0;

	cout<<"downsampling the pointset for ds2..."<<endl;
	start=clock();
	//ds2.Downsample();
	ds2.UniformSample();
	ds2.ComputeMatrices();
	end=clock();
	cout<<"time: "<<end-start<<"ms, there are "<<ds2.SizeDS()<<" points in ds2 after downsampling."<<endl;

	//ofstream fp1;
	//fp1.open("fp2.txt",ios::out);
	//for(int k=0;k<ds2.SizeDS();k++)
	//{
	//	fp1<<ds2[k].index<<endl;
	//	fp1<<ds2[k].x<<" "<<ds2[k].y<<" "<<ds2[k].y<<endl;
	//	fp1<<(int)ds2[k].r<<" "<<(int)ds2[k].g<<" "<<(int)ds2[k].b<<endl;
	//	fp1<<ds2[k].nx<<" "<<ds2[k].ny<<" "<<ds2[k].nz<<endl;
	//	for(int i=0;i<3;i++)
	//	{
	//		for(int j=0;j<3;j++)
	//		{
	//			fp1<<ds2[k].C[i][j]<<" ";
	//		}
	//		fp1<<endl;
	//	}
	//	fp1<<endl;
	//}
	//fp1.close();

	//*************************************************************pose estimate
	cout<<"Test pose estimate library"<<endl;
	start=clock();
	PoseEsti pe;
	pe.SetDebug(debug);
	pe.LoadImg(img1,img2);
	pe.FeatureMatch();
	pe.RansacEsti(&ds1,&ds2);
	Mat Rf(pe.GetRf());
	Mat Tf(pe.GetTf());
	end=clock();
	cout<<"pose estimate time:"<<end-start<<"ms"<<endl;
	//cout<<"Rf:"<<endl;
	//for(int i=0;i<3;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		cout<<Rf.at<float>(i,j)<<" ";
	//	}
	//	cout<<endl;
	//}
	//cout<<"Tf:"<<endl;
	//for(int i=0;i<3;i++)
	//{
	//	cout<<Tf.at<float>(i)<<" ";
	//}
	//cout<<endl;
	//cvWaitKey(0);

	//*************************************************************scanmatching
	cout<<"Test scanmatching library"<<endl;
	start=clock();
	GBPScanMatching sm;
	dgc_transform_t t_base, t0, t1;
	dgc_transform_identity(t_base);
	dgc_transform_identity(t0);
	dgc_transform_identity(t1);

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			t_base[i][j]=Rf.at<float>(i,j);
		}
	}
	for(int i=0;i<3;i++)
	{
		t_base[i][3]=Tf.at<float>(i);
	}

	//apply the transform from PoseEsti to ds1(pointset1)


	sm.SetDebug(debug);
	sm.LoadOptData(&ds1,&ds2,&pe);
	cout<<"data loaded"<<endl;
	sm.AlignScan();

	end=clock();
	cout<<"scan matching time:"<<end-start<<"ms"<<endl;

	dgc_transform_t trans;
	sm.GetTrans(trans);

	PointRansac pt1,pt2;
	double dms=0,dms1=0;
	for(int i=0;i<pe.GetSize();i++)
	{
		pt1=pe.GetPoint1(i);
		pt2=pe.GetPoint2(i);
		dgc_transform_point(&pt1.x,&pt1.y,&pt1.z,t_base);
		dms1+=sqrt((pt1.x-pt2.x)*(pt1.x-pt2.x)+(pt1.y-pt2.y)*(pt1.y-pt2.y)+(pt1.z-pt2.z)*(pt1.z-pt2.z));
		dgc_transform_point(&pt1.x,&pt1.y,&pt1.z,trans);
		dms+=sqrt((pt1.x-pt2.x)*(pt1.x-pt2.x)+(pt1.y-pt2.y)*(pt1.y-pt2.y)+(pt1.z-pt2.z)*(pt1.z-pt2.z));
	}
	dms1/=pe.GetSize();
	dms/=pe.GetSize();
	cout<<"*********************dms1="<<dms1<<endl;
	cout<<"*********************dms="<<dms<<endl;

	for(int i=0;i<ds1.SizeDS();i++)
	{
		dgc_transform_point(&ds1[i].x,&ds1[i].y,&ds1[i].z,t_base);
		dgc_transform_point(&ds1[i].x,&ds1[i].y,&ds1[i].z,trans);
	}
	for(int i=0;i<480;i++)
	{
		for(int j=0;j<640;j++)
		{
			int k=i*640+j;
			dgc_transform_point(&ds1(k).x,&ds1(k).y,&ds1(k).z,t_base);
			dgc_transform_point(&ds1(k).x,&ds1(k).y,&ds1(k).z,trans);
		}
	}


	/*
	//**************************************************************gicp
	cout << "Test program for the gicp library." << endl;
	GICPPointSet p1, p2;
	//dgc_transform_t t_base, t0, t1;
	//// set up the transformations
	//dgc_transform_identity(t_base);
	//dgc_transform_identity(t0);
	//dgc_transform_identity(t1);

	//for(int i=0;i<3;i++)
	//{
	//	for(int j=0;j<3;j++)
	//	{
	//		t_base[i][j]=Rf.at<float>(i,j);
	//	}
	//}
	//for(int i=0;i<3;i++)
	//{
	//	t_base[i][3]=Tf.at<float>(i);
	//}

	////apply the transform from PoseEsti to ds1(pointset1)
	//for(int i=0;i<ds1.SizeDS();i++)
	//{
	//	dgc_transform_point(&ds1[i].x,&ds1[i].y,&ds1[i].z,t_base);
	//}

	//cout<<"t_base"<<endl;
	//for(int i=0;i<4;i++)
	//{
	//	for(int j=0;j<4;j++)
	//	{
	//		cout<<t_base[i][j]<<" ";
	//	}
	//	cout<<endl;
	//}

	// read points clouds
	cout << "Setting up pointclouds..." << endl;
	error = load_points(&p1, ds1);
	if(error) 
		return 1;
	cout << "Loaded " << p1.Size() << " points into GICPPointSet 1." << endl;

	error = load_points(&p2, ds1);
	if(error)
		return 1;
	cout << "Loaded " << p2.Size() << " points into GICPPointSet 2." << endl; 
  
	// build kdtrees and normal matrices
	cout << "Building KDTree and computing surface normals/matrices..." << endl;
	start=clock();
	p1.SetGICPEpsilon(gicp_epsilon);
	p2.SetGICPEpsilon(gicp_epsilon);  
	p1.BuildKDTree();
	p1.ComputeMatrices();
	p2.BuildKDTree();
	p2.ComputeMatrices();
	end=clock();
	cout<<"time="<<end-start<<endl;

	if(debug)
	{
		// save data for debug/visualizations
		p1.SavePoints("pts1.dat");
		p1.SaveMatrices("mats1.dat");
		p2.SavePoints("pts2.dat");
		p2.SaveMatrices("mats2.dat");
	}

	// align the point clouds
	cout << "Aligning point cloud..." << endl;
	start=clock();
	dgc_transform_copy(t1, t0);
	p2.SetDebug(debug);
	p2.SetMaxIterationInner(8);
	p2.SetMaxIteration(100);
	int iterations = p2.AlignScan(&p1, t_base, t1, max_distance);
	end=clock();

	// print the result
	cout << "Converged in " <<end-start<<" ms"<< endl;
	dgc_transform_print(t_base, "t_base");
	dgc_transform_print(t0, "t0");  
	dgc_transform_print(t1, "t1");

	if(debug) 
	{
		ofstream fout("iterations.txt");
		if(!fout) {
			return 0;
	}
	fout << "Converged in " << iterations << " iterations." << endl;
	fout.close();
	}
	*/


	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(10, 200);
	glutInitWindowSize(640, 480);
	glutCreateWindow("3D image");

	glutDisplayFunc(renderScene);
	while (true)
	{
		if(flag_endloop)
		{
			flag_endloop=false;
			break;
		}
		glutReshapeFunc(reshape);
		glutMouseFunc(mouse);
		glutSpecialFunc(key);
		glutMotionFunc(motion);
		glutPostRedisplay();
		glutMainLoopEvent();
	}
	return 0;
}