/*************************************************************
  Generalized-ICP Copyright (c) 2009 Aleksandr Segal.
  All rights reserved.

  Redistribution and use in source and binary forms, with 
  or without modification, are permitted provided that the 
  following conditions are met:

* Redistributions of source code must retain the above 
  copyright notice, this list of conditions and the 
  following disclaimer.
* Redistributions in binary form must reproduce the above
  copyright notice, this list of conditions and the 
  following disclaimer in the documentation and/or other
  materials provided with the distribution.
* The names of the contributors may not be used to endorse
  or promote products derived from this software
  without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
  OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGE.
*************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>

// program options
//#include <boost/program_options.hpp>
//#include <boost/filesystem/path.hpp>
//#include <boost/filesystem/fstream.hpp>
//#include <boost/tokenizer.hpp> 

#include "gicp.h"

using namespace std;
using namespace dgc::gicp;
//namespace po = boost::program_options;


static string filename1="0_pointCloud.txt";
static string filename2="1_pointCloud.txt";
static string filename_t_base;
static bool debug = true;
static double gicp_epsilon = 1e-3;
static double max_distance = 5.;


bool load_points(GICPPointSet *set, const char* filename) {
  bool error = false;
  
  ifstream in(filename);
  if(!in) {
    cout << "Could not open '" << filename << "'." << endl;
    return true;
  }
  string line;
  GICPPoint pt;
  pt.range = -1;
  for(int k = 0; k < 3; k++) {
    for(int l = 0; l < 3; l++) {
      pt.C[k][l] = (k == l)?1:0;
    }
  }    
  while(getline(in, line)) {
    istringstream sin(line);
    sin >> pt.x >> pt.y >> pt.z;   
	if(pt.x==0&&pt.y==0&&pt.z==0)
		continue;
	else
		set->AppendPoint(pt);    
  }
  in.close();

  return false;  
};


int main(int argc, char** argv) {


	clock_t start,finish;

  cout << "Test program for the gicp library." << endl;

  bool error;

  GICPPointSet p1, p2;
  dgc_transform_t t_base, t0, t1;

  // set up the transformations
  dgc_transform_identity(t_base);
  dgc_transform_identity(t0);
  dgc_transform_identity(t1);
  
  // read base transform from file if one is specified
  //if(!filename_t_base.empty()) {
  //  int status = dgc_transform_read(t_base, filename_t_base.c_str());
  //  if(status != 0) {
  //    return 1;
  //  }
  //}
  
  // read points clouds
  cout << "Setting up pointclouds..." << endl;
  error = load_points(&p1, filename1.c_str());
  if(error) {
    return 1;
  }
  cout << "Loaded " << p1.Size() << " points into GICPPointSet 1." << endl;
  error = load_points(&p2, filename2.c_str());
  if(error) {
    return 1;
  }
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
  finish=clock();
  cout<<"time="<<finish-start<<endl;

  if(debug) {
    // save data for debug/visualizations
    p1.SavePoints("pts1.dat");
    p1.SaveMatrices("mats1.dat");
    p2.SavePoints("pts2.dat");
    p2.SaveMatrices("mats2.dat");
  }
  
  // align the point clouds
  cout << "Aligning point cloud..." << endl;
  dgc_transform_copy(t1, t0);
  p2.SetDebug(debug);
  p2.SetMaxIterationInner(8);
  p2.SetMaxIteration(100);
  int iterations = p2.AlignScan(&p1, t_base, t1, max_distance);
  
  // print the result
  cout << "Converged: " << endl;
  dgc_transform_print(t_base, "t_base");
  dgc_transform_print(t0, "t0");  
  dgc_transform_print(t1, "t1");

  if(debug) {
    ofstream fout("iterations.txt");
    if(!fout) {
      return 0;
    }
    fout << "Converged in " << iterations << " iterations." << endl;
    fout.close();
  }

  return 0;
}
