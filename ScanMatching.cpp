#include "ScanMatching.h"

namespace dgc
{
	namespace matching
	{
		GBPScanMatching::GBPScanMatching() : 
			T_min(gsl_multimin_fdfminimizer_vector_bfgs)
		{
			kdtree_done=false;
			kdtree = NULL;
			max_d_sq=25;
			debug_=false;
			max_iter = 100;
			max_iteration=100;
			iter = -1;
			status = GSL_CONTINUE;
			data_loaded=false;
			dgc_transform_identity(t);
			epsilon=5e-3;
			epsilon_rot=2e-2;

			//T_min=gsl_multimin_fdfminimizer_vector_bfgs;
			gsl_minimizer = gsl_multimin_fdfminimizer_alloc(T_min, N);
			if(gsl_minimizer == NULL)
				return;
			x = gsl_vector_alloc(N);
			if(x == NULL)
				return;
		}

		GBPScanMatching::~GBPScanMatching()
		{
		}

		void GBPScanMatching::Clear()
		{
			kdtree_done = false;
			if (kdtree != NULL) {
				delete kdtree;
				kdtree = NULL;
			}
		}

		void GBPScanMatching::GetTrans(dgc_transform_t trans)
		{
			for(int i=0;i<4;i++)
			{
				for(int j=0;j<4;j++)
				{
					trans[i][j]=t[i][j];
				}
			}
		}

		void GBPScanMatching::LoadOptData(GBPPointset *ds1, GBPPointset *ds2, PoseEsti *pe)
		{
			if(data_loaded)
				return;
			data_loaded=true;

			opt_data.ps1=ds1;
			opt_data.ps2=ds2;
			dgc_transform_identity(opt_data.base_t);
			for(int i=0;i<3;i++)
			{
				for(int j=0;j<3;j++)
				{
					opt_data.base_t[i][j]=pe->GetRf().at<float>(i,j);
				}
			}
			for(int i=0;i<3;i++)
			{
				opt_data.base_t[i][3]=pe->GetTf().at<float>(i);
			}
		}

		int GBPScanMatching::AlignScan()
		{
			if(!data_loaded)
			{
				cout<<"data not loaded!"<<endl;
				return 0;
			}

			if(!kdtree_done)
			{
				BuildKDTree();
			}

			ofstream fp_corres;
			if(debug_)
			{
				fp_corres.open("correspondence.txt",ios::out);
			}

			int n=opt_data.ps1->SizeDS();
			ANNpoint query_point=annAllocPt(3);
			ANNdist nn_dist_sq;
			opt_data.nn_indecies=new ANNidx[n];
			int num_matches;
			gsl_matrix *gsl_temp=gsl_matrix_alloc(3,3);
			if(gsl_temp==NULL)
				return -1;
			gsl_matrix *gsl_R = gsl_matrix_alloc(3, 3);
			if (gsl_R == NULL)
				return -1;
			opt_data.M=new mat_t[n];
			bool converged=false;
			int iteration=0;
			bool opt_status=false;
			dgc_transform_t t_last,transform_R;
			double delta=0;

			for(int i=0;i<n;i++)
				for(int k=0;k<3;k++)
					for(int l=0;l<3;l++)
						opt_data.M[i][k][l]=(k==l)?1:0;

			if (debug_) 
			{
				dgc_transform_write(opt_data.base_t, "t_base.tfm");
				dgc_transform_write(t, "t_0.tfm");
			}

			while(!converged)
			{
				dgc_transform_copy(transform_R,opt_data.base_t);
				dgc_transform_left_multiply(transform_R,t);
				for(int i=0;i<3;i++)
					for(int j=0;j<3;j++)
						gsl_matrix_set(gsl_R,i,j,transform_R[i][j]);

				num_matches=0;
				for(int i=0;i<n;i++)
				{
					query_point[0]=(*opt_data.ps1)[i].x;
					query_point[1]=(*opt_data.ps1)[i].y;
					query_point[2]=(*opt_data.ps1)[i].z;

					dgc_transform_point(&query_point[0], &query_point[1], &query_point[2], opt_data.base_t);
					dgc_transform_point(&query_point[0], &query_point[1], &query_point[2], t);

					kdtree->annkSearch(query_point, 1, &opt_data.nn_indecies[i], &nn_dist_sq, 0.0);

					if(nn_dist_sq<10000)
					{
						if(debug_)
						{
							fp_corres<<i<<"\t"<<opt_data.nn_indecies[i]<<endl;
						}
						//**********************************                                                                     
						//gsl_matrix_view C1=gsl_matrix_view_array(&(*opt_data.ps1)[i].C[0][0],3,3);
						//gsl_matrix_view C2=gsl_matrix_view_array(&(*opt_data.ps2)[opt_data.nn_indecies[i]].C[0][0],3,3);
						//gsl_matrix_view M=gsl_matrix_view_array(&opt_data.M[i][0][0],3,3);

						//gsl_matrix_set_zero(&M.matrix);
						//gsl_matrix_set_zero(gsl_temp);

						////M=R*C1
						//gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,gsl_R, &C1.matrix,1.,&M.matrix);
						////temp=R*R'
						//gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., &M.matrix, gsl_R, 0., gsl_temp);
						////temp+=C2
						//gsl_matrix_add(gsl_temp, &C2.matrix);
						//gsl_matrix_set_identity(&M.matrix);

						//gsl_linalg_cholesky_decomp(gsl_temp);
						//for (int k = 0; k < 3; k++) 
						//{
						//	gsl_vector_view row_view = gsl_matrix_row(&M.matrix, k);
						//	gsl_linalg_cholesky_svx(gsl_temp, &row_view.vector);
						//}
						//****************************
						num_matches++;
					}
					else
					{
						opt_data.nn_indecies[i]=-1;
					}
				}

				if (debug_) 
				{
					// save the current M matrices to file for debugging
					ofstream out("mahalanobis.txt");
					if (out) 
					{
						for (int i = 0; i < n; i++) 
						{
							for (int k = 0; k < 3; k++) 
							{
								for (int l = 0; l < 3; l++) 
								{
									out << opt_data.M[i][k][l] << "\t";
								}
								out<<endl;
							}
							out << endl;
						}
					}
					out.close();
				}

				opt_data.num_matches=num_matches;
				dgc_transform_copy(t_last,t);
				opt_status=Optimize();//t

				double x,y,z,rsme=0;
				int num=0;
				for(int i=0;i<n;i++)
				{
					int j = opt_data.nn_indecies[i];
					if(j != -1) 
					{
						x=(*opt_data.ps1)[i].x;
						y=(*opt_data.ps1)[i].y;
						z=(*opt_data.ps1)[i].z;
						dgc_transform_point(&x,&y,&z,opt_data.base_t);
						dgc_transform_point(&x,&y,&z,t);
						rsme+=sqrt((x-(*opt_data.ps2)[j].x)*(x-(*opt_data.ps2)[j].x)
							+(y-(*opt_data.ps2)[j].y)*(y-(*opt_data.ps2)[j].y)
							+(z-(*opt_data.ps2)[j].z)*(z-(*opt_data.ps2)[j].z));
						num++;
					}
				}
				rsme/=num;
				

				if (debug_) 
				{
					cout << "Optimizer converged in " << iter << " iterations." << endl;
					cout<<"num_matches="<<num_matches<<endl;
					cout << "Status: " << status << endl;
					cout<< "rsme="<<rsme<<endl;

					std::ostringstream filename;
					filename << "t_" << iteration + 1 << ".tfm";
					dgc_transform_write(t, filename.str().c_str());
				}

				delta = 0.;
				for (int k = 0; k < 4; k++) 
				{
					for (int l = 0; l < 4; l++) 
					{
						double ratio = 1;
						//if (k < 3 && l < 3) 
						//{ // rotation part of the transform
						//	ratio = 1. / epsilon_rot;//2e-3;
						//}
						//else 
						//{
						//	ratio = 1. / epsilon;//5e-4;
						//}
						double c_delta = ratio*fabs(t_last[k][l] - t[k][l]);

						if (c_delta > delta) 
						{
							delta = c_delta;
						}
					}
				}

				if (debug_) {
					cout << "delta = " << delta << endl;
				}

				iteration++;
				if (iteration >= max_iteration || delta < 0.5) {
					converged = true;
				}
			}

			if(debug_)
			{
				cout << "Converged in " << iteration << " iterations." << endl;
				//fp_corres<<"num_matches="<<num_matches<<endl;
				fp_corres.close();
			}

			if (gsl_R != NULL) {
				gsl_matrix_free(gsl_R);
			}
			if (gsl_temp != NULL) {
				gsl_matrix_free(gsl_temp);
			}
			annDeallocPt(query_point);

			return iteration;
		}

		void GBPScanMatching::BuildKDTree()
		{
			if (kdtree_done)
			{
				cout<<"kdtree error"<<endl;
				return;
			}
			kdtree_done = true;

			int n = opt_data.ps2->SizeDS();

			if (n == 0)
			{
				cout<<"no point in ps2"<<endl;
				return;
			}

			ANNpointArray kdtree_points = annAllocPts(n, 3);
			for (int i = 0; i < n; i++)
			{
				kdtree_points[i][0] = (*opt_data.ps2)[i].x;
				kdtree_points[i][1] = (*opt_data.ps2)[i].y;
				kdtree_points[i][2] = (*opt_data.ps2)[i].z;
			}
			kdtree = new ANNkd_tree(kdtree_points, n, 3);
			
			if(debug_)
				cout<<"there are "<<kdtree->nPoints()<<"points in kdtree"<<endl;
		}

		bool GBPScanMatching::Optimize()
		{
			ofstream fp;
			if(debug_)
			{
				fp.open("optimize.txt",ios::app);
			}

			double line_search_tol = 0.5;
			double gradient_tol = 10.;
			double step_size = 20.;

			gsl_multimin_function_fdf func;
			func.f = f;
			func.df = df;
			func.fdf = fdf;
			func.n = N;
			func.params = &opt_data;

			//dgc_transform_copy(t,opt_data.base_t);
			//
			//if(debug_)
			//{
			//	cout<<"t before optimize:"<<endl;
			//	for(int i=0;i<4;i++)
			//	{
			//		for(int j=0;j<4;j++)
			//		{
			//			cout<<t[i][j]<<" ";
			//		}
			//		cout<<endl;
			//	}
			//	cout<<endl;
			//}

			double tx,ty,tz,rx,ry,rz;
			dgc_transform_get_translation(t, &tx, &ty, &tz);
			dgc_transform_get_rotation(t, &rx, &ry, &rz);
			gsl_vector_set(x, 0, tx);
			gsl_vector_set(x, 1, ty);
			gsl_vector_set(x, 2, tz);
			gsl_vector_set(x, 3, rx);
			gsl_vector_set(x, 4, ry);
			gsl_vector_set(x, 5, rz);
      
			// initialize the minimizer
			gsl_multimin_fdfminimizer_set(gsl_minimizer, &func, x, step_size, line_search_tol);
      
			//iterate the minimization algorithm using gsl primatives
			status = GSL_CONTINUE;
			iter = 0;

			if(debug_)
			{
				fp<<endl<<"iter\tf\ttx\tty\ttz\trx\try\trz\tgtx\tgty\tgtz\tgrx\tgry\tgrz"<<endl;
			}

			while(status == GSL_CONTINUE && iter < max_iter) 
			{
				iter++;
				status = gsl_multimin_fdfminimizer_iterate(gsl_minimizer);
				if(status) 
				{
					break;
				}
				status = gsl_multimin_test_gradient (gsl_minimizer->gradient, gradient_tol);

				if(debug_)
				{
					fp<<iter<<" "<<gsl_minimizer->f<<" "
						<<gsl_vector_get(gsl_minimizer->x,0)<<" "
						<<gsl_vector_get(gsl_minimizer->x,1)<<" "
						<<gsl_vector_get(gsl_minimizer->x,2)<<" "
						<<gsl_vector_get(gsl_minimizer->x,3)<<" "
						<<gsl_vector_get(gsl_minimizer->x,4)<<" "
						<<gsl_vector_get(gsl_minimizer->x,5)<<" "
						<<gsl_vector_get(gsl_minimizer->gradient,0)<<" "
						<<gsl_vector_get(gsl_minimizer->gradient,1)<<" "
						<<gsl_vector_get(gsl_minimizer->gradient,2)<<" "
						<<gsl_vector_get(gsl_minimizer->gradient,3)<<" "
						<<gsl_vector_get(gsl_minimizer->gradient,4)<<" "
						<<gsl_vector_get(gsl_minimizer->gradient,5)<<" "<<endl;
				}
			}
      
			if(status == GSL_SUCCESS || iter == max_iter) 
			{
				//set t to the converged solution
	
				dgc_transform_identity(t);
				// apply the current state to the base
				apply_state(t, gsl_minimizer->x);
				//dgc_transform_print(t, "converged to:");	

				return true;
			}
			else 
			{
				// the algorithm failed to converge
				return false;
			}

			if(debug_)
			{
				fp.close();
			}


			//if(debug_)
			//{
			//	cout<<"t after optimize:"<<endl;
			//	for(int i=0;i<4;i++)
			//	{
			//		for(int j=0;j<4;j++)
			//		{
			//			cout<<t[i][j]<<" ";
			//		}
			//		cout<<endl;
			//	}
			//	cout<<endl;
			//}

		}

		double GBPScanMatching::mat_inner_prod(gsl_matrix const* mat1, gsl_matrix const* mat2)
		{
			double r = 0.;
			int n = mat1->size1;
      
			for(int i = 0; i < n; i++) 
			{
				for(int j = 0; j < n; j++) 
				{ // tr(mat1^t.mat2)
					r += gsl_matrix_get(mat1, j, i)*gsl_matrix_get(mat2, i, j);
				}
			}
      
			return r;
		}
    
		void GBPScanMatching::apply_state(dgc_transform_t t, gsl_vector const* x) 
		{
			double tx, ty, tz, rx, ry, rz;
			tx = gsl_vector_get(x, 0);
			ty = gsl_vector_get(x, 1);
			tz = gsl_vector_get(x, 2);
			rx = gsl_vector_get(x, 3);
			ry = gsl_vector_get(x, 4);
			rz = gsl_vector_get(x, 5);

			dgc_transform_rotate_x(t, rx);
			dgc_transform_rotate_y(t, ry);
			dgc_transform_rotate_z(t, rz);      
			dgc_transform_translate(t, tx, ty, tz);
		}
    
		void GBPScanMatching::compute_dr(gsl_vector const* x, gsl_matrix const* gsl_temp_mat_r, gsl_vector *g) 
		{
			double dR_dPhi[3][3];
			double dR_dTheta[3][3];
			double dR_dPsi[3][3];
			gsl_matrix_view gsl_d_rx = gsl_matrix_view_array(&dR_dPhi[0][0],3, 3);
			gsl_matrix_view gsl_d_ry = gsl_matrix_view_array(&dR_dTheta[0][0],3, 3);
			gsl_matrix_view gsl_d_rz = gsl_matrix_view_array(&dR_dPsi[0][0],3, 3);

			double phi = gsl_vector_get(x ,3);
			double theta = gsl_vector_get(x ,4);
			double psi = gsl_vector_get(x ,5);  
      
			double cphi = cos(phi), sphi = sin(phi);
			double ctheta = cos(theta), stheta = sin(theta);
			double cpsi = cos(psi), spsi = sin(psi);
      
			dR_dPhi[0][0] = 0.;
			dR_dPhi[1][0] = 0.;
			dR_dPhi[2][0] = 0.;
      
			dR_dPhi[0][1] = sphi*spsi + cphi*cpsi*stheta;
			dR_dPhi[1][1] = -cpsi*sphi + cphi*spsi*stheta;
			dR_dPhi[2][1] = cphi*ctheta;
      
			dR_dPhi[0][2] = cphi*spsi - cpsi*sphi*stheta;
			dR_dPhi[1][2] = -cphi*cpsi - sphi*spsi*stheta;
			dR_dPhi[2][2] = -ctheta*sphi;
      
			dR_dTheta[0][0] = -cpsi*stheta;
			dR_dTheta[1][0] = -spsi*stheta;
			dR_dTheta[2][0] = -ctheta;
      
			dR_dTheta[0][1] = cpsi*ctheta*sphi;
			dR_dTheta[1][1] = ctheta*sphi*spsi;
			dR_dTheta[2][1] = -sphi*stheta;
	
			dR_dTheta[0][2] = cphi*cpsi*ctheta;
			dR_dTheta[1][2] = cphi*ctheta*spsi;
			dR_dTheta[2][2] = -cphi*stheta;
      
			dR_dPsi[0][0] = -ctheta*spsi;
			dR_dPsi[1][0] = cpsi*ctheta;
			dR_dPsi[2][0] = 0.;
      
			dR_dPsi[0][1] = -cphi*cpsi - sphi*spsi*stheta;
			dR_dPsi[1][1] = -cphi*spsi + cpsi*sphi*stheta;
			dR_dPsi[2][1] = 0.;
      
			dR_dPsi[0][2] = cpsi*sphi - cphi*spsi*stheta;
			dR_dPsi[1][2] = sphi*spsi + cphi*cpsi*stheta;
			dR_dPsi[2][2] = 0.;
      
			// set d/d_rx = tr(dR_dPhi'*gsl_temp_mat_r) [= <dR_dPhi, gsl_temp_mat_r>]
			gsl_vector_set(g, 3, mat_inner_prod(&gsl_d_rx.matrix, gsl_temp_mat_r));
			// set d/d_ry = tr(dR_dTheta'*gsl_temp_mat_r) = [<dR_dTheta, gsl_temp_mat_r>]
			gsl_vector_set(g, 4, mat_inner_prod(&gsl_d_ry.matrix, gsl_temp_mat_r));
			// set d/d_rz = tr(dR_dPsi'*gsl_temp_mat_r) = [<dR_dPsi, gsl_temp_mat_r>]
			gsl_vector_set(g, 5, mat_inner_prod(&gsl_d_rz.matrix, gsl_temp_mat_r));      

		}

	}
}