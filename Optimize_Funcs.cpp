#include "ScanMatching.h"

namespace dgc
{
	namespace matching
	{
		// GICP cost function
		double GBPScanMatching::f(const gsl_vector *x, void *params) 
		{
			GBPOptData *opt_data = (GBPOptData *)params;
			double pt1[3];
			double pt2[3]; 
			double res[3]; // residual
			double temp[3];
			gsl_vector_view gsl_pt1 = gsl_vector_view_array(pt1, 3);
			gsl_vector_view gsl_pt2 = gsl_vector_view_array(pt2, 3);
			gsl_vector_view gsl_res = gsl_vector_view_array(res, 3);
			gsl_vector_view gsl_temp = gsl_vector_view_array(temp, 3);
			gsl_matrix_view gsl_M;
			dgc_transform_t t;

			// initialize the temp variable; if it happens to be NaN at start, bad things will happen in blas routines below
			temp[0] = 0;
			temp[1] = 0;
			temp[2] = 0;
      

			// take the base transformation
			dgc_transform_copy(t, opt_data->base_t); 
			// apply the current state
			apply_state(t, x);
            
			double f = 0;
			double temp_double = 0;
			int N = opt_data->ps1->SizeDS();
			for(int i = 0; i < N; i++) 
			{
				int j = opt_data->nn_indecies[i];	
				if(j != -1) 
				{
					// get point 1
					pt1[0] = (*opt_data->ps1)[i].x;
					pt1[1] = (*opt_data->ps1)[i].y;
					pt1[2] = (*opt_data->ps1)[i].z;
					// get point 2
					pt2[0] = (*opt_data->ps2)[j].x;
					pt2[1] = (*opt_data->ps2)[j].y;
					pt2[2] = (*opt_data->ps2)[j].z;
					//get M-matrix
					gsl_M = gsl_matrix_view_array(&opt_data->M[i][0][0], 3, 3);
	  
	  
					//transform point 1
					dgc_transform_point(&pt1[0], &pt1[1], &pt1[2], t);
					res[0] = pt1[0] - pt2[0];
					res[1] = pt1[1] - pt2[1];
					res[2] = pt1[2] - pt2[2];

					// temp := M*res
					gsl_blas_dsymv(CblasLower, 1., &gsl_M.matrix, &gsl_res.vector, 0., &gsl_temp.vector);
					// temp_double := res'*temp = res'*M*res
					gsl_blas_ddot(&gsl_res.vector, &gsl_temp.vector, &temp_double);
	  
					f += temp_double/(double)opt_data->num_matches;	  
				}
			}
      
			return f;
		}

		void GBPScanMatching::df(const gsl_vector *x, void *params, gsl_vector *g) 
		{
			GBPOptData *opt_data = (GBPOptData *)params;
			double pt1[3];
			double pt2[3]; 
			double res[3]; // residual
			double temp[3]; // temp local vector
			double temp_mat[9]; // temp matrix used for accumulating the rotation gradient
			gsl_vector_view gsl_pt1 = gsl_vector_view_array(pt1, 3);
			gsl_vector_view gsl_pt2 = gsl_vector_view_array(pt2, 3);
			gsl_vector_view gsl_res = gsl_vector_view_array(res, 3);
			gsl_vector_view gsl_temp = gsl_vector_view_array(temp, 3);
			gsl_vector_view gsl_gradient_t = gsl_vector_subvector(g, 0, 3); // translation comp. of gradient
			gsl_matrix_view gsl_temp_mat_r = gsl_matrix_view_array(temp_mat, 3, 3);
			gsl_matrix_view gsl_M;
			dgc_transform_t t;
			double temp_double;
      
			// take the base transformation
			dgc_transform_copy(t, opt_data->base_t); 
			// apply the current state
			apply_state(t, x);
      
			// zero all accumulator variables
			gsl_vector_set_zero(g);
			gsl_vector_set_zero(&gsl_temp.vector);
			gsl_matrix_set_zero(&gsl_temp_mat_r.matrix);
            
			for(int i = 0; i < opt_data->ps1->SizeDS(); i++) 
			{
				int j = opt_data->nn_indecies[i];	
				if(j != -1) 
				{
					// get point 1
					pt1[0] = (*opt_data->ps1)[i].x;
					pt1[1] = (*opt_data->ps1)[i].y;
					pt1[2] = (*opt_data->ps1)[i].z;
	  
					// get point 2
					pt2[0] = (*opt_data->ps2)[j].x;
					pt2[1] = (*opt_data->ps2)[j].y;
					pt2[2] = (*opt_data->ps2)[j].z;
	  
					//get M-matrix
					gsl_M = gsl_matrix_view_array(&opt_data->M[i][0][0], 3, 3);	  

					//transform point 1
					dgc_transform_point(&pt1[0], &pt1[1], &pt1[2], t);
					res[0] = pt1[0] - pt2[0];
					res[1] = pt1[1] - pt2[1];
					res[2] = pt1[2] - pt2[2];
	  
					// temp := M*res
					gsl_blas_dsymv(CblasLower, 1., &gsl_M.matrix, &gsl_res.vector, 0., &gsl_temp.vector);
					// temp_double := res'*temp = res'*M*res
					gsl_blas_ddot(&gsl_res.vector, &gsl_temp.vector, &temp_double);

					// increment total translation gradient:
					// gsl_gradient_t += 2*M*res/num_matches
					gsl_blas_dsymv(CblasLower, 2./(double)opt_data->num_matches, &gsl_M.matrix, &gsl_res.vector, 1., &gsl_gradient_t.vector);	  

					// compute rotation gradient here
					// get back the original untransformed point to compute the rotation gradient
					pt1[0] = (*opt_data->ps1)[i].x;
					pt1[1] = (*opt_data->ps1)[i].y;
					pt1[2] = (*opt_data->ps1)[i].z;
					dgc_transform_point(&pt1[0], &pt1[1], &pt1[2], opt_data->base_t);
					//gsl_temp_mat_r += 2*pt1*temp'/num_matches
					gsl_blas_dger(2./(double)opt_data->num_matches, &gsl_pt1.vector, &gsl_temp.vector, &gsl_temp_mat_r.matrix);
				
				}
			}
			compute_dr(x, &gsl_temp_mat_r.matrix, g);
		}
    
		void GBPScanMatching::fdf(const gsl_vector *x, void *params, double * f, gsl_vector *g) 
		{
			GBPOptData *opt_data = (GBPOptData *)params;
			double pt1[3];
			double pt2[3]; 
			double res[3]; // residual
			double temp[3]; // temp local vector
			double temp_mat[9]; // temp matrix used for accumulating the rotation gradient
			gsl_vector_view gsl_pt1 = gsl_vector_view_array(pt1, 3);
			gsl_vector_view gsl_pt2 = gsl_vector_view_array(pt2, 3);
			gsl_vector_view gsl_res = gsl_vector_view_array(res, 3);
			gsl_vector_view gsl_temp = gsl_vector_view_array(temp, 3);
			gsl_vector_view gsl_gradient_t = gsl_vector_subvector(g, 0, 3); // translation comp. of gradient
			gsl_vector_view gsl_gradient_r = gsl_vector_subvector(g, 3, 3); // rotation comp. of gradient
			gsl_matrix_view gsl_temp_mat_r = gsl_matrix_view_array(temp_mat, 3, 3);
			gsl_matrix_view gsl_M;
			dgc_transform_t t;
			double temp_double;

			// take the base transformation
			dgc_transform_copy(t, opt_data->base_t); 
			// apply the current state      
			apply_state(t, x);
            
			// zero all accumulator variables
			*f = 0;
			gsl_vector_set_zero(g);
			gsl_vector_set_zero(&gsl_temp.vector);
			gsl_matrix_set_zero(&gsl_temp_mat_r.matrix);

			ofstream fp;
			fp.open("fdf.txt",ios::out);

      
			for(int i = 0; i < opt_data->ps1->SizeDS(); i++) 
			{
				int j = opt_data->nn_indecies[i];	
				if(j != -1) 
				{
					// get point 1
					pt1[0] = (*opt_data->ps1)[i].x;
					pt1[1] = (*opt_data->ps1)[i].y;
					pt1[2] = (*opt_data->ps1)[i].z;
					// get point 2
					pt2[0] = (*opt_data->ps2)[j].x;
					pt2[1] = (*opt_data->ps2)[j].y;
					pt2[2] = (*opt_data->ps2)[j].z;
	  
					//get M-matrix
					gsl_M = gsl_matrix_view_array(&opt_data->M[i][0][0], 3, 3);	  
	  
					//transform point 1
					dgc_transform_point(&pt1[0], &pt1[1], &pt1[2], t);
					res[0] = pt1[0] - pt2[0];
					res[1] = pt1[1] - pt2[1];
					res[2] = pt1[2] - pt2[2];

					fp<<pt1[0]<<" "<<pt1[1]<<" "<<pt1[2]<<endl;
					fp<<pt2[0]<<" "<<pt2[1]<<" "<<pt2[2]<<endl;
					fp<<res[0]<<" "<<res[1]<<" "<<res[2]<<endl<<endl;
	  
					// compute the transformed residual
					// temp := M*res
					gsl_blas_dsymv(CblasLower, 1., &gsl_M.matrix, &gsl_res.vector, 0., &gsl_temp.vector);
	  
					// compute M-norm of the residual
					// temp_double := res'*temp = res'*M*res
					gsl_blas_ddot(&gsl_res.vector, &gsl_temp.vector, &temp_double);
	  
					// accumulate total error: f += res'*M*res
					*f += temp_double/(double)opt_data->num_matches;
	  
					// accumulate translation gradient:
					// gsl_gradient_t += 2*M*res
					gsl_blas_dsymv(CblasLower, 2./(double)opt_data->num_matches, &gsl_M.matrix, &gsl_res.vector, 1., &gsl_gradient_t.vector);	  
					
					// accumulate the rotation gradient matrix
					// get back the original untransformed point to compute the rotation gradient
					pt1[0] = (*opt_data->ps1)[i].x;
					pt1[1] = (*opt_data->ps1)[i].y;
					pt1[2] = (*opt_data->ps1)[i].z;
					dgc_transform_point(&pt1[0], &pt1[1], &pt1[2], opt_data->base_t);
					// gsl_temp_mat_r += 2*(gsl_temp).(gsl_pt1)' [ = (2*M*residual).(gsl_pt1)' ]	  
					gsl_blas_dger(2./(double)opt_data->num_matches, &gsl_pt1.vector, &gsl_temp.vector, &gsl_temp_mat_r.matrix); 
				}
			}      
			compute_dr(x, &gsl_temp_mat_r.matrix, g);

			fp.close();

		}
	}
}
