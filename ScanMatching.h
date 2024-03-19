#ifndef SCANMATCH_H_
#define SCANMATCH_H_

#include "ANN.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "Pointset.h"
#include "transform.h"
#include "PoseEsti.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit_nlin.h>

using namespace std;

namespace dgc
{
	namespace matching
	{
		typedef struct GBPOptData 
		{
			GBPPointset *ps1;
			GBPPointset *ps2;
			//ANNidx *nn_indecies; // nearest point indecies
			ANNidxArray nn_indecies;
			mat_t *M;      // mahalanobis matrices for each pair
			dgc_transform_t base_t;
			int num_matches;
		}GBPOptData;

		class GBPScanMatching
		{
		public:
			GBPScanMatching();
			~GBPScanMatching();
			void SetDebug(bool d) {debug_ = d;}
			void SetMax_d_sq(double max_match_dist) {max_d_sq=pow(max_match_dist,2);}
			void GetTrans(dgc_transform_t trans);

			void Clear();
			void LoadOptData(GBPPointset *ps1, GBPPointset *ps2, PoseEsti *pe);
			void BuildKDTree();//for ps2
			bool Optimize();

			int AlignScan();

		private:
			static double f(const gsl_vector * x, void * params);
			static void df(const gsl_vector * x, void * params, gsl_vector * g);
			static void fdf(const gsl_vector * x, void * params, double * f, gsl_vector * g);

			static void compute_dr(gsl_vector const* x, gsl_matrix const* gsl_temp_mat_r, gsl_vector *g);
			static double mat_inner_prod(gsl_matrix const* mat1, gsl_matrix const* mat2);
			static void apply_state(dgc_transform_t t, gsl_vector const* x);

			GBPOptData opt_data;
			dgc_transform_t t;
			bool debug_;
			bool kdtree_done;
			bool data_loaded;
			ANNkd_tree *kdtree;
			double max_d_sq;
			double epsilon,epsilon_rot;
			int max_iter,max_iteration;
			int status;
			int iter;

			gsl_multimin_fdfminimizer *gsl_minimizer;
			gsl_vector *x;
			const static int N = 6;
			const gsl_multimin_fdfminimizer_type *T_min;
		};
	}
}

#endif