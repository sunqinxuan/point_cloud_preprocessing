//#ifndef DGC_TRANSFORM_H
//#define DGC_TRANSFORM_H
//
//#ifdef __cplusplus
////extern "C" {
//#endif

#define M_PI 3.1415926

  /** 3D transformation - can represent a combined translation and rotation */

typedef double dgc_transform_t[4][4];

  /** Print a 3D transform as a 4x4 matrix 
      @param t - transform to be printed 
      &param str - name of transform, which will be printed above matrix */

void dgc_transform_print(dgc_transform_t t, const char *str);

  /** Initializes a transform as the identity transform */

void dgc_transform_identity(dgc_transform_t t);

  /** Left multiply transform t1 by transform t2 */

void dgc_transform_left_multiply(dgc_transform_t t1, dgc_transform_t t2);

  /** Rotate tranform around global x axis */

void dgc_transform_rotate_x(dgc_transform_t t, double theta);

  /** Rotate transform around global y axis */

void dgc_transform_rotate_y(dgc_transform_t t, double theta);

  /** Rotate transform around global z axis */

void dgc_transform_rotate_z(dgc_transform_t t, double theta);

  /** Add translation to transform */

void dgc_transform_translate(dgc_transform_t t, double x, double y, double z);

  /** Apply transform to 3D point */

//extern inline 
void dgc_transform_point(double *x, double *y, double *z, dgc_transform_t t);

void dgc_transform_point(float *x, float *y, float *z, dgc_transform_t t);

  /** Copy transform src to dest */

void dgc_transform_copy(dgc_transform_t dest, dgc_transform_t src);

  /** Read a transform from file */

//int dgc_transform_read(dgc_transform_t t, const char *filename);

  /** Write a transform from file */

int dgc_transform_write(dgc_transform_t t, const char *filename);

  /** Read a transform from a string */

//int dgc_transform_read_string(dgc_transform_t t, char *str);

  /** Get translation */

void dgc_transform_get_translation(dgc_transform_t t, double *x, double *y,
				   double *z);

  /** Get rotation */

void dgc_transform_get_rotation(dgc_transform_t t, double *x, double *y,
				double *z);

void dgc_transform_rpy(dgc_transform_t dest, dgc_transform_t src, double roll,
		       double pitch, double yaw);

void dgc_transform_inverse(dgc_transform_t in, dgc_transform_t out);

//#ifdef __cplusplus
////}
//#endif
//
//#endif