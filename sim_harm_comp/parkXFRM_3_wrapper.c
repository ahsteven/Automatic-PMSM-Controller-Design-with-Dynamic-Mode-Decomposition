
/*
 * Include Files
 *
 */
#if defined(MATLAB_MEX_FILE)
#include "tmwtypes.h"
#include "simstruc_types.h"
#else
#include "rtwtypes.h"
#endif



/* %%%-SFUNWIZ_wrapper_includes_Changes_BEGIN --- EDIT HERE TO _END */

#include <math.h>
/* %%%-SFUNWIZ_wrapper_includes_Changes_END --- EDIT HERE TO _BEGIN */
#define u_width 1
#define y_width 1

/*
 * Create external references here.  
 *
 */
/* %%%-SFUNWIZ_wrapper_externs_Changes_BEGIN --- EDIT HERE TO _END */

/* extern double func(double a); */
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Output function
 *
 */
void parkXFRM_3_Outputs_wrapper(const real_T *fa,
			const real_T *fb,
			const real_T *fc,
			const real_T *theta,
			real_T *fd,
			real_T *fq)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */

/* 
 *Theta between d axes and alpha
 *fs = vq-jvd

*/
#include <math.h>

double Alpha, Beta;
//Clark
// Alpha = fa[0];												
// Beta = (fa[0] +2*fb[0])/sqrt(3.0);
Alpha = .666667*fa[0]-0.333333*fb[0]-0.333333*fc[0];												
Beta = 0.57735*(fb[0]-fc[0]);

//Park
fd[0] = Alpha*cos(theta[0]) + Beta*sin(theta[0]);	
fq[0] = Beta*cos(theta[0]) - Alpha*sin(theta[0]);
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}


