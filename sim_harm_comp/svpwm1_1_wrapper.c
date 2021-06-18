
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
void svpwm1_1_Outputs_wrapper(const real_T *i_alpha,
			const real_T *i_beta,
			real_T *Ta,
			real_T *Tb,
			real_T *Tc)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */

float vtmp1, vtmp2, vtmp3, vVecSector, vTa, vTb, vTc;
\
	vtmp1= i_beta[0];															\
	vtmp2= (i_beta[0]/2) + ((0.866)*i_alpha[0]);					\
    vtmp3= vtmp2 - vtmp1;													\
																				\
	vVecSector=3;																\
	vVecSector=(vtmp2> 0)?( vVecSector-1):vVecSector;						\
	vVecSector=(vtmp3> 0)?( vVecSector-1):vVecSector;						\
	vVecSector=(vtmp1< 0)?(7-vVecSector) :vVecSector;						\
																				\
	if     (vVecSector==1 || vVecSector==4)                                   \
      {     vTa= vtmp2; 														\
      		vTb= vtmp1-vtmp3; 												\
      		vTc=-vtmp2;														\
      }								    										\
   																				\
    else if(vVecSector==2 || vVecSector==5)                                   \
      {     vTa= vtmp3+vtmp2; 												\
      		vTb= vtmp1; 														\
      		vTc=-vtmp1;														\
      }																	   		\
   																				\
    else                                                                        \
      {     vTa= vtmp3; 														\
      		vTb=-vtmp3; 														\
      		vTc=-(vtmp1+vtmp2);												\
      }																	   		\
																				\
      Ta[0]=vTa;
      Tb[0]=vTb;
      Tc[0]=vTc;
              
              
 // __SVGEN_H__
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}


