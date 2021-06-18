
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
void speedb4_Outputs_wrapper(const real_T *oldTheta_in,
			const real_T *theta_in,
			const real_T *baseFreq_in,
			const real_T *T_in,
			const real_T *P_in,
			const real_T *speed_in,
			real_T *theta_pu,
			real_T *speed_pu)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */

//Calculates speed from theta using low pass filter
#include <math.h>;

#define pi 3.14159265358979323846
//For Saturation of Speed
#define DERIVE_MAX (1.0)	
#define SAT(x) ( ((x) > DERIVE_MAX) ? DERIVE_MAX : ( (-(x) > DERIVE_MAX) ? (-DERIVE_MAX) : (x) ) )
        
// double Tmp,theta,oldTheta,K1,K2,K3,speed_pu,baseRpm,SpeedRpm;
double Tmp,theta,oldTheta,K1,K2,K3,Speed;
double P,T,baseFreq;
//Normalize angle between zero and one
oldTheta = oldTheta_in[0];
// theta = theta_in[0]/(2*pi);
theta = theta_in[0];
baseFreq = baseFreq_in[0];
T = T_in[0];
P = P_in[0];
Speed = speed_in[0];

// Initialize the Speed module for speed calculation from QEP/RESOLVER
    K1 = 1/(baseFreq*T);
    K2 = 1/(1+T*2*pi*5);      // Low-pass cut-off frequency rad/s
    K3 = 1-K2;
//     baseRpm = 120*(baseFreq/P);
    
//Get difference between current theta and old theta
Tmp = theta - oldTheta;		                    
//Normalize between -.5 and.5
    if (Tmp < -0.5)			                                
     Tmp = Tmp + 1.0;                                      
   else if (Tmp > 0.5)			                           
     Tmp = Tmp - 1.0;                                      
   
    Tmp = K1*Tmp;		                                
/* Low-pass filter*/																			
   	Tmp = K2*Speed + K3*Tmp;		
/* Saturate the output */    
	//Tmp=_IQsat(v.Tmp,_IQ21(1),_IQ21(-1));	
    Tmp = SAT(Tmp);

Speed = Tmp;
    						
/* Change motor speed from pu value to rpm value */				
//     SpeedRpm = baseRpm*Speed;
    
/* Update the electrical angle */
    theta_pu[0] = theta;
    speed_pu[0] = Speed;
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}


