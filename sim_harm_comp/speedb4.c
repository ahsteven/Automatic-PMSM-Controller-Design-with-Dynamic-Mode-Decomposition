/*
 * File: speedb4.c
 *
 *
 *   --- THIS FILE GENERATED BY S-FUNCTION BUILDER: 3.0 ---
 *
 *   This file is an S-function produced by the S-Function
 *   Builder which only recognizes certain fields.  Changes made
 *   outside these fields will be lost the next time the block is
 *   used to load, edit, and resave this file. This file will be overwritten
 *   by the S-function Builder block. If you want to edit this file by hand, 
 *   you must change it only in the area defined as:  
 *
 *        %%%-SFUNWIZ_defines_Changes_BEGIN
 *        #define NAME 'replacement text' 
 *        %%% SFUNWIZ_defines_Changes_END
 *
 *   DO NOT change NAME--Change the 'replacement text' only.
 *
 *   For better compatibility with the Simulink Coder, the
 *   "wrapper" S-function technique is used.  This is discussed
 *   in the Simulink Coder's Manual in the Chapter titled,
 *   "Wrapper S-functions".
 *
 *  -------------------------------------------------------------------------
 * | See matlabroot/simulink/src/sfuntmpl_doc.c for a more detailed template |
 *  ------------------------------------------------------------------------- 
 *
 * Created: Thu Jun 17 09:57:22 2021
 */

#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME speedb4
/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
/* %%%-SFUNWIZ_defines_Changes_BEGIN --- EDIT HERE TO _END */
#define NUM_INPUTS            6
/* Input Port  0 */
#define IN_PORT_0_NAME        oldTheta_in
#define INPUT_0_WIDTH         1
#define INPUT_DIMS_0_COL      1
#define INPUT_0_DTYPE         real_T
#define INPUT_0_COMPLEX       COMPLEX_NO
#define IN_0_FRAME_BASED      FRAME_NO
#define IN_0_BUS_BASED        0
#define IN_0_BUS_NAME         
#define IN_0_DIMS             1-D
#define INPUT_0_FEEDTHROUGH   1
#define IN_0_ISSIGNED         0
#define IN_0_WORDLENGTH       8
#define IN_0_FIXPOINTSCALING  1
#define IN_0_FRACTIONLENGTH   9
#define IN_0_BIAS             0
#define IN_0_SLOPE            0.125
/* Input Port  1 */
#define IN_PORT_1_NAME        theta_in
#define INPUT_1_WIDTH         1
#define INPUT_DIMS_1_COL      1
#define INPUT_1_DTYPE         real_T
#define INPUT_1_COMPLEX       COMPLEX_NO
#define IN_1_FRAME_BASED      FRAME_NO
#define IN_1_BUS_BASED        0
#define IN_1_BUS_NAME         
#define IN_1_DIMS             1-D
#define INPUT_1_FEEDTHROUGH   1
#define IN_1_ISSIGNED         0
#define IN_1_WORDLENGTH       8
#define IN_1_FIXPOINTSCALING  1
#define IN_1_FRACTIONLENGTH   9
#define IN_1_BIAS             0
#define IN_1_SLOPE            0.125
/* Input Port  2 */
#define IN_PORT_2_NAME        baseFreq_in
#define INPUT_2_WIDTH         1
#define INPUT_DIMS_2_COL      1
#define INPUT_2_DTYPE         real_T
#define INPUT_2_COMPLEX       COMPLEX_NO
#define IN_2_FRAME_BASED      FRAME_NO
#define IN_2_BUS_BASED        0
#define IN_2_BUS_NAME         
#define IN_2_DIMS             1-D
#define INPUT_2_FEEDTHROUGH   1
#define IN_2_ISSIGNED         0
#define IN_2_WORDLENGTH       8
#define IN_2_FIXPOINTSCALING  1
#define IN_2_FRACTIONLENGTH   9
#define IN_2_BIAS             0
#define IN_2_SLOPE            0.125
/* Input Port  3 */
#define IN_PORT_3_NAME        T_in
#define INPUT_3_WIDTH         1
#define INPUT_DIMS_3_COL      1
#define INPUT_3_DTYPE         real_T
#define INPUT_3_COMPLEX       COMPLEX_NO
#define IN_3_FRAME_BASED      FRAME_NO
#define IN_3_BUS_BASED        0
#define IN_3_BUS_NAME         
#define IN_3_DIMS             1-D
#define INPUT_3_FEEDTHROUGH   1
#define IN_3_ISSIGNED         0
#define IN_3_WORDLENGTH       8
#define IN_3_FIXPOINTSCALING  1
#define IN_3_FRACTIONLENGTH   9
#define IN_3_BIAS             0
#define IN_3_SLOPE            0.125
/* Input Port  4 */
#define IN_PORT_4_NAME        P_in
#define INPUT_4_WIDTH         1
#define INPUT_DIMS_4_COL      1
#define INPUT_4_DTYPE         real_T
#define INPUT_4_COMPLEX       COMPLEX_NO
#define IN_4_FRAME_BASED      FRAME_NO
#define IN_4_BUS_BASED        0
#define IN_4_BUS_NAME         
#define IN_4_DIMS             1-D
#define INPUT_4_FEEDTHROUGH   1
#define IN_4_ISSIGNED         0
#define IN_4_WORDLENGTH       8
#define IN_4_FIXPOINTSCALING  1
#define IN_4_FRACTIONLENGTH   9
#define IN_4_BIAS             0
#define IN_4_SLOPE            0.125
/* Input Port  5 */
#define IN_PORT_5_NAME        speed_in
#define INPUT_5_WIDTH         1
#define INPUT_DIMS_5_COL      1
#define INPUT_5_DTYPE         real_T
#define INPUT_5_COMPLEX       COMPLEX_NO
#define IN_5_FRAME_BASED      FRAME_NO
#define IN_5_BUS_BASED        0
#define IN_5_BUS_NAME         
#define IN_5_DIMS             1-D
#define INPUT_5_FEEDTHROUGH   1
#define IN_5_ISSIGNED         0
#define IN_5_WORDLENGTH       8
#define IN_5_FIXPOINTSCALING  1
#define IN_5_FRACTIONLENGTH   9
#define IN_5_BIAS             0
#define IN_5_SLOPE            0.125

#define NUM_OUTPUTS           2
/* Output Port  0 */
#define OUT_PORT_0_NAME       theta_pu
#define OUTPUT_0_WIDTH        1
#define OUTPUT_DIMS_0_COL     1
#define OUTPUT_0_DTYPE        real_T
#define OUTPUT_0_COMPLEX      COMPLEX_NO
#define OUT_0_FRAME_BASED     FRAME_NO
#define OUT_0_BUS_BASED       0
#define OUT_0_BUS_NAME        
#define OUT_0_DIMS            1-D
#define OUT_0_ISSIGNED        1
#define OUT_0_WORDLENGTH      8
#define OUT_0_FIXPOINTSCALING 1
#define OUT_0_FRACTIONLENGTH  3
#define OUT_0_BIAS            0
#define OUT_0_SLOPE           0.125
/* Output Port  1 */
#define OUT_PORT_1_NAME       speed_pu
#define OUTPUT_1_WIDTH        1
#define OUTPUT_DIMS_1_COL     1
#define OUTPUT_1_DTYPE        real_T
#define OUTPUT_1_COMPLEX      COMPLEX_NO
#define OUT_1_FRAME_BASED     FRAME_NO
#define OUT_1_BUS_BASED       0
#define OUT_1_BUS_NAME        
#define OUT_1_DIMS            1-D
#define OUT_1_ISSIGNED        1
#define OUT_1_WORDLENGTH      8
#define OUT_1_FIXPOINTSCALING 1
#define OUT_1_FRACTIONLENGTH  3
#define OUT_1_BIAS            0
#define OUT_1_SLOPE           0.125

#define NPARAMS               0

#define SAMPLE_TIME_0         INHERITED_SAMPLE_TIME
#define NUM_DISC_STATES       0
#define DISC_STATES_IC        [0]
#define NUM_CONT_STATES       0
#define CONT_STATES_IC        [0]

#define SFUNWIZ_GENERATE_TLC  1
#define SOURCEFILES           "__SFB____SFB__"
#define PANELINDEX            0
#define USE_SIMSTRUCT         0
#define SHOW_COMPILE_STEPS    0
#define CREATE_DEBUG_MEXFILE  0
#define SAVE_CODE_ONLY        0
#define SFUNWIZ_REVISION      3.0
/* %%%-SFUNWIZ_defines_Changes_END --- EDIT HERE TO _BEGIN */
/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
#include "simstruc.h"


extern void speedb4_Outputs_wrapper(const real_T *oldTheta_in,
			const real_T *theta_in,
			const real_T *baseFreq_in,
			const real_T *T_in,
			const real_T *P_in,
			const real_T *speed_in,
			real_T *theta_pu,
			real_T *speed_pu);
/*====================*
 * S-function methods *
 *====================*/
/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *   Setup sizes of the various vectors.
 */
static void mdlInitializeSizes(SimStruct *S)
{

    DECL_AND_INIT_DIMSINFO(inputDimsInfo);
    DECL_AND_INIT_DIMSINFO(outputDimsInfo);
    ssSetNumSFcnParams(S, NPARAMS);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; /* Parameter mismatch will be reported by Simulink */
    }

    ssSetArrayLayoutForCodeGen(S, SS_COLUMN_MAJOR);

    ssSetOperatingPointCompliance(S, USE_DEFAULT_OPERATING_POINT);

    ssSetNumContStates(S, NUM_CONT_STATES);
    ssSetNumDiscStates(S, NUM_DISC_STATES);


    if (!ssSetNumInputPorts(S, NUM_INPUTS)) return;
    /* Input Port 0 */
    ssSetInputPortWidth(S, 0, INPUT_0_WIDTH);
    ssSetInputPortDataType(S, 0, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 0, INPUT_0_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 0, INPUT_0_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 0, 1); /*direct input signal access*/

    /* Input Port 1 */
    ssSetInputPortWidth(S, 1, INPUT_1_WIDTH);
    ssSetInputPortDataType(S, 1, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 1, INPUT_1_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 1, INPUT_1_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 1, 1); /*direct input signal access*/

    /* Input Port 2 */
    ssSetInputPortWidth(S, 2, INPUT_2_WIDTH);
    ssSetInputPortDataType(S, 2, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 2, INPUT_2_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 2, INPUT_2_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 2, 1); /*direct input signal access*/

    /* Input Port 3 */
    ssSetInputPortWidth(S, 3, INPUT_3_WIDTH);
    ssSetInputPortDataType(S, 3, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 3, INPUT_3_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 3, INPUT_3_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 3, 1); /*direct input signal access*/

    /* Input Port 4 */
    ssSetInputPortWidth(S, 4, INPUT_4_WIDTH);
    ssSetInputPortDataType(S, 4, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 4, INPUT_4_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 4, INPUT_4_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 4, 1); /*direct input signal access*/

    /* Input Port 5 */
    ssSetInputPortWidth(S, 5, INPUT_5_WIDTH);
    ssSetInputPortDataType(S, 5, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 5, INPUT_5_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 5, INPUT_5_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 5, 1); /*direct input signal access*/


    if (!ssSetNumOutputPorts(S, NUM_OUTPUTS)) return;
    /* Output Port 0 */
    ssSetOutputPortWidth(S, 0, OUTPUT_0_WIDTH);
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);
    ssSetOutputPortComplexSignal(S, 0, OUTPUT_0_COMPLEX);
    /* Output Port 1 */
    ssSetOutputPortWidth(S, 1, OUTPUT_1_WIDTH);
    ssSetOutputPortDataType(S, 1, SS_DOUBLE);
    ssSetOutputPortComplexSignal(S, 1, OUTPUT_1_COMPLEX);
    ssSetNumPWork(S, 0);

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    ssSetSimulinkVersionGeneratedIn(S, "10.1");

    /* Take care when specifying exception free code - see sfuntmpl_doc.c */
    ssSetOptions(S, (SS_OPTION_EXCEPTION_FREE_CODE |
                     SS_OPTION_USE_TLC_WITH_ACCELERATOR |
                     SS_OPTION_WORKS_WITH_CODE_REUSE));
}

/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    Specifiy  the sample time.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, SAMPLE_TIME_0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_SET_INPUT_PORT_DATA_TYPE
static void mdlSetInputPortDataType(SimStruct *S, int port, DTypeId dType)
{
    ssSetInputPortDataType(S, 0, dType);
}

#define MDL_SET_OUTPUT_PORT_DATA_TYPE
static void mdlSetOutputPortDataType(SimStruct *S, int port, DTypeId dType)
{
    ssSetOutputPortDataType(S, 0, dType);
}

#define MDL_SET_DEFAULT_PORT_DATA_TYPES
static void mdlSetDefaultPortDataTypes(SimStruct *S)
{
    ssSetInputPortDataType(S, 0, SS_DOUBLE);
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);
}

#define MDL_START  /* Change to #undef to remove function */
#if defined(MDL_START)
/* Function: mdlStart =======================================================
 * Abstract:
 *    This function is called once at start of model execution. If you
 *    have states that should be initialized once, this is the place
 *    to do it.
 */
static void mdlStart(SimStruct *S)
{
}
#endif /*  MDL_START */

/* Function: mdlOutputs =======================================================
 *
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
    const real_T *oldTheta_in = (real_T *) ssGetInputPortRealSignal(S, 0);
    const real_T *theta_in = (real_T *) ssGetInputPortRealSignal(S, 1);
    const real_T *baseFreq_in = (real_T *) ssGetInputPortRealSignal(S, 2);
    const real_T *T_in = (real_T *) ssGetInputPortRealSignal(S, 3);
    const real_T *P_in = (real_T *) ssGetInputPortRealSignal(S, 4);
    const real_T *speed_in = (real_T *) ssGetInputPortRealSignal(S, 5);
    real_T *theta_pu = (real_T *) ssGetOutputPortRealSignal(S, 0);
    real_T *speed_pu = (real_T *) ssGetOutputPortRealSignal(S, 1);

    speedb4_Outputs_wrapper(oldTheta_in, theta_in, baseFreq_in, T_in, P_in, speed_in, theta_pu, speed_pu);

}

/* Function: mdlTerminate =====================================================
 * Abstract:
 *    In this function, you should perform any actions that are necessary
 *    at the termination of a simulation.  For example, if memory was
 *    allocated in mdlStart, this is the place to free it.
 */
static void mdlTerminate(SimStruct *S)
{

}


#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif



