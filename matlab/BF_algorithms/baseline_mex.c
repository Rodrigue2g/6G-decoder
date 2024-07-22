#include "mex.h"
#include "math.h"
#include "float.h"
// Baseline decoding function. 

/* 

This doesn't flip 

    1: takes in the input's y from the AWGN channel (y=c+z where c is the bipolar codeword and z is the noise)
    2: outputs y<0 as 1 and y>=0 as 0

It is used to compare with bitflipping algorithms to see if they actually improve the decoding process.

*/


void baseline_decoder(double *input, double *output, int L, int w, int Z, int nb, int mb);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *param_Hfull;
    double *param_input;
    int param_L;
    int param_w;
    int param_Z;
    int param_nb;
    int param_mb;
    
    param_Hfull = mxGetPr(prhs[0]);
    param_input = mxGetPr(prhs[1]);
    param_L = mxGetScalar(prhs[2]);
    param_w = mxGetScalar(prhs[3]);
    param_Z = mxGetScalar(prhs[4]);
    param_nb = mxGetScalar(prhs[5]);
    param_mb = mxGetScalar(prhs[6]);
    
    double *param_output;
    plhs[0] = mxCreateDoubleMatrix(param_L * param_nb * param_Z, 1, mxREAL);
    param_output = mxGetPr(plhs[0]);

    baseline_decoder(param_input, param_output, param_L, param_w, param_Z, param_nb, param_mb);
    return;
}

void baseline_decoder(double *input, double *output, int L, int w, int Z, int nb, int mb)
{
    for(int i=0; i<L*nb*Z; i++)
    {
        output[i] = (input[i] < 0 ? 1 : 0);
    }
}