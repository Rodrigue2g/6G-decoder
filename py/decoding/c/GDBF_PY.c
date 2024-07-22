//#include "mex.h"
#include <Python.h>
#include "math.h"
#include "float.h"
// Single bit gradient descent decoding function.
/* 
    Decoding based on the paper 6 "Gradient Descent Bit Flipping Algorithms for Decoding LDPC Codes"
    It follows the basic single bit Gradient descent algorithm but must have some changes to account for SC-LDPC codes.
    Uses (H2 H1) matrix as the LDPC matrix to decode the codewords as they come in.
    A codeword can start it's decoding process only once the next one is received, the decoding latency is w (2 in our case).
*/


void GDBF(double *H_full, double *input, double *output, int L, int w, int Z, int nb, int mb);

/*void GDBF_PY(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
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

}*/
static PyObject * GDBF_PY(double *H_full, double *input, double *output, int L, int w, int Z, int nb, int mb)
{
    double *param_Hfull;
    double *param_input;
    int param_L;
    int param_w;
    int param_Z;
    int param_nb;
    int param_mb;
    
    param_Hfull = H_full;
    param_input = input;
    param_L = L;
    param_w = w;
    param_Z = Z;
    param_nb = nb;
    param_mb = mb;
    
    double *param_output;
    param_output = malloc(param_L * param_nb * param_Z, 1);

    GDBF(param_Hfull, param_input, param_output, param_L, param_w, param_Z, param_nb, param_mb);
    return Py_BuildValue("i", param_output);

}


void GDBF(double *H_full, double *input, double *output, int L, int w, int Z, int nb, int mb)
{   

    int itera=45;

    // SC-LDPC full matrix dimensions
    int tail_length = 128; // ! TODO: Check this
    int M_tot = mb * Z * L + tail_length;
    int N_tot = nb * Z * L + tail_length;

    // LDPC H matrix dimensions
    int M = mb * Z;
    int N = w* nb * Z; 

    int code_length = nb*Z;

    // Allocate memory for H
    int** H = (int**)malloc(M * sizeof(int*)); 
    for (int i = 0; i < M; i++) {
        H[i] = (int*)malloc(N * sizeof(int));
    }

    // Copy H from MATLAB to C
    for(int i=0; i<M; i++)
    {
        for(int j=0; j<N; j++)
        {
            H[i][j] = (int)H_full[(M+i)*N_tot+j];
        }
    }
    
    // estimated bits
    int *xhat_tot=calloc(N_tot+code_length,sizeof(int));  

    for(int i=0; i<N_tot; i++)
    {
        xhat_tot[i] = (input[i]>=0 ? 1 : -1);
    }



    // iterate over all the codewords (t denotes the spatial position of the codeword = time index)
    for(int t=0; t<L-1; t++)
    {   
        // the new y (input) and out (output) arrays are shifted with the spatial position
        double* y       = &input[t*code_length];
        double* out     = &output[t*code_length];
        int *xhat       = &xhat_tot[t*code_length];

        for (int i_itera = 0; i_itera < itera; i_itera++){
            // Check parity check matrix //* OK
            bool holds = true;
            for (int i = 0; i < M; i++)
            {
                int accum = 1;
                for (int j = 0; j < N; j++)
                {   
                    if (H[i][j] == 1)
                        accum *= xhat[j];
                }
                if (accum != 1)
                {
                    holds = false;
                    break;
                }
            }

            if (holds)
            {
                printf("Converged after %d iterations\n", i_itera);
                break;
            }
            else if (i_itera == itera - 1)
            {
                printf("Did not converge after %d iterations\n", i_itera);
            }
            
            int lowestIdx = 0;
            int lowestValue = 0;
            double inverFunc;
            for (int k = 0; k < N; k++)
            {
                inverFunc = xhat[k] * y[k];
                for (int mk = 0; mk < M; mk++)
                {
                    if (H[mk][k] == 1)
                    {
                        int accum = 1;
                        for (int nk = 0; nk < N; nk++)
                        {   
                            if(H[mk][nk]==1)
                                accum *= xhat[nk];
                        }
                        inverFunc += accum;
                    }
                }

                // find the smallest value in inverFunc
                if(inverFunc < lowestValue || k==0)
                {
                    lowestIdx = k;
                    lowestValue=inverFunc;
                }
            }

            // flip single bit
            xhat[lowestIdx] = -xhat[lowestIdx];
        }

        // output decoded bits
        for (int i = 0; i < code_length; i++)
        {   
            out[i] = (xhat[i] < 0);
        }
    }


    free(xhat_tot);

    // Free Memory
    for (int i = 0; i < M; i++) {
        free(H[i]);
    }

    free(H);
}