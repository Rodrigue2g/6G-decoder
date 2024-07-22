#include "mex.h"
#include "math.h"
#include "float.h"
// Single bit gradient descent decoding function with small tweaks.
/* 
    Decoding based on the paper 6 "Gradient Descent Bit Flipping Algorithms for Decoding LDPC Codes"
    It follows the basic single bit Gradient descent algorithm but must have some changes to account for SC-LDPC codes (to be implemented).
    Uses (H2 H1) matrix as the LDPC matrix to decode the codewords as they come in.

    The small tweak is that we keep the best codeword found so far and we return it if the last codeword is not better than the best one.
*/
#define ABS(x) ((x) > 0 ? (x) : -(x))

#define MAX_ITER 100
#define EPSILON 10
#define DEBUG 1


void SGDBF_full(double *H_full, double *input, double *output, int L, int w, int Z, int nb, int mb);
void SGDBF(int **H, int *x, double *y, int M, int N);
int ParityCheck(int **H, int *x, int M, int N);
double ObjectiveFunction(int **H, int *x, double *y, int M, int N);

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

    SGDBF_full(param_Hfull, param_input, param_output, param_L, param_w, param_Z, param_nb, param_mb);
    return;
}

void SGDBF_full(double *H_full, double *input, double *output, int L, int w, int Z, int nb, int mb)
{
    // SC-LDPC full matrix dimensions
    int tail_length = 128; // ! TODO: Check this
    int M_tot = mb * Z * L + tail_length;
    int N_tot = nb * Z * L + tail_length;

    // LDPC H matrix dimensions
    int M = mb * Z;
    int N = w * nb * Z;

    // length of one codeword
    int code_length = nb * Z;

    // Allocate memory for H
    int **H = (int **)mxMalloc(M * sizeof(int *));
    for (int i = 0; i < M; i++)
    {
        H[i] = (int *)mxMalloc(N * sizeof(int));
    }

    // Copy H from MATLAB to C to create the LDPC matrix
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            H[i][j] = (int)H_full[(M + i) * N_tot + j];
        }
    }

    // estimated bits
    int *xhat_tot = mxCalloc(N_tot + code_length, sizeof(int));

    for (int i = 0; i < N_tot; i++)
    {
        xhat_tot[i] = (input[i] >= 0 ? 1 : -1);
    }

    // iterate over all the codewords (t denotes the spatial position of the codeword = time index)
    for (int t = 0; t < L - 1; t++)
    {
        // offset the input and output pointers to the current codeword
        double *y = &input[t * code_length];
        double *out = &output[t * code_length];
        int    *xhat = &xhat_tot[t * code_length];

        // decode the current codeword
        SGDBF(H, xhat, y, M, N);
        
        // output decoded bits
        for (int i = 0; i < code_length; i++)
        {
            out[i] = (xhat[i] < 0);
        }
    }

    mxFree(xhat_tot);

    // Free Memory
    for (int i = 0; i < M; i++)
    {
        mxFree(H[i]);
    }

    mxFree(H);
}

void SGDBF(int **H, int *x, double *y, int M, int N)
{   

    double f_first = ObjectiveFunction(H, x, y, M, N);
    double f_last = f_first;

    if(DEBUG) printf("SGDBF V2 ");
    for (int i_itera = 0; i_itera < MAX_ITER; i_itera++)
    {

        // break if parity check is successful
        if(ParityCheck(H, x, M, N)==1)
        {   
            if(DEBUG) printf("Breaking (Par Check) at iter %d and ", i_itera);
            break;
        }

        int lowestIdx = 0;
        double lowestValue = 0;
        double inverFunc;
        for (int k = 0; k < N; k++)
        {
            inverFunc = x[k] * y[k];
            for (int mk = 0; mk < M; mk++)
            {
                if (H[mk][k] == 1)
                {
                    int accum = 1;
                    for (int nk = 0; nk < N; nk++)
                    {
                        if (H[mk][nk] == 1)
                            accum *= x[nk];
                    }
                    inverFunc += accum;
                }
            }

            // find the smallest value in inverFunc
            if (inverFunc < lowestValue || k == 0)
            {
                lowestIdx = k;
                lowestValue = inverFunc;
            }
        }

        // flip single bit
        x[lowestIdx] = -x[lowestIdx];

        double tmp = ObjectiveFunction(H, x, y, M, N);
        /*
        if(tmp<f_last)
        {   
            x[lowestIdx] = -x[lowestIdx];
            if(DEBUG) printf("Breaking (Conv) at iter %d and ", i_itera);
            break;
        }
        f_last = tmp;
        */
        if (tmp < f_last)
        {   
            x[lowestIdx] = -x[lowestIdx];
            if(DEBUG) printf("Breaking (Conv) at iter %d and ", i_itera);
            break;
        }
        if (ABS(f_last - tmp) < EPSILON)
        {
            if(DEBUG) printf("Breaking (Epsilon) at iter %d and ", i_itera);
            break;
        }
        f_last = tmp;
    }

    if(DEBUG) printf("f first %f and f last %f \n", f_first, f_last);

}

int ParityCheck(int **H, int *x, int M, int N)
{
    for (int i = 0; i < M; i++)
    {
        int accum = 1;
        for (int j = 0; j < N; j++)
        {
            if (H[i][j] == 1)
                accum *= x[j];
        }
        if (accum != 1)
        {
            return 0;
        }
    }
    return 1;
}


double ObjectiveFunction(int **H, int *x, double *y, int M, int N)
{
    double obj = 0;
    for (int i = 0; i < N; i++)
    {
        obj += x[i] * y[i];
    }
    for (int i = 0; i < M; i++)
    {
        int accum = 1;
        for (int j = 0; j < N; j++)
        {
            if (H[i][j] == 1)
                accum *= x[j];
        }
        obj += accum;
    }
    return obj;
}