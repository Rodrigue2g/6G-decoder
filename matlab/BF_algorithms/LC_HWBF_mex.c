#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

// Single bit gradient descent decoding function with small tweaks.
/* 
    Decoding based on the paper 6 "Gradient Descent Bit Flipping Algorithms for Decoding LDPC Codes"
    It follows the basic single bit Gradient descent algorithm but must have some changes to account for SC-LDPC codes (to be implemented).
    Uses (H2 H1) matrix as the LDPC matrix to decode the codewords as they come in.

    The small tweak is that we keep the best codeword found so far and we return it if the last codeword is not better than the best one.
*/
#define ABS(x) ((x) > 0 ? (x) : -(x))

#define MAX_ITER 100
#define EPSILON 0.001
#define DEBUG 1
#define THETA_ATTN 0.9

void LC_HWBF_full(double *H_full, double *input, double *output, double *counter_flip, int L, int w, int Z, int nb, int mb);
void LC_HWBF(int **H, int *x, double *y, int M, int N, double *counter_flip);
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
    plhs[1] = mxCreateDoubleMatrix(1, MAX_ITER, mxREAL);
    double *num_flips = mxGetPr(plhs[1]);

    LC_HWBF_full(param_Hfull, param_input, param_output, num_flips, param_L, param_w, param_Z, param_nb, param_mb);
    return;
}

void LC_HWBF_full(double *H_full, double *input, double *output, double *counter_flip, int L, int w, int Z, int nb, int mb)
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

    //*counter_flip = 0;

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
        LC_HWBF(H, xhat, y, M, N, counter_flip);
        
        // output decoded bits
        for (int i = 0; i < code_length; i++)
        {
            //printf("xhat[i] %d\n", xhat[i]);
            out[i] = xhat[i]; //(xhat[i] < 0);
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

void LC_HWBF(int **H, int *x, double *y, int M, int N, double *counter_flip) {
    int k = 0;
    const int max_iterations = 100;
    double *e = (double *)mxCalloc(N, sizeof(double)); // Extrinsic information
    int *Z = (int *)mxCalloc(N, sizeof(int)); // Hard decisions on bits
    double **w_mn = (double **)mxCalloc(M, sizeof(double *));
    for (int i = 0; i < M; i++)
    {
        w_mn[i] = (double *)mxCalloc(N, sizeof(double));
    }
    int *S = (int *)mxCalloc(M, sizeof(int)); // Syndrome bits

    // Step 1: Initialization
    for (int n = 0; n < N; n++) {
        Z[n] = (y[n] > 0) ? 1 : 0;
    }

    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            w_mn[m][n] = DBL_MAX;
            if (H[m][n] == 1) {
                for(int n2=0; n2<N; n2++){
                    if (H[m][n2] == 1 && ABS(y[n2]) < w_mn[m][n] && n2!=n) {
                        w_mn[m][n] = ABS(y[n2]);   
                    }
                }
            }
        }
    }

    while (k < max_iterations) {
        
        for (int n = 0; n < N; n++) {
            e[n] = 0.0;
        }
        
        // Step 3: Check node update process
        // Compute syndrome for each check node
        for (int m = 0; m < M; m++) {
            S[m] = 0;
            for (int n = 0; n < N; n++) {
                S[m] += Z[n] * H[m][n];
            }
            S[m] = S[m] % 2;
        }
        int all_zero = 1;
        for(int m = 0; m < M; m++) {
            if (S[m] != 0) {
                all_zero = 0;
                break;
            }
        }
        if(all_zero) {
            break;
        }

        // Step 4: Variable node update process
        for (int n = 0; n < N; n++) {
            // printf("Z[n]: %d\n", Z[n]);
            // printf("y[n]: %f , &y[n]: %f , ABS(y[n]): %f\n", y[n], &y[n], ABS(y[n]));
            for (int m = 0; m < M; m++) {
                //printf("S[m]: %d , w_mn[m]: %d\n", S[m], w_mn[m]);
                //e[n] += 1/ABS(y[n]) * (2 * S[m] - 1) * w_mn[m] * THETA_ATTN;                
                if(H[m][n]) {
                    e[n] += 1/ABS(y[n]) * (2 * S[m] - 1) * w_mn[m][n] * THETA_ATTN;
                }
            }
            // printf("e[n]: %f\n", e[n]);
        }

        // Step 5: Update hard decisions based on extrinsic information
        for (int n = 0; n < N; n++) {
            //printf("e[n]: %f\n", e[n]);
            //Z[n] = (e[n] <= 0) ? 0 : 1;
        }

        // Step 6: Flip the bit with the highest uncertainty
        int max_e_index = 0;
        double max_e_value = ABS(e[0]);
        for (int n = 1; n < N; n++) {
            if (ABS(e[n]) > max_e_value) {
                max_e_value = ABS(e[n]);
                max_e_index = n;
            }
        }
        Z[max_e_index] = 1 - Z[max_e_index];
        //(*counter_flip)++;

        // Check if all syndromes are zero
        all_zero = 1;
        for (int m = 0; m < M; m++) {
            if (S[m] != 0) {
                all_zero = 0;
                break;
            }
        }
        if (all_zero) {
            break; // Stop if decoding is successful
        }

        counter_flip[k] = 1;

        k++;
    }

    // Step 7: Output the final decoded codeword
    for (int n = 0; n < N; n++) {
        x[n] = Z[n];
    }

    mxFree(e);
    mxFree(Z);
    for (int i = 0; i < M; i++)
    {
        mxFree(w_mn[i]);
    }
    mxFree(w_mn);
    mxFree(S);
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