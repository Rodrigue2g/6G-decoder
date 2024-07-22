#include "mex.h"
#include "math.h"
#include "float.h"
// Multi bit gradient descent decoding function.
/*
    Decoding based on the paper 6 "Gradient Descent Bit Flipping Algorithms for Decoding LDPC Codes"
    It follows the basic multi bit Gradient descent algorithm but must have some changes to account for SC-LDPC codes.
    Uses (H2 H1) matrix as the LDPC matrix to decode the codewords as they come in.
    A codeword can start it's decoding process only once the next one is received, the decoding latency is w (2 in our case).

    Small tweaks were made
*/


#define ABS(x) ((x) > 0 ? (x) : -(x))

#define MAX_ITER 100
#define THETA -0.5
#define ALPHA 1
#define BETA 1

#define EPSILON 0.001
#define DEBUG 0

void MGDBF_full(double *H_full, double *input, double *output, int L, int w, int Z, int nb, int mb);
void MGDBF(int **H, int *x, double *y, int M, int N);
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

    MGDBF_full(param_Hfull, param_input, param_output, param_L, param_w, param_Z, param_nb, param_mb);
    return;
}

void MGDBF_full(double *H_full, double *input, double *output, int L, int w, int Z, int nb, int mb)
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
        MGDBF(H, xhat, y, M, N);
        
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

void MGDBF(int **H, int *x, double *y, int M, int N)
{   

    int *flip_pos = mxMalloc(N * sizeof(int));
    int mode = 0; // 0 multi bit, 1 single bit

    int counter = 0;

    double f_first = ObjectiveFunction(H, x, y, M, N);

    if(DEBUG) printf("MGDBF ");
    for (int i_itera = 0; i_itera < MAX_ITER; i_itera++)
    {

        // break if parity check is successful
        if(ParityCheck(H, x, M, N)==1)
        {   
            if(DEBUG) printf("Breaking (Par Check) at iter %d and ", i_itera);
            break;
        }

        int flip_pos_length = 0;
        int lowestIdx = 0;
        double lowestValue = 0;
        double inverFunc;
        for (int k = 0; k < N; k++)
        {
            inverFunc = ALPHA * x[k] * y[k]; //alpha
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
                    inverFunc += BETA * accum; //beta
                }
            }

            // find the smallest value in inverFunc
            if (inverFunc <lowestValue || k == 0)
            {
                lowestIdx = k;
                lowestValue = inverFunc;
            }

            // find the values that are smaller than theta
            if (inverFunc < THETA)
            {
                flip_pos[flip_pos_length] = k;
                flip_pos_length++;
            }
        }

        if (mode == 0)
        {   
            counter+=flip_pos_length;
            double f1 = ObjectiveFunction(H, x, y, M, N);
            for(int j=0; j<flip_pos_length; j++)
            {
                x[flip_pos[j]] = -x[flip_pos[j]];
            }
            double f2 = ObjectiveFunction(H, x, y, M, N);

            if (f2 < f1)
            {       
                for(int j=0; j<flip_pos_length; j++)
                {
                    x[flip_pos[j]] = -x[flip_pos[j]];
                }
                mode = 1;
                //break;
            }
        }
        else
        {   
            counter += 1;
            double f1 = ObjectiveFunction(H, x, y, M, N);
            // flip single bit
            x[lowestIdx] = -x[lowestIdx];
            double f2 = ObjectiveFunction(H, x, y, M, N);

            if (f2 < f1)
            {
                x[lowestIdx] = -x[lowestIdx];
                break;
            }
        }
    }

    double f_last = ObjectiveFunction(H, x, y, M, N);

    if(DEBUG) printf("f first %f and f last %f, flipped: %d times\n", f_first, f_last, counter);


    mxFree(flip_pos);
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
        obj += ALPHA * x[i] * y[i]; //alpha
    }
    for (int i = 0; i < M; i++)
    {
        int accum = 1;
        for (int j = 0; j < N; j++)
        {
            if (H[i][j] == 1)
                accum *= x[j];
        }
        obj += BETA * accum; //beta
    }
    return obj;
}