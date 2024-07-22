#include "mex.h"
#include "math.h"
#include "float.h"
// REVISIT: This window decoding refers to TIT 2012 'Window decoding of LDPC-CC on erasure channels'. However, since we still instantiate memories for the decoded ms-1 region (the memory of the code), it is slightly resource-wasting

// prototypes
// BF Algorithms
void   ldpc_weighted_BF_decoding(double *H, int H_rows, int H_cols, int itera, double *LLRIn, int LLRIn_rows, int frames, double *decode_bits, double *decode_itera);
void   ldpc_multiple_gradient_descent_BF_decoding(double **H, int itera, double *LLRIn, int frames, int delta, double *debits, double *deiter);

void flip_bits_based_on_weights(double *H, int H_rows, int H_cols, double *decode_bits, double *weights);
void ldpc_weighted_CBF_decoding(double *H, int H_rows, int H_cols, int itera, double *LLRIn, int LLRIn_rows, int frames, double *decode_bits, double *decode_itera);
// Helpers
void   Rotation_mex (float *xIn, int shiftNum, int Z, float *yOut);
int    signOP_mex (float x);
int    HD_mex (float x);
float  sum_mex (float *x, int size);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *param_LLRIn;
    double *param_H;
    double *param_debits;
    double *param_deiter;
    int param_itera;
    int param_frames;
    int param_H_rows;
    int param_H_cols;

    // Input arguments
    param_H = mxGetPr(prhs[0]);
    param_H_rows = mxGetM(prhs[0]);
    param_H_cols = mxGetN(prhs[0]);
    param_itera = (int)mxGetScalar(prhs[1]);
    param_LLRIn = mxGetPr(prhs[2]);
    int param_LLRIn_rows = mxGetM(prhs[2]);
    param_frames = (int)mxGetScalar(prhs[3]);

    // Output matrices
    plhs[0] = mxCreateDoubleMatrix(param_H_cols, param_frames, mxREAL);
    param_debits = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1, param_frames, mxREAL);
    param_deiter = mxGetPr(plhs[1]);

    ldpc_weighted_BF_decoding(param_H, param_H_rows, param_H_cols, param_itera, param_LLRIn, param_LLRIn_rows, param_frames, param_debits, param_deiter);
    return;
}

/**
 * LDPC Wighted BF
 * 
 * check parity matrix
 * ~ check if err rate > treshold // later
 * if < --> return decoded bits directly to the ouput
 * elif >= --> BF
 * 
 * @returns{ debits, deiter } 
 * debits = decode_bits; - return the LLRs to the aposteriori memory of the windowed decoder + decoded bits to the ouput
 * deiter = decode_itera;
*/
void ldpc_weighted_BF_decoding(double *H, int H_rows, int H_cols, int itera, double *LLRIn, int LLRIn_rows, int frames, double *decode_bits, double *decode_itera) { //Use Bgs instead of Hs
    double *weight =    (double *)mxMalloc(H_rows * sizeof(double));
    double *inverFunc = (double *)mxMalloc(H_cols * sizeof(double));
    double *xhat =      (double *)mxMalloc(H_cols * sizeof(double));
    double *syndrome =  (double *)mxMalloc(H_rows * sizeof(double));
    
    for (int i_frames = 0; i_frames < frames; i_frames++) {
        // Syndrome Calculation
        for (int i_rows = 0; i_rows < H_rows; i_rows++) {
            weight[i_rows] = 0;
            for (int i_cols = 0; i_cols < H_cols; i_cols++) {
                weight[i_rows] += H[i_cols * H_rows + i_rows] * fabs(LLRIn[i_cols * LLRIn_rows + i_frames]);
            }
        }
        
        for (int i_itera = 0; i_itera < itera; i_itera++) {
            for (int i_cols = 0; i_cols < H_cols; i_cols++) {
                inverFunc[i_cols] = 0;
                for (int i_rows = 0; i_rows < H_rows; i_rows++) {
                    inverFunc[i_cols] += H[i_cols * H_rows + i_rows] * weight[i_rows] * (2 * syndrome[i_rows] - 1);
                }
            }
            
            int flipIdx = 0;
            for (int i_cols = 0; i_cols < H_cols; i_cols++) {
                if (inverFunc[i_cols] > inverFunc[flipIdx])
                    flipIdx = i_cols;
            }
            
            if (flipIdx < H_cols) {
                xhat[flipIdx] = fmod(xhat[flipIdx] + 1, 2);
                for (int i_rows = 0; i_rows < H_rows; i_rows++) {
                    syndrome[i_rows] = fmod(syndrome[i_rows] + H[i_rows * H_cols + flipIdx], 2);
                }
            }
            
            int syndrome_sum = 0;
            for (int i_rows = 0; i_rows < H_rows; i_rows++) {
                syndrome_sum += syndrome[i_rows];
            }
            if (syndrome_sum == 0) {
                break;
            }
        }
        
        for (int i_cols = 0; i_cols < H_cols; i_cols++) {
            decode_bits[i_cols * LLRIn_rows + i_frames] = xhat[i_cols];
        }
        decode_itera[i_frames] = i_itera;
    }

    mxFree(weight);
    mxFree(inverFunc);
    mxFree(xhat);
    mxFree(syndrome);
}


// Example function to flip bits based on the weighted criteria.
void flip_bits_based_on_weights(double *H, int H_rows, int H_cols, double *decode_bits, double *weights) {
    for (int j = 0; j < H_cols; j++) {
        double weight_sum = 0.0;
        for (int i = 0; i < H_rows; i++) {
            if (H[i * H_cols + j] != 0) {
                weight_sum += weights[i];
            }
        }
        // flip if the weighted sum is greater than a threshold
        if (weight_sum > H_rows / 2) { 
            decode_bits[j] =  ~decode_bits[j]; //1 - decode_bits[j];
        }
    }
}

void ldpc_weighted_CBF_decoding(double *H, int H_rows, int H_cols, int itera, double *LLRIn, int LLRIn_rows, int frames, double *decode_bits, double *decode_itera) {
    double *weights = (double *)mxMalloc(H_rows * sizeof(double));
    for (int frame = 0; frame < frames; frame++) {
        for (int j = 0; j < H_cols; j++) {
            decode_bits[j] = LLRIn[j] > 0 ? 0 : 1;
        }
        for (int iter = 0; iter < itera; iter++) {
            // Reset weights each iteration (simplified approach)
            for (int i = 0; i < H_rows; i++) {
                weights[i] = 0;  // Reset weights
            }
            // Compute syndrome and adjust weights
            for (int i = 0; i < H_rows; i++) {
                double syndrome = 0;
                for (int j = 0; j < H_cols; j++) {
                    if (H[i * H_cols + j] != 0) {
                        syndrome += decode_bits[j];
                        weights[i] += 1;  // Increase weight if the bit is involved in the check
                    }
                }
                syndrome = fmod(syndrome, 2.0);  // Compute syndrome for each row
                if (syndrome != 0) { 
                    // Increase the weights of bits involved in unsatisfied checks
                }
            }
            // Flip bits based on weights
            flip_bits_based_on_weights(H, H_rows, H_cols, decode_bits, weights);
            // check for decoding success and break out early (not implemented)
        }
        // store the number of iterations actually performed (not implemented)
    }
    mxFree(weights);
}

void Rotation_mex(float *xIn, int shiftNum, int Z, float *yOut)
{

    int i;
    for (i = 0; i < Z; i++)
    {
        if (i < (Z - shiftNum))
        {
            yOut[i] = xIn[i + shiftNum];
        }
        else 
        {
            yOut[i] = xIn[i - (Z - shiftNum)];
        }
    }
}

int signOP_mex (float x)
{
    return ((x >= 0) - (x < 0));
}

int HD_mex (float x)
{
    return (int)(x < 0);
}

float sum_mex (float *x, int size)
{   
    float sum = 0;
    int i;
    for (i = 0; i < size; i++)
    {
        sum = sum + x[i];
    }
    return sum;
}


void ldpc_multiple_gradient_descent_BF_decoding(double **H, int itera, double *LLRIn, int frames, int delta, double *debits, double *deiter) // function [decode_bits, decode_itera] =
{
    // debits = decode_bits;
    // deiter = decode_itera;
}
