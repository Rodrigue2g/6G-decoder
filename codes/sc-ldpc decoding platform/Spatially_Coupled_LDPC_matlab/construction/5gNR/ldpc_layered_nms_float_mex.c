#include "mex.h"
#include "math.h"
#include "float.h"

// prototypes
void   ldpc_layered_nms_float_mex(double *BG, double norm, int Z, int itera, int H_rows, int H_cols, int Frames, double *LLRIn, double *layer_order, double *debits, double *deitera);
void   Rotation_mex (float *xIn, int shiftNum, int Z, float *yOut);
int    signOP_mex (float x);
int    HD_mex (float x);
float  sum_mex (float *x, int size);
float  max_mex (float x, float y);

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *param_BG;
    double  param_norm;
    int     param_Z;
    int     param_itera;
    int     param_H_rows;
    int     param_H_cols;
    int     param_Frames;
    double *param_LLRIn;
    double *param_layer_order;
    double *param_debits;
    double *param_deitera;

    param_BG          = mxGetPr(prhs[0]);
    param_norm        = mxGetScalar(prhs[1]);
    param_Z           = (int)mxGetScalar(prhs[2]);
    param_itera       = (int)mxGetScalar(prhs[3]);
    param_H_rows      = (int)mxGetScalar(prhs[4]);
    param_H_cols      = (int)mxGetScalar(prhs[5]);
    param_Frames      = (int)mxGetScalar(prhs[6]);
    param_LLRIn       = mxGetPr(prhs[7]);
    param_layer_order = mxGetPr(prhs[8]);

    // Note that any output must have an array expression
    // It can be optimized later, but we leave it now
    plhs[0]       = mxCreateDoubleMatrix(1, param_Z*param_H_cols*param_Frames, mxREAL);
    param_debits  = mxGetPr(plhs[0]);
    plhs[1]       = mxCreateDoubleMatrix(1, 1, mxREAL);
    param_deitera = mxGetPr(plhs[1]);

    ldpc_layered_nms_float_mex(param_BG, param_norm, param_Z, param_itera, param_H_rows, param_H_cols, param_Frames, param_LLRIn, param_layer_order, param_debits, param_deitera);
    return;
}

void ldpc_layered_nms_float_mex(double *BG, double norm, int Z, int itera, int H_rows, int H_cols, int Frames, double *LLRIn, double *layer_order, double *debits, double *deitera)
{
    // parameter and memory allocations
    int i_frame, i_itera, i_rows, i_cols, i_z;

    int TDFlag;
    int shiftNum;

    // 1D or 2D matrix
    float **QMem;
    float **RMem;
    float **TMem;
    int   **QSign;

    float *min1;
    float *min2;
    int   *sign;
    int   *minIdx;
    float *minVec;
    int   *SynCheck;
    float *Qc;
    float *Qtmp;
    float *TMemAbs;

    QMem     = (float**)mxMalloc(sizeof(float*)*H_cols);
    RMem     = (float**)mxMalloc(sizeof(float*)*H_rows*H_cols);
    TMem     = (float**)mxMalloc(sizeof(float*)*H_cols);
    QSign    = (int**)mxMalloc(sizeof(int*)*H_cols);

    min1     = (float*)mxMalloc(sizeof(float)*Z);
    min2     = (float*)mxMalloc(sizeof(float)*Z);
    sign     = (int*)mxMalloc(sizeof(int)*Z);
    minIdx   = (int*)mxMalloc(sizeof(int)*Z);
    minVec   = (float*)mxMalloc(sizeof(float)*Z);
    SynCheck = (int*)mxMalloc(sizeof(int)*Z);
    Qc       = (float*)mxMalloc(sizeof(float)*Z);
    Qtmp     = (float*)mxMalloc(sizeof(float)*Z);
    TMemAbs  = (float*)mxMalloc(sizeof(float)*Z);
    
    // open 2D matrix
    for (i_rows = 0; i_rows < H_cols; i_rows++)
    {
        QMem[i_rows]  = (float*)mxMalloc(sizeof(float)*Z);
        TMem[i_rows]  = (float*)mxMalloc(sizeof(float)*Z);
        QSign[i_rows] = (int*)mxMalloc(sizeof(int)*Z);
    }

    for (i_rows = 0; i_rows < H_rows*H_cols; i_rows++)
    {
        RMem[i_rows]  = (float*)mxMalloc(sizeof(float)*Z);
    }
    
    for (i_frame = 0; i_frame < Frames; i_frame++)
    {
        // initilization
        for (i_rows = 0; i_rows < H_cols; i_rows++)
        {
            for (i_cols = 0; i_cols < Z; i_cols++)
            {
                QMem[i_rows][i_cols]  = (float)LLRIn[i_frame*H_cols*Z + i_rows*Z + i_cols];
                TMem[i_rows][i_cols]  = 0;
                QSign[i_rows][i_cols] = 0;
            }
        }
        for (i_rows = 0; i_rows < H_rows*H_cols; i_rows++)
        {
            for (i_cols = 0; i_cols < Z; i_cols++)
            {
                RMem[i_rows][i_cols]  = 0;
            }
        }

        // start decoding
        for (i_itera = 0; i_itera < itera; i_itera++)
        {
            TDFlag = 1;

            for (i_rows = 0; i_rows < H_rows; i_rows++)
            {
                // initilization
                for (i_z = 0; i_z < Z; i_z++)
                {
                    min1[i_z] = FLT_MAX; // head file "float.h"
                    min2[i_z] = FLT_MAX;
                    sign[i_z] = 1;

                    minIdx[i_z] = 0;
                    minVec[i_z] = 0;

                    SynCheck[i_z] = 0;
                }

                // PHASE 1: MIN
                for (i_cols = 0; i_cols < H_cols; i_cols++)
                {
                    // checkout whether there is -1 in base matrix, if yes, directly skip
                    if (BG[(int)(layer_order[i_rows]-1)*H_cols + i_cols] == -1)
                    {
                        continue;
                    }
                    shiftNum = (int)BG[(int)(layer_order[i_rows]-1)*H_cols + i_cols];

                    // calculate the Qc, TMem, and its sign and magnitude
                    Rotation_mex(QMem[i_cols], shiftNum, Z, Qc);

                    for (i_z = 0; i_z < Z; i_z++)
                    {
                        TMem[i_cols][i_z]  = Qc[i_z] - RMem[(int)(layer_order[i_rows]-1)*H_cols + i_cols][i_z];
                        sign[i_z]          = sign[i_z]*signOP_mex(TMem[i_cols][i_z]);
                        TMemAbs[i_z]       = fabs(TMem[i_cols][i_z]); // abs() only for int
                        QSign[i_cols][i_z] = signOP_mex(Qc[i_z]);
                    }

                    // search the min1, min2, and the column index corresponding to min1
                    for (i_z = 0; i_z < Z; i_z++)
                    {
                        if (TMemAbs[i_z] < min1[i_z])
                        {
                            min2[i_z]   = min1[i_z];
                            min1[i_z]   = TMemAbs[i_z];
                            minIdx[i_z] = i_cols;
                        }
                        else if (TMemAbs[i_z] < min2[i_z])
                        {
                            min2[i_z]   = TMemAbs[i_z];
                        }
                    }
                }

                // PHASE 2: Q- and R- messages update
                for (i_cols = 0; i_cols < H_cols; i_cols++)
                {
                    // checkout whether there is -1 in base matrix, if yes, directly skip
                    if (BG[(int)(layer_order[i_rows]-1)*H_cols + i_cols] == -1)
                    {
                        continue;
                    }
                    shiftNum = (int)BG[(int)(layer_order[i_rows]-1)*H_cols + i_cols];

                    // deal with the special cases that the minina belongs to the current column
                    for (i_z = 0; i_z < Z; i_z++)
                    {
                        if (minIdx[i_z] == i_cols)
                        {
                            minVec[i_z] = min2[i_z];
                        }
                        else
                        {
                            minVec[i_z] = min1[i_z];
                        }
                    }

                    for (i_z = 0; i_z < Z; i_z++)
                    {
                        RMem[(int)(layer_order[i_rows]-1)*H_cols + i_cols][i_z] = sign[i_z]*signOP_mex(TMem[i_cols][i_z])*minVec[i_z]*norm;
                        Qtmp[i_z] = TMem[i_cols][i_z] + RMem[(int)(layer_order[i_rows]-1)*H_cols + i_cols][i_z];

                        if (HD_mex(Qtmp[i_z]) != HD_mex((float)QSign[i_cols][i_z]))
                        {
                            TDFlag = 0;
                        }

                        SynCheck[i_z] = SynCheck[i_z]^HD_mex(Qtmp[i_z]);
                    }

                    // inverse rotation
                    Rotation_mex(Qtmp, Z - shiftNum, Z, QMem[i_cols]);
                }

                // parity check
                if (sum_mex((float *)SynCheck, Z) > 0) // change the type of matrix pointer
                {
                    TDFlag = 0;
                }
            }

            if (TDFlag == 1)
            {   
                break;
            }
        }

        // ouput the results
        for (i_rows = 0; i_rows < H_cols; i_rows++)
        {
            for (i_cols = 0; i_cols < Z; i_cols++)
            {
                debits[i_frame*H_cols*Z + i_rows*Z + i_cols] = (double)HD_mex(QMem[i_rows][i_cols]);
            }
        }
        deitera[0] = deitera[0] + (double)i_itera;
    }

    mxFree(QMem);
    mxFree(RMem);
    mxFree(TMem);
    mxFree(QSign);

    mxFree(min1);
    mxFree(min2);
    mxFree(sign);
    mxFree(minIdx);
    mxFree(minVec);
    mxFree(SynCheck);
    mxFree(Qc);
    mxFree(TMemAbs);
    mxFree(Qtmp);
}

void Rotation_mex (float *xIn, int shiftNum, int Z, float *yOut)
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

float max_mex (float x, float y)
{
    if (x >= y)
    {
        return x;
    }
    else
    {
        return y;
    }
}
