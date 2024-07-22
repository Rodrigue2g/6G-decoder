#include "mex.h"
#include "math.h"
#include "float.h"
// REVISIT: This window decoding refers to TIT 2012 'Window decoding of LDPC-CC on erasure channels'. However, since we still instantiate memories for the decoded ms-1 region (the memory of the code), it is slightly resource-wasting

// prototypes
void   sc_ldpc_layered_nms_float_decoding_mex(double *LLRIn, double *BGWin, double *BGTerm, int winSize, int cpd_L, int cpd_w, int Z, int itera, int nb, int mb, int Nb, int Mb, double alpha, int Frames, double *debits, double *deiter);
void   Rotation_mex (float *xIn, int shiftNum, int Z, float *yOut);
int    signOP_mex (float x);
int    HD_mex (float x);
float  sum_mex (float *x, int size);

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double  *param_LLRIn;
    double  *param_BGWin;
    double  *param_BGTerm;
    int      param_winSize;
    int      param_cpd_L;
    int      param_cpd_w;
    int      param_Z;
    int      param_itera;
    int      param_nb;
    int      param_mb;
    int      param_Nb;
    int      param_Mb;
    double   param_alpha;
    int      param_Frames;
    double  *param_debits;
    double  *param_deiter;

    param_LLRIn  = mxGetPr(prhs[0]);
    param_BGWin  = mxGetPr(prhs[1]);
    param_BGTerm = mxGetPr(prhs[2]);
    param_winSize = (int)mxGetScalar(prhs[3]);
    param_cpd_L   = (int)mxGetScalar(prhs[4]);
    param_cpd_w   = (int)mxGetScalar(prhs[5]);
    param_Z       = (int)mxGetScalar(prhs[6]);
    param_itera   = (int)mxGetScalar(prhs[7]);
    param_nb      = (int)mxGetScalar(prhs[8]);
    param_mb      = (int)mxGetScalar(prhs[9]);
    param_Nb      = (int)mxGetScalar(prhs[10]);
    param_Mb      = (int)mxGetScalar(prhs[11]);
    param_alpha   = mxGetScalar(prhs[12]);
    param_Frames  = (int)mxGetScalar(prhs[13]);

    // mwSize m = mxGetM(prhs[0]);
    // mwSize n = mxGetN(prhs[0]);
    // mexPrintf("rows = %d\n", m);
    // mexPrintf("columns = %d\n", n);

    // Note that any output must have an array expression
    plhs[0]       = mxCreateDoubleMatrix(param_Frames*param_cpd_L*param_Nb, 1, mxREAL);
    param_debits  = mxGetPr(plhs[0]);
    plhs[1]       = mxCreateDoubleMatrix(1, 1, mxREAL);
    param_deiter  = mxGetPr(plhs[1]);

    sc_ldpc_layered_nms_float_decoding_mex(param_LLRIn, param_BGWin, param_BGTerm, param_winSize, param_cpd_L, param_cpd_w, param_Z, param_itera, param_nb, param_mb, param_Nb, param_Mb, param_alpha, param_Frames, param_debits, param_deiter);
    return;
}

void sc_ldpc_layered_nms_float_decoding_mex(double *LLRIn, double *BGWin, double *BGTerm, int winSize, int cpd_L, int cpd_w, int Z, int itera, int nb, int mb, int Nb, int Mb, double alpha, int Frames, double *debits, double *deiter)
{
    // parameter and memory allocations
    int i_frame, i_run, i_itera, i_rows, i_cols, i_z;
    int TermBGCnt, shiftNum, TDFlag;

    int BGWin_cols  = nb*(winSize + cpd_w - 1); // window length
    int BGWin_rows  = mb*winSize;               // window width
    int winDelay    = winSize - 1;              // delay of sliding in and sliding out
    int run_len     = cpd_L + winDelay;         // running steps
    int code_len    = Nb*run_len;               // single block length
    int BGTerm_cols = nb*(winSize + cpd_w - 1 + winDelay - 1); // length of auxilary BG for termination

    // mexPrintf("BGTerm_cols = %d\n", BGTerm_cols);

    // 1D, 2D, or 3D matrix
    float  **PosteriorMem;
    float ***ExtrinsicMem;

    float  **QMem;
    float ***RMem;
    float  **TMem;
    int    **QSign;

    int     *BGWD;

    float   *min1;
    float   *min2;
    int     *sign;
    int     *minIdx;
    float   *minVec;
    int     *SynCheck;
    float   *Qc;
    float   *Qtmp;
    float   *TMemAbs;
 
    PosteriorMem = (float**)mxMalloc(sizeof(float*)*BGWin_cols);
    ExtrinsicMem = (float***)mxMalloc(sizeof(float**)*BGWin_rows);

    QMem         = (float**)mxMalloc(sizeof(float*)*BGWin_cols);
    RMem         = (float***)mxMalloc(sizeof(float**)*BGWin_rows);
    TMem         = (float**)mxMalloc(sizeof(float*)*BGWin_cols);
    QSign        = (int**)mxMalloc(sizeof(int*)*BGWin_cols);

    BGWD         = (int*)mxMalloc(sizeof(int)*BGWin_rows*BGWin_cols);

    min1         = (float*)mxMalloc(sizeof(float)*Z);
    min2         = (float*)mxMalloc(sizeof(float)*Z);
    sign         = (int*)mxMalloc(sizeof(int)*Z);
    minIdx       = (int*)mxMalloc(sizeof(int)*Z);
    minVec       = (float*)mxMalloc(sizeof(float)*Z);
    SynCheck     = (int*)mxMalloc(sizeof(int)*Z);
    Qc           = (float*)mxMalloc(sizeof(float)*Z);
    Qtmp         = (float*)mxMalloc(sizeof(float)*Z);
    TMemAbs      = (float*)mxMalloc(sizeof(float)*Z);
    
    // open multiple dimensional matrix
    for (i_rows = 0; i_rows < BGWin_cols; i_rows++)
    {
        PosteriorMem[i_rows] = (float*)mxMalloc(sizeof(float)*Z);
        QMem[i_rows]         = (float*)mxMalloc(sizeof(float)*Z);
        TMem[i_rows]         = (float*)mxMalloc(sizeof(float)*Z);
        QSign[i_rows]        = (int*)mxMalloc(sizeof(int)*Z);
    }

    for (i_rows = 0; i_rows < BGWin_rows; i_rows++)
    {
        ExtrinsicMem[i_rows] = (float**)mxMalloc(sizeof(float*)*BGWin_cols);
        RMem[i_rows]         = (float**)mxMalloc(sizeof(float*)*BGWin_cols);
        for (i_cols = 0; i_cols < BGWin_cols; i_cols++)
        {
            ExtrinsicMem[i_rows][i_cols] = (float*)mxMalloc(sizeof(float)*Z);
            RMem[i_rows][i_cols]         = (float*)mxMalloc(sizeof(float)*Z);
        }
    }

    // Level 1 preparation: frame
    for (i_frame = 0; i_frame < Frames; i_frame++)
    {
        TermBGCnt = 0;

        // initilization
        for (i_rows = 0; i_rows < BGWin_cols; i_rows++)
        {
            for (i_z = 0; i_z < Z; i_z++)
            {
                PosteriorMem[i_rows][i_z]  = 1000;
            }
        }

        for (i_rows = 0; i_rows < BGWin_rows; i_rows++)
        {
            for (i_cols = 0; i_cols < BGWin_cols; i_cols++)
            {
                for (i_z = 0; i_z < Z; i_z++)
                {
                    ExtrinsicMem[i_rows][i_cols][i_z] = 0;
                }
            }
        }

        // Level 2 preparation: window
        for (i_run = 0; i_run < run_len; i_run++)
        {
            // shift two memories, load channel LLRs, and clean
            for (i_rows = 0; i_rows < BGWin_cols; i_rows++)
            {
                for (i_z = 0; i_z < Z; i_z++)
                {
                    if (i_rows < (BGWin_cols - nb))
                    {
                        PosteriorMem[i_rows][i_z] = PosteriorMem[i_rows + nb][i_z];
                    }
                    else
                    {
                        PosteriorMem[i_rows][i_z] = (float)LLRIn[i_frame*code_len + i_run*Nb + (i_rows-(BGWin_cols - nb))*Z + i_z];
                    }

                }
            }

            for (i_rows = 0; i_rows < BGWin_rows; i_rows++)
            {
                for (i_cols = 0; i_cols < BGWin_cols; i_cols++)
                {
                    for (i_z = 0; i_z < Z; i_z++)
                    {
                        if (i_rows < (BGWin_rows - mb) && i_cols < (BGWin_cols - nb))
                        {
                            ExtrinsicMem[i_rows][i_cols][i_z] = ExtrinsicMem[i_rows + mb][i_cols + nb][i_z];
                        }
                        else
                        {
                            ExtrinsicMem[i_rows][i_cols][i_z] = 0;
                        }
                    }
                }
            }

            // load PCM of sliding window
            if (i_run >= cpd_L)
            {
                for (i_rows = 0; i_rows < BGWin_rows; i_rows++)
                {
                    for (i_cols = 0; i_cols < BGWin_cols; i_cols++)
                    {
                        BGWD[i_rows*BGWin_cols + i_cols] = (int)BGTerm[(TermBGCnt*mb + i_rows)*BGTerm_cols + (TermBGCnt*nb + i_cols)];
                    }
                }
                TermBGCnt = TermBGCnt + 1;
            }
            else
            {
                for (i_rows = 0; i_rows < BGWin_rows; i_rows++)
                {
                    for (i_cols = 0; i_cols < BGWin_cols; i_cols++)
                    {
                        BGWD[i_rows*BGWin_cols + i_cols] = (int)BGWin[i_rows*BGWin_cols + i_cols];
                    }
                }
            }

            // initilization for common LDPC decoding
            for (i_rows = 0; i_rows < BGWin_cols; i_rows++)
            {
                for (i_z = 0; i_z < Z; i_z++)
                {
                    QMem[i_rows][i_z]  = PosteriorMem[i_rows][i_z];
                    TMem[i_rows][i_z]  = 0;
                    QSign[i_rows][i_z] = 0;
                }
            }

            for (i_rows = 0; i_rows < BGWin_rows; i_rows++)
            {
                for (i_cols = 0; i_cols < BGWin_cols; i_cols++)
                {
                    for (i_z = 0; i_z < Z; i_z++)
                    {
                        RMem[i_rows][i_cols][i_z] = ExtrinsicMem[i_rows][i_cols][i_z];
                    }
                }
            }

            // Level 3 preparation: iteration
            for (i_itera = 0; i_itera < itera; i_itera++)
            {
                TDFlag = 1;

                for (i_rows = 0; i_rows < BGWin_rows; i_rows++)
                {
                    // initilization
                    for (i_z = 0; i_z < Z; i_z++)
                    {
                        min1[i_z]   = FLT_MAX; // head file "float.h"
                        min2[i_z]   = FLT_MAX;
                        sign[i_z]   = 1;

                        minIdx[i_z] = 0;
                        minVec[i_z] = 0;

                        SynCheck[i_z] = 0;
                    }

                    // PHASE 1: MIN
                    for (i_cols = 0; i_cols < BGWin_cols; i_cols++)
                    {
                        // checkout whether there is -1 in base matrix, if yes, directly skip
                        if (BGWD[i_rows*BGWin_cols + i_cols] == -1)
                        {
                            continue;
                        }
                        shiftNum = BGWD[i_rows*BGWin_cols + i_cols];

                        // calculate the Qc, TMem, and its sign and magnitude
                        Rotation_mex(QMem[i_cols], shiftNum, Z, Qc);

                        for (i_z = 0; i_z < Z; i_z++)
                        {
                            TMem[i_cols][i_z]  = Qc[i_z] - RMem[i_rows][i_cols][i_z];
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
                    for (i_cols = 0; i_cols < BGWin_cols; i_cols++)
                    {
                        // checkout whether there is -1 in base matrix, if yes, directly skip
                        if (BGWD[i_rows*BGWin_cols + i_cols] == -1)
                        {
                            continue;
                        }
                        shiftNum = BGWD[i_rows*BGWin_cols + i_cols];

                        if (i_cols < (cpd_w-1)*nb)
                        {
                            Rotation_mex(QMem[i_cols], shiftNum, Z, Qtmp);

                            for (i_z = 0; i_z < Z; i_z++)
                            {
                                SynCheck[i_z] = SynCheck[i_z]^HD_mex(Qtmp[i_z]);
                            }
                        }
                        else
                        {
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
                                RMem[i_rows][i_cols][i_z] = (float)sign[i_z]*signOP_mex(TMem[i_cols][i_z])*minVec[i_z]*alpha;
                                Qtmp[i_z] = TMem[i_cols][i_z] + RMem[i_rows][i_cols][i_z];
                            
                                if (HD_mex(Qtmp[i_z]) != HD_mex((float)QSign[i_cols][i_z]))
                                {
                                    TDFlag = 0;
                                }

                                SynCheck[i_z] = SynCheck[i_z]^HD_mex(Qtmp[i_z]);
                            }

                            // inverse rotation
                            Rotation_mex(Qtmp, Z - shiftNum, Z, QMem[i_cols]);
                        }
                    }

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

            // update two memories
            for (i_rows = 0; i_rows < BGWin_cols; i_rows++)
            {
                for (i_z = 0; i_z < Z; i_z++)
                {
                    PosteriorMem[i_rows][i_z] = QMem[i_rows][i_z];
                }
            }

            for (i_rows = 0; i_rows < BGWin_rows; i_rows++)
            {
                for (i_cols = 0; i_cols < BGWin_cols; i_cols++)
                {
                    for (i_z = 0; i_z < Z; i_z++)
                    {
                        ExtrinsicMem[i_rows][i_cols][i_z] = RMem[i_rows][i_cols][i_z];
                    }
                }
            }

            // output decoded results
            if (i_run >= winDelay)
            {
                for (i_rows = 0; i_rows < nb; i_rows++)
                {
                    for (i_z = 0; i_z < Z; i_z++)
                    {
                        debits[i_frame*cpd_L*Nb + (i_run - winDelay)*Nb + i_rows*Z + i_z] = (double)HD_mex(PosteriorMem[i_rows + (cpd_w - 1)*nb][i_z]);
                    }
                }

                deiter[0] = deiter[0] + i_itera;
            }
        }
    }

    mxFree(PosteriorMem);
    mxFree(ExtrinsicMem);
    mxFree(QMem);
    mxFree(RMem);
    mxFree(TMem);
    mxFree(QSign);

    mxFree(BGWD);

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