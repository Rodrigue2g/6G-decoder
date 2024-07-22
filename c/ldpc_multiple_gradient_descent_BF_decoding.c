#include "mex.h"
#include "math.h"
#include "float.h"
#define H_COLS

// prototypes
//void   sc_ldpc_layered_nms_float_decoding_mex(double *LLRIn, double *BGWin, double *BGTerm, int winSize, int cpd_L, int cpd_w, int Z, int itera, int nb, int mb, int Nb, int Mb, double alpha, int Frames, double *debits, double *deiter);
void    ldpc_multiple_gradient_descent_BF_decoding(int **H, int itera, double **LLRIn, int frames, double delta, int **decode_bits, int *decode_itera);
void    Rotation_mex (float *xIn, int shiftNum, int Z, float *yOut);
int     signOP_mex (float x);
int     HD_mex (float x);
float   sum_mex (float *x, int size);

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double  *param_H;
    double  *param_LLRIn;
    int      param_itera;
    int      param_frames;
    double   param_delta;
    double  *param_debits;
    double  *param_deiter;

    param_H         = mxGetPr(prhs[0]);
    param_LLRIn     = mxGetPr(prhs[1]);
    param_itera     = (int)mxGetScalar(prhs[2]);
    param_frames    = (int)mxGetScalar(prhs[3]);
    param_delta     = mxGetScalar(prhs[4]);

    // mwSize m = mxGetM(prhs[0]);
    // mwSize n = mxGetN(prhs[0]);
    // mexPrintf("rows = %d\n", m);
    // mexPrintf("columns = %d\n", n);

    // Note that any output must have an array expression
    plhs[0]       = mxCreateDoubleMatrix(param_frames*param_delta, 1, mxREAL);
    param_debits  = mxGetPr(plhs[0]);
    plhs[1]       = mxCreateDoubleMatrix(1, 1, mxREAL);
    param_deiter  = mxGetPr(plhs[1]);

    ldpc_multiple_gradient_descent_BF_decoding(param_H, param_itera, param_LLRIn, param_frames, param_delta, param_debits, param_deiter);
    return;
}

void ldpc_multiple_gradient_descent_BF_decoding(int **H, int itera, double **LLRIn, int frames, double delta, int **decode_bits, int *decode_itera)
{
     // parameter and memory allocations
    int i_frame, i_run, i_itera, i_rows, i_cols, i_z;
    int TermBGCnt, shiftNum, TDFlag;

    int BGWin_cols = nb * (winSize + cpd_w - 1);                 // window length
    int BGWin_rows = mb * winSize;                               // window width
    int winDelay = winSize - 1;                                  // delay of sliding in and sliding out
    int run_len = cpd_L + winDelay;                              // running steps
    int code_len = Nb * run_len;                                 // single block length
    int BGTerm_cols = nb * (winSize + cpd_w - 1 + winDelay - 1); // length of auxilary BG for termination

    // mexPrintf("BGTerm_cols = %d\n", BGTerm_cols);

    // 1D, 2D, or 3D matrix
    float **PosteriorMem;
    float ***ExtrinsicMem;


    float **QMem;
    float ***RMem;
    float **TMem;
    int **QSign;

    int *BGWD;
    int **H;

    float *min1;
    float *min2;
    int *sign;
    int *minIdx;
    float *minVec;
    int *SynCheck;
    float *Qc;
    float *Qtmp;
    float *TMemAbs;

    PosteriorMem = (float **)mxMalloc(sizeof(float *) * BGWin_cols);
    ExtrinsicMem = (float ***)mxMalloc(sizeof(float **) * BGWin_rows);
    //H = (int **)mxMalloc(sizeof(int *) * mb*(winSize+cpd_w-1)*Z);

    H = (int **)mxMalloc(sizeof(int *) * mb*Z*(cpd_L+cpd_w-1));

    QMem = (float **)mxMalloc(sizeof(float *) * BGWin_cols);
    RMem = (float ***)mxMalloc(sizeof(float **) * BGWin_rows);
    TMem = (float **)mxMalloc(sizeof(float *) * BGWin_cols);
    QSign = (int **)mxMalloc(sizeof(int *) * BGWin_cols);

    BGWD = (int *)mxMalloc(sizeof(int) * BGWin_rows * BGWin_cols);

    min1 = (float *)mxMalloc(sizeof(float) * Z);
    min2 = (float *)mxMalloc(sizeof(float) * Z);
    sign = (int *)mxMalloc(sizeof(int) * Z);
    minIdx = (int *)mxMalloc(sizeof(int) * Z);
    minVec = (float *)mxMalloc(sizeof(float) * Z);
    SynCheck = (int *)mxMalloc(sizeof(int) * Z);
    Qc = (float *)mxMalloc(sizeof(float) * Z);
    Qtmp = (float *)mxMalloc(sizeof(float) * Z);
    TMemAbs = (float *)mxMalloc(sizeof(float) * Z);

    // open multiple dimensional matrix
    for (i_rows = 0; i_rows < BGWin_cols; i_rows++)
    {
        PosteriorMem[i_rows] = (float *)mxMalloc(sizeof(float) * Z);
        QMem[i_rows] = (float *)mxMalloc(sizeof(float) * Z);
        TMem[i_rows] = (float *)mxMalloc(sizeof(float) * Z);
        QSign[i_rows] = (int *)mxMalloc(sizeof(int) * Z);
    }


    for(i_rows=0; i_rows < mb*(cpd_L+cpd_w-1)*Z; i_rows++){
        H[i_rows] = (int *)mxMalloc(sizeof(int) * (nb*cpd_L*Z+128));
    }

    for (i_rows = 0; i_rows < BGWin_rows; i_rows++)
    {
        ExtrinsicMem[i_rows] = (float **)mxMalloc(sizeof(float *) * BGWin_cols);
        RMem[i_rows] = (float **)mxMalloc(sizeof(float *) * BGWin_cols);
        for (i_cols = 0; i_cols < BGWin_cols; i_cols++)
        {
            ExtrinsicMem[i_rows][i_cols] = (float *)mxMalloc(sizeof(float) * Z);
            RMem[i_rows][i_cols] = (float *)mxMalloc(sizeof(float) * Z);
        }
    }

    
    // copy the H matrix
    for (i_rows = 0; i_rows < mb*(cpd_L+cpd_w-1)*Z; i_rows++)
    {
        for (i_cols = 0; i_cols < nb*cpd_L*Z+128; i_cols++)
        {
            H[i_rows][i_cols] = Hwin[i_rows * (nb*cpd_L*Z+128) + i_cols];
        }
    }
        int total_itera = 0;
    int correction = 0;
    
    for(int i_frame=0; i_frame<Frames; i_frame++){
        if(!correction){
            for(int i_frame=0; i_frame<Frames; i_frame++){
                for(int i=0; i<32000; i++){
                    debits[i+32000*i_frame] = LLRIn[i+i_frame*34560]<0;
                }
            }
        }
        else{
                for(int i=0; i<49; i++){
                    total_itera+= test(&debits[i*640+i_frame*32000], H, 10000, &LLRIn[i*640+i_frame*34560], mb, nb, Z, winSize, cpd_w, (i==-1));
                    
                }
            }
        correction = 1;
    }

    printf("Total iterations: %d\n", total_itera);  
  

    mxFree(PosteriorMem);
    mxFree(ExtrinsicMem);
    mxFree(QMem);
    mxFree(RMem);
    mxFree(TMem);
    mxFree(QSign);


    mxFree(BGWD);

    for(i_rows=0; i_rows < mb*(cpd_L+cpd_w-1)*Z; i_rows++){
        mxFree(H[i_rows]);
    }

    mxFree(H);

    mxFree(min1);
    mxFree(min2);
    mxFree(sign);
    mxFree(minIdx);
    mxFree(minVec);
    mxFree(SynCheck);
    mxFree(Qc);
    mxFree(TMemAbs);
    mxFree(Qtmp);

/*
    int H_rows = sizeof(H) / sizeof(H[0]);
    int H_cols = sizeof(H[0]) / sizeof(H[0][0]);
    int i_frames, i_itera, i_cols;

    double *inverFunc = (double *)mxMalloc(H_cols * sizeof(double));

    for (i_frames = 0; i_frames < frames; i_frames++) 
    {
        // Synchrome/Objective function Calculation
        double *xhat = (double *)mxMalloc(H_cols * sizeof(double));
        int *syndrome = (int *)mxMalloc(H_cols * sizeof(int));

        for (i_cols = 0; i_cols < H_cols; i_cols++) 
        {
            xhat[i_cols] = (LLRIn[i_frames][i_cols] < 0) ? 1.0 : 0.0; // hard decision
            syndrome[i_cols] = 0;
        }

        double objectFunc = 0.0;
        int modeFlag = 1;

        for (i_itera = 0; i_itera < itera; i_itera++) 
        {
            for (i_cols = 0; i_cols < H_cols; i_cols++) 
            {
                double sum_val = 0.0;
                for (int row = 0; row < H_rows; row++) 
                {
                    sum_val += H[row][i_cols] * (1 - 2 * syndrome[row]);
                }
                inverFunc[i_cols] = sum_val + (1 - 2 * xhat[i_cols]) * LLRIn[i_frames][i_cols];
            }

            // Sort inverFunc in ascending order -- sorting algorithm?

            int flipIdx;
            if (modeFlag == 1) {
                // delta check + BF
            } else {
                // BF without delta check
            }

            // Update objective function
            double objectLast = objectFunc;
            objectFunc = 0.0;
            for (i_cols = 0; i_cols < H_cols; i_cols++) 
            {
                objectFunc += (1 - 2 * xhat[i_cols]) * LLRIn[i_frames][i_cols];
            }
            for (i_cols = 0; i_cols < H_cols; i_cols++) 
            {
                objectFunc += (1 - 2 * syndrome[i_cols]);
            }

            if (objectFunc < objectLast) {
                modeFlag = 0;
            }

            // Check if syndrome is all zeros
            int syndrome_sum = 0;
            for (i_cols = 0; i_cols < H_cols; i_cols++) 
            {
                syndrome_sum += syndrome[i_cols];
            }
            if (syndrome_sum == 0) 
                break;
        }

        // Copy xhat to decode_bits and update decode_itera
        for (i_cols = 0; i_cols < H_cols; i_cols++) 
        {
            decode_bits[i_frames][i_cols] = (int)xhat[i_cols];
        }
        decode_itera[i_frames] = i_itera;

        // return the LLrs as well?

        mxFree(xhat);
        mxFree(syndrome);
    }

    mxFree(inverFunc);
    */

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

int signOP_mex(float x)
{
    return ((x >= 0) - (x < 0));
}

int HD_mex(float x)
{
    return (int)(x < 0);
}

float sum_mex(float *x, int size)
{   
    float sum = 0;
    int i;
    for (i = 0; i < size; i++)
    {
        sum = sum + x[i];
    }
    return sum;
}