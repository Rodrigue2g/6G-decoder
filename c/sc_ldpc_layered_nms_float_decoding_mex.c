#include "mex.h"
#include "math.h"
#include "float.h"
// REVISIT: This window decoding refers to TIT 2012 'Window decoding of LDPC-CC on erasure channels'. However, since we still instantiate memories for the decoded ms-1 region (the memory of the code), it is slightly resource-wasting

#define D(i) printf("debug: %d \n", i);

// prototypes
void sc_ldpc_layered_nms_float_decoding_mex(double *LLRIn, double *BGWin, double *BGTerm, int winSize, int cpd_L, int cpd_w, int Z, int itera, int nb, int mb, int Nb, int Mb, double alpha, int Frames, double *debits, double *deiter, double *Hwin);
void Rotation_mex(float *xIn, int shiftNum, int Z, float *yOut);
int signOP_mex(float x);
int HD_mex(float x);
float sum_mex(float *x, int size);
int ldpc_gradient_BF_decoding_mex(double *debits, int **H, int itera, double *LLRIn, int mb, int nb, int Z, int WinSize, int w);
int test(double *debits, int **H, int itera, double *LLRIn, int mb, int nb, int Z, int WinSize, int w, int firstRun);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *param_LLRIn;
    double *param_BGWin;
    double *param_BGTerm;
    double *param_Hwin;
    int param_winSize;
    int param_cpd_L;
    int param_cpd_w;
    int param_Z;
    int param_itera;
    int param_nb;
    int param_mb;
    int param_Nb;
    int param_Mb;
    double param_alpha;
    int param_Frames;
    double *param_debits;
    double *param_deiter;
    

    param_LLRIn = mxGetPr(prhs[0]);
    param_BGWin = mxGetPr(prhs[1]);
    param_BGTerm = mxGetPr(prhs[2]);
    param_Hwin = mxGetPr(prhs[3]);
    param_winSize = (int)mxGetScalar(prhs[4]);
    param_cpd_L = (int)mxGetScalar(prhs[5]);
    param_cpd_w = (int)mxGetScalar(prhs[6]);
    param_Z = (int)mxGetScalar(prhs[7]);
    param_itera = (int)mxGetScalar(prhs[8]);
    param_nb = (int)mxGetScalar(prhs[9]);
    param_mb = (int)mxGetScalar(prhs[10]);
    param_Nb = (int)mxGetScalar(prhs[11]);
    param_Mb = (int)mxGetScalar(prhs[12]);
    param_alpha = mxGetScalar(prhs[13]);
    param_Frames = (int)mxGetScalar(prhs[14]);
    


    // mwSize m = mxGetM(prhs[0]);
    // mwSize n = mxGetN(prhs[0]);
    // mexPrintf("rows = %d\n", m);
    // mexPrintf("columns = %d\n", n);

    // Note that any output must have an array expression
    plhs[0] = mxCreateDoubleMatrix(param_Frames * param_cpd_L * param_Nb, 1, mxREAL);
    param_debits = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    param_deiter = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

    sc_ldpc_layered_nms_float_decoding_mex(param_LLRIn, param_BGWin, param_BGTerm, param_winSize, param_cpd_L, param_cpd_w, param_Z, param_itera, param_nb, param_mb, param_Nb, param_Mb, param_alpha, param_Frames, param_debits, param_deiter, param_Hwin);
    return;
}

void sc_ldpc_layered_nms_float_decoding_mex(double *LLRIn, double *BGWin, double *BGTerm, int winSize, int cpd_L, int cpd_w, int Z, int itera, int nb, int mb, int Nb, int Mb, double alpha, int Frames, double *debits, double *deiter, double *Hwin)
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

    // Level 1 preparation: frame
    /*for (i_frame = 0; i_frame < Frames; i_frame++)
    {
        TermBGCnt = 0;

        // initilization
        for (i_rows = 0; i_rows < BGWin_cols; i_rows++)
        {
            for (i_z = 0; i_z < Z; i_z++)
            {
                PosteriorMem[i_rows][i_z] = 1000;
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
            // shift two memories, load channel LLRs, and clean   #### 1. (Receive and shift)
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
                        PosteriorMem[i_rows][i_z] = (float)LLRIn[i_frame * code_len + i_run * Nb + (i_rows - (BGWin_cols - nb)) * Z + i_z];
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
                        BGWD[i_rows * BGWin_cols + i_cols] = (int)BGTerm[(TermBGCnt * mb + i_rows) * BGTerm_cols + (TermBGCnt * nb + i_cols)];
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
                        BGWD[i_rows * BGWin_cols + i_cols] = (int)BGWin[i_rows * BGWin_cols + i_cols];
                    }
                }
            }

            // initilization for common LDPC decoding
            for (i_rows = 0; i_rows < BGWin_cols; i_rows++)
            {
                for (i_z = 0; i_z < Z; i_z++)
                {
                    QMem[i_rows][i_z] = PosteriorMem[i_rows][i_z];
                    TMem[i_rows][i_z] = 0;
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
                        min1[i_z] = FLT_MAX; // head file "float.h"
                        min2[i_z] = FLT_MAX;
                        sign[i_z] = 1;

                        minIdx[i_z] = 0;
                        minVec[i_z] = 0;

                        SynCheck[i_z] = 0;
                    }

                    // PHASE 1: MIN
                    for (i_cols = 0; i_cols < BGWin_cols; i_cols++)
                    {
                        // checkout whether there is -1 in base matrix, if yes, directly skip
                        if (BGWD[i_rows * BGWin_cols + i_cols] == -1)
                        {
                            continue;
                        }
                        shiftNum = BGWD[i_rows * BGWin_cols + i_cols];

                        // calculate the Qc, TMem, and its sign and magnitude
                        Rotation_mex(QMem[i_cols], shiftNum, Z, Qc);

                        for (i_z = 0; i_z < Z; i_z++)
                        {
                            TMem[i_cols][i_z] = Qc[i_z] - RMem[i_rows][i_cols][i_z];
                            sign[i_z] = sign[i_z] * signOP_mex(TMem[i_cols][i_z]);
                            TMemAbs[i_z] = fabs(TMem[i_cols][i_z]); // abs() only for int
                            QSign[i_cols][i_z] = signOP_mex(Qc[i_z]);
                        }

                        // search the min1, min2, and the column index corresponding to min1
                        for (i_z = 0; i_z < Z; i_z++)
                        {
                            if (TMemAbs[i_z] < min1[i_z])
                            {
                                min2[i_z] = min1[i_z];
                                min1[i_z] = TMemAbs[i_z];
                                minIdx[i_z] = i_cols;
                            }
                            else if (TMemAbs[i_z] < min2[i_z])
                            {
                                min2[i_z] = TMemAbs[i_z];
                            }
                        }
                    }

                    // PHASE 2: Q- and R- messages update
                    for (i_cols = 0; i_cols < BGWin_cols; i_cols++)
                    {
                        // checkout whether there is -1 in base matrix, if yes, directly skip
                        if (BGWD[i_rows * BGWin_cols + i_cols] == -1)
                        {
                            continue;
                        }
                        shiftNum = BGWD[i_rows * BGWin_cols + i_cols];

                        if (i_cols < (cpd_w - 1) * nb)
                        {
                            Rotation_mex(QMem[i_cols], shiftNum, Z, Qtmp);

                            for (i_z = 0; i_z < Z; i_z++)
                            {
                                SynCheck[i_z] = SynCheck[i_z] ^ HD_mex(Qtmp[i_z]);
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
                                RMem[i_rows][i_cols][i_z] = (float)sign[i_z] * signOP_mex(TMem[i_cols][i_z]) * minVec[i_z] * alpha;
                                Qtmp[i_z] = TMem[i_cols][i_z] + RMem[i_rows][i_cols][i_z];

                                if (HD_mex(Qtmp[i_z]) != HD_mex((float)QSign[i_cols][i_z]))
                                {
                                    TDFlag = 0;
                                }

                                SynCheck[i_z] = SynCheck[i_z] ^ HD_mex(Qtmp[i_z]);
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
            //if (i_run<cpd_L+winDelay) // if (i_run >= winDelay)
            //{
                if(TDFlag==1){
                    for (i_rows = 0; i_rows < nb; i_rows++)
                    {
                        for (i_z = 0; i_z < Z; i_z++)
                        {
                            //if(TDFlag==1)
                            //printf("PosteriorMem[%d + (%d - %d)*%d][%d] = %f\n", i_rows, cpd_w, 1, nb, i_z, PosteriorMem[i_rows + (cpd_w - 1)*nb][i_z]);
                            debits[i_frame*cpd_L*Nb + (i_run - winDelay)*Nb + i_rows*Z + i_z] = (double)HD_mex(PosteriorMem[i_rows + (cpd_w - 1)*nb][i_z]);
                                // the bits respect the parity check matrix so we output them to the outside
                            //else ;
                                //LLRin=PosteriorMem

                                //ldpc_gradient_BF_decoding_mex(H, itera, (double)HD_mex(PosteriorMem[i_rows + (cpd_w - 1)*nb][i_z]));
                                // reinject the corrected llr's into the windowed decoder
                                // output the corrected bits to the outside

                        }
                    }
                }
                else{


                //ldpc_gradient_BF_decoding_mex(&debits[i_frame*cpd_L*Nb + (i_run-winDelay)*Nb], H, 2, &LLRIn[(i_run-winDelay)*512], mb, nb, Z, winSize, cpd_w);
                
                
                
                // Check improvement relative to no correction
                for (int i=0; i<Nb; i++){
                    debits[i_frame*cpd_L*Nb + (i_run-winDelay)*Nb + i] = LLRIn[i_frame*cpd_L*Nb + (i_run-winDelay)*Nb + i]<0;
                }



                deiter[0] = deiter[0] + i_itera;

    
            //}
        
        }
    }
    */

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
                    total_itera+= test(&debits[i*640+i_frame*32000], H, 100, &LLRIn[i*640+i_frame*34560], mb, nb, Z, winSize, cpd_w, (i==-1));
                    
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

int test(double *debits, int **H, int itera, double *LLRIn, int mb, int nb, int Z, int WinSize, int w, int firstRun){

    double delta = -0.6;

    // window size
    int M = mb * Z;
    int N = w * nb * Z; 
    int offset_x = 0;
    int offset_y = 128;

    if(firstRun){
        N=nb*Z;
        offset_x=640;
    }

    // LLR
    double *LLR = (double *)mxMalloc(N * sizeof(double));
    for (int i = 0; i < N; i++)
    {  // N
        LLR[i] = LLRIn[i]; //LLRIn[i / Z][i % Z];
    }


    // inverFunc
    double *inverFunc = (double *)mxMalloc(sizeof(double) * N);

    // xhat
    int *xhat = mxMalloc(sizeof(int) * N);

    // STEP 1: hard decision //* OK
    for (int i = 0; i < N; i++)
    {
        xhat[i] = (LLR[i] >= 0) ? 1 : -1;
    }

    bool flag = 0;

    int *Idx = mxMalloc(sizeof(int) * N);

    int i_itera;
    for (i_itera = 0; i_itera < itera; i_itera++){

        // Check parity check matrix //* OK
        bool holds = true;
        for (int i = 0; i < M; i++)
        {
            int accum = 1;
            for (int j = 0; j < N; j++)
            {   
                if (H[i+offset_y][j+offset_x] == 1)
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

        // INVERSION FUNCTION
        int lowestIdx = 0;
        for (int k = 0; k < N; k++)
        {
            inverFunc[k] = xhat[k] * LLR[k];
            for (int mk = 0; mk < M; mk++)
            {
                if (H[mk+offset_y][k+offset_x] == 1)
                {
                    double accum = 1;
                    for (int nk = 0; nk < N; nk++)
                    {   
                        if(H[mk+offset_y][nk+offset_x]==1)
                            accum *= xhat[nk];
                    }
                    inverFunc[k] += accum;
                }
            }
        
            // find the indices that are smaller than delta
            if (inverFunc[k] < delta)
            {
                Idx[k] = 1;
            }
            else
            {
                Idx[k] = 0;
            }
            // find the smallest value in inverFunc
            if(inverFunc[k] < inverFunc[lowestIdx])
            {
                lowestIdx = k;
            }
        }

        // flip single bit
        xhat[lowestIdx] = -xhat[lowestIdx];


        // STEP 3: objective function calculation
        /*if (flag == 0) 
        {       
            double objectFunc1 = 0;
            for (int j = 0; j < N; j++)
            {
                objectFunc1 += xhat[j] * LLR[j];
            }
            for (int i = 0; i < M; i++)
            {
                int accum = 1;
                for (int j = 0; j < N; j++)
                {   
                    if(H[i+offset_y][j+offset_x]==1)
                        accum *= xhat[j];
                }
                objectFunc1 += accum;
            }

            // flip the bits
            for (int k = 0; k < N; k++)
            {
                if (Idx[k] == 1)
                {  
                    xhat[k] = -xhat[k];
                }
            }

            // object function calculation
            double objectFunc2 = 0;
            for (int j = 0; j < N; j++)
            {
                objectFunc2 += xhat[j] * LLR[j];
            }
            for (int i = 0; i < M; i++)
            {
                int accum = 1;
                for (int j = 0; j < N; j++)
                {   
                    if(H[i+offset_y][j+offset_x]==1)
                        accum *= xhat[j];
                }
                objectFunc2 += accum;
            }

            if (objectFunc2 < objectFunc1)
            {
                flag = 1;
            }
        }
        else
        {   
            // flip the bit
            xhat[lowestIdx] = -xhat[lowestIdx];
        }*/





    }


    for (int i = 0; i < nb*Z; i++)
    {
        debits[i] = (xhat[i] < 0);
    }
    


    mxFree(LLR);
    mxFree(xhat);
    mxFree(inverFunc);
    mxFree(Idx);

    return i_itera;
}

int ldpc_gradient_BF_decoding_mex(double *debits, int **H, int itera, double *LLRIn, int mb, int nb, int Z, int WinSize, int w)
{
    // H: parity check matrix of the window
    // itera: number of iterations
    // LLRIn: input LLR (output of the window)

    /*
    The idea is to:
    1) Check if the decoded bits from the windowed decoder respect the parity check matrix
    2) If not, we need to correct the bits:
        2.1) We need to find the bits that are the most likely to be wrong
        2.2) We need to correct these bits until the parity check matrix is respected
    3) We output the corrected bits to the outside of the window and we reinject the most likely llr's in the windowed decoder
    */

    // mb=8 nb=40 Z=16 WinSize=5 w=2

    // M=5*8*16=640 N=6*40*16=3840
    double delta = -0.6;
    int M = WinSize * mb * Z;
    int N = (WinSize+w-1)  * nb * Z;
    
    int offset_x = 3200;
    int offset_y = 768;

    // Flatten LLRIn
    // create 3840x1 array of LLR's (y)
    double *LLR = (double *)mxMalloc(N * sizeof(double));
    for (int i = 0; i < N; i++)
    { // N
        LLR[i] = -LLRIn[i]; //LLRIn[i / Z][i % Z];
    }

    /*
    printf("LLRIn[%d] = %f\n", 1,   LLRIn[0]);
    printf("LLRIn[%d] = %f\n", 1+1, LLRIn[1]);
    printf("LLRIn[%d] = %f\n", 2+1, LLRIn[2]);
    printf("LLRIn[%d] = %f\n", 3+1, LLRIn[3]);
    printf("LLRIn[%d] = %f\n", 4+1, LLRIn[4]);
    printf("LLRIn[%d] = %f\n", 5+1, LLRIn[5]);
    */

    // correct window shape, we need to recontruct the windowed H matrix
    // H is a 640 x 3840 matrix

    /*int **H = (int **)mxCalloc(M, sizeof(int *));

    for (int i = 0; i < M; i++)
    {
        H[i] = (int *)mxCalloc(N, sizeof(int));
    }

    //* OK
    for (int i = 0; i < M / Z; i++)
    {
        for (int j = 0; j < N / Z; j++)
        {

            int a = Window[i * N / Z + j];
            for (int z = 0; z < Z && a != -1; z++)
            {
                H[z + Z * i][(z + a) % Z + Z * j] = 1;
            }
        }
    }*/

    // inverFunc: 3840x1 array
    double *inverFunc = (double *)mxMalloc(sizeof(double) * N);

    // xhat: 3840x1 array of hard decision bits
    int *xhat = mxMalloc(sizeof(int) * N);

    // STEP 1: hard decision //* OK
    for (int i = 0; i < N; i++)
    {
        xhat[i] = (LLR[i] > 0) ? 1 : -1;
    }

    static int hasPrinted = 1;
    static int cntr = 0;


    if (hasPrinted>0)
    {

        hasPrinted -= 1;

        if (hasPrinted == 0)
        {
            //for (int i = 0; i < 200; i++)
            //{
            //    printf("%d: xhat = %d and of H = %d %d\n", i+1,xhat[i], H[0][i], H[1][i]);
            //}

            /*for(int i=0; i<M; i++){
                for(int j=0; j<N; j++){
                    printf("%d ", H[i][j]);
                }
                printf("\n");
            }*/

        }
    }

    bool flag = 0;

    int *Idx = mxMalloc(sizeof(int) * N);

    bool holds = true;
    for (int i = 0; i <640; i++)
    {
        int accum = 0;
        for (int j = 0; j < N; j++)
        {
            
            accum += H[i+offset_y][j+offset_x] * (LLR[j] > 0);
            //printf("H[%d][%d] = %d, xhat[%d] = %d, accum = %d\n", i, j, H[i][j], j, xhat[j], accum);
            
        }
        accum = accum % 2;
        if (accum != 0)
        {   
            holds = false;
            break;
        }
    }
    if(holds){
        printf("%d: holds\n", cntr);
        
    }
    else{
        printf("%d: does not hold\n", cntr);
    }
    cntr += 1;


    for (int i_itera = 0; i_itera < 0; i_itera++)
    {   

        // STEP 2: syndrome calculation
        
        // * OK
        bool holds = true;
        for (int i = 0; i <M; i++)
        {
            int accum = 0;
            for (int j = 0; j < N; j++)
            {
                

                accum += H[i][j] * ((xhat[j] + 1) / 2);
                //printf("H[%d][%d] = %d, xhat[%d] = %d, accum = %d\n", i, j, H[i][j], j, xhat[j], accum);
                
            }

            accum = accum % 2;
            if (accum != 0)
            {   
                holds = false;
                break;
            }
            
        }
        if(holds){
            printf("exited after %d iterations with flag %d\n", i_itera, flag);
            break;
        }
        else if(i_itera == itera-1){
            printf("DID NOT CONVERGE after %d iterations with flag %d \n", i_itera, flag);

        }


        // INVERSION FUNCTION
        int lowestIdx = 0;
        double lowestValue;
        for (int k = 0; k < N; k++)
        {
            inverFunc[k] = xhat[k] * LLR[k];
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
                    inverFunc[k] += accum;
                }
            }

            // find the indices that are smaller than delta
            if (inverFunc[k] < delta)
            {
                Idx[k] = 1;
            }
            else
            {
                Idx[k] = 0;
            }

            // find the smallest value in inverFunc
            if (k == 0)
            {
                lowestValue = inverFunc[k];
            }
            else if(inverFunc[k] < lowestValue)
            {
                lowestValue = inverFunc[k];
                lowestIdx = k;
            }
        }

        // STEP 3: objective function calculation
        if (flag == 0) 
        {       
            double objectFunc1 = 0;
            for (int j = 0; j < N; j++)
            {
                objectFunc1 += xhat[j] * LLR[j];
            }
            for (int i = 0; i < M; i++)
            {
                int accum = 1;
                for (int j = 0; j < N; j++)
                {   
                    if(H[i][j]==1)
                        accum *= xhat[j];
                }
                objectFunc1 += accum;
            }

            if(hasPrinted)
                printf("objectFunc1 = %f\n", objectFunc1);

            // flip the bits
            for (int k = 0; k < N; k++)
            {
                if (Idx[k] == 1)
                {  
                    xhat[k] = -xhat[k];
                }
            }

            // object function calculation
            double objectFunc2 = 0;
            for (int j = 0; j < N; j++)
            {
                objectFunc2 += xhat[j] * LLR[j];
            }
            for (int i = 0; i < M; i++)
            {
                int accum = 1;
                for (int j = 0; j < N; j++)
                {   
                    if(H[i][j]==1)
                        accum *= xhat[j];
                }
                objectFunc2 += accum;
            }

            if (objectFunc2 < objectFunc1)
            {
                flag = 1;
            }
        }
        else
        {   
            // flip the bit
            xhat[lowestIdx] = -xhat[lowestIdx];
        }
    }

    hasPrinted = 0;

    // output the corrected bits to the outside
    for (int i = 0; i < nb*Z; i++)
    {
        debits[i] = (xhat[i] > 0);
    }

    // Free each sub-array
    /*for (int i = 0; i < M; i++)
    {

        mxFree(H[i]);
    }*/



    // Finally, free the top-level array
    //mxFree(H);
    mxFree(LLR);
    mxFree(xhat);
    mxFree(inverFunc);
    mxFree(Idx);

    return 1;
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

/*
    // mb=8 nb=40 Z=16 WinSize=5 w=2


    // M=5*8*16=640 N=6*40*16=3840
    double delta=0.1;
    int     M=WinSize*mb*Z;
    int     N=(WinSize+w-1)*nb*Z;


    // Flatten LLRIn
    // create 3840x1 array of LLR's (y)
    double* LLR=(double*)mxMalloc(N*sizeof(double));
    for(int i=0; i<N; i++){ //N
        LLR[i]=LLRIn[i/Z][i%Z];
    }



    // correct window shape, we need to recontruct the windowed H matrix
    // H is a 640 x 3840 matrix
    int**    H=(int**)mxCalloc(M, sizeof(int*));

    for(int i = 0; i < M; i++) {
        H[i] = (int*)mxCalloc(N, sizeof(int));
    }

    for(int i=0; i<M/Z; i++){
        for(int j=0; j<N/Z; j++){

            int a=Window[i*N/Z+j];
            for(int z=0; z<Z && a!=-1; z++){
                H[z+Z*i][(z+a)%Z+Z*j]=1;
            }

        }
    }

    // inverFunc: 3840x1 array
    double* inverFunc=(double*)mxMalloc(sizeof(double)*N);

    // xhat: 3840x1 array of hard decision bits
    int*   xhat=mxMalloc(sizeof(int)*N);

    // syndrome: 640x1 array
    int*   syndrome=mxMalloc(sizeof(int)*M);

    // hard decision
    for (int i=0;i<N;i++){
        xhat[i]=(LLRIn[i]>0);
    }


    // syndrome calculation: H*xhat=scalar
    for (int i=0;i<M;i++){
        syndrome[i]=0;
        for (int j=0;j<N;j++){
            syndrome[i]=(syndrome[i]+xhat[j]*H[i][j])%2;
        }
    }

    // object function: scalar
    double objectFunc = 0;

    // object function calculation
    for (int i=0;i<N;i++){
        objectFunc+=(1-2*xhat[i])*LLR[i];
    }
    for (int i=0;i<M;i++){
        objectFunc+=1-2*syndrome[i];
    }


    int modeFlag=1;

    int* Idx=mxMalloc(sizeof(int)*N);

    for(int i_itera=0; i_itera<itera; i_itera++){

        for(int i_cols=0; i_cols<N; i_cols++){
            inverFunc[i_cols]=(1-2*xhat[i_cols])*LLR[i_cols];
            for(int i_rows=0; i_rows<M; i_rows++){
               inverFunc[i_cols]+=H[i_rows][i_cols]*(1-2*syndrome[i_rows]);
            }
        }

        //inverFunc -> values, Idx
        double* values=inverFunc;


        // fill Idx
        for(int i=0; i<N; i++){
            Idx[i]=i;
        }

        //sort inverse function (ascend)
        for(int i=0; i<N; i++){
            for(int j=i+1; j<N; j++){
                if(values[i]>values[j]){
                    double temp=values[i];
                    values[i]=values[j];
                    values[j]=temp;

                    int tempIdx=Idx[i];
                    Idx[i]=Idx[j];
                    Idx[j]=tempIdx;
                }
            }
        }

        int IdxNb;

        // multi bit flipping
        if(modeFlag==1){

            // find the amount of indices that are smaller than delta
            IdxNb=0;
            for(int i=0; i<N; i++){
                if(values[i]>=delta) break;
                IdxNb++;
            }

            if(IdxNb==0){
                IdxNb=1;
                modeFlag=0;
            }

            // flip the bits
            for(int i=0; i<IdxNb; i++){
                xhat[Idx[i]]=1-xhat[Idx[i]];
            }

            // syndrome calculation
            for(int i=0; i<M; i++){
                for (int j=0;j<IdxNb;j++){
                    syndrome[i]+=H[i][Idx[j]];
                }
                syndrome[i]=syndrome[i]%2;
            }

        }
        // single bit flipping
        else{

            IdxNb=1;
            xhat[Idx[0]]=1-xhat[Idx[0]]; // flip the bit

            // syndrome calculation
            for(int i=0; i<M; i++){
                syndrome[i]+=H[i][Idx[0]];
                syndrome[i]=syndrome[i]%2;
            }
        }


        // object function calculation

        double objectFuncLast=objectFunc;
        objectFunc=0;
        for (int i=0;i<N;i++){
            objectFunc+=(1-2*xhat[i])*LLR[i];
        }
        for (int i=0;i<M;i++){
            objectFunc+=1-2*syndrome[i];
        }

        if(objectFunc<objectFuncLast){
            modeFlag=0;
        }


        int sum_syndrome=0;
        for(int i=0; i<M; i++){
            sum_syndrome+=syndrome[i];
        }

        if(sum_syndrome==0){
            break;
        }



    }



    // output the corrected bits to the outside
    for (int i=0;i<N;i++){
        debits[i]=(double)xhat[i];
    }

    // Free each sub-array
    for(int i = 0; i < M; i++) {

        mxFree(H[i]);
    }

    // Finally, free the top-level array
    mxFree(H);
    mxFree(LLR);
    mxFree(xhat);
    mxFree(syndrome);
    mxFree(inverFunc);
    mxFree(Idx);


    return 1;





*/