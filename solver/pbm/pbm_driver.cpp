
#include "mex.h"
#include "objfunc.h"
#include "float.h"

#include <string.h>

#ifdef _MSC_VER     // if using Visual Studio compiler
    #define mpbngc MPBNGC
#else
    #define mpbngc mpbngc_
#endif

extern "C" int mpbngc(int *n, double *x, int *IX, double *BL, double *BU,
           int *M, int *MG, int *MC, int *IC, double *CL, double *CU, double *CG,
           double *F, double *lsParam, int *LMAX, double *GAM, double *EPS, double *FEAS,
           int *JMAX, int *NITER, int *NFASG, int *NOUT, int *verbose, int *IERR,
           int *IWORK, int *LIWORK, double *WORK, int *LWORK, int *IUSER, double *USER);

inline int round(double x) { 
    return static_cast<int>(x + (x > 0.0 ? +0.5 : -0.5)); 
}


static void printDoubleVector(double *x, int n)
{
  int i;
  printf("[ ");
  for(i = 0; i < n; i++)
  {
    if(x[i] >= 0.01 || x[i] < DBL_EPSILON)
      printf("%.2f ", x[i]);
    else
      printf("%.2e ", x[i]);
  }
  printf("]");
}

static void printIntVector(int *x, int n)
{
  int i;
  printf("[ ");
  for(i = 0; i < n; i++)
    printf("%d ", x[i]);
  printf("]");
}

void checkScalar(const mxArray *p, const char *name)
{
  char *errorMsg;
  errorMsg = new char [1000];
  int error = 0;
  
  errorMsg[0] = '\0';
  strcat(errorMsg, name);
  
  if(!mxIsDouble(p) || !(mxGetNumberOfElements(p) == 1))
  {
    strcat(errorMsg, " must be a scalar.");
    error = 1;
  }
  
  if(error)
    mexErrMsgTxt(errorMsg);
  delete []errorMsg;
}

void checkIntScalar(const mxArray *p, const char *name)
{
  char errorMsg[1000];
  int i, error = 0;
  double d;
  
  
  errorMsg[0] = '\0';
  strcat(errorMsg, name);
  
  if(!mxIsDouble(p) || !(mxGetNumberOfElements(p) == 1))
  {
    strcat(errorMsg, " must be a scalar.");
    error = 1;
  }
  
  d = mxGetScalar(p);
  i = (int)mxGetScalar(p);
  if(d != i)
  {
    strcat(errorMsg, " must be an integer.");
    error = 1;
  }
  
  if(error)
    mexErrMsgTxt(errorMsg);
}

void checkArray(const mxArray *p, const char *name, int n)
{
  char errorMsg[1000];
  int error = 0;
  int ne;
  
  errorMsg[0] = '\0';
  strcat(errorMsg, name);
  
  if(!mxIsDouble(p))
  {
    strcat(errorMsg, " must be a scalar array.");
    error = 1;
  }
  
  ne = mxGetNumberOfElements(p);
  if(ne != n)
  {
    sprintf(errorMsg, "The length of %s is %d != %d.", name, ne, n);
    error = 1;
  }
  
  if(error)
    mexErrMsgTxt(errorMsg);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char *funcname;
    int i, j;
    bool flag1, flag2;
    char argBuf[2047];
    bool debug = false;

    int n;
    int *ix;    
    int m;
    int ic[1];
    int mc;
    int mg;
    int mp1;
    int maxIter;
    int maxLsIter;
    int maxBdl;
    int maxFnCall;
    int nout;
    int verbose;
    int ierr;
    int liwork;
    int *iwork;
    int lwork;
    int iuser[1];
    
    double *x;
    double *bl, *bu;
    double cl[1], cu[1], cg[1];
    double opt_f;
    double lsParam;
    double tolFun;
    double tolCon;
    double *gam;
    double *work;
    double user[1];
    double *ptr1, *ptr2;

    
    //////////////////////////////////////////////////////////////////////////
    // Default values
    //////////////////////////////////////////////////////////////////////////
    
    m = 1;
    mg = 0;
    mc = 0;
    ic[0] = 0;
    nout = 0;
    maxIter = 1000;
    maxLsIter = 50;
    maxBdl = 100;
    maxFnCall = 1000;
    verbose = -1;

    cl[0] = cu[0] = cg[0] = 0.0;
    lsParam = 0.01;
    tolFun = 1e-5;
    tolCon = 1e-5;

    if (nlhs < 1)
        mexErrMsgTxt("pbm requires at least one input\n");


    //////////////////////////////////////////////////////////////////////////
    // Read in the command line arguments
    //////////////////////////////////////////////////////////////////////////

    func_arg = (mxArray *) prhs[0];

    n = mxGetM(prhs[1]) * mxGetN(prhs[1]);

    x = new double [n];
    memcpy(x, mxGetPr(prhs[1]), sizeof(double) * n);


    //////////////////////////////////////////////////////////////////////////
    // Some optimization parameters to set
    //////////////////////////////////////////////////////////////////////////

    ix = new int [n];
    bl = new double [n];
    bu = new double [n];

    ptr1 = mxGetPr(prhs[2]);    // lower bound
    ptr2 = mxGetPr(prhs[3]);    // upper bound
    for (i = 0; i < n; i ++)    {
        flag1 = mxIsFinite(ptr1[i]);
        flag2 = mxIsFinite(ptr2[i]);
        if (flag1 && flag2)        {
            ix[i] = 3;  bl[i] = ptr1[i];  bu[i] = ptr2[i];
            if (x[i] > bu[i] || x[i] < bl[i])   { 
                sprintf(argBuf, "Initial x[%d] violates box constraints.\n", i);
                mexErrMsgTxt(argBuf);
            }
        }
        else if (!flag1 && flag2)        {
            ix[i] = 2;  bl[i] = 0.0;  bu[i] = ptr2[i];
            if (x[i] > bu[i]){ 
                sprintf(argBuf, "Initial x[%d] violates box constraints.\n", i);
                mexErrMsgTxt(argBuf);
            }
        }
        else if (flag1 && !flag2)        {
            ix[i] = 1;  bl[i] = ptr1[i];  bu[i] = 0.0;
            if (x[i] < bl[i]){ 
                sprintf(argBuf, "Initial x[%d] violates box constraints.\n", i);
                mexErrMsgTxt(argBuf);
            }
        }
        else        {
            ix[i] = 0;  bl[i] = 0.0;  bu[i] = 0.0;
        }
    }

    if (nrhs > 4)
    {
        for(i = 0; i < mxGetNumberOfFields(prhs[4]); i++) 
        {
           const char * fname = mxGetFieldNameByNumber(prhs[4], i); // field name
           mxArray *tmp = mxGetFieldByNumber(prhs[4], 0, i);        // field value

           if (strcmp(fname, "maxIter") == 0)   {
                maxIter = round(*((double*)mxGetData(tmp)));
                if (maxIter < 0)
                    mexErrMsgTxt("maxIter must be >= 0.\n");
                if(debug) mexPrintf("Got maxIter = %d\n", maxIter);
           }
           else if (strcmp(fname, "maxLsIter") == 0)   {
               maxLsIter = round(*((double*)mxGetData(tmp)));
               if (maxLsIter < 0)
                   mexErrMsgTxt("maxLsIter must be >= 0.\n");
               if(debug) mexPrintf("Got maxLsIter = %d\n", maxLsIter);
           }
           else if (strcmp(fname, "maxBdl") == 0)   {
               maxBdl = round(*((double*)mxGetData(tmp)));
               if (maxBdl < 0)
                   mexErrMsgTxt("maxBdl must be >= 0.\n");
              if(debug) mexPrintf("Got maxBdl = %d\n", maxBdl);
           }
           else if (strcmp(fname, "maxFnCall") == 0)   {
               maxFnCall = round(*((double*)mxGetData(tmp)));
               if (maxFnCall < 0)
                   mexErrMsgTxt("maxFnCall must be >= 0.\n");
              if(debug) mexPrintf("Got maxFnCall = %d\n", maxFnCall);
           }
           else if (strcmp(fname, "verbose") == 0)   {
               verbose = round(*((double*)mxGetData(tmp)));
               if (verbose > 0)  mexErrMsgTxt("verbose must be <= 0.\n");
               if(debug) mexPrintf("Got verbose = %d\n", verbose);
           }
           else if (strcmp(fname, "lsParam") == 0)   {
               lsParam = *((double*)mxGetData(tmp));
               if (lsParam <= 0 || lsParam >= 0.5) mexErrMsgTxt("lsParam must be in (0, 0.5)\n");
               if(debug) mexPrintf("Got lsParam = %f\n", lsParam);
           }
           else if (strcmp(fname, "tolFun") == 0)   {
               tolFun = *((double*)mxGetData(tmp));
               if (tolFun < 1e-8) mexErrMsgTxt("tolFun must be at least 1e-8\n");
               if(debug) mexPrintf("Got tolFun = %f\n", tolFun);
           }
           else if (strcmp(fname, "tolCon") == 0)   {
               tolCon = *((double*)mxGetData(tmp));
               if (tolCon < 1e-8) mexErrMsgTxt("tolCon must be at least 1e-8\n");
               if(debug) mexPrintf("Got tolCon = %f\n", tolCon);
           }
           else {
               if(mxIsChar(tmp))    {
                   int buflen = (mxGetM(tmp) * mxGetN(tmp)) + 1;
                   mxGetString(tmp, argBuf, buflen);
                   if(debug) mexPrintf("Unrecognized param: %s = %s\n", fname, argBuf);
               }
               else     {
                   if(debug) mexPrintf("Unrecognized param: %s = %f\n", fname, 
                                        *((double*)mxGetData(tmp)));
               }
           }
        }
    }


    //////////////////////////////////////////////////////////////////////////
    // Allocate work space using the parameters
    //////////////////////////////////////////////////////////////////////////

    mp1 = mg > 0 ? (m + 1) : m;

    liwork = 2 * (mp1 * (maxBdl+1) + mc + n);
    iwork = new int [liwork];

    lwork = n * (n+2*mp1*maxBdl+2*mp1+2*mc+2*mg+2*m+31)/2
                +mp1*(6*maxBdl+10)+maxBdl+5*mc+mg+m+18;
    work = new double [lwork];
    gam = new double [m+1];
    for (j = 0; j <= m; j ++)
        gam[j] = 0.0;

    // Call the solver
    mpbngc(&n, x, ix, bl, bu, &m, &mg, &mc, ic, cl, cu, cg,
            &opt_f, &lsParam, &maxLsIter, gam, &tolFun, &tolCon,
            &maxBdl, &maxIter, &maxFnCall, &nout, &verbose, &ierr,
            iwork, &liwork, work, &lwork, iuser, user);

    // Pass back the solution
    if(nlhs >= 1)
    {
        plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
        memcpy(mxGetPr(plhs[0]), x, n*sizeof(double));
    }
    if(nlhs >= 2)
    {
        plhs[1] = mxCreateDoubleScalar(0.0);
        *mxGetPr(plhs[1]) = opt_f;
    }
    if(nlhs >= 3)
    {
        plhs[2] = mxCreateDoubleScalar(0.0);
        *mxGetPr(plhs[2]) = (double) maxIter;
    }
    if(nlhs >= 4)
    {
        plhs[3] = mxCreateDoubleScalar(0.0);
        *mxGetPr(plhs[3]) = (double) maxFnCall;
    }
    if(nlhs >= 5)
    {
        plhs[4] = mxCreateDoubleScalar(0.0);
        *mxGetPr(plhs[4]) = (double) ierr;
    }

    delete [] x;
    delete [] work;
    delete [] iwork;
    delete [] gam;    
    delete [] bl;
    delete [] bu;
    delete [] ix;
}
