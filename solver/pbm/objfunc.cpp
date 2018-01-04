
#include "mex.h"
#include "objfunc.h"

#include <string.h>

mxArray * func_arg = NULL;

void objfunc(int *n, double *x, double *f, double *g)
{
  int i;
  double *arr;

  mxArray *lhs[2];
  mxArray *rhs[2];

  rhs[0] = func_arg;
  rhs[1] = mxCreateDoubleMatrix(*n, 1, mxREAL);
  
  arr = mxGetPr(rhs[1]);
  for(i = 0; i < *n; i++)
    arr[i] = x[i];
  
  mexCallMATLAB(2, lhs, 2, rhs, "feval");
  
  *f = (double)(mxGetScalar(lhs[0]));
  arr = mxGetPr(lhs[1]);
  if(g != NULL)
  {
    for(i = 0; i < *n; i++)
      g[i] = arr[i];
  }

  mxDestroyArray(rhs[1]);
  mxDestroyArray(lhs[0]);
  mxDestroyArray(lhs[1]);
}

