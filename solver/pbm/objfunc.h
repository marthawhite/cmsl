
#ifndef OBJFUNC_H
#define OBJFUNC_H

#ifdef _MSC_VER     // if using Visual Studio compiler
    #define objfunc OBJFUNC
#else
    #define objfunc objfunc_
#endif

#ifndef mex_h
    class mxArray;
#endif

extern "C" 
{
    void objfunc(int *n, double *x, double *f, double *g);
}

extern mxArray * func_arg;

#endif
