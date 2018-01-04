#include "mex.h"
#include "math.h"


inline double _max(register double x, register double y)
{
    return x >= y ? x : y;
}

inline double _min(register double x, register double y)
{
    return x <= y ? x : y;
}


// Function definitions. 
// downsample(img, mx, my)
// Will return a vector;
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[]) 
{
    double *ori_vec, *down_vec;

    double **ori_mat, *p, value;
    int ori_sx, ori_sy, new_sx, new_sy;
    int step_x, step_y;
    int ori_x_from, ori_y_from, ori_x_to, ori_y_to, num_y;
    int new_px, new_py, new_idx;
    int i, j, denom;
    double x_frac, y_frac, step_x_frac, step_y_frac;

    ori_vec = mxGetPr(prhs[0]);

    ori_sx = mxGetM(prhs[0]);
    ori_sy = mxGetN(prhs[0]);

    new_sx = (int)(*mxGetPr(prhs[1]) + 0.5);
    new_sy = (int)(*mxGetPr(prhs[2]) + 0.5);

    plhs[0] = mxCreateDoubleMatrix(new_sx*new_sy, 1, mxREAL);
    down_vec = mxGetPr(plhs[0]);

    ori_mat = new double *[ori_sy];    
    for (i = 0; i < ori_sy; i ++)
    {
      ori_mat[i] = ori_vec;
      ori_vec += ori_sx;
    }

    step_x = ori_sx / new_sx;
    step_y = ori_sy / new_sy;
    step_x_frac = (ori_sx-1) / (double)(new_sx+1);
    step_y_frac = (ori_sy-1) / (double)(new_sy+1);
    new_idx = 0;  ori_y_to = 0;    

//     mexPrintf("ori_sx = %d, ori_sy = %d, new_sx = %d, new_sy = %d, step_x = %d, step_y = %d\n", 
//       ori_sx, ori_sy, new_sx, new_sy, step_x, step_y);

    for (new_py = 0; new_py < new_sy; new_py ++)
    {

      ori_y_from = floor(step_y_frac * new_py);
      ori_y_to = ceil(step_y_frac * (new_py+1));
      
//       mexPrintf("ori_y_from = %d, ori_y_to = %d\n", ori_y_from, ori_y_to);
      for (new_px = 0; new_px < new_sx; new_px ++)
      {        
        value = 0.0;
        denom = 0;
        ori_x_from = floor(step_x_frac * new_px);
        ori_x_to = ceil(step_x_frac * (new_px+1));

        for (i = ori_y_from; i < ori_y_to; i++)
        {
          p = ori_mat[i] + ori_x_from;
          for (j = ori_x_from; j <= ori_x_to; j ++)
          {
            value += *(p++);    
//             mexPrintf("%f\n", *(p-1));
          }
          denom += ori_x_to - ori_x_from + 1;
        }
//         ori_x_from = ori_x_from + step_x;
        
        down_vec[new_idx ++] = value / denom;
      }      
    }
    
    delete [] ori_mat;
}

