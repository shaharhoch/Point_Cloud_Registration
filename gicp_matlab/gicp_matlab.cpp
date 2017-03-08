/*
 * gpic_matlab.c - Wrapper for gicp algorithm to run from Matlab
 *
 * Runs the icp registration algorithm 
 *
 * The calling syntax is:
 *
 *		reg_mtx = arrayProduct(local_mtx, global_mtx, epsilon, d_max)
 *
 * This is a MEX file for MATLAB.
*/

#include "mex.h"
#include "gicp.h"
#include <string.h>

using namespace dgc::gicp;

//static function declearion:
static void verifyInput(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static void getPointCloud(const double* in_mtx, int num_of_points, GICPPointSet *out_set);


/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{  
    verifyInput(nlhs, plhs, nrhs, prhs);
    
    double epsilon = mxGetScalar(prhs[2]);
    double d_max = mxGetScalar(prhs[3]);
    
    // Initialzie point clouds
    GICPPointSet local_cloud, global_cloud;
    getPointCloud(mxGetPr(prhs[0]), mxGetM(prhs[0]), &local_cloud);
    getPointCloud(mxGetPr(prhs[1]), mxGetM(prhs[1]), &global_cloud);
    
    local_cloud.SetGICPEpsilon(epsilon);
    global_cloud.SetGICPEpsilon(epsilon);  
    local_cloud.BuildKDTree();
    local_cloud.ComputeMatrices();
    global_cloud.BuildKDTree();
    global_cloud.ComputeMatrices();
    
    global_cloud.SetDebug(false);
    global_cloud.SetMaxIterationInner(8);
    global_cloud.SetMaxIteration(100);
    
    // Initialzie transformation to identity transformations     
    dgc_transform_t t_base, t1; 
    dgc_transform_identity(t_base);
    dgc_transform_identity(t1);
    
    int iterations = global_cloud.AlignScan(&local_cloud, t_base, t1, d_max);
	mexPrintf("GICP algorithm converged after %d iterations\n", iterations);
    
    // Prepate output argument
    plhs[0] = mxCreateDoubleMatrix(4,4,mxREAL);
    memcpy(mxGetPr(plhs[0]), &t1, 4*4*sizeof(double));
    
}

// This function verifies the correctness of the inputs from matlab
static void verifyInput(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Verify number of inputs
    if(nrhs != 4) 
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Four inputs required.");
    }
    
    // Verify input 1's type
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) 
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble", "First input matrix must be type double.");
    }
    
    if(mxGetN(prhs[0]) != 3) 
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:InvalidMtxSize", "First input matrix has to have 3 columns");
    }
    
    // Verify input 2's type
    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) 
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble", "Second input matrix must be type double.");
    }
    
    if(mxGetN(prhs[1]) != 3) 
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:InvalidMtxSize", "Second input matrix has to have 3 columns");
    }
    
    // Verify input 3's type
    if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1 ) 
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar", "Input epsilon must be a scalar.");
    }
    
    // Verify input 3's type
    if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxGetNumberOfElements(prhs[3]) != 1 ) 
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar", "Input d_max must be a scalar.");
    }
}


/*
    This function converts a matrix of size (num_of_points X 3) to a GICPPointSet.  
    
    in_mtx in a 2-D array of size (num_of_points X 3). 
    output: out_set
*/
static void getPointCloud(const double* in_mtx, int num_of_points, GICPPointSet *out_set)
{
    int ind; 
    for(ind=0; ind<num_of_points; ind++)
    {
        GICPPoint pt;
        
        pt.range = -1; 
        
        int k, l; 
        for(k = 0; k < 3; k++) 
        {
            for(l = 0; l < 3; l++) 
            {
                pt.C[k][l] = (k == l)?1:0;
            }
        }
        
        pt.x = *(in_mtx + (ind+0*num_of_points));
        pt.y = *(in_mtx + (ind+1*num_of_points));
        pt.z = *(in_mtx + (ind+2*num_of_points));
        
        out_set->AppendPoint(pt);
    }
}
