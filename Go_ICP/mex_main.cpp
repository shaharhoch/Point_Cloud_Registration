/*
 * mex_main.c - Wrapper for GoICP algorithm to run from Matlab
 *
 * Runs the GoICP registration algorithm 
 *
 * The calling syntax is:
 *
 *		[reg_mtx, iter_out] = GoICP(local_mtx, global_mtx)
 *
 * This is a MEX file for MATLAB.
*/

#include "jly_goicp.h"
#include "ConfigMap.hpp"
#include <mex.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>

//Defines: 
#define CONFIG_FILE_MAME "GoICP_config.cfg"
#define TEXT_SEPERATOR "----------------------------------------------------------------------\n\n\n"

//Function prototypes: 
static void loadPointCloud(const double* in_mtx, int num_of_points, POINT3D **pp_point);
static void verifyInput(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static void goicpConfig(string FName, GoICP & goicp);
static void printGoICP(GoICP & goicp);

//Static variables:
static std::ofstream Debug_Out_Stream;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Debug_Out_Stream.open("debug.txt");
    
	verifyInput(nlhs, plhs, nrhs, prhs);

	POINT3D *p_global, *p_local;
	GoICP goicp;

	goicpConfig(CONFIG_FILE_MAME, goicp);

	loadPointCloud(mxGetPr(prhs[0]), mxGetM(prhs[0]), &p_local);
	loadPointCloud(mxGetPr(prhs[1]), mxGetM(prhs[1]), &p_global);

	goicp.pModel = p_global;
	goicp.Nm = mxGetM(prhs[1]);
	goicp.pData = p_local;
	goicp.Nd = mxGetM(prhs[0]);

	//Build distance transform
	goicp.BuildDT();

	//Run GoICP algorithm
	float error; 
	error = goicp.Register();
	Debug_Out_Stream << "Registration Error: " << error << "\n\n";

	printGoICP(goicp);

	//Prepare output arguments
	//Currently I set number of iterations to always be 0, 
	//because getting the number of iterations is not supported in GoICP registration
	plhs[1] = mxCreateDoubleScalar(0);
	double* iter_out = mxGetPr(plhs[1]);
	*iter_out = 0;

	plhs[0] = mxCreateDoubleMatrix(4, 4, mxREAL);
	double* p_out_mtx4 = mxGetPr(plhs[0]);
	
	//Init matrix
	p_out_mtx4[0 * 4 + 3] = 0.0;
	p_out_mtx4[1 * 4 + 3] = 0.0;
	p_out_mtx4[2 * 4 + 3] = 0.0;
	p_out_mtx4[3 * 4 + 3] = 1.0;

	//Fill in translation
	p_out_mtx4[3 * 4 + 0] = goicp.optT.val[0][0];
	p_out_mtx4[3 * 4 + 1] = goicp.optT.val[1][0];
	p_out_mtx4[3 * 4 + 2] = goicp.optT.val[2][0];

	//Fill in rotation
	p_out_mtx4[0 * 4 + 0] = goicp.optR.val[0][0];
	p_out_mtx4[1 * 4 + 1] = goicp.optR.val[1][1];
	p_out_mtx4[2 * 4 + 2] = goicp.optR.val[2][2];

	p_out_mtx4[0 * 4 + 1] = goicp.optR.val[1][0];
	p_out_mtx4[0 * 4 + 2] = goicp.optR.val[2][0];

	p_out_mtx4[1 * 4 + 0] = goicp.optR.val[0][1];
	p_out_mtx4[1 * 4 + 2] = goicp.optR.val[2][1];

	p_out_mtx4[2 * 4 + 0] = goicp.optR.val[0][2];
	p_out_mtx4[2 * 4 + 1] = goicp.optR.val[1][2];

	Debug_Out_Stream.close();
	delete(p_global);
	delete(p_local);
}

/*
This function converts a matrix of size (num_of_points X 3) to an array of POINT3D elements.
This function also allocates room for that array.

in_mtx in a 2-D array of size (num_of_points X 3).
output: out_set
*/
static void loadPointCloud(const double* in_mtx, int num_of_points, POINT3D **pp_point)
{
	*pp_point = (POINT3D *)malloc(sizeof(POINT3D) * num_of_points);

	int ind;
	for (ind = 0; ind<num_of_points; ind++)
	{
		(*pp_point)[ind].x = *(in_mtx + (ind + 0 * num_of_points));
		(*pp_point)[ind].y = *(in_mtx + (ind + 1 * num_of_points));
		(*pp_point)[ind].z = *(in_mtx + (ind + 2 * num_of_points));
	}
}

// This function verifies the correctness of the inputs from matlab
static void verifyInput(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Verify number of inputs
	if (nrhs != 2)
	{
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Four inputs required.");
	}

	// Verify input 1's type
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
	{
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble", "First input matrix must be type double.");
	}

	if (mxGetN(prhs[0]) != 3)
	{
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:InvalidMtxSize", "First input matrix has to have 3 columns");
	}

	// Verify input 2's type
	if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]))
	{
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble", "Second input matrix must be type double.");
	}

	if (mxGetN(prhs[1]) != 3)
	{
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:InvalidMtxSize", "Second input matrix has to have 3 columns");
	}
}

static void goicpConfig(string FName, GoICP & goicp)
{
	// Open and parse the associated config file
	ConfigMap config(FName.c_str());

	goicp.MSEThresh = config.getF("MSEThresh");
	goicp.initNodeRot.a = config.getF("rotMinX");
	goicp.initNodeRot.b = config.getF("rotMinY");
	goicp.initNodeRot.c = config.getF("rotMinZ");
	goicp.initNodeRot.w = config.getF("rotWidth");
	goicp.initNodeTrans.x = config.getF("transMinX");
	goicp.initNodeTrans.y = config.getF("transMinY");
	goicp.initNodeTrans.z = config.getF("transMinZ");
	goicp.initNodeTrans.w = config.getF("transWidth");
	goicp.trimFraction = config.getF("trimFraction");
	// If < 0.1% trimming specified, do no trimming
	if (goicp.trimFraction < 0.001)
	{
		goicp.doTrim = false;
	}
	goicp.dt.SIZE = config.getI("distTransSize");
	goicp.dt.expandFactor = config.getF("distTransExpandFactor");
}

static void printGoICP(GoICP & goicp)
{	
	Debug_Out_Stream << "GoICP Model Details: \n"; 
	Debug_Out_Stream << "MSEThresh: " << goicp.MSEThresh << "\n";
	Debug_Out_Stream << "initNodeRot.a: " << goicp.initNodeRot.a << "\n";
	Debug_Out_Stream << "initNodeRot.b: " << goicp.initNodeRot.b << "\n";
	Debug_Out_Stream << "initNodeRot.c: " << goicp.initNodeRot.c << "\n";
	Debug_Out_Stream << "initNodeRot.w: " << goicp.initNodeRot.w << "\n";
	Debug_Out_Stream << "initNodeTrans.x: " << goicp.initNodeTrans.x << "\n";
	Debug_Out_Stream << "initNodeTrans.y: " << goicp.initNodeTrans.y << "\n";
	Debug_Out_Stream << "initNodeTrans.z: " << goicp.initNodeTrans.z << "\n";
	Debug_Out_Stream << "initNodeTrans.w: " << goicp.initNodeTrans.w << "\n";
	Debug_Out_Stream << "trimFraction: " << goicp.trimFraction << "\n";
	Debug_Out_Stream << "doTrim: " << goicp.doTrim << "\n";
	Debug_Out_Stream << "dt.SIZE: " << goicp.dt.SIZE << "\n";
	Debug_Out_Stream << "dt.expandFactor: " << goicp.dt.expandFactor << "\n";
	Debug_Out_Stream << TEXT_SEPERATOR;

	Debug_Out_Stream << "Registration Results:\n";
	Debug_Out_Stream << "Translation: \n" << goicp.optT << "\n";
	Debug_Out_Stream << "Rotation: \n" << goicp.optR << "\n";
	Debug_Out_Stream << TEXT_SEPERATOR;

	Debug_Out_Stream << "Model points: \n";
	for (int i = 0; i < goicp.Nm; i++)
	{
		Debug_Out_Stream << goicp.pModel[i].x << ", " 
			<< goicp.pModel[i].y << ", " << goicp.pModel[i].z << "\n"; 
	}
	Debug_Out_Stream << TEXT_SEPERATOR;


	Debug_Out_Stream << "Data points: \n";
	for (int i = 0; i < goicp.Nd; i++)
	{
		Debug_Out_Stream << goicp.pData[i].x << ", "
			<< goicp.pData[i].y << ", " << goicp.pData[i].z << "\n";
	}
	Debug_Out_Stream << TEXT_SEPERATOR;

	Debug_Out_Stream.flush();
}