#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/registration/gicp.h>
#include "mex.h"

//static function declearion:
static void verifyInput(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static void getPointCloud(const double* in_mtx, int num_of_points, pcl::PointCloud<pcl::PointXYZ>::Ptr p_point_cloud);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mexPrintf("Bla\n");
	/*pcl::PointCloud<pcl::PointXYZ>::Ptr global_cloud (new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr local_cloud(new pcl::PointCloud<pcl::PointXYZ>);

	getPointCloud(mxGetPr(prhs[0]), mxGetM(prhs[0]), local_cloud);
	getPointCloud(mxGetPr(prhs[1]), mxGetM(prhs[1]), global_cloud);*/
}

// This function verifies the correctness of the inputs from matlab
static void verifyInput(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Verify number of inputs
	if (nrhs != 4)
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

	// Verify input 3's type
	if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1)
	{
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar", "Input epsilon must be a scalar.");
	}

	// Verify input 3's type
	if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxGetNumberOfElements(prhs[3]) != 1)
	{
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar", "Input d_max must be a scalar.");
	}
}

/*
This function converts a matrix of size (num_of_points X 3) to a GICPPointSet.

in_mtx in a 2-D array of size (num_of_points X 3).
output: out_set
*/
static void getPointCloud(const double* in_mtx, int num_of_points, pcl::PointCloud<pcl::PointXYZ>::Ptr p_point_cloud)
{
	int ind;
	for (ind = 0; ind<num_of_points; ind++)
	{
		pcl::PointXYZ pt;

		pt.x = *(in_mtx + (ind + 0 * num_of_points));
		pt.y = *(in_mtx + (ind + 1 * num_of_points));
		pt.z = *(in_mtx + (ind + 2 * num_of_points));

		p_point_cloud->push_back(pt); 
	}
}

int old(int argc, char** argv)
{
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_in(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_out(new pcl::PointCloud<pcl::PointXYZ>);

	// Fill in the CloudIn data
	cloud_in->width = 5;
	cloud_in->height = 1;
	cloud_in->is_dense = false;
	cloud_in->points.resize(cloud_in->width * cloud_in->height);
	for (size_t i = 0; i < cloud_in->points.size(); ++i)
	{
		cloud_in->points[i].x = 1024 * rand() / (RAND_MAX + 1.0f);
		cloud_in->points[i].y = 1024 * rand() / (RAND_MAX + 1.0f);
		cloud_in->points[i].z = 1024 * rand() / (RAND_MAX + 1.0f);
	}
	std::cout << "Saved " << cloud_in->points.size() << " data points to input:" << std::endl;

	for (size_t i = 0; i < cloud_in->points.size(); ++i)
	{
		std::cout << "    " <<
			cloud_in->points[i].x << " " << cloud_in->points[i].y << " " <<
			cloud_in->points[i].z << std::endl;
	}

	*cloud_out = *cloud_in;
	std::cout << "size:" << cloud_out->points.size() << std::endl;

	for (size_t i = 0; i < cloud_in->points.size(); ++i)
	{
		cloud_out->points[i].x = cloud_in->points[i].x + 0.7f;
	}
	std::cout << "Transformed " << cloud_in->points.size() << " data points:" << std::endl;

	for (size_t i = 0; i < cloud_out->points.size(); ++i)
	{
		std::cout << "    " << cloud_out->points[i].x << " " <<
			cloud_out->points[i].y << " " << cloud_out->points[i].z << std::endl;
	}

	pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
	icp.setInputCloud(cloud_in);
	icp.setInputTarget(cloud_out);
	pcl::PointCloud<pcl::PointXYZ> Final;
	icp.align(Final);

	std::cout << "has converged:" << icp.hasConverged() << " score: " <<
		icp.getFitnessScore() << std::endl;
	std::cout << icp.getFinalTransformation() << std::endl;

	return (0);
}