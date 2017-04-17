#include <iostream>
#include <fstream>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/registration/gicp.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/features/normal_3d.h>
#include <Eigen/SVD>
#include "mex.h"

static const bool DEBUG = true; 
static const double GLOBAL_CLOUD_EPSILON = 1e-6; 
static std::ofstream Debug_Out_Stream;

//static function declearion:
static void verifyInput(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static void getPointCloud(const double* in_mtx, int num_of_points, pcl::PointCloud<pcl::PointXYZ>::Ptr p_point_cloud);
static void savePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr p_point_cloud, const char* p_file_name);
static void showPointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr p_point_cloud, const char* p_cloud_name);
static void computeCovariances(const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
	pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::MatricesVectorPtr cloud_covariances, double gicp_epsilon);
static void computeCovariancesPt2Pl(const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
	pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::MatricesVectorPtr cloud_covariances);
static void computeCovariancesZeros(const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
	pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::MatricesVectorPtr cloud_covariances);
static void computeCovariancesOnes(const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
	pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::MatricesVectorPtr cloud_covariances);
static void calcMatrixPsudoInverse(Eigen::Matrix3d &mat_in, Eigen::Matrix3d &mat_out);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Debug_Out_Stream.open("debug.txt");

	verifyInput(nlhs, plhs, nrhs, prhs);

	double epsilon = mxGetScalar(prhs[2]);
	double d_max = mxGetScalar(prhs[3]);
	
	pcl::PointCloud<pcl::PointXYZ>::Ptr p_global_cloud(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr p_local_cloud(new pcl::PointCloud<pcl::PointXYZ>);

	getPointCloud(mxGetPr(prhs[0]), mxGetM(prhs[0]), p_local_cloud);
	getPointCloud(mxGetPr(prhs[1]), mxGetM(prhs[1]), p_global_cloud);

	if (DEBUG == true)
	{
		//savePointCloud(p_local_cloud, "Local_Cloud.txt");
		//savePointCloud(p_global_cloud, "Global_Cloud.txt");
		//showPointCloud(p_global_cloud, "Global Cloud");
		//showPointCloud(p_local_cloud, "Local Cloud");
	}
	
	pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> gicp;
	gicp.setMaximumOptimizerIterations(100);
	gicp.setRotationEpsilon(2e-5);
	gicp.setTransformationEpsilon(1e-5);
	gicp.setInputCloud(p_local_cloud);
	gicp.setInputTarget(p_global_cloud);
	gicp.setMaxCorrespondenceDistance(d_max);

	//Compute covariance
	pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::MatricesVectorPtr p_global_cov(new pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::MatricesVector);
	pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::MatricesVectorPtr p_local_cov(new pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::MatricesVector);

	computeCovariances(p_global_cloud, p_global_cov, GLOBAL_CLOUD_EPSILON);
	gicp.setTargetCovariances(p_global_cov);

	computeCovariances(p_local_cloud, p_local_cov, std::max(epsilon, GLOBAL_CLOUD_EPSILON));
	gicp.setSourceCovariances(p_local_cov);

	pcl::PointCloud<pcl::PointXYZ> final;
	pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::Matrix4 transform_matrix; 
	gicp.align(final);
	transform_matrix = gicp.getFinalTransformation();

	if (DEBUG == true)
	{
		mexPrintf("GICP algorithm converged? %s\n", gicp.hasConverged() ? "true" : "false");
		mexPrintf("GICP algorithm fitness score: %lf\n", gicp.getFitnessScore());
	}

	// Prepare output arguments
	//Currently I set number of iterations to always be 0, 
	//because getting the number of iterations is not supported in PCL registration
	plhs[1] = mxCreateDoubleScalar(0);
	double* iter_out = mxGetPr(plhs[1]);
	*iter_out = 0;

	plhs[0] = mxCreateDoubleMatrix(4, 4, mxREAL);
	double* p_out_mtx4 = mxGetPr(plhs[0]); 

	//Copy transform matrix
	for (int i = 0; i<4; i++)
	{
		for (int k = 0; k < 4; k++)
		{
			p_out_mtx4[i * 4 + k] = transform_matrix(k, i); 
		}
	}

	Debug_Out_Stream.close();
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

	p_point_cloud->is_dense = true;
}

static void savePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr p_point_cloud, const char* p_file_name)
{
	std::ofstream out_file; 
	out_file.open(p_file_name);

	for (size_t i = 0; i < p_point_cloud->size(); i++)
	{
		out_file << "(" << p_point_cloud->points[i].x << "," <<
			p_point_cloud->points[i].y << "," << p_point_cloud->points[i].z
			<< ")\n"; 
	}
	out_file.close();
}

static void showPointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr p_point_cloud, const char* p_cloud_name)
{
	pcl::visualization::CloudViewer viewer(p_cloud_name);
	viewer.showCloud(p_point_cloud); 
	while (viewer.wasStopped() == false); 
}

static void computeCovariances(const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,  
	pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::MatricesVectorPtr cloud_covariances, double gicp_epsilon)
{
	static const int k_correspondences_ = 20; 

	pcl::search::KdTree<pcl::PointXYZ> kdtree;
	kdtree.setInputCloud(cloud);

	if (k_correspondences_ > int(cloud->size()))
	{
		PCL_ERROR("[pcl::GeneralizedIterativeClosestPoint::computeCovariances] Number or points in cloud (%lu) is less than k_correspondences_ (%lu)!\n", cloud->size(), k_correspondences_);
		return;
	}

	Eigen::Vector3d mean;
	std::vector<int> nn_indecies; nn_indecies.reserve(k_correspondences_);
	std::vector<float> nn_dist_sq; nn_dist_sq.reserve(k_correspondences_);

	if (cloud_covariances->size() < cloud->size())
		cloud_covariances->resize(cloud->size());

	pcl::PointCloud<pcl::PointXYZ>::const_iterator points_iterator = cloud->begin();
	pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::MatricesVector::iterator matrices_iterator = cloud_covariances->begin();
	for (;
		points_iterator != cloud->end();
		++points_iterator, ++matrices_iterator)
	{
		const pcl::PointXYZ &query_point = *points_iterator;
		Eigen::Matrix3d &cov = *matrices_iterator;
		// Zero out the cov and mean
		cov.setZero();
		mean.setZero();

		// Search for the K nearest neighbours
		kdtree.nearestKSearch(query_point, k_correspondences_, nn_indecies, nn_dist_sq);

		// Find the covariance matrix
		for (int j = 0; j < k_correspondences_; j++) {
			const pcl::PointXYZ &pt = (*cloud)[nn_indecies[j]];

			mean[0] += pt.x;
			mean[1] += pt.y;
			mean[2] += pt.z;

			cov(0, 0) += pt.x*pt.x;

			cov(1, 0) += pt.y*pt.x;
			cov(1, 1) += pt.y*pt.y;

			cov(2, 0) += pt.z*pt.x;
			cov(2, 1) += pt.z*pt.y;
			cov(2, 2) += pt.z*pt.z;
		}

		mean /= static_cast<double> (k_correspondences_);
		// Get the actual covariance
		for (int k = 0; k < 3; k++)
			for (int l = 0; l <= k; l++)
			{
				cov(k, l) /= static_cast<double> (k_correspondences_);
				cov(k, l) -= mean[k] * mean[l];
				cov(l, k) = cov(k, l);
			}

		// Compute the SVD (covariance matrix is symmetric so U = V')
		Eigen::JacobiSVD<Eigen::Matrix3d> svd(cov, Eigen::ComputeFullU);
		cov.setZero();
		Eigen::Matrix3d U = svd.matrixU();
		// Reconstitute the covariance matrix with modified singular values using the column     // vectors in V.
		for (int k = 0; k < 3; k++) {
			Eigen::Vector3d col = U.col(k);
			double v = 1.; // biggest 2 singular values replaced by 1
			if (k == 2)   // smallest singular value replaced by gicp_epsilon
				v = gicp_epsilon;
			cov += v * col * col.transpose();
		}
	}
}

// This function probably crashes Matlab.
static void computeCovariancesPt2Pl(const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
	pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::MatricesVectorPtr cloud_covariances)
{
	static const int k_correspondences_ = 20;

	if (cloud_covariances->size() < cloud->size())
		cloud_covariances->resize(cloud->size());

	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
	ne.setInputCloud(cloud);

	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>());
	ne.setSearchMethod(tree);

	pcl::PointCloud<pcl::Normal>::Ptr cloud_normals(new pcl::PointCloud<pcl::Normal>);

	ne.setKSearch(k_correspondences_);
	ne.compute(*cloud_normals);

	pcl::PointCloud<pcl::Normal>::const_iterator points_iterator = cloud_normals->begin();
	pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::MatricesVector::iterator matrices_iterator = cloud_covariances->begin();
	for (;
		points_iterator != cloud_normals->end();
		++points_iterator, ++matrices_iterator)
	{
		const pcl::Normal &query_point = *points_iterator;
		Eigen::Matrix3d &cov = *matrices_iterator;
		Eigen::Matrix3d cov_inverse; 
		double dot_product = 0; 

		for (int i=0; i<3; i++)
		{
			dot_product += query_point.normal[i] * query_point.normal[i];
			for (int k=0; k<3; k++)
			{
				cov_inverse(i, k) = query_point.normal[i] * query_point.normal[k];
			}
		}
		cov_inverse /= dot_product;
		
		calcMatrixPsudoInverse(cov_inverse, cov);
	}
}

static void computeCovariancesZeros(const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
	pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::MatricesVectorPtr cloud_covariances)
{
	if (cloud_covariances->size() < cloud->size())
		cloud_covariances->resize(cloud->size());

	pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::MatricesVector::iterator matrices_iterator = cloud_covariances->begin();
	for (;
		matrices_iterator != cloud_covariances->end();
		++matrices_iterator)
	{
		Eigen::Matrix3d &cov = *matrices_iterator;
		cov = Eigen::Matrix3d::Zero();
	}
}

static void computeCovariancesOnes(const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
	pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::MatricesVectorPtr cloud_covariances)
{
	if (cloud_covariances->size() < cloud->size())
		cloud_covariances->resize(cloud->size());

	pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>::MatricesVector::iterator matrices_iterator = cloud_covariances->begin();
	for (;
		matrices_iterator != cloud_covariances->end();
		++matrices_iterator)
	{
		Eigen::Matrix3d &cov = *matrices_iterator;
		cov.setIdentity();
	}
}

static void calcMatrixPsudoInverse(Eigen::Matrix3d &mat_in, Eigen::Matrix3d &mat_out)
{
	static const double  P_INV_TOLERANCE = 1.e-6;

	Eigen::JacobiSVD<Eigen::Matrix3d> svd(mat_in, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::JacobiSVD<Eigen::Matrix3d>::SingularValuesType singularValues_inv, singularValues;
	singularValues = svd.singularValues();
	singularValues_inv.resizeLike(singularValues);
 
	for (long i = 0; i<mat_in.cols(); ++i) {
		singularValues_inv(i) = 0;
		if (singularValues(i) > P_INV_TOLERANCE)
		{
			singularValues_inv(i) = 1.0 / singularValues(i);
		}
	}
	mat_out = (svd.matrixV()*singularValues_inv.asDiagonal()*svd.matrixU().transpose());
}