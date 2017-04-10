close all
clear all
clc 

% Parameters
NUM_OF_CLOUDS = 24;
global_cloud = pcread('Database\Global_Cloud\Global_Cloud.ply'); 
COARS_REG_DB = csvread('Database\Local2Global_Coarse.csv'); 
GROUND_TRUTH_DB = csvread('Database\Local2Global_GT.csv');

GICP_EPSILON = 0.0004;
GICP_D_MAX = 15;
NUM_OF_POINTS_FOR_NORMAL_CALC = 6;

LOCAL_CLOUD_RANGE = [2, 60]; %[min, max]

ICP_EXTRAPOLATE = false;
ICP_MAX_ITERATIONS = 100;
MAX_INLINER_DISTANCE = 7;
MAX_INLINER_RATIO = MAX_INLINER_DISTANCE/100;
% (!!!) Effective value is 100*maxInlierDistance {For example, maxInlierDistance=0.07  --->  effectively 7m}
ICP_TOLERANCE = [0.1,0.1*pi/180];

% ICP Methods
no_icp.name = sprintf('No ICP');
no_icp.errors = []; 
no_icp.function_handle = @(local, global_c) eye(4,4);

icp_point_to_point.name = 'ICP Point to Point' ; 
icp_point_to_point.errors = []; 
icp_point_to_point.function_handle = @(local, global_c) ...
    transpose(getfield(pcregrigid_v3(local,global_c, ...
    'Metric','pointToPoint',...
    'Extrapolate',ICP_EXTRAPOLATE,...
    'MaxIterations',ICP_MAX_ITERATIONS,...
    'InlierRatio',MAX_INLINER_RATIO,... % value not actually used as inlier ratio - it is used as maxInlierDistance!
    'Tolerance',ICP_TOLERANCE), 'T'));

icp_point_to_plane.name = 'ICP Point to Plane' ; 
icp_point_to_plane.errors = []; 
icp_point_to_plane.function_handle = @(local, global_c) ...
    transpose(getfield(pcregrigid_v3(local,global_c, ...
    'Metric','pointToPlane',...
    'Extrapolate',ICP_EXTRAPOLATE,...
    'MaxIterations',ICP_MAX_ITERATIONS,...
    'InlierRatio',MAX_INLINER_RATIO,... % value not actually used as inlier ratio - it is used as maxInlierDistance!
    'Tolerance',ICP_TOLERANCE), 'T'));

gicp.name = 'Generalized ICP' ; 
gicp.errors = []; 
gicp.function_handle = @(local, global_c) transpose(gicp_alg(double(getfield(local, 'Location')),...
    double(getfield(global_c, 'Location')), GICP_EPSILON, GICP_D_MAX));

ICP_METHODS = {no_icp, icp_point_to_point, icp_point_to_plane, gicp}; 

for i=1:NUM_OF_CLOUDS
    % Get cloud name
    assert(NUM_OF_CLOUDS < 100 && NUM_OF_CLOUDS > 0)
    if(i < 10)
        cloud_name = ['00', num2str(i)];
    else
        cloud_name = ['0', num2str(i)];
    end
    cloud_name = [cloud_name, '.ply'];
    
    % Load local cloud
    local_cloud = pcread(['Database\Local_Clouds\', cloud_name]); 
    
    %Remove all points in the local cloud that are out of range
    local_cloud_range = sqrt(sum(local_cloud.Location.^2, 2)); 
    valid_idxs = find((local_cloud_range >= LOCAL_CLOUD_RANGE(1)).* ...
        (local_cloud_range <= LOCAL_CLOUD_RANGE(2)));
    local_cloud = select(local_cloud, valid_idxs); 
    
    % Get coarse registration
    coarse_reg_mat = COARS_REG_DB(i,:); 
    coarse_reg_mat = transpose(reshape(coarse_reg_mat, [4,4]));
    
    % Get ground truth registration
    gt_reg_mat = GROUND_TRUTH_DB(i,:);
    gt_reg_mat = transpose(reshape(gt_reg_mat, [4,4]));
    
    %Make sure this is a rotation matrix
    if(det(coarse_reg_mat(1:3,1:3)) ~= 1)
        assert(abs(det(coarse_reg_mat(1:3,1:3))-1) < 1e-4);
        [U,~,V] = svd(coarse_reg_mat(1:3,1:3));
        S = eye(size(U)); 
        coarse_reg_mat(1:3,1:3) = U*S*V';
    end
    
    % Apply coarse transform on the local cloud
    local_cloud_coarse = pctransform(local_cloud, affine3d(transpose(coarse_reg_mat)));
    
    % Get relevant part of global cloud
    global_cloud_relevant = extractRelevantPartFromGlobalCloud(global_cloud,...
        local_cloud, coarse_reg_mat);
    
    % Calculate normals for the point cloud (required for Point-to-Plane)
    local_cloud.Normal = pcnormals(local_cloud, ...
        NUM_OF_POINTS_FOR_NORMAL_CALC);
    global_cloud_relevant.Normal = pcnormals(global_cloud_relevant,...
        NUM_OF_POINTS_FOR_NORMAL_CALC);
    
    % Apply ICP registration
    for ind=1:length(ICP_METHODS)
        icp_method = ICP_METHODS{ind}; 
        
        reg_trans = icp_method.function_handle(local_cloud_coarse,...
            global_cloud_relevant);
        
        % Get the total registration transformation
        total_reg_mtx = reg_trans*coarse_reg_mat;
        
        %Make sure this is a rotation matrix
        if(det(total_reg_mtx(1:3,1:3)) ~= 1)
            assert(abs(det(total_reg_mtx(1:3,1:3))-1) < 1e-3);
            [U,~,V] = svd(total_reg_mtx(1:3,1:3));
            S = eye(size(U)); 
            total_reg_mtx(1:3,1:3) = U*S*V';
        end
        
        icp_method.errors = [icp_method.errors; ...
            getTransformationDiff(total_reg_mtx, gt_reg_mat)];
        
        ICP_METHODS{ind} = icp_method;
    end    
end

% Display average errors 
fprintf('\n\n');
fprintf('\t\t\t\t\t\t Localization Error[m] \t\t Yaw Absolute Error[deg] \t\t Pitch Absolute Error[deg] \t\t Roll Absolute Error[deg]\n');
for ind=1:length(ICP_METHODS)
    icp_method = ICP_METHODS{ind};
    abs_avg_errs = mean(abs(icp_method.errors), 1);
    abs_std_errs = std(abs(icp_method.errors), 0, 1);
    
    % Print Mean
    fprintf('%-18s Mean \t\t %f \t\t\t\t\t %f \t\t\t\t\t\t %f \t\t\t\t\t\t %f\n', ...
        icp_method.name, abs_avg_errs(1), abs_avg_errs(2),...
        abs_avg_errs(3), abs_avg_errs(4));
    
    % Print STD
    fprintf('%-18s STD \t\t\t %f \t\t\t\t\t %f \t\t\t\t\t\t %f \t\t\t\t\t\t %f\n', ...
        icp_method.name, abs_std_errs(1), abs_std_errs(2),...
        abs_std_errs(3), abs_std_errs(4));
end
fprintf('\n\n');

% Display errors 
fprintf('\t\t\t\t\t\t Localization Error[m] \t\t Yaw Absolute Error[deg] \t\t Pitch Absolute Error[deg] \t\t Roll Absolute Error[deg]\n');
for cloud_ind=1:NUM_OF_CLOUDS
    for ind=1:length(ICP_METHODS)
        icp_method = ICP_METHODS{ind};
        abs_errs = abs(icp_method.errors(cloud_ind,:));

        fprintf('Cloud %-2d %-18s \t\t %f \t\t\t\t\t %f \t\t\t\t\t\t %f \t\t\t\t\t\t %f\n', ...
            cloud_ind, icp_method.name, abs_errs(1), abs_errs(2),...
            abs_errs(3), abs_errs(4));
    end
end
fprintf('\n\n');