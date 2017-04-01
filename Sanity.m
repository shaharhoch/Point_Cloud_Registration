close all
clear all
clc 

% Parameters
NUM_OF_CLOUDS = 20;
global_cloud = pcread('Database\Global_Cloud\Global_Cloud.ply'); 
COARS_REG_DB = csvread('Database\Local2Global_Coarse.csv'); 
GROUND_TRUTH_DB = csvread('Database\Local2Global_GT.csv');

GICP_EPSILON = 1e-6;
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
no_icp.num_iter = [];
no_icp.function_handle = @(local, global_c) unitTransform();

icp_point_to_point.name = 'ICP Point to Point' ; 
icp_point_to_point.errors = []; 
icp_point_to_point.num_iter = [];
icp_point_to_point.function_handle = @(local, global_c) ...
    runMatlabICP(local,global_c, ...
    'pointToPoint',...
    ICP_EXTRAPOLATE,...
    ICP_MAX_ITERATIONS,...
    MAX_INLINER_RATIO,... % value not actually used as inlier ratio - it is used as maxInlierDistance!
    ICP_TOLERANCE);

icp_point_to_plane.name = 'ICP Point to Plane' ; 
icp_point_to_plane.errors = []; 
icp_point_to_plane.num_iter = [];
icp_point_to_plane.function_handle = @(local, global_c) ...
    runMatlabICP(local,global_c, ...
    'pointToPlane',...
    ICP_EXTRAPOLATE,...
    ICP_MAX_ITERATIONS,...
    MAX_INLINER_RATIO,... % value not actually used as inlier ratio - it is used as maxInlierDistance!
    ICP_TOLERANCE);

gicp.name = 'Generalized ICP' ; 
gicp.errors = []; 
gicp.num_iter = [];
gicp.function_handle = @(local, global_c) gicpWrapper(local, global_c,...
    MAX_INLINER_DISTANCE, GICP_EPSILON);

ICP_METHODS = {no_icp, icp_point_to_point, icp_point_to_plane, gicp}; 
%ICP_METHODS = {gicp, no_icp};

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
    
    % Get ground truth registration
    gt_reg_mat = [0.9999    0.0000   -0.0169   -0.6197; 0.0002    0.9999    0.0123   -0.8551; 0.0169   -0.0123    0.9998   -0.0062;   0         0         0    1.0000];
    
    %Make sure this is a rotation matrix
    if(det(gt_reg_mat(1:3,1:3)) ~= 1)
        [U,~,V] = svd(gt_reg_mat(1:3,1:3));
        S = eye(size(U)); 
        gt_reg_mat(1:3,1:3) = U*S*V';
    end
    
    global_cloud_relevant = pctransform(local_cloud, affine3d(transpose(gt_reg_mat)));
    coarse_reg_mat = eye(4);
    
    % Calculate normals for the point cloud (required for Point-to-Plane)
    local_cloud.Normal = pcnormals(local_cloud, ...
        NUM_OF_POINTS_FOR_NORMAL_CALC);
    global_cloud_relevant.Normal = pcnormals(global_cloud_relevant,...
        NUM_OF_POINTS_FOR_NORMAL_CALC);
    
    % Apply ICP registration
    for ind=1:length(ICP_METHODS)
        icp_method = ICP_METHODS{ind}; 
        
        [reg_trans, num_iter] = icp_method.function_handle(local_cloud,...
            global_cloud_relevant);
        
        % Get the total registration transformation
        total_reg_mtx = reg_trans*coarse_reg_mat;
        
        %Make sure this is a rotation matrix
        if(det(total_reg_mtx(1:3,1:3)) ~= 1)
            [U,~,V] = svd(total_reg_mtx(1:3,1:3));
            S = eye(size(U)); 
            total_reg_mtx(1:3,1:3) = U*S*V';
        end
        
        icp_method.errors = [icp_method.errors; ...
            getTransformationDiff(total_reg_mtx, gt_reg_mat)];
        
        icp_method.num_iter = [icp_method.num_iter; ...
            num_iter];
        
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
for ind=1:length(ICP_METHODS)
    icp_method = ICP_METHODS{ind};
    for cloud_ind=1:NUM_OF_CLOUDS
        abs_errs = abs(icp_method.errors(cloud_ind,:));

        fprintf('Cloud %-2d %-18s \t\t %f \t\t\t\t\t %f \t\t\t\t\t\t %f \t\t\t\t\t\t %f\n', ...
            cloud_ind, icp_method.name, abs_errs(1), abs_errs(2),...
            abs_errs(3), abs_errs(4));
    end
end
fprintf('\n\n');

fprintf('\t\t\t\t\t\t Num of iterations (Mean) \t\t\t Num of iterations (STD)\n');
for ind=1:length(ICP_METHODS)
    icp_method = ICP_METHODS{ind};
    avg_num_iter = mean(icp_method.num_iter, 1);
    std_num_iter = std(icp_method.num_iter, 0, 1);
    
    % Print
    fprintf('%-18s \t\t\t %f \t\t\t\t\t\t\t %f \n', ...
        icp_method.name, avg_num_iter, std_num_iter);
end
fprintf('\n\n');

% Plot error histograms
histogram_types = {'Localization Error', 'Realtive Rotation Error'};
units = {'[m]', '[deg]'};

% Print histograms
for i=1:length(ICP_METHODS)
    figure;
    error_vec = [ICP_METHODS{i}.errors(:,1), ...
        sum(abs(ICP_METHODS{i}.errors(:,2:4)),2)];
    for j=1:length(histogram_types)
        errors = abs(ICP_METHODS{i}.errors(:,j));
        subplot(1,2,j);
        histogram(errors, linspace(0,2,6))
        title_str = sprintf('%s for %s \n Mean: %.2f, STD: %.2f',...
            histogram_types{j}, ICP_METHODS{i}.name, mean(errors), std(errors));
        title(title_str)
        xlabel([histogram_types{j}, units{j}])
        ylabel('Num of Occurences')
    end
end

% Print graphs of all 24 errors
figure;
subplot(1,2,1); 
title('Cloud Localization Error')
xlabel('Cloud Index')
ylabel('Localization Error[m]')
hold on

subplot(1,2,2);
title('Cloud Realtive Rotation Error')
xlabel('Cloud Index')
ylabel('Realtive Rotation Error[deg]')
hold on

plot_types = {'-x', '-o', '-*', '-+'};

name_cell = {};
for i=1:length(ICP_METHODS)
    error_vec = [ICP_METHODS{i}.errors(:,1), ...
        sum(abs(ICP_METHODS{i}.errors(:,2:4)),2)];
    
    subplot(1,2,1)
    plot(error_vec(:,1), plot_types{i})
    hold on
    
    subplot(1,2,2)
    plot(error_vec(:,2), plot_types{i})
    hold on
    
    name_cell{i} = ICP_METHODS{i}.name;
end
subplot(1,2,1)
legend(name_cell)

subplot(1,2,2)
legend(name_cell)
    