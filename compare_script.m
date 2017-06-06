close all
clear all
clc 

% Open log file
f_id = fopen('Results\log.txt', 'w'); 

% Parameters
NUM_OF_CLOUDS = 24;
global_cloud = pcread('Database\Global_Cloud\Global_Cloud.ply'); 
COARS_REG_DB = csvread('Database\Local2Global_Coarse.csv'); 
GROUND_TRUTH_DB = csvread('Database\Local2Global_GT.csv');

GICP_EPSILON = 1e-6;
NUM_OF_POINTS_FOR_NORMAL_CALC = 6;

LOCAL_CLOUD_RANGE = [2, 60]; %[min, max]

ICP_EXTRAPOLATE = false;
ICP_MAX_ITERATIONS = 100;
GICP_MAX_ITERATIONS = 1200;
MAX_INLINER_DISTANCE = 7;
MAX_INLINER_RATIO = MAX_INLINER_DISTANCE/100;
% (!!!) Effective value is 100*maxInlierDistance {For example, maxInlierDistance=0.07  --->  effectively 7m}
ICP_TOLERANCE = [0.1,0.1*pi/180];

% ICP Methods
no_icp.name = sprintf('No ICP');
no_icp.errors = []; 
no_icp.num_iter = [];
no_icp.function_handle = @(local, global_c) unitTransform();
no_icp.run_time = [];

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
icp_point_to_point.run_time = [];

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
icp_point_to_plane.run_time = [];

gicp.name = 'Generalized ICP' ; 
gicp.errors = []; 
gicp.num_iter = [];
gicp.function_handle = @(local, global_c) gicpWrapper(local, global_c,...
    GICP_EPSILON, MAX_INLINER_DISTANCE, GICP_MAX_ITERATIONS);
gicp.run_time = []; 

go_icp.name = 'Go-ICP' ; 
go_icp.errors = []; 
go_icp.num_iter = [];
go_icp.function_handle = @(local, global_c) GoICPWrapper(local, global_c);
go_icp.run_time = []; 

ICP_METHODS = {no_icp, icp_point_to_point, icp_point_to_plane, gicp, go_icp}; 
%ICP_METHODS = {go_icp}; 

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
        
        [reg_trans, num_iter, run_time] = icp_method.function_handle(local_cloud_coarse,...
            global_cloud_relevant);
        
        % Get the total registration transformation
        total_reg_mtx = reg_trans*coarse_reg_mat;
        
        %Make sure this is a rotation matrix
        if(det(total_reg_mtx(1:3,1:3)) ~= 1)
            [U,~,V] = svd(total_reg_mtx(1:3,1:3));
            S = eye(size(U)); 
            total_reg_mtx(1:3,1:3) = U*S*V';
        end
        
        local_shifted = pctransform(local_cloud,...
            affine3d(transpose(total_reg_mtx)));
        
        icp_method.errors = [icp_method.errors; ...
            getTransformationDiff(total_reg_mtx, gt_reg_mat),...
            getCloudsRMSE( global_cloud_relevant, local_shifted, MAX_INLINER_DISTANCE )];
        
        icp_method.num_iter = [icp_method.num_iter; ...
            num_iter];
        
        icp_method.run_time = [icp_method.run_time, run_time]; 
        
        ICP_METHODS{ind} = icp_method;
    end    
end

% Display average errors 
general_print(f_id,'\n\n');
general_print(f_id,'\t\t\t\t\t\t Localization Error[m] \t\t Yaw Absolute Error[deg] \t\t Pitch Absolute Error[deg] \t\t Roll Absolute Error[deg] \t\t RMSE[m]\n');
for ind=1:length(ICP_METHODS)
    icp_method = ICP_METHODS{ind};
    abs_avg_errs = mean(abs(icp_method.errors), 1);
    abs_std_errs = std(abs(icp_method.errors), 0, 1);
    
    % Print Mean
    general_print(f_id,'%-18s Mean \t\t %f \t\t\t\t\t %f \t\t\t\t\t\t %f \t\t\t\t\t\t %f \t\t\t\t %f\n', ...
        icp_method.name, abs_avg_errs(1), abs_avg_errs(2),...
        abs_avg_errs(3), abs_avg_errs(4), abs_avg_errs(5));
    
    % Print STD
    general_print(f_id,'%-18s STD \t\t\t %f \t\t\t\t\t %f \t\t\t\t\t\t %f \t\t\t\t\t\t %f \t\t\t\t %f\n', ...
        icp_method.name, abs_std_errs(1), abs_std_errs(2),...
        abs_std_errs(3), abs_std_errs(4), abs_std_errs(5));
end
general_print(f_id,'\n\n');

% Display errors 
general_print(f_id,'\t\t\t\t\t\t Localization Error[m] \t\t Yaw Absolute Error[deg] \t\t Pitch Absolute Error[deg] \t\t Roll Absolute Error[deg] \t\t RMSE[m]\n');
for ind=1:length(ICP_METHODS)
    icp_method = ICP_METHODS{ind};
    for cloud_ind=1:NUM_OF_CLOUDS
        abs_errs = abs(icp_method.errors(cloud_ind,:));

        general_print(f_id,'Cloud %-2d %-18s \t %f \t\t\t\t\t %f \t\t\t\t\t\t %f \t\t\t\t\t\t %f \t\t\t\t %f\n', ...
            cloud_ind, icp_method.name, abs_errs(1), abs_errs(2),...
            abs_errs(3), abs_errs(4), abs_errs(5));
    end
end
general_print(f_id,'\n\n');

general_print(f_id,'\t\t\t\t\t\t Num of iterations (Mean) \t\t\t Num of iterations (STD)\n');
for ind=1:length(ICP_METHODS)
    icp_method = ICP_METHODS{ind};
    avg_num_iter = mean(icp_method.num_iter, 1);
    std_num_iter = std(icp_method.num_iter, 0, 1);
    
    % Print
    general_print(f_id,'%-18s \t\t\t %f \t\t\t\t\t\t\t %f \n', ...
        icp_method.name, avg_num_iter, std_num_iter);
end
general_print(f_id,'\n\n');

% Plot error histograms
histogram_types = {'Localization Error', 'Realtive Rotation Error', 'RMSE'};
units = {'[m]', '[deg]', '[m]'};

% Print histograms
for i=1:length(ICP_METHODS)
    figure;
    error_vec = [ICP_METHODS{i}.errors(:,1), ...
        sum(abs(ICP_METHODS{i}.errors(:,2:4)),2), ICP_METHODS{i}.errors(:,5)];
    for j=1:length(histogram_types)
        errors = abs(ICP_METHODS{i}.errors(:,j));
        subplot(3,1,j);
        histogram(errors, linspace(0,2,6))
        title_str = sprintf('%s for %s \n Mean: %.2f, STD: %.2f',...
            histogram_types{j}, ICP_METHODS{i}.name, mean(errors), std(errors));
        title(title_str)
        xlabel([histogram_types{j}, units{j}])
        ylabel('Num of Occurences')
    end
    savefig(sprintf('Results\\%s_Histogram', ICP_METHODS{i}.name))
    saveas(gcf, sprintf('Results\\%s_Histogram.jpg', ICP_METHODS{i}.name))
end

% Print graphs of all 24 errors
figure;
subplot(3,1,1); 
title('Cloud Localization Error')
xlabel('Cloud Index')
ylabel('Localization Error[m]')
hold on

subplot(3,1,2);
title('Cloud Realtive Rotation Error')
xlabel('Cloud Index')
ylabel('Realtive Rotation Error[deg]')
hold on

subplot(3,1,3);
title('RMSE')
xlabel('Cloud Index')
ylabel('RMSE[m]')
hold on

plot_types = {'-x', '-o', '-*', '-+', '-d'};

name_cell = {};
for i=1:length(ICP_METHODS)
    error_vec = [ICP_METHODS{i}.errors(:,1), ...
        sum(abs(ICP_METHODS{i}.errors(:,2:4)),2), ICP_METHODS{i}.errors(:,5)];
    
    subplot(3,1,1)
    plot(error_vec(:,1), plot_types{i})
    hold on
    
    subplot(3,1,2)
    plot(error_vec(:,2), plot_types{i})
    hold on
    
    subplot(3,1,3)
    plot(error_vec(:,3), plot_types{i})
    hold on
    
    name_cell{i} = ICP_METHODS{i}.name;
end
subplot(3,1,1)
legend(name_cell)

subplot(3,1,2)
legend(name_cell)

subplot(3,1,3)
legend(name_cell)

savefig('Results\Error_Graph')
saveas(gcf, 'Results\Error_Graph.jpg')

% Print graphs of all 24 run_times
figure;
plot_types = {'-x', '-o', '-*', '-+', '-d'};

name_cell = {};
for i=1:length(ICP_METHODS)
    plot(ICP_METHODS{i}.run_time, plot_types{i})
    hold on
    
    name_cell{i} = ICP_METHODS{i}.name;
end

title('Registration time')
xlabel('Cloud Index')
ylabel('Run Time[sec]')
legend(name_cell)   
savefig('Results\Run_Time')
saveas(gcf, 'Results\Run_Time.jpg')

% Save entire workspce 
save('Results\Workspae_Save')
fclose(f_id);