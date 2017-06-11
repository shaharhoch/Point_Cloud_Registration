COARS_REG_DB = csvread('Database\Local2Global_Coarse.csv'); 
GROUND_TRUTH_DB = csvread('Database\Local2Global_GT.csv');
global_cloud = pcread('Database\Global_Cloud\Global_Cloud.ply'); 

LOCAL_CLOUD_RANGE = [2, 60]; %[min, max]
NUM_OF_POINTS_FOR_NORMAL_CALC = 6;
MAX_INLINER_DISTANCE = 7;

LOCAL_CLOUD_IDX = 1; 

RMSE_TH_VEC = logspace(1,-3, 10);
NUM_OF_TEST_POINTS = length(RMSE_TH_VEC); 

% Load local cloud
local_cloud = pcread(sprintf('Database\\Local_Clouds\\00%d.ply', LOCAL_CLOUD_IDX)); 

%Remove all points in the local cloud that are out of range
local_cloud_range = sqrt(sum(local_cloud.Location.^2, 2)); 
valid_idxs = find((local_cloud_range >= LOCAL_CLOUD_RANGE(1)).* ...
    (local_cloud_range <= LOCAL_CLOUD_RANGE(2)));
local_cloud = select(local_cloud, valid_idxs); 

% Get coarse registration
coarse_reg_mat = COARS_REG_DB(LOCAL_CLOUD_IDX,:); 
coarse_reg_mat = transpose(reshape(coarse_reg_mat, [4,4]));

% Get ground truth registration
gt_reg_mat = GROUND_TRUTH_DB(LOCAL_CLOUD_IDX,:);
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

rmse_vec = zeros(1, NUM_OF_TEST_POINTS);
run_time_vec = zeros(1, NUM_OF_TEST_POINTS);
% Apply ICP registration
for ind=1:NUM_OF_TEST_POINTS
    writeGoICPCfgFile( RMSE_TH_VEC(ind) );
    [reg_trans, num_iter, run_time] = GoICPWrapper(local_cloud_coarse,...
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
    
    rmse = getCloudsRMSE(global_cloud_relevant, local_shifted, MAX_INLINER_DISTANCE);
    
    rmse_vec(ind) = rmse;
    run_time_vec(ind) = run_time; 
end

figure; 
plot(run_time_vec, rmse_vec, 'x');

title('GO-ICP Performacne VS Run Time')
ylabel('RMSE[m]')
xlabel('Run Time[s]')

savefig('Results\Run_Time_VS_Performance')
saveas(gcf, 'Results\Run_Time_VS_Performance.jpg')