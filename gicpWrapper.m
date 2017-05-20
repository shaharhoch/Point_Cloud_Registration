function [ transform, num_iter, run_time ] = gicpWrapper( local_cloud,  global_cloud, ...
    gicp_epsilon, max_inliner_distance, max_iter)

local_c = removeInvalidPoints(local_cloud);
global_c = removeInvalidPoints(global_cloud);

tic;
[transform, num_iter] = gicp_alg(double(local_c.Location),...
    double(global_c.Location), gicp_epsilon, max_inliner_distance, max_iter);
run_time = toc; 

transform = transpose(transform);
 end

