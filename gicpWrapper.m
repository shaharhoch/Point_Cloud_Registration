function [ transform, max_iter ] = gicpWrapper( local_cloud,  global_cloud, ...
    gicp_epsilon, max_inliner_distance)

local_c = removeInvalidPoints(local_cloud);
global_c = removeInvalidPoints(global_cloud);

tic
[transform, max_iter] = gicp_alg(double(local_c.Location),...
    double(global_c.Location), gicp_epsilon, max_inliner_distance);
toc

transform = transpose(transform);
 end

