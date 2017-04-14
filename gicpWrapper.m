function [ transform, max_iter ] = gicpWrapper( local_cloud,  global_cloud, ...
    gicp_epsilon, max_inliner_distance)

[transform, max_iter] = gicp_alg(double(local_cloud.Location),...
    double(global_cloud.Location), gicp_epsilon, max_inliner_distance);
    
clear mex;

transform = transpose(transform);
end

