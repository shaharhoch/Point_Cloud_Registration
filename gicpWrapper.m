function [ transform, max_iter ] = gicpWrapper( local_cloud,  global_cloud, ...
    max_inliner_distance, global_epsilon, local_epsilon_0, local_epsilon_1, pov)

if (nargin == 4)
    local_epsilon_0 = global_epsilon;
    local_epsilon_1 = 0; 
    pov = [0; 0; 0];
end

[transform, max_iter] = gicp_alg(double(local_cloud.Location),...
    double(global_cloud.Location), max_inliner_distance, global_epsilon,...
    local_epsilon_0, local_epsilon_1, pov);

transform = transpose(transform);
end

