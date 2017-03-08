function [ rmse ] = getCloudsRMSE( global_cloud, local_cloud, d_max )
% This function returns the RMSE of the clouds.
% Global cloud is the one that is supposed to contain the local cloud. 
% d_max is the maximum distance allowed for a point before it counts as an
% outlier. This is to prevent large errors in case a point the the local
% cloud doesn't exist in the global cloud. 

T = delaunayn(double(global_cloud.Location));
[~, dist] = dsearchn(double(global_cloud.Location),T,double(local_cloud.Location));

% Remove distances who are bigger than d_max
filtered_dist = dist(dist < d_max);

rmse = sqrt(mean(filtered_dist.^2)); 


end

