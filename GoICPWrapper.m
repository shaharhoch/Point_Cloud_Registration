function [ transform,  num_iter, run_time] = GoICPWrapper(local_cloud, global_cloud)

% Make both clouds be in [-0.5,0.5]^3 cube. 
% Centralize global cloud
global_cloud_trans = mean(global_cloud.Location, 1); 
global_cloud_shifted = global_cloud.Location - global_cloud_trans;

% Centralize local cloud
% Currently I centralize the local cloud with the same translation as the
% glonal cloud. This doesn't have to be very good, but becasue the input to
% this function is the local cloud that is already roughly shifted to the
% correct location, I do this in order to not ruin that. 
%local_cloud_trans = mean(local_cloud.Location, 1); 
local_cloud_shifted = local_cloud.Location - global_cloud_trans;
%assert(abs(mean(local_cloud_shifted(:))) < 5);

%Scale both clouds
scaling_factor = max([global_cloud_shifted(:); local_cloud_shifted(:)]); 
local_cloud_shifted = local_cloud_shifted/scaling_factor;
global_cloud_shifted = global_cloud_shifted/scaling_factor;

tic; 
[transform_orig, num_iter] = GoICP(double(local_cloud_shifted),...
    double(global_cloud_shifted));
run_time = toc; 

% Shift transform matrix to fit original clouds
transform_scaled = transform_orig;
transform_scaled(1:3,end) = transform_orig(1:3,end)*scaling_factor + ...
    global_cloud_trans.' - transform_scaled(1:3,1:3)*(global_cloud_trans.');

transform = transform_scaled;
end