function [ transform,  num_iter, run_time] = GoICPWrapper(local_cloud, global_cloud)

max_dist1 = max(local_cloud.Location(:))-min(global_cloud.Location(:));
max_dist2 = max(global_cloud.Location(:))-min(local_cloud.Location(:));
max_dist = max([max_dist1, max_dist2]);

min_translation = double(-1*abs(max_dist));
translation_width = double(abs(max_dist) * 2);

% min_translation = -0.5;
% translation_width = 1;

tic; 
[transform, num_iter] = GoICP(double(local_cloud.Location),...
    double(global_cloud.Location), min_translation, translation_width);
run_time = toc; 
end