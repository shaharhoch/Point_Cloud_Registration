function [ relevant_point_cloud ] = extractRelevantPartFromGlobalCloud( global_cloud, local_cloud, transform_matrix )
% This function extracts a relevant point cloud out of the global cloud.
% The relevant part starts at the translation of the transform and is at
% the size of the local cloud, with a security margin of MARGIN_SIZE_METERS

MARGIN_SIZE_METERS = 10; 

translation_vec = transform_matrix(1:3, 4); 

% Get global cloud relevant are
xlim = local_cloud.XLimits + translation_vec(1) + [-MARGIN_SIZE_METERS, MARGIN_SIZE_METERS];
ylim = local_cloud.YLimits + translation_vec(2) + [-MARGIN_SIZE_METERS, MARGIN_SIZE_METERS];
zlim = local_cloud.ZLimits + translation_vec(3) + [-MARGIN_SIZE_METERS, MARGIN_SIZE_METERS];
roi_matrix = [xlim; ylim; zlim]; 

relevant_ind = findPointsInROI(global_cloud, roi_matrix);
relevant_point_cloud = select(global_cloud, relevant_ind);
end

