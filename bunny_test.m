close all
clear all
clc 

bunny1 = pcread('Database\Bunny\bun000.ply');
bunny1.Color = uint8(repmat([255, 0, 0],size(bunny1.Location,1),1));

bunny2 = pcread('Database\Bunny\bun045.ply');
bunny2.Color = uint8(repmat([0, 0, 255],size(bunny2.Location,1),1));

GT_REG_MATRIX = [quat2rotm([0.955586 0.00548449 -0.294635 -0.0038555]).',...
    [-0.0520211; -0.000383981; -0.0109223]; [0, 0, 0, 1]];

GICP_EPSILON = 1e-6; 
D_MAX = 1; 

ICP_EXTRAPOLATE = false;
ICP_MAX_ITERATIONS = 100;
GICP_MAX_ITERATIONS = 2000;
MAX_INLINER_DISTANCE = 7;
MAX_INLINER_RATIO = MAX_INLINER_DISTANCE/100;
% (!!!) Effective value is 100*maxInlierDistance {For example, maxInlierDistance=0.07  --->  effectively 7m}
ICP_TOLERANCE = [0.1,0.1*pi/180];

% ICP Methods 
no_icp.name = sprintf('No ICP');
no_icp.errors = []; 
no_icp.num_iter = [];
no_icp.function_handle = @(local, global_c) unitTransform();

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

gicp.name = 'Generalized ICP' ; 
gicp.errors = []; 
gicp.num_iter = [];
gicp.function_handle = @(local, global_c) gicpWrapper(local, global_c,...
    GICP_EPSILON, D_MAX, GICP_MAX_ITERATIONS);

go_icp.name = 'Go-ICP' ; 
go_icp.errors = []; 
go_icp.num_iter = [];
go_icp.function_handle = @(local, global_c) GoICPWrapper(local, global_c);

%ICP_METHODS = {no_icp, icp_point_to_point, icp_point_to_plane, gicp, go_icp}; 
ICP_METHODS = {go_icp}; 

% Apply ICP registration
for ind=1:length(ICP_METHODS)
    icp_method = ICP_METHODS{ind}; 

    reg_trans = icp_method.function_handle(bunny2, bunny1);

    %Make sure this is a rotation matrix
    if(det(reg_trans(1:3,1:3)) ~= 1)
        [U,~,V] = svd(reg_trans(1:3,1:3));
        S = eye(size(U)); 
        reg_trans(1:3,1:3) = U*S*V';
    end
    
    bunny2_shifted = pctransform(bunny2, affine3d(transpose(reg_trans)));
    
    error_vec = [getTransformationDiff(reg_trans, GT_REG_MATRIX),...
        getCloudsRMSE( bunny1, bunny2_shifted, D_MAX )];
    icp_method.errors = [icp_method.errors; error_vec]; 
    
    ICP_METHODS{ind} = icp_method;
    
    % Display result
    subplot(2,3,ind);
    pcshow(bunny2_shifted);
    hold on
    pcshow(bunny1);
    title(icp_method.name)
    view([0,90])
end

% Display average errors 
fprintf('\t\t\t\t\t Localization Error[m] \t\t Yaw Absolute Error[deg] \t\t Pitch Absolute Error[deg] \t\t Roll Absolute Error[deg] \t\t RMSE[m]\n');
for ind=1:length(ICP_METHODS)
    icp_method = ICP_METHODS{ind};
    abs_errs = abs(icp_method.errors);
    
    %Print
    fprintf('%-18s \t\t %f \t\t\t\t\t\t %f \t\t\t\t\t %f \t\t\t\t\t\t %f \t\t\t\t\t %f\n', ...
        icp_method.name, abs_errs(1), abs_errs(2),...
        abs_errs(3), abs_errs(4), abs_errs(5));
end