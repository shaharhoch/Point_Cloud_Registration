close all
clear all
clc 

bunny1 = pcread('Database\Bunny\bun000.ply');
bunny1.Color = uint8(repmat([255, 0, 0],size(bunny1.Location,1),1));

bunny2_orig = pcread('Database\Bunny\bun045.ply');
bunny2_orig.Color = uint8(repmat([0, 0, 255],size(bunny2_orig.Location,1),1));

GT_REG_MATRIX = [quat2rotm([0.955586 0.00548449 -0.294635 -0.0038555]).',...
    [-0.0520211; -0.000383981; -0.0109223]; [0, 0, 0, 1]];

GICP_EPSILON = 1e-6; 
D_MAX = 1; 

ICP_EXTRAPOLATE = false;
ICP_MAX_ITERATIONS = 100;
MAX_INLINER_DISTANCE = D_MAX;
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
    GICP_EPSILON, D_MAX);

ICP_METHODS = {no_icp, icp_point_to_point, icp_point_to_plane, gicp}; 
%ICP_METHODS = {gicp}; 

EPSILON_VALS = logspace(-8, -3, 10); 
rmse_bunny1 = zeros(size(EPSILON_VALS));
rmse_bunny2_orig = zeros(size(EPSILON_VALS));
rre_mtx = zeros(length(EPSILON_VALS), length(ICP_METHODS));
location_err_mtx = zeros(length(EPSILON_VALS), length(ICP_METHODS));

for eps_ind=1:length(EPSILON_VALS)
    % Modify GICP according to noise epsilon
    gicp.function_handle = @(local, global_c) gicpWrapper(local, global_c,...
    EPSILON_VALS(eps_ind)*10, D_MAX);

    ICP_METHODS = {no_icp, icp_point_to_point, icp_point_to_plane, gicp}; 
    
    % Add noise to bunny2 
    normals = pcnormals(bunny2_orig);
    noise_size = sqrt(EPSILON_VALS(eps_ind))*normrnd(0, 1,...
        [size(normals, 1), 1]); 
    location = double(bunny2_orig.Location)+noise_size.*normals;
    bunny2 = pointCloud(location);
    bunny2.Color = uint8(repmat([0, 0, 255],size(bunny2.Location,1),1));
    
    % Calculate RMSE
    rmse_bunny1(eps_ind) = getCloudsRMSE(bunny1, bunny2, D_MAX);
    rmse_bunny2_orig(eps_ind) = getCloudsRMSE(bunny2_orig, bunny2, D_MAX);
    
    % Apply ICP registration
    figure('Visible','off');
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
        
        rre_mtx(eps_ind, ind) = sum(abs(error_vec(1:3)));
        location_err_mtx(eps_ind, ind) = abs(error_vec(4));

        ICP_METHODS{ind} = icp_method;

        % Display result 
        subplot(2, 2, ind)
        pcshow(bunny2_shifted);
        hold on
        pcshow(bunny1);
        title(icp_method.name)
        view([0,90])
    end
    pic_name = sprintf('Results/Epsilon_%.2e', EPSILON_VALS(eps_ind));
    pic_name = strrep(pic_name,'.','_');
    saveas(gcf, pic_name, 'jpg');
end

% Plot RMSE
figure; 
plot(EPSILON_VALS, rmse_bunny1);
title('Bunny1 RMSE VS Epsilon')
xlabel('Epsilon[A.U]')
ylabel('RMSE[A.U]')
saveas(gcf, 'Results/RMSE_Bunny1', 'jpg');

figure; 
loglog(EPSILON_VALS, rmse_bunny2_orig);
title('Bunny2 RMSE VS Epsilon')
xlabel('Epsilon[A.U]')
ylabel('RMSE[A.U]')
saveas(gcf, 'Results/RMSE_Bunny2', 'jpg');

% Plot RRE
figure; 
legend_list = {}; 
for i=1:length(ICP_METHODS)
    semilogx(EPSILON_VALS, rre_mtx(:, i));
    hold on
    legend_list{i} = ICP_METHODS{i}.name; 
end
title('RRE VS Epsilon')
xlabel('Epsilon[A.U]')
ylabel('RRE[A.U]')
legend(legend_list);
saveas(gcf, 'Results/RRE', 'jpg');
savefig('Results/RRE')

% Plot Localization Error
figure; 
legend_list = {}; 
for i=1:length(ICP_METHODS)
    semilogx(EPSILON_VALS, location_err_mtx(:, i));
    hold on
    legend_list{i} = ICP_METHODS{i}.name; 
end
title('Localization Error VS Epsilon')
xlabel('Epsilon[A.U]')
ylabel('Localization Error[A.U]')
legend(legend_list);
saveas(gcf, 'Results/Localization_Error', 'jpg');
savefig('Results/Localization_Error')

close all
clear mex

