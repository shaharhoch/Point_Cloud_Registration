close all
clear all
clc 

bunny1 = pcread('Database\Bunny\bun000.ply');
bunny1.Color = uint8(repmat([255, 0, 0],size(bunny1.Location,1),1));

bunny2 = pcread('Database\Bunny\bun045.ply');
bunny2.Color = uint8(repmat([0, 0, 255],size(bunny2.Location,1),1));

POV = [0; 0; 0]; 
LOCAL_EPS_0 = logspace(-10, -1, 20); 
LOCAL_EPS_1 = zeros(size(LOCAL_EPS_0)); 
GLOBAL_EPS = 1e-6; 

GT_REG_MATRIX = [quat2rotm([0.955586 0.00548449 -0.294635 -0.0038555]).',...
    [-0.0520211; -0.000383981; -0.0109223]; [0, 0, 0, 1]];

GICP_EPSILON = 1e-4; 
D_MAX = 4e-2; 

ICP_EXTRAPOLATE = false;
ICP_MAX_ITERATIONS = 100;
MAX_INLINER_DISTANCE = 0.1;
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
    D_MAX, GICP_EPSILON);

gicp_change.name = 'Generalized ICP CE' ; 
gicp_change.errors = []; 
gicp_change.num_iter = [];
gicp_change.function_handle = @(local, global_c) gicpWrapper(local, global_c,...
        D_MAX, GLOBAL_EPS);

ICP_METHODS = {no_icp, icp_point_to_point, icp_point_to_plane, gicp_change}; 
%ICP_METHODS = {gicp_change};

location_err_vec = zeros(length(ICP_METHODS), length(LOCAL_EPS_0)); 
rre_vec = zeros(length(ICP_METHODS), length(LOCAL_EPS_0));
for run_ind=1:length(LOCAL_EPS_0)
    % Add noise to bunny2 (local cloud)
    local_normals = double(pcnormals(bunny2));
    dist = sum((double(bunny2.Location)-POV.').^2, 2);
    noise_size2 = sqrt((LOCAL_EPS_0(run_ind)+LOCAL_EPS_1(run_ind)*dist)).*normrnd(0,1,[size(local_normals,1),1]);
    bunny2 = pointCloud(double(bunny2.Location)+noise_size2.*local_normals);
    bunny2.Color = uint8(repmat([0, 0, 255],size(bunny2.Location,1),1));

    % Apply ICP registration
    figure('Visible','Off');
    for ind=1:length(ICP_METHODS)
        icp_method = ICP_METHODS{ind}; 

        reg_trans = icp_method.function_handle(bunny2, bunny1);

        %Make sure this is a rotation matrix
        rot_det = det(reg_trans(1:3,1:3));
        if(rot_det ~= 1)
            assert(abs(rot_det-1) < 1e-4, num2str(rot_det));
            [U,~,V] = svd(reg_trans(1:3,1:3));
            S = eye(size(U)); 
            reg_trans(1:3,1:3) = U*S*V';
        end

        bunny2_shifted = pctransform(bunny2, affine3d(transpose(reg_trans)));

        error_vec = [getTransformationDiff(reg_trans, GT_REG_MATRIX),...
            getCloudsRMSE( bunny1, bunny2_shifted, D_MAX )];
        icp_method.errors = [icp_method.errors; error_vec]; 

        location_err = error_vec(1); 
        rre = sum(abs(error_vec(2:4)));
        location_err_vec(ind, run_ind) = location_err; 
        rre_vec(ind, run_ind) = rre;

        ICP_METHODS{ind} = icp_method;

        % Display result
        subplot(2,2,ind);
        pcshow(bunny2_shifted);
        hold on
        pcshow(bunny1);
        title(icp_method.name)
        view([0,90])
    end
    suptitle(['ICP Performance for Epsilon = ', num2str(LOCAL_EPS_0(run_ind))])
    saveas(gcf, ['Results/epsilon_', num2str(LOCAL_EPS_0(run_ind)), '.jpg'])
    close all

    % Display average errors 
%     fprintf('\t\t\t\t\t Localization Error[m] \t\t Yaw Absolute Error[deg] \t\t Pitch Absolute Error[deg] \t\t Roll Absolute Error[deg] \t\t RMSE[m]\n');
%     for ind=1:length(ICP_METHODS)
%         icp_method = ICP_METHODS{ind};
%         abs_errs = abs(icp_method.errors);
% 
%         %Print
%         fprintf('%s \t\t %f \t\t\t\t\t\t %f \t\t\t\t\t %f \t\t\t\t\t\t %f \t\t\t\t\t %f\n', ...
%             icp_method.name, abs_errs(1), abs_errs(2),...
%             abs_errs(3), abs_errs(4), abs_errs(5));
%     end

end

% Plot RRE graph
figure('Visible','Off'); 
icp_legend = {};
for i=1:length(ICP_METHODS)
    semilogx(LOCAL_EPS_0, rre_vec(i,:))
    hold on 
    icp_legend{i} = ICP_METHODS{i}.name;
end
title('RRE vs Epsilon')
xlabel('Epsilon[AU]')
ylabel('RRE[deg]')
legend(icp_legend)
saveas(gcf, 'Results/RRE.jpg')

% Plot Localization Error graph
figure('Visible','Off'); 
icp_legend = {};
for i=1:length(ICP_METHODS)
    semilogx(LOCAL_EPS_0, location_err_vec(i,:))
    hold on 
    icp_legend{i} = ICP_METHODS{i}.name;
end
title('Localization Error vs Epsilon')
xlabel('Epsilon[AU]')
ylabel('Localization Error[m]')
legend(icp_legend)
saveas(gcf, 'Results/Localization.jpg')
    
    