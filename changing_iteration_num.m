close all
clear all
clc 

NUM_OF_GICP_METHODS = 4;
NOISE_STD = 1e-3;

bunny1 = pcread('Database\Bunny\bun000.ply');
bunny1.Color = uint8(repmat([255, 0, 0],size(bunny1.Location,1),1));

bunny2_orig = pcread('Database\Bunny\bun045.ply');
bunny2_orig.Color = uint8(repmat([0, 0, 255],size(bunny2_orig.Location,1),1));

 % Add noise to bunny2 
normals = pcnormals(bunny2_orig);
noise_size = NOISE_STD*normrnd(0, 1, [size(normals, 1), 1]); 
location = double(bunny2_orig.Location)+noise_size.*normals;
bunny2 = pointCloud(location);
bunny2.Color = uint8(repmat([0, 0, 255],size(bunny2.Location,1),1));

GT_REG_MATRIX = [quat2rotm([0.955586 0.00548449 -0.294635 -0.0038555]).',...
    [-0.0520211; -0.000383981; -0.0109223]; [0, 0, 0, 1]];

GICP_EPSILON = 1e-6; 
D_MAX = 1; 

ICP_EXTRAPOLATE = false;
MAX_INLINER_DISTANCE = D_MAX;
MAX_INLINER_RATIO = MAX_INLINER_DISTANCE/100;
% (!!!) Effective value is 100*maxInlierDistance {For example, maxInlierDistance=0.07  --->  effectively 7m}
ICP_TOLERANCE = [1e-6,1e-6];%[0.1,0.1*pi/180];

MAX_ITER_VEC = floor(linspace(10, 500, 30)); 

rre_mtx = zeros(length(MAX_ITER_VEC), NUM_OF_GICP_METHODS);
location_err_mtx = zeros(length(MAX_ITER_VEC), NUM_OF_GICP_METHODS);
running_time_vec = zeros(length(MAX_ITER_VEC), NUM_OF_GICP_METHODS);

for iter_ind=1:length(MAX_ITER_VEC)
    max_iter = MAX_ITER_VEC(iter_ind);
    
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
        max_iter,...
        MAX_INLINER_RATIO,... % value not actually used as inlier ratio - it is used as maxInlierDistance!
        ICP_TOLERANCE);

    icp_point_to_plane.name = 'ICP Point to Plane' ; 
    icp_point_to_plane.errors = []; 
    icp_point_to_plane.num_iter = [];
    icp_point_to_plane.function_handle = @(local, global_c) ...
        runMatlabICP(local,global_c, ...
        'pointToPlane',...
        ICP_EXTRAPOLATE,...
        max_iter,...
        MAX_INLINER_RATIO,... % value not actually used as inlier ratio - it is used as maxInlierDistance!
        ICP_TOLERANCE);

    gicp.name = 'Generalized ICP' ; 
    gicp.errors = []; 
    gicp.num_iter = [];
    gicp.function_handle = @(local, global_c) gicpWrapper(local, global_c,...
        GICP_EPSILON, D_MAX, max_iter);

    ICP_METHODS = {no_icp, icp_point_to_point, icp_point_to_plane, gicp}; 
    %ICP_METHODS = {gicp}; 
    assert(NUM_OF_GICP_METHODS == length(ICP_METHODS))

    % Apply ICP registration
    figure('Visible','off');
    for ind=1:length(ICP_METHODS)
        icp_method = ICP_METHODS{ind}; 
        
        tic; 
        reg_trans = icp_method.function_handle(bunny2, bunny1);
        running_time = toc; 
        
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

        rre_mtx(iter_ind, ind) = sum(abs(error_vec(1:3)));
        location_err_mtx(iter_ind, ind) = abs(error_vec(4));
        running_time_vec(iter_ind, ind) = running_time;

        ICP_METHODS{ind} = icp_method;

        % Display result 
        subplot(2, 2, ind)
        pcshow(bunny2_shifted);
        hold on
        pcshow(bunny1);
        title(icp_method.name)
        view([0,90])
    end
    suptitle(sprintf('Registration for %d Max Iterations', max_iter))
    pic_name = sprintf('Results/Max_Iter_%d', max_iter);
    pic_name = strrep(pic_name,'.','_');
    saveas(gcf, pic_name, 'jpg');
    close all
end

% Plot RRE
figure; 
legend_list = {}; 
for i=1:length(ICP_METHODS)
    plot(MAX_ITER_VEC, rre_mtx(:, i));
    hold on
    legend_list{i} = ICP_METHODS{i}.name; 
end
title('RRE VS Max Iterations')
xlabel('Max Iter[A.U]')
ylabel('RRE[A.U]')
legend(legend_list);
saveas(gcf, 'Results/RRE_Max_Iter', 'jpg');
savefig('Results/RRE_Max_Iter')

% Plot Localization Error
figure; 
legend_list = {}; 
for i=1:length(ICP_METHODS)
    plot(MAX_ITER_VEC, location_err_mtx(:, i));
    hold on
    legend_list{i} = ICP_METHODS{i}.name; 
end
title('Localization Error VS Max Iterations')
xlabel('Max Iter[A.U]')
ylabel('Localization Error[A.U]')
legend(legend_list);
saveas(gcf, 'Results/Localization_Error_Max_Iter', 'jpg');
savefig('Results/Localization_Error_Max_Iter')

% Plot Running Time
figure; 
legend_list = {}; 
for i=1:length(ICP_METHODS)
    plot(MAX_ITER_VEC, running_time_vec(:, i));
    hold on
    legend_list{i} = ICP_METHODS{i}.name; 
end
title('Running Time VS Max Iterations')
xlabel('Max Iter[A.U]')
ylabel('Running Time[sec]')
legend(legend_list);
saveas(gcf, 'Results/Running_Time', 'jpg');
savefig('Results/Running_Time')

close all
clear mex

