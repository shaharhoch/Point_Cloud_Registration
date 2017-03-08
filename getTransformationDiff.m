function [ diff_vec] = getTransformationDiff(mtx1, mtx2)
% This function returns the difference between transformations mtx1 and
% mtx2. The difference is given in [translation, yaw, pitch roll] vector.
% The angles are given in degrees.

diff_trans = mtx1(1:3,1:3)*transpose(mtx2(1:3,1:3)); 

diff_vec = zeros(1,4); 
diff_vec(1) = norm(mtx1(1:3, 4)-mtx2(1:3, 4));
diff_vec(2:4) = rad2deg(rotm2eul(diff_trans(1:3,1:3)));

end

