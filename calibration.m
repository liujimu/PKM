% 机构标定
clear;
pkm = PKM();
target_poses = dlmread('.\calibration_20180101\calibration_target_poses_20180101.txt',',');
target_poses = target_poses(:,1:6);
measured_poses = dlmread('.\calibration_20180101\calibration_measured_poses_20180101.txt',',');
target_poses = target_poses';
measured_poses = measured_poses';

n = length(target_poses);
target_input = zeros(6,n);
for i = 1:n
    pkm.setPose(target_poses(:,i));
    target_input(:,i) = pkm.q;
end

updated_poses = target_poses;

param_errors = zeros(54,1);
accumulated_param_errors = param_errors;

delta = 1;
eps = 10e-3;
k = 0;%迭代次数
while delta > eps
    W = zeros(6*n,54);
    for i = 1:n
%         pkm.setPose(updated_poses(:,i));
%         if k == 0
%             qin(:,i) = pkm.q;
%         end
%         pkm.calVelJac();

%         Jx = [pkm.l_dir; cross(pkm.rotm*pkm.S_init,pkm.l_dir)]';
%         Jm = zeros(6,54);
%         for j = 1:6
%             % Jm(j,(j*9-8):j*9) = [-pkm.l_dir(:,j)', pkm.l_dir(:,j)'*pkm.rotm, -pkm.q(j)*pkm.l_dir(:,j)'];
%             Jm(j,(3*j-2):(3*j)) = -pkm.l_dir(:,j)';
%             Jm(j,(18+3*j-2):(18+3*j)) = pkm.l_dir(:,j)'*pkm.rotm;
%             Jm(j,(36+3*j-2):(36+3*j)) = -pkm.q(j)*pkm.l_dir(:,j)';
%         end
%         J = [inv(Jx), -inv(Jx)*Jm];
%         W((6*i-5):6*i,:) = J;

%         % 曹睿的标定公式
%         % 杆长误差
%         E_l = [pkm.l_dir; cross(pkm.rotm * pkm.S_init, (pkm.S_cur - pkm.U_cur))]';
%         % U副初始位置误差
%         E_u = zeros(6,18);
%         for j = 1:6
%             E_u(j,(3*j-2):3*j) = pkm.l_dir(:,j)' / (pkm.l_dir(:,j)' * pkm.P_dir(:,j));
%         end
%         E_u = pkm.jac * E_u;
%         % S副初始位置误差
%         E_s = zeros(6,18);
%         for j = 1:6
%             E_s(j,(3*j-2):3*j) = pkm.l_dir(:,j)' * pkm.rotm / (pkm.l_dir(:,j)' * pkm.P_dir(:,j));
%         end
%         E_s = pkm.jac * E_s;
%         % P副方向误差
%         E_p = zeros(6,18);
%         for j = 1:6
%             H = pkm.pose(1:3) + pkm.rotm * pkm.S_init(:,j) - pkm.U_init(:,j);
%             E_p(j,(3*j-2):3*j) = pkm.q(j) * H' / (pkm.l(j) * pkm.l_dir(:,j)' * pkm.P_dir(:,j));
%         end
%         E_p = -pkm.jac * E_p;
%         
%         E = [inv(E_l) E_u E_s E_p];

        %数值法求误差雅可比
        E = cal_numerical_jacobian(@(e)cal_pose_error(e,target_input(:,i),target_poses(:,i)),accumulated_param_errors);
        W((6*i-5):6*i,:) = E;
    end
%     disp(det(W'*W));
    error_poses = measured_poses - updated_poses;
    dX = reshape(error_poses,[],1);
    param_errors = W\dX;
    accumulated_param_errors = accumulated_param_errors + param_errors;
    
    % 将误差参数带入机构，用正解计算位姿
    pkm = PKM(accumulated_param_errors);
    for i = 1:n
        pkm.forKin(target_input(:,i), target_poses(:,i));
        updated_poses(:,i) = pkm.pose;
    end

    % 迭代终止条件
    delta_position = norm(error_poses(1:3,:),1);
    delta_orientation = norm(error_poses(4:6,:),1);
    delta = delta_position + 50 * delta_orientation;
    k = k + 1;
    fprintf('In the %dth iteration, delta = %f\n', k, delta);
    if k > 5
        break;
    end
end

% 将误差参数带入机构，用正解计算位姿
pkm_new = PKM(accumulated_param_errors);
updated_poses = zeros(size(target_poses));
updated_input = zeros(size(target_input));
for i = 1:n
    pkm_new.forKin(target_input(:,i), target_poses(:,i));
    updated_poses(:,i) = pkm_new.pose;
    pkm_new.setPose(measured_poses(:,i));
    updated_input(:,i) = pkm_new.q;
end
pose_errors = updated_poses - measured_poses;
input_errors = updated_input - target_input;

max_pose_errors = max(abs(pose_errors),[],2);
max_input_errors = max(abs(input_errors),[],2);

function pose_error = cal_pose_error(param_error, input, pose)
    pkm = PKM(param_error);
    pkm.forKin(input,pose);
    pose_error = pkm.pose - pose;
end
