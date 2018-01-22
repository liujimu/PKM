% 机构标定
clear;
pkm = PKM();
target_poses = dlmread('.\calibration_20180101\calibration_target_poses_20180101.txt',',');
target_poses = target_poses(:,1:6);
measured_poses = dlmread('.\calibration_20180101\calibration_measured_poses_20180101.txt',',');
target_poses = target_poses';
measured_poses = measured_poses';

n = length(target_poses);
updated_poses = target_poses;
qin = zeros(6,n);
load('calibration_PSO_result.mat', 'param_errors')
param_errors = param_errors';
accumulated_param_errors = param_errors;
delta = 1;
eps = 10e-3;
k = 0;%迭代次数
while delta > eps
    W = zeros(6*n,60);
    for i = 1:n
        pkm.setPose(updated_poses(:,i));
        if k == 0
            qin(:,i) = pkm.q;
        end
        pkm.calVelJac();
        % 判断测试点是否在工作空间内，是否是奇异点
        %     if ~pkm.isInWorkspace()
        %         disp(i);
        %         disp(pkm.q);
        %         msg = 'Target pose is out of the workspace.';
        %         error(msg);
        %     end
        %     if ~pkm.isSingular()
        %         disp(i);
        %         disp(det(pkm.jac));
        %         msg = 'Target pose is a singular pose.';
        %         error(msg);
        %     end
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

        % 曹睿的标定公式
        % 杆长误差
        E_l = [pkm.l_dir; cross(pkm.rotm * pkm.S_init, (pkm.S_cur - pkm.U_cur))]';
        % U副初始位置误差
        E_u = zeros(6,18);
        for j = 1:6
            E_u(j,(3*j-2):3*j) = pkm.l_dir(:,j)' / (pkm.l_dir(:,j)' * pkm.P_dir(:,j));
        end
        E_u = pkm.jac * E_u;
        % S副初始位置误差
        E_s = zeros(6,18);
        for j = 1:6
            E_s(j,(3*j-2):3*j) = pkm.l_dir(:,j)' * pkm.rotm / (pkm.l_dir(:,j)' * pkm.P_dir(:,j));
        end
        E_s = pkm.jac * E_s;
        % P副方向误差
        E_p = zeros(6,18);
        for j = 1:6
            H = pkm.pose(1:3) + pkm.rotm * pkm.S_init(:,j) - pkm.U_init(:,j);
            E_p(j,(3*j-2):3*j) = pkm.q(j) * H' / (pkm.l(j) * pkm.l_dir(:,j)' * pkm.P_dir(:,j));
        end
        E_p = -pkm.jac * E_p;
        
        E = [inv(E_l) E_u E_s E_p];
        W((6*i-5):6*i,:) = E;
    end
%     disp(det(W'*W));
    error_poses = measured_poses - updated_poses;
    dX = reshape(error_poses,[],1);
    delta_position = norm(error_poses(1:3,:),1);
    delta_orientation = norm(error_poses(4:6,:),1);
    delta = delta_position + 50 * delta_orientation;
    param_errors = W\dX;
    accumulated_param_errors = accumulated_param_errors + param_errors;
    
    % 将误差参数带入机构，用正解计算位姿
    pkm = PKM(accumulated_param_errors);
    for i = 1:n
        pkm.forKin(qin(:,i), target_poses(:,i));
        updated_poses(:,i) = pkm.pose;
    end

    k = k + 1;
    fprintf('In the %dth iteration, delta = %f\n', k, delta);
    if k > 10
        break;
    end
end
