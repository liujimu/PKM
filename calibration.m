% 机构标定
clear;
pkm = PKM();
target_poses = xlsread('target_poses.xlsx');
measured_poses = xlsread('measured_poses.xlsx');
% target_poses(:,6) = [];
% measured_poses(:,6) = [];

n = length(target_poses);
dX = reshape(measured_poses - target_poses,[],1);
updated_poses = zeros(6,n);
qin = zeros(6,n);
param_errors = zeros(60,1);
accumulated_param_errors = zeros(60,1);
delta = 1;
alw = 10e-5;
k = 0;%迭代次数
while delta > alw
    W = zeros(6*n,60);
    for i = 1:n
        if k == 0
            pkm.pose = target_poses(:,i);
            qin(:,i) = pkm.q;
        else
            pkm.pose = updated_poses(:,i);
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
        %             Jm(j,(j*9-8):j*9) = [-pkm.l_dir(:,j)', pkm.l_dir(:,j)'*pkm.rotm, -pkm.q(j)*pkm.l_dir(:,j)'];
        %         end
        %         J = [inv(Jx), -inv(Jx)*Jm];

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
    param_errors = (W'*W)\W'*dX;
    accumulated_param_errors = accumulated_param_errors + param_errors;
    pkm = PKM(accumulated_param_errors);
    for i = 1:n
        pkm.forKin(qin(:,i), target_poses(:,i));
        updated_poses(:,i) = pkm.pose;
    end
    dX = reshape(measured_poses - updated_poses,[],1);
    delta = max(abs(dX));

    k = k + 1;
    fprintf('In the %dth iteration, delta = %f\n', k, delta);
    if k > 10
        break;
    end
end
