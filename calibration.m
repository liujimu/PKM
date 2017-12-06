% 机构标定
clear;
pkm = PKM();
target_poses = xlsread('target_poses.xlsx');
measured_poses = xlsread('measured_poses.xlsx');
% target_poses(:,6) = [];
% measured_poses(:,6) = [];

dX = reshape(measured_poses - target_poses,[],1);
n = length(target_poses);
updated_poses = zeros(6,n);
qin = zeros(6,n);
param_errors = zeros(60,1);
accumulated_param_errors = zeros(60,1);
delta = 1;
alw = 10e-5;
m = 0;%迭代次数
while delta > alw
    %W = zeros(6*n,60);
    W = [];
    for i = 1:n
        pkm.pose = target_poses(:,i);
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
        if m == 0
            qin(:,i) = pkm.q;
        end
        Jx = [pkm.l_dir; cross(pkm.rotm*pkm.S_init,pkm.l_dir)]';
        Jm = zeros(6,54);
        for j = 1:6
            Jm(j,(j*9-8):j*9) = [-pkm.l_dir(:,j)', pkm.l_dir(:,j)'*pkm.rotm, -pkm.q(j)*pkm.l_dir(:,j)'];
        end
        J = [inv(Jx), -inv(Jx)*Jm];
        
        W = [W; J];
    end
    disp(det(W'*W));
    param_errors = (W'*W)\W'*dX;
    accumulated_param_errors = accumulated_param_errors + param_errors;
    delta = max(abs(param_errors));
    pkm = PKM(accumulated_param_errors);
    for i = 1:n
        pkm.forKin(qin(:,i), target_poses(:,i));
        updated_poses(:,i) = pkm.pose;
    end
    dX = reshape(updated_poses - target_poses,[],1);

    m = m + 1;
    disp(m);
    if m > 10
        break;
    end
end
