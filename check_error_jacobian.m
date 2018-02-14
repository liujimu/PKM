% 对比数值法求出的误差传递矩阵与曹睿公式求出的误差传递矩阵
clear

pkm = PKM();
% global target_pose target_input
target_poses = dlmread('.\calibration_20180101\calibration_target_poses_20180101.txt',',');
target_poses = target_poses(:,1:6);
n = length(target_poses);

W = [];
Wn = [];

for i = 1:n
    target_pose = target_poses(i,:)';
    pkm.setPose(target_pose);
    target_input = pkm.q;
    pkm.calVelJac();
    
    %% 曹睿的标定公式
%     % 杆长误差
%     E_l = [pkm.l_dir; cross(pkm.rotm * pkm.S_init, (pkm.S_cur - pkm.U_cur))]';
%     % U副初始位置误差
%     E_u = zeros(6,18);
%     for j = 1:6
%         E_u(j,(3*j-2):3*j) = pkm.l_dir(:,j)' / (pkm.l_dir(:,j)' * pkm.P_dir(:,j));
%     end
%     E_u = pkm.jac * E_u;
%     % S副初始位置误差
%     E_s = zeros(6,18);
%     for j = 1:6
%         E_s(j,(3*j-2):3*j) = pkm.l_dir(:,j)' * pkm.rotm / (pkm.l_dir(:,j)' * pkm.P_dir(:,j));
%     end
%     E_s = pkm.jac * E_s;
%     % P副方向误差
%     E_p = zeros(6,18);
%     for j = 1:6
%         H = pkm.pose(1:3) + pkm.rotm * pkm.S_init(:,j) - pkm.U_init(:,j);
%         E_p(j,(3*j-2):3*j) = pkm.q(j) * H' / (pkm.l(j) * pkm.l_dir(:,j)' * pkm.P_dir(:,j));
%     end
%     E_p = -pkm.jac * E_p;
%     
%     E = [inv(E_l) E_u E_s E_p];
%     W = [W;E];
    
    %% 数值法求出的误差传递矩阵
    param_error = zeros(54,1);
    error_jac = cal_numerical_jacobian(@(e)cal_pose_error(e,target_input,target_pose),param_error);
    Wn = [Wn; error_jac];
end

function pose_error = cal_pose_error(param_error, input, pose)
    pkm = PKM(param_error);
    pkm.forKin(input,pose);
    pose_error = pkm.pose - pose;
end