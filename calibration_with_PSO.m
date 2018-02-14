% PSO机构标定
clear;
global real_qin measured_poses
target_poses = dlmread('.\calibration_20180101\calibration_target_poses_20180101.txt',',');
target_poses = target_poses(:,1:6);
measured_poses = dlmread('.\calibration_20180101\calibration_measured_poses_20180101.txt',',');
target_poses = target_poses';
measured_poses = measured_poses';

errors = zeros(54,1);
% ru = 82.1 / 1000;
% rl = 240 / 1000;
% thetau = 86.0151 / 180 * pi;
% thetal = 26.9868 / 180 * pi;
% l = 260 / 1000;
% qmin = -0.005;
% qmax = 0.26;
pkm = PKM();
n = length(target_poses);
real_qin = zeros(6,n);
for i = 1:n
    pkm.setPose(target_poses(:,i));
    real_qin(:,i) = pkm.q;
end
% 计算测量位姿对应的适应度值
calc_qin = zeros(6,n);
for i = 1:n
    pkm.setPose(measured_poses(:,i));
    calc_qin(:,i) = pkm.q;
end
dq = reshape(calc_qin - real_qin,[],1);
initial_fitness = norm(dq);

% 粒子群优化
nvars = length(errors);
ub = [1*ones(1,6), 0.01*ones(1,12), 2*ones(1,18), 0.5*ones(1,18)];
lb = -ub;
options = optimoptions('particleswarm','Display','iter','SwarmSize',100,'HybridFcn',@fmincon);
[param_errors,fval,exitflag,output] = particleswarm(@calibration_fitness,nvars,lb,ub,options);

% 将误差参数带入机构，用正解计算位姿
pkm_new = PKM(param_errors');
updated_poses = zeros(size(target_poses));
updated_qin = zeros(size(real_qin));
for i = 1:n
    pkm_new.forKin(real_qin(:,i), target_poses(:,i));
    updated_poses(:,i) = pkm_new.pose;
    pkm_new.setPose(measured_poses(:,i));
    updated_qin(:,i) = pkm_new.q;
end
pose_errors = updated_poses - measured_poses;
q_errors = updated_qin - real_qin;

max_pose_errors = max(abs(pose_errors),[],2);
max_q_errors = max(abs(q_errors),[],2);

pose_errors_init = target_poses - measured_poses;
max_error_init = max(pose_errors_init,[],2) - min(pose_errors_init,[],2);

