% PSO�����궨
clear;
global real_qin measured_poses
target_poses = dlmread('.\calibration_20180101\calibration_target_poses_20180101.txt',',');
target_poses = target_poses(:,1:6);
measured_poses = dlmread('.\calibration_20180101\calibration_measured_poses_20180101.txt',',');
target_poses = target_poses';
measured_poses = measured_poses';

errors = zeros(60,1);
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
% �������λ�˶�Ӧ����Ӧ��ֵ
calc_qin = zeros(6,n);
for i = 1:n
    pkm.setPose(measured_poses(:,i));
    calc_qin(:,i) = pkm.q;
end
dq = reshape(calc_qin - real_qin,[],1);
initial_fitness = norm(dq);

% ����Ⱥ�Ż�
nvars = 60;
ub = [2*ones(1,6), 4*ones(1,18), 1*ones(1,18), 0.02*ones(1,18)];
lb = -ub;
options = optimoptions('particleswarm','Display','iter','SwarmSize',100,'HybridFcn',@fmincon);
[param_errors,fval,exitflag,output] = particleswarm(@calibration_fitness,nvars,lb,ub,options);

% ��������������������������λ��
pkm_new = PKM(param_errors');
updated_poses = zeros(size(target_poses));
for i = 1:n
    pkm_new.forKin(real_qin(:,i), target_poses(:,i));
    updated_poses(:,i) = pkm_new.pose;
end
pose_errors = updated_poses - measured_poses;
for i=i:6
    max_error(i) = max(abs(pose_errors(i,:)));
end

pose_errors_init = target_poses - measured_poses;
max_error_init = max(pose_errors_init') - min(pose_errors_init');

