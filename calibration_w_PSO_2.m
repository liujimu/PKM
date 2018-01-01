% PSO机构标定
clear;
global real_qin measured_poses
target_poses = xlsread('target_poses.xlsx');
measured_poses = xlsread('measured_poses.xlsx');

errors = zeros(60,1);
ru = 82.1 / 1000;
rl = 240 / 1000;
thetau = 86.0151 / 180 * pi;
thetal = 26.9868 / 180 * pi;
l = 260 / 1000;
qmin = -0.005;
qmax = 0.2;
pkm = PKM(errors, ru, rl, thetau, thetal, l, qmin, qmax);
n = length(target_poses);
real_qin = zeros(6,n);
for i = 1:n
    pkm.pose = target_poses(:,i);
    real_qin(:,i) = pkm.q;
end

nvars = 60;
ub = [0.002*ones(1,6), 0.02*ones(1,18), 0.001*ones(1,18), 0.02*ones(1,18)];
lb = -ub;
options = optimoptions('particleswarm','Display','iter','SwarmSize',100,'HybridFcn',@fmincon);
[param_errors,fval,exitflag,output] = particleswarm(@fitness_calibration,nvars,lb,ub,options);
