% PSO机构标定
clear;
global target_poses measured_poses
target_poses = xlsread('target_poses.xlsx');
measured_poses = xlsread('measured_poses.xlsx');
nvars = 60;
ub = [0.002*ones(1,6), 0.02*ones(1,18), 0.001*ones(1,18), 0.02*ones(1,18)];
lb = -ub;
options = optimoptions('particleswarm','Display','iter','SwarmSize',100,'HybridFcn',@fmincon);
[param_errors,fval,exitflag,output] = particleswarm(@fitness_calibration,nvars,lb,ub,options);
