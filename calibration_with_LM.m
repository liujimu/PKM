% 机构标定
clear;
pkm = PKM();
target_poses = dlmread('.\calibration_20180227\calibration_target_poses_20180227.txt',',');
target_poses = target_poses(:,1:6);
measured_poses = dlmread('.\calibration_20180227\calibration_measured_poses_20180227.txt',',');
target_poses = target_poses';
measured_poses = measured_poses';

n = length(target_poses);
target_input = zeros(6,n);
updated_input = zeros(6,n);

for i = 1:n
    pkm.setPose(target_poses(:,i));
    target_input(:,i) = pkm.q;
    pkm.setPose(measured_poses(:,i));
    updated_input(:,i) = pkm.q;
end

param_errors = zeros(54,1);

fun = @(e)cal_input_error(e,measured_poses,target_input,50);
options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','iter');
accumulated_param_errors = lsqnonlin(fun,param_errors,[],[],options);

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

function dq = cal_input_error(param_error, pose, input, n)
    %n = length(pose);
    input_error = zeros(6,n);
    for i = 1:n
        pkm = PKM(param_error);
        pkm.setPose(pose(:,i));
        input_error(:,i) = pkm.q - input(:,i);
    end
    dq = reshape(input_error,[],1);
end
