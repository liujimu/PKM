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
updated_input = zeros(6,n);

for i = 1:n
    pkm.setPose(target_poses(:,i));
    target_input(:,i) = pkm.q;
    pkm.setPose(measured_poses(:,i));
    updated_input(:,i) = pkm.q;
end

param_errors = zeros(54,1);
accumulated_param_errors = param_errors;

delta = 1;
eps = 10e-3;
k = 0;%迭代次数
while delta > eps
    W = zeros(6*n,54);
    for i = 1:n
        %数值法求误差雅可比
        E = cal_numerical_jacobian(@(e)cal_input_error(e,measured_poses(:,i),target_input(:,i)),accumulated_param_errors);
        W((6*i-5):6*i,:) = E;
    end
    
    input_error = updated_input - target_input;
    dq = reshape(input_error,[],1);
    param_errors = W\dq;
    accumulated_param_errors = accumulated_param_errors - param_errors;
    
    % 将误差参数带入机构，用正解计算位姿
    pkm = PKM(accumulated_param_errors);
    for i = 1:n
        pkm.setPose(measured_poses(:,i));
        updated_input(:,i) = pkm.q;
    end

    % 迭代终止条件
    delta = norm(dq,Inf);
    k = k + 1;
    fprintf('In the %dth iteration, delta = %f\n', k, delta);
    if k > 10
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

function dq = cal_input_error(param_error, pose, input)
    pkm = PKM(param_error);
    pkm.setPose(pose);
    dq = pkm.q - input;
end
