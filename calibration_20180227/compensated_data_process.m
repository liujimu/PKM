%% 将标定测量点转化为动平台位姿
clear;
load('compensated_data.mat');
n = length(center1);
xAxis = (right1 - center1)'; %转成列向量
yAxis = (up1 - center1)';


dx = zeros(1,n);
dy = zeros(1,n);
theta = zeros(1,n);
for i = 1:n
    dx(i) = norm(xAxis(:,i));
    dy(i) = norm(yAxis(:,i));
    theta(i) = acos(xAxis(:,i)'*yAxis(:,i)/dx(i)/dy(i));
end
fprintf('右靶球与中间靶球的距离偏差：%f\n',max(dx)-min(dx));
fprintf('上靶球与中间靶球的距离偏差：%f\n',max(dy)-min(dy));
fprintf('x轴和y轴夹角的平均值：%f\n',mean(theta)/pi*180);
fprintf('x轴和y轴夹角的偏差：%f\n',(max(theta)-min(theta))/pi*180);

% 靶球坐标系到动平台坐标系的变换矩阵，m为测量点构建的坐标系，p为动平台坐标系
Tm2p = eye(4);
Tm2p(3,4) = -91.65;

x_dir = zeros(3,n);
y_dir = zeros(3,n);
z_dir = zeros(3,n);
measured_poses = zeros(6,n);
for i = 1:n
    x_dir(:,i) = xAxis(:,i) / norm(xAxis(:,i));
    y_dir(:,i) = yAxis(:,i) / norm(yAxis(:,i));
    z_dir(:,i) = cross(x_dir(:,i),y_dir(:,i));
    z_dir(:,i) = z_dir(:,i) / norm(z_dir(:,i));
    y_dir(:,i) = cross(z_dir(:,i),x_dir(:,i));
    rotm = [x_dir(:,i) y_dir(:,i) z_dir(:,i)];
    p = [center1(i,:) 1]';
    Tm = [[rotm; zeros(1,3)], p]; %机架坐标系到测量点坐标系的变换矩阵
    Tp = Tm*Tm2p;
    beta = atan2(Tp(1,3), sqrt(Tp(2,3)^2 + Tp(3,3)^2));
    alpha = atan2(-Tp(2,3)/cos(beta), Tp(3,3)/cos(beta)); 
    gamma = atan2(-Tp(1,2)/cos(beta), Tp(1,1)/cos(beta)); 
    measured_poses(:,i) = [Tp(1:3,4); alpha; beta; gamma];
end    

% 以动平台零位为原点建立机器人坐标系，姿态由机架确定
offset_pose = [measured_poses(1:3,1);zeros(3,1)];
output_poses = measured_poses - repmat(offset_pose,1,n);
output_poses = output_poses(:,2:end)';
dlmwrite('compensated_measured_poses_20180227.txt',output_poses);
% 长度单位mm
target_poses = dlmread('calibration_target_poses_20180227.txt',',');
target_poses = target_poses(:,1:6);
pose_errors = output_poses - target_poses;
disp(max(abs(pose_errors)));

