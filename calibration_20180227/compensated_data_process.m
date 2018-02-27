%% ���궨������ת��Ϊ��ƽ̨λ��
clear;
load('compensated_data.mat');
n = length(center1);
xAxis = (right1 - center1)'; %ת��������
yAxis = (up1 - center1)';


dx = zeros(1,n);
dy = zeros(1,n);
theta = zeros(1,n);
for i = 1:n
    dx(i) = norm(xAxis(:,i));
    dy(i) = norm(yAxis(:,i));
    theta(i) = acos(xAxis(:,i)'*yAxis(:,i)/dx(i)/dy(i));
end
fprintf('�Ұ������м����ľ���ƫ�%f\n',max(dx)-min(dx));
fprintf('�ϰ������м����ľ���ƫ�%f\n',max(dy)-min(dy));
fprintf('x���y��нǵ�ƽ��ֵ��%f\n',mean(theta)/pi*180);
fprintf('x���y��нǵ�ƫ�%f\n',(max(theta)-min(theta))/pi*180);

% ��������ϵ����ƽ̨����ϵ�ı任����mΪ�����㹹��������ϵ��pΪ��ƽ̨����ϵ
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
    Tm = [[rotm; zeros(1,3)], p]; %��������ϵ������������ϵ�ı任����
    Tp = Tm*Tm2p;
    beta = atan2(Tp(1,3), sqrt(Tp(2,3)^2 + Tp(3,3)^2));
    alpha = atan2(-Tp(2,3)/cos(beta), Tp(3,3)/cos(beta)); 
    gamma = atan2(-Tp(1,2)/cos(beta), Tp(1,1)/cos(beta)); 
    measured_poses(:,i) = [Tp(1:3,4); alpha; beta; gamma];
end    

% �Զ�ƽ̨��λΪԭ�㽨������������ϵ����̬�ɻ���ȷ��
offset_pose = [measured_poses(1:3,1);zeros(3,1)];
output_poses = measured_poses - repmat(offset_pose,1,n);
output_poses = output_poses(:,2:end)';
dlmwrite('compensated_measured_poses_20180227.txt',output_poses);
% ���ȵ�λmm
target_poses = dlmread('calibration_target_poses_20180227.txt',',');
target_poses = target_poses(:,1:6);
pose_errors = output_poses - target_poses;
disp(max(abs(pose_errors)));

