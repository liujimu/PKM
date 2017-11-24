% 机构标定
clear;
pkm = PKM();
target_poses = zeros(6,20);
target_poses(:,1) = [0 0 0.05 0 0 0]';
target_poses(:,2) = [0.05 0 0 0.1 0.12 0.1]';
target_poses(:,3) = [0 0.05 0 0.05 0.08 0.09]';
target_poses(:,4) = [0 0 0.05 0.14 0.1 0.07]';
target_poses(:,5) = [0.02 0.03 0 0.11 0.06 0.04]';
target_poses(:,6) = [-0.06 0.03 0 -0.03 -0.06 -0.1]';
target_poses(:,7) = [0 0 0.05 0.15 -0.1 0.05]';
target_poses(:,8) = [0 -0.05 0 -0.05 0.1 0.05]';
target_poses(:,9) = [-0.05 0 0 0.1 0.15 -0.15]';
target_poses(:,10) = [0 0 -0.05 0 0 0]';
%measured_poses = load('measured_poses.txt');
%pose_errors = measured_poses - target_poses;
param_errors = zeros(60,1);

W = [];
for i = 1:10
    pkm.pose = target_poses(:,i);
    pkm.pose = [0.05 0 0 5/180*pi 10/180*pi 15/180*pi]';
    Jx = [pkm.l_dir; cross(pkm.rotm*pkm.S_init,pkm.l_dir)]';
    Jm = zeros(6,54);
    for j = 1:6
        Jm(j,(j*9-8):j*9) = [-pkm.l_dir(:,j)', pkm.l_dir(:,j)'*pkm.rotm, -pkm.q(j)*pkm.l_dir(:,j)'];
    end
    J = [inv(Jx), -inv(Jx)*Jm];
    
    W = [W; J];
end
B=W'*W;
det(W'*W)