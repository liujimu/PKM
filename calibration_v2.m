% �����궨
clear;
pkm = PKM();
target_poses = xlsread('target_poses.xlsx');
measured_poses = load('measured_poses.txt');
dX = reshape(measured_poses - target_poses,[],1);

n = length(target_poses)+1;
section_poses = [zeros(6,1) target_poses zeros(6,1)];

%W = zeros(6*n,60);
T = 100;
count = 0:n*T;
interpolated_pose = zeros(6,length(count));
for i = 1:length(count)-1
    j = floor((i - 1)/T) + 1;
    s = 1/2 - cos(mod(count(i),T) / T * pi) / 2;
    interpolated_pose(:,i) = section_poses(:,j).*(1-s) + section_poses(:,j+1).*s;
end
interpolated_pose(:,end) = section_poses(:,end);
q = [];
for i = 1:length(count)
    pkm.pose = interpolated_pose(:,i);
    if ~pkm.isInWorkspace()
        disp(i);
        disp(pkm.q);
        msg = 'Target pose is out of the workspace.';
        error(msg);
    end
    if ~pkm.isSingular()
        disp(i);
        disp(det(pkm.jac));
        msg = 'Target pose is a singular pose.';
        error(msg);
    end
    q = [q pkm.q];
%     Jx = [pkm.l_dir; cross(pkm.rotm*pkm.S_init,pkm.l_dir)]';
%     Jm = zeros(6,54);
%     for j = 1:6
%         Jm(j,(j*9-8):j*9) = [-pkm.l_dir(:,j)', pkm.l_dir(:,j)'*pkm.rotm, -pkm.q(j)*pkm.l_dir(:,j)'];
%     end
%     J = [inv(Jx), -inv(Jx)*Jm];
%     
%     W = [W; J];
end
% det(W'*W)
% param_errors = (W'*W)\W'*dX;

% target_poses = target_poses + 0.1*repmat([0 0 1 0 0 0]',[1,n]);
% target_poses_out = [1000 .* target_poses(1:3,:); (180/pi) .* target_poses(4:6,:)];
% target_poses_out = target_poses_out';

q_for_sim = [count./100;q]';