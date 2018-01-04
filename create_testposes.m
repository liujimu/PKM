clear
%% 查找满足条件的标定测试点
pkm = PKM();
tic
poses_num = 50;
target_poses = zeros(6,poses_num);
target_inputs = zeros(6,poses_num);
i_pose = 0;
x = -0.06:0.01:0.06;
y = -0.06:0.01:0.06;
z = 0.04:0.01:0.2;
a = -0.3:0.05:0.3;
b = -0.3:0.05:0.3;
c = -0.5:0.05:0.5;
while i_pose < poses_num
    ix = ceil(length(x)*rand);
    iy = ceil(length(y)*rand);
    iz = ceil(length(z)*rand);
    ia = ceil(length(a)*rand);
    ib = ceil(length(b)*rand);
    ic = ceil(length(c)*rand);
    pose = [x(ix) y(iy) z(iz) a(ia) b(ib) c(ic)]';
    if ~ismember(pose',target_poses','rows')
        pkm.setPose(pose);
        if pkm.isInWorkspace(0.02) && pkm.checkSingularity() < 70
            i_pose = i_pose + 1;
            target_poses(:,i_pose) = pose;
            target_inputs(:,i_pose) = pkm.q;
        end
    end
end
toc

output_poses = target_poses';
output_poses(:,1:3) = 1000.*output_poses(:,1:3);
output_poses = sortrows(output_poses,3,'descend'); %按z坐标降序排序
output_poses(1,1:2) = zeros(1,2);
output_poses(1,4:6) = zeros(1,3);
dlmwrite('calibration_poses.txt',output_poses,'precision','%7.2f');


