clear;
target_poses = xlsread('testdata.xlsx');
T = 1;
t = 0.05:0.05:T;
nt = length(t);
s = (1 - cos(pi * t / T))/2;
target_poses(:,1:3) = 0.001.*target_poses(:,1:3);
np = length(target_poses);
target_poses = [zeros(1,6);target_poses];
spline_poses = zeros(nt*np+1,6);
for i = 1:np
    spline_poses(2+(i-1)*nt:1+i*nt,:) = (1 - s)' * target_poses(i,:) + s' * target_poses(i+1,:);
end
spline_data = [(0:0.05:T*np)', spline_poses];
dlmwrite('spline_poses.txt',spline_data,'precision','%8.5f');
name
output_table = table(spline_data)