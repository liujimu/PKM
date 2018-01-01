clear
%% 查找满足条件的标定测试点
test_poses = [];
test_inputs = [];
pkm = PKM();
tic
for x = -0.06:0.02:0.06
    for y = -0.06:0.02:0.06
        for z = 0.04:0.02:0.2
            for a = -0.4:0.05:0.4
                for b = -0.4:0.05:0.4
                    for c = -0.5:0.1:0.5
                        pose = [x y z a b c]';
                        pkm.setPose(pose);
                        if pkm.isInWorkspace(0.02) && pkm.checkSingularity() < 70
                            test_poses = [test_poses, pose];
                            test_inputs = [test_inputs, pkm.q];
                        end
                    end
                end
            end
        end
    end
    fprintf('x = %f\n',x);
end
toc
n = length(test_poses);
index = round(n*rand(1,50));
target_poses = test_poses(:,index);
target_inputs = test_inputs(:,index);
output_poses = target_poses';
output_poses(:,1:3) = 1000.*output_poses(:,1:3);
sortrows(output_poses,3,'descend'); %按z坐标降序排序
output_poses(1,4:6) = zeros(1,3);
dlmwrite('calibration_poses.txt',output_poses,'precision','%7.2f');


