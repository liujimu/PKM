function fitness = fitness_PSO(x)
%fitness_PSO 粒子群算法的适应度函数
%   此处显示详细说明
    global target_poses measured_poses qin
    param_errors = x';
    pkm = PKM(param_errors);
    calculated_poses = zeros(size(target_poses));
    n = length(target_poses);
    for j = 1:n
        pkm.forKin(qin(:,j), target_poses(:,j)); 
        calculated_poses(:,j) = pkm.pose;
    end
    dX = measured_poses - calculated_poses;
    dX1 = dX(1:3,:);
    dX2 = dX(4:6,:);
    fitness = max(max(abs(dX1))) + 10*max(max(abs(dX2)));
end

