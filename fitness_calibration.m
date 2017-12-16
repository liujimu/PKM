function fitness_value = fitness_calibration( input )
%fitness_calibration 用于标定的适应度函数
%   此处显示详细说明
    pkm = PKM();
    target_poses = xlsread('target_poses.xlsx');
    measured_poses = xlsread('measured_poses.xlsx');
    n = length(target_poses);
    real_qin = zeros(6,n);
    for i = 1:n
        pkm.pose = target_poses(:,i);
        real_qin(:,i) = pkm.q;
    end
    param_errors = reshape(input,[],1);
    pkm_ = PKM(param_errors);
    calc_qin = zeros(6,n);
    for i = 1:n
        pkm_.pose = measured_poses(:,i);
        calc_qin(:,i) = pkm_.q;
    end
    dq = reshape(calc_qin - real_qin,[],1);
    fitness_value = norm(dq);
end

