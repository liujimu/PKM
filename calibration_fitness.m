function fitness_value = calibration_fitness( input )
%fitness_calibration 用于标定的适应度函数
%   此处显示详细说明
    global real_qin measured_poses

    param_errors = reshape(input,[],1);
    pkm_ = PKM(param_errors);
    n = length(measured_poses);
    calc_qin = zeros(6,n);
    for i = 1:n
        pkm_.setPose(measured_poses(:,i));
        if isreal(pkm_.q)
            calc_qin(:,i) = pkm_.q;
        else
            calc_qin(:,i) = zeros(size(pkm_.q));
        end
    end
    dq = reshape(calc_qin - real_qin,[],1);
    %fitness_value = max(abs(dq));
    fitness_value = norm(dq);
end

