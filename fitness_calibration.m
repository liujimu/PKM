function fitness_value = fitness_calibration( input )
%fitness_calibration ���ڱ궨����Ӧ�Ⱥ���
%   �˴���ʾ��ϸ˵��
    global real_qin measured_poses

    param_errors = reshape(input,[],1);
    pkm_ = PKM(param_errors);
    n = length(measured_poses);
    calc_qin = zeros(6,n);
    for i = 1:n
        pkm_.pose = measured_poses(:,i);
        if isreal(pkm_.q)
            calc_qin(:,i) = pkm_.q;
        else
            calc_qin(:,i) = ones(size(pkm_.q));
        end
    end
    dq = reshape(calc_qin - real_qin,[],1);
    fitness_value = norm(dq);
end

