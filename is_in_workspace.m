function boolOutput = is_in_workspace( p, rotm )
%�жϸ���λ���Ƿ��ڻ��������ռ���
%   ���Ϊ����ֵ
%   pΪĩ��λ��
%   rotmΪĩ����ת����

    boolOutput = true;
    param = default_param();
    q_min = param.q_min;
    q_max = param.q_max;
    q = inverse_kinematics(p,rotm);
    for i=1:length(q)
        if(q(i) < q_min || q(i) > q_max)
            boolOutput = false;
            break
        end
    end
end

