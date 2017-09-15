function boolOutput = is_in_workspace( p, rotm )
%判断给定位姿是否在机构工作空间内
%   输出为布尔值
%   p为末端位置
%   rotm为末端旋转矩阵

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

