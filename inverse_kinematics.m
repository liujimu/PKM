function [q, J] = inverse_kinematics(p, rotm)
% Inverse kinematics solution of a new 6-PUS platform, designed by Liu Jimu
    p = [p(1);p(2);p(3)];
    
    %%定义机构尺寸
    r1 = 100/1e3;       % 动平台竖直三角形分布半径
    r2 = 140/1e3;       % 动平台倾斜三角形分布半径
    h1 = -6/1e3;        % 动平台竖直三角形沿Z方向的偏置
    limb = 154/1e3;     % 连杆长度
    alpha = 30;         % 动平台竖直三角形与倾斜三角形的旋转角
    beta = 15;          % 倾斜三角形P副方向与xy平面的夹角
    d_l = 0.02;         % 连杆直径
    d_s = 0.02;         % 丝杠小径
    q_min = -0.1;       % 滑块的最小行程
    q_max = 0.1         % 滑块的最大行程
    %旋转矩阵
    RzB1  = [cosd(alpha) -sind(alpha) 0;
             sind(alpha) cosd(alpha)  0;
             0           0            1];
    Rz120 = [cosd(120) -sind(120) 0;
             sind(120) cosd(120)  0;
             0         0          1];
    RyP4 = [cosd(beta)  0 sind(beta);
            0           1          0;
            -sind(beta) 0 cosd(beta)];
    % 动平台S副位置
    B1_= [0 r1 h1]';
    B1 = RzB1 * B1_;
    B2 = Rz120 * B1;
    B3 = Rz120 * B2;
    B4 = [0 r2 0]';
    B5 = Rz120 * B4;
    B6 = Rz120 * B5;
    B = [B1, B2, B3, B4, B5, B6];
    % P副方向
    P1_dir = [0 0 1]';
    P2_dir = P1_dir;
    P3_dir = P1_dir;
    P4_dir_ = [-1 0 0]';
    P4_dir = RyP4 * P4_dir_;
    P5_dir = Rz120 * P4_dir;
    P6_dir = Rz120 * P5_dir;
    P_dir = [P1_dir, P2_dir, P3_dir, P4_dir, P5_dir, P6_dir];
    % U副初始位置
    A = -limb * P_dir + B;   
    % 运动后S副位置
    B_ = rotm * B + p;
    AB_ = B_ - A;
    A_ = zeros(size(A));
    for i = 1:6
        q(i) = AB_(:,i)'*P_dir(:,i) - sqrt( limb^2 - (AB_(:,i)'*AB_(:,i) - (AB_(:,i)'*P_dir(:,i))^2));
        A_(:,i) = A(:,i) + q(i) * P_dir(:,i);
    end
    
    %% 雅可比矩阵
    % 运动后连杆的向量
    l = B_ - A_;
    % 6-PUS的力雅可比
    l_dot_p = zeros(1,6);
    G_p = zeros(6);
    for i = 1:6
        l_dot_p(i) = l(:,i)'*P_dir(:,i);
        G_p(:,i) = [l(:,i);cross(rotm*B(:,i),l(:,i))]./l_dot_p(i);
    end
    % 6-UPS的力雅可比（连杆变形对末端位姿的影响）
    G_l = zeros(6);
    for i = 1:6
        G_l(:,i) = [l(:,i);cross(rotm*B(:,i),l(:,i))]./limb;
    end
    
    %% 刚度矩阵
    % 丝杠刚度
    K_s = zeros(6);
    for i = 1:6
        K_s(i,i) = limb_stiffness(d_s, q-q_min);
    end
    K_p = G_p*K_s*G_p';
    % 连杆变形对末端输出刚度的影响
    K_l = limb_stiffness(d_l,limb).*eye(6);
    K_d = G_l*K_l*G_l';
    % 总刚度
    K_k = inv(inv(K_p)+inv(K_d));
    
    %% 画图验证杆件位置
    for i=1:6
        plot3([B(1,i) A(1,i)],[B(2,i) A(2,i)],[B(3,i) A(3,i)]);
        hold on
    end
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal
        
end