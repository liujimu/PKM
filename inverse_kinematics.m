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