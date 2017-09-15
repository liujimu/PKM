function jac = cal_vel_jacobian( p, rotm )
% 计算速度雅可比
    p = [p(1);p(2);p(3)];
    %%定义机构尺寸
    param = default_param();
    r1 = param.r1;          % 动平台竖直三角形分布半径
    r2 = param.r2;          % 动平台倾斜三角形分布半径
    h1 = param.h1;          % 动平台竖直三角形沿Z方向的偏置
    limb = param.limb;      % 连杆长度
    alpha = param.alpha;    % 动平台竖直三角形与倾斜三角形的旋转角
    beta = param.beta;      % 倾斜三角形P副方向与xy平面的夹角
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

    jac = zeros(6, 6);
    
end