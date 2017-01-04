function [PAP PB0 Dr lr U_Pb U_Bb] = Structure_6PSS(a,b,lr,angx)
%Structure_6PSS Summary of this function goes here
%   返回6PSS机构的各个铰点参数 驱动方向参数等
%   PAP - 动平台铰点布置
%   PB0 - 滑块铰点初始位置
%   Dr - 滑块驱动方向
%   lr - 杆长
%   U_Pb - 动平台坐标系下，动平台虎克铰基座固定轴轴线方向
%   U_Bb - 滑块虎克铰基座固定轴轴线方向
    
    %机构基本参数定义，单位取为m
    if nargin < 4 %当输入参数小于4个时，返回默认给定的参数
        a=450/1000;%数值杆铰点圆半径
        b=450/1000;%水平杆铰点圆半径
        lr=1200/1000;%连杆长度
        ang=150/180*pi;%水平杆铰点 布置偏角
        ang2=-0/180*pi;
    else
        ang=120/180*pi;
        ang2=-angx/180*pi;
    end

    
    c=0;%ang_0=0/180*pi;%初始模型中的尺寸参数，已废弃，置0
    h_ver = -19/1000;%动平台竖直三杆的虎克铰平面相对于动平台参考平面的位置
    
    %动平台上下三角形夹角，上三角形绕x轴转动ang角
    R_X=[1 0 0;0 cos(ang) -sin(ang);0 sin(ang) cos(ang);];%绕X轴旋转矩阵
    
    R_X2=[1 0 0;0 cos(ang2) -sin(ang2);0 sin(ang2) cos(ang2);];%绕X轴旋转矩阵
    %上三角形0位铰点布置
    PAPu=[0, b*cos(pi/3*2) - c*sin(pi/3*2), c*cos(pi/3*2) + b*sin(pi/3*2);
        0,b,c;
        0, b*cos(-pi/3*2) - c*sin(-pi/3*2), c*cos(-pi/3*2) + b*sin(-pi/3*2);
        ].';
    %下三角形铰点布置
    PAPd=[0 + h_ver,-a/2,-a*3^0.5/2;
        0 + h_ver,-a/2,a*3^0.5/2;
        0 + h_ver,a,0].';
    %动平台铰点布置
    PAP=[R_X*PAPu R_X2*PAPd];

    

    %滑块铰点初始位置（0姿态时）
    PBu=[0, b*cos(pi/3*2) - (c+lr)*sin(pi/3*2), (c+lr)*cos(pi/3*2) + b*sin(pi/3*2);
        0,b,c+lr;
        0, b*cos(-pi/3*2) - (c+lr)*sin(-pi/3*2), (c+lr)*cos(-pi/3*2) + b*sin(-pi/3*2);
        ].';
    PBd=[
        lr+ h_ver,-a/2,-a*3^0.5/2
        lr+ h_ver,-a/2,a*3^0.5/2
        lr+ h_ver,a,0].';
    PB0=[R_X*PBu R_X2*PBd];

    %滑块驱动方向  远离动平台中心为正
    vec_p=[0 -sin(pi/3*2) cos(pi/3*2)
    0 0 1 
    0 -sin(-pi/3*2) cos(-pi/3*2)]';
    vec_p=R_X*vec_p;
    Dr=[vec_p.';[1 0 0];[1 0 0];[1 0 0];]';

    ang_U=ang2+30/180*pi;
    R_XU=[1 0 0;0 cos(ang_U) -sin(ang_U);0 sin(ang_U) cos(ang_U);];%绕X轴旋转矩阵
    %动平台坐标系下，动平台虎克铰基座固定轴轴线方向
    U_Pb_u=[ 
    1 0 0;
    1 0 0;
    1 0 0;]';
    U_Pb_d=[ 
    0 sin(pi/6) cos(pi/6);
    0 sin(-pi/6) cos(-pi/6);
    0 -1 0;]';
    U_Pb=[U_Pb_u R_XU*U_Pb_d];
    %滑块虎克铰基座固定轴轴线方向
    U_Bb_u=[ 
    1 0 0;
    1 0 0;
    1 0 0;]';
    U_Bb_d=[ 
    0 sin(pi/6) cos(pi/6);
    0 sin(-pi/6) cos(-pi/6);
    0 -1 0;]';
    U_Bb=[U_Bb_u R_XU*U_Bb_d];
end

