function [q, J] = InverseSolution_MOD(pu, angu)
%% Mod for newly designed 805 project 6-PSS robot
   %% geometry parameter initialization

    pu = [pu(1);pu(2);pu(3)];
    a   = 100/1e3;    % vertical limbs' upper joints radius
    b   = 140/1e3;    % horizontal limbs' upper joints radius
    hv  = 6/1e3;      % vertical limbs' upper joints offset along X axis
    lr  = 154/1e3;    % limb length
    ang =150/180*pi;  % angle between horizontal limb's upper joint and vertical limb's upper joint
    l   = [lr lr lr lr lr lr];
    
    Rx120=[1 0         0;
           0 cosd(120) -sind(120);
           0 sind(120) cosd(120)];
    Rx240=[1 0         0;
           0 cosd(240) -sind(240);
           0 sind(240) cosd(240)];
    
    %动平台水平/垂直三角形夹角，水平三角形绕x轴转动ang角
    RX =[1        0         0;
         0 cos(ang) -sin(ang);
         0 sin(ang) cos(ang);];%绕X轴旋转矩阵
    
    %水平杆倾斜角
    beta=pi/180*45;
    RY =[cos(beta)  0 sin(beta);
         0          1         0;
         -sin(beta) 0 cos(beta);];%绕Y轴旋转矩阵
     
    %水平三角形0位铰点布置
     A2=[0 b 0]';
     Auh = [Rx120*A2, A2, Rx240*A2];
        
    %垂直三角形铰点布置
    A6=[hv a 0]';
    Auv=[Rx240*A6, Rx120*A6, A6];
      
    %动平台铰点布置
    Au = [RX * Auh Auv];

    %滑块铰点初始位置
    B2=[0 b lr]';
    B2_=RY*B2;
    Buh = [Rx120*B2_, B2_, Rx240*B2_];
     
    Buv = [hv + lr,  a*cos(-pi/3*2),  a*sin(-pi/3*2);
           hv + lr,  a*cos( pi/3*2),  a*sin( pi/3*2);
           hv + lr,               a,               0;]';

    Bu = [RX * Buh Buv];

    %滑块驱动方向
    Dr2 = RY*[0; 0; 1];
    Drh = [Rx120*Dr2, Dr2, Rx240*Dr2];
    Drh = RX * Drh;

    Dr=[ Drh';
        [1 0 0];
        [1 0 0];
        [1 0 0];]';



    %% configuration

    q = zeros(1, 6);

    %% position and orientation of upper platform

    alpha = angu(3); beta = angu(2); gamma = angu(1);

    %  RPY(XYZ) fixed angle, gamma on x, then beta on y, alpha on z in the end;
    Ru = [ cos(alpha)*cos(beta), cos(alpha)*sin(beta)*sin(gamma) - cos(gamma)*sin(alpha), sin(alpha)*sin(gamma) + cos(alpha)*cos(gamma)*sin(beta);
           cos(beta)*sin(alpha), cos(alpha)*cos(gamma) + sin(alpha)*sin(beta)*sin(gamma), cos(gamma)*sin(alpha)*sin(beta) - cos(alpha)*sin(gamma);
                     -sin(beta),                                    cos(beta)*sin(gamma),                                    cos(beta)*cos(gamma)];



    %% position of A w.r.t. upper platform base and lower platform base
    Al = zeros(3, 6);
    RrA = Ru * Au;
    for i = 1:6
        Al(:, i) = RrA(:, i) + pu;
    end

    %% calculate inverse solution, update Bl
    H = Al - Bu;
    Bl = zeros(3, 6);
    for i = 1:6
        q(i) = H(:,i)'*Dr(:,i)+ sqrt( (H(:,i)'*Dr(:,i))^2 - H(:,i)'*H(:,i) + lr^2 );
        Bl(:,i) = Bu(:,i) + q(i) * Dr(:,i);
    end

    %% calculate limb unit vector, n
    n = Al - Bl;
    for i = 1:6
        n(:, i) = n(:, i) / norm(n(:, i));   % normalize
    end
    
    
    %% calculate limb kinetics
    J = zeros(6, 6);

    
    for i = 1:6
        
        ni = n(:, i);                               %limb unit vector
        nte = ni' * Dr(:, i);
        qA = RrA(:, i);                             % Represent Ru * Au(i)
        
        J(i, :) = [ni'/nte (cross(qA, ni)/nte)'];   % jacobian
    end
end