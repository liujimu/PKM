function [q, J] = inverse_kinematics(pu, rotm)
% Inverse kinematics solution of a new 6-PUS platform, designed by Liu Jimu
    pu = [pu(1);pu(2);pu(3)];
    
    %%��������ߴ�
    r1 = 100/1e3;       % ��ƽ̨��ֱ�����ηֲ��뾶
    r2 = 140/1e3;       % ��ƽ̨��б�����ηֲ��뾶
    h1 = -6/1e3;        % ��ƽ̨��ֱ��������Z�����ƫ��
    lr = 154/1e3;       % ���˳���
    alpha = 30;         % ��ƽ̨��ֱ����������б�����ε���ת��
    %��ת����
    RzB1  = [cosd(alpha) -sind(alpha) 0;
             sind(alpha) cosd(alpha)  0;
             0           0            1];
    Rz120 = [cosd(120) -sind(120) 0;
             sind(120) cosd(120)  0;
             0         0          1];
    B1_= [0 r1 h1]';
    B1 = RzB1*B1_;
    B2 = Rz120*B1;
    B3 = Rz120*B2;
    B4 = [0 r2 0]';
    B5 = Rz120*B4;
    B6 = Rz120*B5;
    B = [B1,B2,B3,B4,B5,B6];
    
    
    
    
    %��ƽ̨ˮƽ/��ֱ�����μнǣ�ˮƽ��������x��ת��ang��
    RX =[1        0         0;
         0 cos(ang) -sin(ang);
         0 sin(ang) cos(ang);];%��X����ת����
    
    %ˮƽ����б��
    beta=pi/180*45;
    RY =[cos(beta)  0 sin(beta);
         0          1         0;
         -sin(beta) 0 cos(beta);];%��Y����ת����
     
    %ˮƽ������0λ�µ㲼��
     A2=[0 r2 0]';
     Auh = [Rz120*A2, A2, Rz240*A2];
        
    %��ֱ�����νµ㲼��
    A6=[h1 r1 0]';
    Auv=[Rz240*A6, Rz120*A6, A6];
      
    %��ƽ̨�µ㲼��
    Au = [RX * Auh Auv];

    %����µ��ʼλ��
    B2=[0 r2 lr]';
    B2_=RY*B2;
    Buh = [Rz120*B2_, B2_, Rz240*B2_];
     
    Buv = [h1 + lr,  r1*cos(-pi/3*2),  r1*sin(-pi/3*2);
           h1 + lr,  r1*cos( pi/3*2),  r1*sin( pi/3*2);
           h1 + lr,               r1,               0;]';

    Bu = [RX * Buh Buv];

    %������������
    Dr2 = RY*[0; 0; 1];
    Drh = [Rz120*Dr2, Dr2, Rz240*Dr2];
    Drh = RX * Drh;

    Dr=[ Drh';
        [1 0 0];
        [1 0 0];
        [1 0 0];]';

    %% configuration

    q = zeros(1, 6);

    %% position of A w.r.t. upper platform base and lower platform base
    Al = zeros(3, 6);
    RrA = rotm * Au;
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