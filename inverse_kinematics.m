function [q, J] = inverse_kinematics(p, rotm)
% Inverse kinematics solution of a new 6-PUS platform, designed by Liu Jimu
    p = [p(1);p(2);p(3)];
    
    %%��������ߴ�
    r1 = 100/1e3;       % ��ƽ̨��ֱ�����ηֲ��뾶
    r2 = 140/1e3;       % ��ƽ̨��б�����ηֲ��뾶
    h1 = -6/1e3;        % ��ƽ̨��ֱ��������Z�����ƫ��
    limb = 154/1e3;     % ���˳���
    alpha = 30;         % ��ƽ̨��ֱ����������б�����ε���ת��
    beta = 15;          % ��б������P��������xyƽ��ļн�
    d_l = 0.02;         % ����ֱ��
    d_s = 0.02;         % ˿��С��
    q_min = -0.1;       % �������С�г�
    q_max = 0.1         % ���������г�
    %��ת����
    RzB1  = [cosd(alpha) -sind(alpha) 0;
             sind(alpha) cosd(alpha)  0;
             0           0            1];
    Rz120 = [cosd(120) -sind(120) 0;
             sind(120) cosd(120)  0;
             0         0          1];
    RyP4 = [cosd(beta)  0 sind(beta);
            0           1          0;
            -sind(beta) 0 cosd(beta)];
    % ��ƽ̨S��λ��
    B1_= [0 r1 h1]';
    B1 = RzB1 * B1_;
    B2 = Rz120 * B1;
    B3 = Rz120 * B2;
    B4 = [0 r2 0]';
    B5 = Rz120 * B4;
    B6 = Rz120 * B5;
    B = [B1, B2, B3, B4, B5, B6];
    % P������
    P1_dir = [0 0 1]';
    P2_dir = P1_dir;
    P3_dir = P1_dir;
    P4_dir_ = [-1 0 0]';
    P4_dir = RyP4 * P4_dir_;
    P5_dir = Rz120 * P4_dir;
    P6_dir = Rz120 * P5_dir;
    P_dir = [P1_dir, P2_dir, P3_dir, P4_dir, P5_dir, P6_dir];
    % U����ʼλ��
    A = -limb * P_dir + B;   
    % �˶���S��λ��
    B_ = rotm * B + p;
    AB_ = B_ - A;
    A_ = zeros(size(A));
    for i = 1:6
        q(i) = AB_(:,i)'*P_dir(:,i) - sqrt( limb^2 - (AB_(:,i)'*AB_(:,i) - (AB_(:,i)'*P_dir(:,i))^2));
        A_(:,i) = A(:,i) + q(i) * P_dir(:,i);
    end
    
    %% �ſɱȾ���
    % �˶������˵�����
    l = B_ - A_;
    % 6-PUS�����ſɱ�
    l_dot_p = zeros(1,6);
    G_p = zeros(6);
    for i = 1:6
        l_dot_p(i) = l(:,i)'*P_dir(:,i);
        G_p(:,i) = [l(:,i);cross(rotm*B(:,i),l(:,i))]./l_dot_p(i);
    end
    % 6-UPS�����ſɱȣ����˱��ζ�ĩ��λ�˵�Ӱ�죩
    G_l = zeros(6);
    for i = 1:6
        G_l(:,i) = [l(:,i);cross(rotm*B(:,i),l(:,i))]./limb;
    end
    
    %% �նȾ���
    % ˿�ܸն�
    K_s = zeros(6);
    for i = 1:6
        K_s(i,i) = limb_stiffness(d_s, q-q_min);
    end
    K_p = G_p*K_s*G_p';
    % ���˱��ζ�ĩ������նȵ�Ӱ��
    K_l = limb_stiffness(d_l,limb).*eye(6);
    K_d = G_l*K_l*G_l';
    % �ܸն�
    K_k = inv(inv(K_p)+inv(K_d));
    
    %% ��ͼ��֤�˼�λ��
    for i=1:6
        plot3([B(1,i) A(1,i)],[B(2,i) A(2,i)],[B(3,i) A(3,i)]);
        hold on
    end
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal
        
end