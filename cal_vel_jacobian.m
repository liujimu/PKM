function jac = cal_vel_jacobian( p, rotm )
% �����ٶ��ſɱ�
    p = [p(1);p(2);p(3)];
    %%��������ߴ�
    param = default_param();
    r1 = param.r1;          % ��ƽ̨��ֱ�����ηֲ��뾶
    r2 = param.r2;          % ��ƽ̨��б�����ηֲ��뾶
    h1 = param.h1;          % ��ƽ̨��ֱ��������Z�����ƫ��
    limb = param.limb;      % ���˳���
    alpha = param.alpha;    % ��ƽ̨��ֱ����������б�����ε���ת��
    beta = param.beta;      % ��б������P��������xyƽ��ļн�
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

    jac = zeros(6, 6);
    
end