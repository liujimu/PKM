function [PAP PB0 Dr lr U_Pb U_Bb] = Structure_6PSS(a,b,lr,angx)
%Structure_6PSS Summary of this function goes here
%   ����6PSS�����ĸ����µ���� �������������
%   PAP - ��ƽ̨�µ㲼��
%   PB0 - ����µ��ʼλ��
%   Dr - ������������
%   lr - �˳�
%   U_Pb - ��ƽ̨����ϵ�£���ƽ̨���˽»����̶������߷���
%   U_Bb - ���黢�˽»����̶������߷���
    
    %���������������壬��λȡΪm
    if nargin < 4 %���������С��4��ʱ������Ĭ�ϸ����Ĳ���
        a=450/1000;%��ֵ�˽µ�Բ�뾶
        b=450/1000;%ˮƽ�˽µ�Բ�뾶
        lr=1200/1000;%���˳���
        ang=150/180*pi;%ˮƽ�˽µ� ����ƫ��
        ang2=-0/180*pi;
    else
        ang=120/180*pi;
        ang2=-angx/180*pi;
    end

    
    c=0;%ang_0=0/180*pi;%��ʼģ���еĳߴ�������ѷ�������0
    h_ver = -19/1000;%��ƽ̨��ֱ���˵Ļ��˽�ƽ������ڶ�ƽ̨�ο�ƽ���λ��
    
    %��ƽ̨���������μнǣ�����������x��ת��ang��
    R_X=[1 0 0;0 cos(ang) -sin(ang);0 sin(ang) cos(ang);];%��X����ת����
    
    R_X2=[1 0 0;0 cos(ang2) -sin(ang2);0 sin(ang2) cos(ang2);];%��X����ת����
    %��������0λ�µ㲼��
    PAPu=[0, b*cos(pi/3*2) - c*sin(pi/3*2), c*cos(pi/3*2) + b*sin(pi/3*2);
        0,b,c;
        0, b*cos(-pi/3*2) - c*sin(-pi/3*2), c*cos(-pi/3*2) + b*sin(-pi/3*2);
        ].';
    %�������νµ㲼��
    PAPd=[0 + h_ver,-a/2,-a*3^0.5/2;
        0 + h_ver,-a/2,a*3^0.5/2;
        0 + h_ver,a,0].';
    %��ƽ̨�µ㲼��
    PAP=[R_X*PAPu R_X2*PAPd];

    

    %����µ��ʼλ�ã�0��̬ʱ��
    PBu=[0, b*cos(pi/3*2) - (c+lr)*sin(pi/3*2), (c+lr)*cos(pi/3*2) + b*sin(pi/3*2);
        0,b,c+lr;
        0, b*cos(-pi/3*2) - (c+lr)*sin(-pi/3*2), (c+lr)*cos(-pi/3*2) + b*sin(-pi/3*2);
        ].';
    PBd=[
        lr+ h_ver,-a/2,-a*3^0.5/2
        lr+ h_ver,-a/2,a*3^0.5/2
        lr+ h_ver,a,0].';
    PB0=[R_X*PBu R_X2*PBd];

    %������������  Զ�붯ƽ̨����Ϊ��
    vec_p=[0 -sin(pi/3*2) cos(pi/3*2)
    0 0 1 
    0 -sin(-pi/3*2) cos(-pi/3*2)]';
    vec_p=R_X*vec_p;
    Dr=[vec_p.';[1 0 0];[1 0 0];[1 0 0];]';

    ang_U=ang2+30/180*pi;
    R_XU=[1 0 0;0 cos(ang_U) -sin(ang_U);0 sin(ang_U) cos(ang_U);];%��X����ת����
    %��ƽ̨����ϵ�£���ƽ̨���˽»����̶������߷���
    U_Pb_u=[ 
    1 0 0;
    1 0 0;
    1 0 0;]';
    U_Pb_d=[ 
    0 sin(pi/6) cos(pi/6);
    0 sin(-pi/6) cos(-pi/6);
    0 -1 0;]';
    U_Pb=[U_Pb_u R_XU*U_Pb_d];
    %���黢�˽»����̶������߷���
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

