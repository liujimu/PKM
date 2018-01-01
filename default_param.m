function param = default_param()
%����6-PUS�����ĳߴ����
r1 = 100/1e3;       % ��ƽ̨��ֱ�����ηֲ��뾶
r2 = 140/1e3;       % ��ƽ̨��б�����ηֲ��뾶
h1 = -6/1e3;        % ��ƽ̨��ֱ��������Z�����ƫ��
limb = 154/1e3;     % ���˳���
alpha = 30;         % ��ƽ̨��ֱ����������б�����ε���ת��
beta = 0;          % ��б������P��������xyƽ��ļн�
d_l = 0.02;         % ����ֱ��
d_s = 0.02;         % ˿��С��
q_min = -0.2;       % �������С�г�
q_max = 0.2;        % ���������г�
param = struct('r1',r1,'r2',r2,'h1',h1,'limb',limb,'alpha',alpha,'beta',beta,'d_l',d_l,'d_s',d_s,'q_min',q_min,'q_max',q_max);
end

