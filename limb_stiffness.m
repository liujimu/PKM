function k = limb_stiffness( d, l )
%�������˸ն�
%   kΪ���˸ն�
%   dΪ����ֱ��
%   lΪ���˳���
    E = 2.1e11; % ��̼ͨ�ֵĵ���ģ��
    A = pi/4*d^2;
    k = E*A/l;
end

