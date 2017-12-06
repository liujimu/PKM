%��֤PKM�е��ſɱȹ�ʽ�Ƿ���ȷ
clear;
pkm = PKM();
pose = [0 0 0 0 0 0]';
pkm.pose = pose;
pkm.calVelJac();
%��ʽ�Ƶ����ſɱ�
jac = pkm.jac;
h = 1e-8;
h12 = 1.2e-7;
dydx = zeros(6);
b = h * eye(6);
for i = 1:6
    br = b(:,i);
    pkm.pose = pose - 2 * br;
    y1 = pkm.q;
    pkm.pose = pose - br;
    y2 = pkm.q;
    pkm.pose = pose + br;
    y4 = pkm.q;
    pkm.pose = pose + 2 * br;
    y5 = pkm.q;
    dydx(:,i) = (8 * y4 - 8 * y2 + y1 - y5)' ./ h12;
end
%��ֵ������õ��ſɱȵ���
inv_jac = dydx;
delta = jac*inv_jac - eye(6)