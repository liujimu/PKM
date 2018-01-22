% �Ա���ֵ����������ݾ�������ʽ��������ݾ���
clear

pkm = PKM();
global target_pose target_input
target_pose = [0 0 100 0 0 0]';
pkm.setPose(target_pose);
target_input = pkm.q;
pkm.calVelJac();

%% ��ı궨��ʽ
% �˳����
E_l = [pkm.l_dir; cross(pkm.rotm * pkm.S_init, (pkm.S_cur - pkm.U_cur))]';
% U����ʼλ�����
E_u = zeros(6,18);
for j = 1:6
    E_u(j,(3*j-2):3*j) = pkm.l_dir(:,j)' / (pkm.l_dir(:,j)' * pkm.P_dir(:,j));
end
E_u = pkm.jac * E_u;
% S����ʼλ�����
E_s = zeros(6,18);
for j = 1:6
    E_s(j,(3*j-2):3*j) = pkm.l_dir(:,j)' * pkm.rotm / (pkm.l_dir(:,j)' * pkm.P_dir(:,j));
end
E_s = pkm.jac * E_s;
% P���������
E_p = zeros(6,18);
for j = 1:6
    H = pkm.pose(1:3) + pkm.rotm * pkm.S_init(:,j) - pkm.U_init(:,j);
    E_p(j,(3*j-2):3*j) = pkm.q(j) * H' / (pkm.l(j) * pkm.l_dir(:,j)' * pkm.P_dir(:,j));
end
E_p = -pkm.jac * E_p;

E = [inv(E_l) E_u E_s E_p];

%% ��ֵ����������ݾ���
param_error = zeros(60,1);
error_jac = cal_numerical_jacobian(@error_transfer_matrix,param_error);

dE = error_jac - E;

function pose_error = error_transfer_matrix(param_error)
    global target_pose target_input
    pkm = PKM(param_error);
    pkm.forKin(target_input,target_pose);
    pose_error = pkm.pose - target_pose;
end