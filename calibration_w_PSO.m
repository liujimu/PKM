% �����궨
clear;
pkm = PKM();
target_poses = xlsread('target_poses.xlsx');
measured_poses = xlsread('measured_poses.xlsx');
% target_poses(:,6) = [];
% measured_poses(:,6) = [];

n = length(target_poses);
dX = reshape(measured_poses - target_poses,[],1);
updated_poses = zeros(6,n);
qin = zeros(6,n);
param_errors = zeros(60,1);
accumulated_param_errors = zeros(60,1);
delta = 1;

%------������ʼ������----------------------------------------------
c1=2;             %ѧϰ����1
c2=2;             %ѧϰ����2
w=0.7;            %����Ȩ��
MaxDT=1000;       %����������
D=60;             %�����ռ�ά����δ֪��������
M=40;             %��ʼ��Ⱥ�������Ŀ
eps=10e-6;      %���þ���(����֪��Сֵʱ����)


%------��ʼ����Ⱥ�ĸ���(�����������޶�λ�ú��ٶȵķ�Χ)------------
for i=1:M
    for j=1:D
        x(i,j)=randn; %�����ʼ��λ��
        v(i,j)=randn; %�����ʼ���ٶ�
    end
end

%------�ȼ���������ӵ���Ӧ�ȣ�����ʼ��p(i)��gbest--------------------
for i=1:M
    p(i)=f1(x(i,:),D);
    y(i,:)=x(i,:);
end
gbest=x(1,:);             %gbestΪȫ������

for i=2:M
    if f1(x(i,:),D) < f1(gbest,D)
        gbest=x(i,:);
    end
end

%------������Ҫѭ�������չ�ʽ���ε�����ֱ�����㾫��Ҫ��------------
for t=1:MaxDT
    for i=1:M
        v(i,:)=w*v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(gbest-x(i,:));
        x(i,:)=x(i,:)+v(i,:);
        if f1(x(i,:),D)<p(i)
            p(i)=f1(x(i,:),D);
            y(i,:)=x(i,:);
        end
        if p(i)<f1(gbest,D)
            gbest=y(i,:);
        end
    end
end


k = 0;%��������
while delta > alw
    W = zeros(6*n,60);
    for i = 1:n
        if k == 0
            pkm.pose = target_poses(:,i);
            qin(:,i) = pkm.q;
        else
            pkm.pose = updated_poses(:,i);
        end
        pkm.calVelJac();
        % �жϲ��Ե��Ƿ��ڹ����ռ��ڣ��Ƿ��������
        %     if ~pkm.isInWorkspace()
        %         disp(i);
        %         disp(pkm.q);
        %         msg = 'Target pose is out of the workspace.';
        %         error(msg);
        %     end
        %     if ~pkm.isSingular()
        %         disp(i);
        %         disp(det(pkm.jac));
        %         msg = 'Target pose is a singular pose.';
        %         error(msg);
        %     end
%         Jx = [pkm.l_dir; cross(pkm.rotm*pkm.S_init,pkm.l_dir)]';
%         Jm = zeros(6,54);
%         for j = 1:6
%             % Jm(j,(j*9-8):j*9) = [-pkm.l_dir(:,j)', pkm.l_dir(:,j)'*pkm.rotm, -pkm.q(j)*pkm.l_dir(:,j)'];
%             Jm(j,(3*j-2):(3*j)) = -pkm.l_dir(:,j)';
%             Jm(j,(18+3*j-2):(18+3*j)) = pkm.l_dir(:,j)'*pkm.rotm;
%             Jm(j,(36+3*j-2):(36+3*j)) = -pkm.q(j)*pkm.l_dir(:,j)';
%         end
%         J = [inv(Jx), -inv(Jx)*Jm];
%         W((6*i-5):6*i,:) = J;

        % ��ı궨��ʽ
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
        W((6*i-5):6*i,:) = E;
    end
%     disp(det(W'*W));
    param_errors = W\dX;
    accumulated_param_errors = accumulated_param_errors + param_errors;
    pkm = PKM(accumulated_param_errors);
    for i = 1:n
        pkm.forKin(qin(:,i), target_poses(:,i));
        updated_poses(:,i) = pkm.pose;
    end
    dX = reshape(measured_poses - updated_poses,[],1);
    delta = max(abs(dX));

    k = k + 1;
    fprintf('In the %dth iteration, delta = %f\n', k, delta);
    if k > 10
        break;
    end
end
