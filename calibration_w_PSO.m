% 机构标定
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

%------给定初始化条件----------------------------------------------
c1=2;             %学习因子1
c2=2;             %学习因子2
w=0.7;            %惯性权重
MaxDT=1000;       %最大迭代次数
D=60;             %搜索空间维数（未知数个数）
M=40;             %初始化群体个体数目
eps=10e-6;      %设置精度(在已知最小值时候用)


%------初始化种群的个体(可以在这里限定位置和速度的范围)------------
for i=1:M
    for j=1:D
        x(i,j)=randn; %随机初始化位置
        v(i,j)=randn; %随机初始化速度
    end
end

%------先计算各个粒子的适应度，并初始化p(i)和gbest--------------------
for i=1:M
    p(i)=f1(x(i,:),D);
    y(i,:)=x(i,:);
end
gbest=x(1,:);             %gbest为全局最优

for i=2:M
    if f1(x(i,:),D) < f1(gbest,D)
        gbest=x(i,:);
    end
end

%------进入主要循环，按照公式依次迭代，直到满足精度要求------------
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


k = 0;%迭代次数
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
        % 判断测试点是否在工作空间内，是否是奇异点
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

        % 曹睿的标定公式
        % 杆长误差
        E_l = [pkm.l_dir; cross(pkm.rotm * pkm.S_init, (pkm.S_cur - pkm.U_cur))]';
        % U副初始位置误差
        E_u = zeros(6,18);
        for j = 1:6
            E_u(j,(3*j-2):3*j) = pkm.l_dir(:,j)' / (pkm.l_dir(:,j)' * pkm.P_dir(:,j));
        end
        E_u = pkm.jac * E_u;
        % S副初始位置误差
        E_s = zeros(6,18);
        for j = 1:6
            E_s(j,(3*j-2):3*j) = pkm.l_dir(:,j)' * pkm.rotm / (pkm.l_dir(:,j)' * pkm.P_dir(:,j));
        end
        E_s = pkm.jac * E_s;
        % P副方向误差
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
