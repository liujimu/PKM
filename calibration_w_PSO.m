% 机构标定
clear;
global target_poses measured_poses qin
pkm = PKM();
target_poses = xlsread('target_poses.xlsx'); %目标位姿
measured_poses = xlsread('measured_poses.xlsx'); %激光跟踪仪实测位姿
n = length(target_poses);
qin = zeros(6,n); %驱动
for i = 1:n
    pkm.pose = target_poses(:,i);
    qin(:,i) = pkm.q;
end
calculated_poses = zeros(6,n); %正解算出的位姿

param_errors = zeros(60,1);
delta = 1;
dX = reshape(measured_poses - target_poses,[],1);

%------给定初始化条件----------------------------------------------
c1=2;             %学习因子1
c2=2;             %学习因子2
w=0.7;            %惯性权重
MaxDT=1000;       %最大迭代次数
D=60;             %搜索空间维数（未知数个数）
M=100;            %初始化群体个体数目
eps=10e-6;        %设置精度(在已知最小值时候用)


%------初始化种群的个体(可以在这里限定位置和速度的范围)------------
xrange = [0.001*ones(1,6) 0.005*ones(1,18) 0.001*ones(1,18) 0.005*ones(1,18)];

for i=1:M
    for j=1:D
        x(i,j)=randn * xrange(j); %随机初始化位置
        v(i,j)=randn * xrange(j); %随机初始化速度
    end
end

%------先计算各个粒子的适应度，并初始化p(i)和gbest--------------------
for i=1:M
    % 对于每个粒子，计算所有更新后的模型的正解位姿与实测位姿的差值，取其最大值作为适应度函数
    p(i) = fitness_PSO(x(i,:));
    pbest(i,:) = x(i,:);
end

gbest=x(1,:);             %gbest为全局最优
for i=2:M
    if fitness_PSO(x(i,:)) < fitness_PSO(gbest)
        gbest=x(i,:);
    end
end

%------进入主要循环，按照公式依次迭代，直到满足精度要求------------
for t=1:MaxDT
    for i=1:M
        v(i,:) = w*v(i,:) + c1*rand*(pbest(i,:) - x(i,:)) + c2*rand*(gbest - x(i,:));
        x(i,:) = x(i,:) + v(i,:);
        if fitness_PSO(x(i,:)) < p(i)
            p(i) = fitness_PSO(x(i,:));
            pbest(i,:) = x(i,:);
        end
        if p(i) < fitness_PSO(gbest)
            gbest = pbest(i,:);
        end
    end
end

