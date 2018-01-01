% �����궨
clear;
global target_poses measured_poses qin
pkm = PKM();
target_poses = xlsread('target_poses.xlsx'); %Ŀ��λ��
measured_poses = xlsread('measured_poses.xlsx'); %���������ʵ��λ��
n = length(target_poses);
qin = zeros(6,n); %����
for i = 1:n
    pkm.pose = target_poses(:,i);
    qin(:,i) = pkm.q;
end
calculated_poses = zeros(6,n); %���������λ��

param_errors = zeros(60,1);
delta = 1;
dX = reshape(measured_poses - target_poses,[],1);

%------������ʼ������----------------------------------------------
c1=2;             %ѧϰ����1
c2=2;             %ѧϰ����2
w=0.7;            %����Ȩ��
MaxDT=1000;       %����������
D=60;             %�����ռ�ά����δ֪��������
M=100;            %��ʼ��Ⱥ�������Ŀ
eps=10e-6;        %���þ���(����֪��Сֵʱ����)


%------��ʼ����Ⱥ�ĸ���(�����������޶�λ�ú��ٶȵķ�Χ)------------
xrange = [0.001*ones(1,6) 0.005*ones(1,18) 0.001*ones(1,18) 0.005*ones(1,18)];

for i=1:M
    for j=1:D
        x(i,j)=randn * xrange(j); %�����ʼ��λ��
        v(i,j)=randn * xrange(j); %�����ʼ���ٶ�
    end
end

%------�ȼ���������ӵ���Ӧ�ȣ�����ʼ��p(i)��gbest--------------------
for i=1:M
    % ����ÿ�����ӣ��������и��º��ģ�͵�����λ����ʵ��λ�˵Ĳ�ֵ��ȡ�����ֵ��Ϊ��Ӧ�Ⱥ���
    p(i) = fitness_PSO(x(i,:));
    pbest(i,:) = x(i,:);
end

gbest=x(1,:);             %gbestΪȫ������
for i=2:M
    if fitness_PSO(x(i,:)) < fitness_PSO(gbest)
        gbest=x(i,:);
    end
end

%------������Ҫѭ�������չ�ʽ���ε�����ֱ�����㾫��Ҫ��------------
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

