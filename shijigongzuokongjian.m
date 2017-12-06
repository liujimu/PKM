function [  ] = shijigongzuokongjian( )
clc;
clear;
%给定碰撞初始条件，求解模拟器实际所需工作空间

P0=[0 0 0 1]';
rad=230; %给定坐标系工作空间半径
in=0;ina=0;
ia=5; ib=5;ic=5;
Rrx=linspace(-6,6,ia)*pi/180;
Rry=linspace(-6,6,ib)*pi/180;
Rrz=linspace(-6,6,ic)*pi/180;
for ang=0:pi/180:2*pi
    for h=0:170/20:170 %给定工作空间高度
        A=[-h;cos(ang)*rad;sin(ang)*rad;];
        ina=ina+1;
        B0(:,ina)=-A;
        for i=1:ia
            for j=1:ib
                for k=1:ic
                    %pth=[6 6 6]*pi/180;
                    %RPY角
                    %pth=[Rrx(i) Rry(j) Rrz(k)];
                    %Rrpy=[cos(pth(3))*cos(pth(2)),cos(pth(3))*sin(pth(2))*sin(pth(1))-sin(pth(3))*cos(pth(1)),cos(pth(3))*sin(pth(2))*cos(pth(1))+sin(pth(3))*sin(pth(1));
                    %      sin(pth(3))*cos(pth(2)),sin(pth(3))*sin(pth(2))*sin(pth(1))+cos(pth(3))*cos(pth(1)),sin(pth(3))*sin(pth(2))*cos(pth(1))-cos(pth(3))*sin(pth(1));
                    %      -sin(pth(2)),cos(pth(2))*sin(pth(1)),cos(pth(2))*cos(pth(1))];
                    
                    %312欧拉角
                    rx=Rrx(i);ry=Rry(j);rz=Rrz(k);
                    Rrpy=[ cos(ry)*cos(rz) - sin(rx)*sin(ry)*sin(rz), -cos(rx)*sin(rz), cos(rz)*sin(ry) + cos(ry)*sin(rx)*sin(rz);
                        cos(ry)*sin(rz) + cos(rz)*sin(rx)*sin(ry),  cos(rx)*cos(rz), sin(ry)*sin(rz) - cos(ry)*cos(rz)*sin(rx);
                        -cos(rx)*sin(ry),          sin(rx),                           cos(rx)*cos(ry)];
                    
                    T=[Rrpy A; 0 0 0 1];
                    in=in+1;
                    B(:,in)=DHinv(T)*P0;
                end
            end
        end
    end
end

max(abs(B([1,2,3],:)'))
min((B([1,2,3],:)'))

figure
view(138,26);
%view(-90,-66);
hold on
%plot3(B0(1,:),B0(2,:),B0(3,:),'k')
plot3(B(2,:),B(3,:),-B(1,:),'.b')
%plot3(B2(1,:),B2(2,:),B2(3,:),'r')
axis equal;



% RPY角度偏差变换
% in=0;
% ia=36; ib=36;ic=36;
% Rrx=linspace(-6,6,ia)*pi/180;
% Rry=linspace(-6,6,ib)*pi/180;
% Rrz=linspace(-6,6,ic)*pi/180;
% for i=1:ia
%     for j=1:ib
%         for k=1:ic
%             %pth=[6 6 6]*pi/180;
%             pth=[Rrx(i) Rry(j) Rrz(k)];
%             Rrpy=[cos(pth(3))*cos(pth(2)),cos(pth(3))*sin(pth(2))*sin(pth(1))-sin(pth(3))*cos(pth(1)),cos(pth(3))*sin(pth(2))*cos(pth(1))+sin(pth(3))*sin(pth(1));
%                sin(pth(3))*cos(pth(2)),sin(pth(3))*sin(pth(2))*sin(pth(1))+cos(pth(3))*cos(pth(1)),sin(pth(3))*sin(pth(2))*cos(pth(1))-cos(pth(3))*sin(pth(1));
%                -sin(pth(2)),cos(pth(2))*sin(pth(1)),cos(pth(2))*cos(pth(1))];
%            in=in+1;
%            Ang(in,:)=mattoRPY( Rrpy' )/pi*180;
%         end
%     end
% end
% max(Ang)
end

function [ Ai ] = DHinv(A)
%DASTENCE Summary of this function goes here
%   Detailed explanation goes here
% 本函数用于计算机械臂坐标变换矩阵的逆矩阵

% 计算两点间距离  取一个矩阵其中的第1,3行B=A([1,3],:);取第1,3行和第1,5列相交的数据:C=A([1,3],[1,5]);
R=A([1,2,3],[1,2,3]);
P=A([1,2,3],4);

Ai=[R.' -R.'*P; %注意，.'表示的转置中不会出项conj(x)复数共轭
    0 0 0 1];
end