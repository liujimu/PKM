clear;
% close all
%%参数初始化
qh_max=0.0435;
qh_min=qh_max-0.13;
qv_max=0.0375;
qv_min=qv_max-0.13;
dx=0.005;
dr=0.005;
da=pi/180*5;
r_max=0.1;
a_max=pi/180*120;
q_min=[qh_min*ones(1,3) qv_min*ones(1,3)];
q_max=[qh_max*ones(1,3) qv_max*ones(1,3)];
map_color=jet(round(r_max/dr)+1);

ws=[];
wsc=[];
rm3=[];
for x=qv_min:dx:qv_max
    rm2=[];
    for a=0:da:a_max
        rm1=[];
        for r=0:dr:r_max
            pu=[x r*cos(a) r*sin(a)];
            angu=[0 0 0];
            [q,J]=InverseSolution_MOD(pu,[0 0 0]);
            in_limit=(q>=q_min)&(q<=q_max);
            if(in_limit==ones(size(in_limit)))
                ws=[ws;[pu r]];
                wsc=[wsc;map_color(round(r/dr)+1,:)];
                rm1=r;
            end
        end
        rm2=[rm2 rm1];
    end
    rm3=[rm3;[x min(rm2)]];
end

%%画图
figure
Rx=[1 0          0;
    0 cos(a_max) -sin(a_max);
    0 sin(a_max) cos(a_max)];
ws1=ws(:,1:3)';
ws2=Rx*ws1;
ws3=Rx*ws2;
ws_all=[ws1 ws2 ws3];
wsc=[wsc;wsc;wsc];
scatter3(ws_all(1,:),ws_all(2,:),ws_all(3,:),20,wsc,'.');
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

%%计算工作空间内切圆柱
% figure
% plot(rm3(:,1),rm3(:,2));
% xlabel('x');
% ylabel('r');
% rm4=rm3(rm3(:,1)>qv_min/2,:);
% rmin=min(rm4(:,2));
% rm5=rm3(rm3(:,2)>rmin,:);
% figure
% plot(rm5(:,1),rm5(:,2));