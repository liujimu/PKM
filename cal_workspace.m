clear;
close all
%%参数初始化
qh_max=0.0435;
qh_min=qh_max-0.13;
qv_max=0.0375;
qv_min=qv_max-0.13;
dx=0.002;
dr=0.005;
da=pi/180*5;
r_max=0.079;
a_max=pi/180*120;
q_min=[qh_min*ones(1,3) qv_min*ones(1,3)];
q_max=[qh_max*ones(1,3) qv_max*ones(1,3)];

ws=[];
wsc=[];
for x=qv_min:dx:qv_max
    for a=0:da:a_max
        for r=0:dr:r_max
            pu=[x r*cos(a) r*sin(a)];
            angu=[0 0 0];
            [q,J]=InverseSolution_MOD(pu,[0 0 0]);
            in_limit=(q>=q_min)&(q<=q_max);
            if(in_limit==ones(size(in_limit)))
                ws=[ws;pu];
%                 wsc=[wsc;map_color(round((x-qv_min)./dx)+1,:)];
            end
        end
    end
end

%%画图
Rx=[1 0          0;
    0 cos(a_max) -sin(a_max);
    0 sin(a_max) cos(a_max)];
ws1=ws';
ws2=Rx*ws1;
ws3=Rx*ws2;
ws_all=[ws1 ws2 ws3];

n_color=size(qv_min:dx:qv_max,2);
map_color=parula(n_color);
c=round((ws_all(1,:)-qv_min)/dx)+1;
wsc=map_color(c,:);

scatter3(ws_all(1,:),ws_all(2,:),ws_all(3,:),20,wsc,'.');
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
