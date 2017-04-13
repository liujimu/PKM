clear;
close all
%%参数初始化
qh_max=0.25;
qh_min=qh_max-0.5;
qv_max=0.25;
qv_min=qv_max-0.5;
dx=0.01;
dr=0.005;
da=pi/180*5;
r_max=0.5;
a_max=pi/180*120;
q_min=[qh_min*ones(1,3) qv_min*ones(1,3)];
q_max=[qh_max*ones(1,3) qv_max*ones(1,3)];
map_color=jet(round(r_max/dr)+1);

ws=[];
wsc=[];
rm3=[];

for x=qv_min:dx:qv_max
    for a=0:da:a_max
        rm1=[];
        r=0:dr:r_max;
        isWithinLimit=zeros(size(r));
        for i=1:length(r)
            pu=[x r(i)*cos(a) r(i)*sin(a)];
            angu=[0 0 0];
            [q,J]=InverseSolution_MOD(pu,[0 0 0]);
            in_limit=(q>=q_min)&(q<=q_max);
            %判断是否超行程
            if(in_limit==ones(size(in_limit)))
                isWithinLimit(i)=1;
            end
        end
        %找边界
        dif=diff(isWithinLimit);
        for j=1:length(dif)
            if(dif(j)==-1)
                ws=[ws;[a,r(j),x]];
            elseif(dif(j)==1)
                ws=[ws;[a,r(j+1),x]];
            else
            end
        end
    end
end

%%画图
figure
ws_all=[ws;ws+repmat([pi/180*120 0 0],size(ws,1),1);ws+repmat([pi/180*240 0 0],size(ws,1),1)];

[Y,Z,X] = pol2cart(ws_all(:,1),ws_all(:,2),ws_all(:,3));

scatter3(X,Y,Z,'.');
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