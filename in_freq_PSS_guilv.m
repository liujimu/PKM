function [ output_args ] = find_qudong_max_new( )
clc;clear;
tic;


%灵巧工作空间
tx=[];ty=[];tz=[];trrx=[];trry=[];trrz=[];
r_step=360/15; %圆周点数
for ang_r=0:360/r_step:360%0:360/60:360       %工作空间圆柱边线扫描
    ia=3+1; ib=3+1;ic=3+1;i_angs=0;
    r0x=linspace(-6,6,ia)*pi/180;
    r0y=linspace(-6,6,ib)*pi/180;
    r0z=linspace(-6,6,ic)*pi/180;
    for i=1:ia      %组合姿态角 灵巧工作空间 扫描
        for j=1:ib
            for k=1:ic
                i_angs=i_angs+1;
                trrx1(i_angs)=r0x(i);
                trry1(i_angs)=r0y(j);
                trrz1(i_angs)=r0z(k);
            end
        end
    end
    angx=ones(1,i_angs)*ang_r/180*pi;
    tx=[tx 0.001*ones(1,i_angs)*-100];
    ty=[ty 0.001*[250*cos(angx)]];
    tz=[tz 0.001*[250*sin(angx)]];
    trrx=[trrx trrx1];trry=[trry trry1];trrz=[trrz trrz1];
end
i_angs=i_angs*(r_step+1)


........在固定坐标系下C点坐标，为方便以下力计算，单位取为m...................
a=450/1000;lra=1200/1000;b=450/1000;ang=30;


lable_x='a';i_a=0;
for a=0.2:0.05:0.7%100:500/30:400
    i_a=i_a+1
    Ap(i_a)=a;
    
    %机构结构尺寸定义，单位取为m...................
    [PAP PB0 Dr lr U_Pb U_Bb] = Structure_6PSS(a,b,lra,ang);
    
    iia=0;
    for i=1:i_angs
        px=tx(i);py=ty(i);pz=tz(i);
        rx=trrx(i);ry=trry(i);rz=trrz(i);
        p_x=px;p_y=py;p_z=pz;
        Pos=[p_x p_y p_z rx ry rz]';

        [In_Fri MKi]= In_freq_6PSSU(PAP,PB0,Dr,lr,Pos );

        if mod(i,5000)==0
            toc
        end
        iia=iia+1;
        MK(:,iia)=MKi;%
        In_Fr(:,iia)=In_Fri;%

    end
    toc

    min(In_Fr');
    fr=min(abs(In_Fr'));
    FR(i_a)=min(fr);


    KKi=min(abs(MK'));
    KK(i_a)=min(KKi);
end

left=0.1;wideth=0.38;hight=0.78;bot=0.13;
figure;
set(gcf,'Position',[300 300 800 350]);%边距 左 下 宽 高
subplot(1,2,1)
pDat=FR
x=Ap;
plot(x,pDat,'-');
%plot(x,pDat(1,:),'-.',x,pDat(2,:),'*-',x,pDat(3,:),'-+',x,pDat(4,:),'-d',x,pDat(5,:),'-s',x,pDat(6,:),'-^');
XLAB=xlabel(lable_x);
ylabel('Frequency');
xlim([ min(x) max(x)]);
set(XLAB,'FontName','Times New Roman','FontAngle','italic')
set(gca,'Position',[left bot wideth hight]);%边距 左 下 宽 高


subplot(1,2,2)
pDat=KK
x=Ap;
plot(x,pDat,'-');
%plot(x,pDat(1,:),'-.',x,pDat(2,:),'*-',x,pDat(3,:),'-+',x,pDat(4,:),'-d',x,pDat(5,:),'-s',x,pDat(6,:),'-^');
%legend('f1','f21','f22','f31','f32','f33');
XLAB=xlabel(lable_x);
ylabel('Stiffness');
xlim([ min(x) max(x)]);
set(XLAB,'FontName','Times New Roman','FontAngle','italic')

left=left+wideth+0.1;
set(gca,'Position',[left bot wideth hight]);%边距 左 下 宽 高
end

