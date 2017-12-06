function [ In_Freq,MK] = In_freq_6PSSU(PAP,PB0,Dr,lr,Pos)
%��6PSS��������ĳһ�̶�λ��ʱ�������˵���ѹ�նȿ������ڣ�����һ��SPS���� ��������ĸն� ����Ƶ��
%   In_Freq - ����Ƶ��
%   MK - �����ն�

%Pos=[0.3 0 0 0 0 0]';
vel_p=[0 0 0 0 0 0]';

I_p=[0.03270264 -0.00000301 0.00000000;
	-0.00000301 0.01953580 -0.00000012;
	0.00000000 -0.00000012 0.01953098];%��ƽ̨��������
plat_cm=[-0.02703714 0.00002610 0.00000000]'; %��ƽ̨��������
m_p=2.53308082;%��ƽ̨����

I_l=[0.00131289 0.00000000 0.00000000;
     0.00000000 0.00132473 0.00000000
     0.00000000 0.00000000 0.00006213];%���˹�������
m_l=0.444*[1 1 1 1 1 1];%��������
m_s=1.2;%��������
%Is=8E-004*1.1;%˿��ת������


%�����ṹ�ߴ綨�壬��λȡΪm
sH=0.054;%�Խ�������ƽ̨�µ����ĵĸ߶�

        vec_e=Dr;
        p_x=Pos(1); p_y=Pos(2); p_z=Pos(3);
        pth=[Pos(4) Pos(5) Pos(6)];
        P=[p_x p_y p_z]';
        %����ת������RPY��................. 
        % T=[cos(pth(3))*cos(pth(2)),cos(pth(3))*sin(pth(2))*sin(pth(1))-sin(pth(3))*cos(pth(1)),cos(pth(3))*sin(pth(2))*cos(pth(1))+sin(pth(3))*sin(pth(1));
        %                    sin(pth(3))*cos(pth(2)),sin(pth(3))*sin(pth(2))*sin(pth(1))+cos(pth(3))*cos(pth(1)),sin(pth(3))*sin(pth(2))*cos(pth(1))-cos(pth(3))*sin(pth(1));
        %                    -sin(pth(2)),cos(pth(2))*sin(pth(1)),cos(pth(2))*cos(pth(1))];
        %312ŷ����
        rx=pth(1);ry=pth(2);rz=pth(3);
        T=[ cos(ry)*cos(rz) - sin(rx)*sin(ry)*sin(rz), -cos(rx)*sin(rz), cos(rz)*sin(ry) + cos(ry)*sin(rx)*sin(rz);
            cos(ry)*sin(rz) + cos(rz)*sin(rx)*sin(ry),  cos(rx)*cos(rz), sin(ry)*sin(rz) - cos(ry)*cos(rz)*sin(rx);
                               -cos(rx)*sin(ry),          sin(rx),                           cos(rx)*cos(ry)];
        P=P+T*[sH 0 0]'-[sH 0 0]';
        PA=T*PAP+repmat(P,1,numel(PAP(1,:)));%��ƽ̨�µ��ڻ�����ϵ�µ�����
        LI=PA-PB0;

        %�˶�ѧ���⣬��������λ����
        q1 = LI(:,1).'*Dr(:,1)+( (LI(:,1).'*Dr(:,1))^2 - LI(:,1).'*LI(:,1) + lr^2 )^0.5;
        q2 = LI(:,2).'*Dr(:,2)+( (LI(:,2).'*Dr(:,2))^2 - LI(:,2).'*LI(:,2) + lr^2 )^0.5;
        q3 = LI(:,3).'*Dr(:,3)+( (LI(:,3).'*Dr(:,3))^2 - LI(:,3).'*LI(:,3) + lr^2 )^0.5;
        q4 = LI(:,4).'*Dr(:,4)+( (LI(:,4).'*Dr(:,4))^2 - LI(:,4).'*LI(:,4) + lr^2 )^0.5;
        q5 = LI(:,5).'*Dr(:,5)+( (LI(:,5).'*Dr(:,5))^2 - LI(:,5).'*LI(:,5) + lr^2 )^0.5;
        q6 = LI(:,6).'*Dr(:,6)+( (LI(:,6).'*Dr(:,6))^2 - LI(:,6).'*LI(:,6) + lr^2 )^0.5;
        L=[q1 q2 q3 q4 q5 q6].';
        
        vec_qi=[Dr(:,1)*q1 Dr(:,2)*q2 Dr(:,3)*q3 Dr(:,4)*q4 Dr(:,5)*q5 Dr(:,6)*q6];
        LL = LI - vec_qi;%���˷�������  ��Bָ��A
        PB = PB0 + vec_qi;%����µ�λ������
        %���ſ˱Ⱦ���������ǵ��������ľ��󣬶������˵�ѹ���� 6SPS�ſ˱�
        RAP=T*PAP;
        G=[LL;cross(RAP,LL,1)]/lr;%�����ˣ����Ͻ���ѧ���壩 ���Դ��ӿ�ֱ��˵ļ����ٶ�
        GSPS=G;
        cang=[-Dr(:,1)'*LL(:,1)/lr
        -Dr(:,2)'*LL(:,2)/lr
        -Dr(:,3)'*LL(:,3)/lr
        -Dr(:,4)'*LL(:,4)/lr
        -Dr(:,5)'*LL(:,5)/lr
        -Dr(:,6)'*LL(:,6)/lr];
        %6PSS�ſ˱�
        G=-[G(:,1)/cang(1) G(:,2)/cang(2) G(:,3)/cang(3) G(:,4)/cang(4) G(:,5)/cang(5) G(:,6)/cang(6)];
 
    %�����ٶ�
    v_p=vel_p([1,2,3],:);%��ƽ̨ƽ��
    w_p=vel_p([4,5,6],:);%��ƽ̨ת��    
    
    %��ƽ̨����
    r_cp=plat_cm;%��ƽ̨����ϵ�£���ƽ̨���ĵ�λ��
    plat_cm_0=T*plat_cm;
    %����һ���Զ�ƽ̨����ϵΪ�ο���
    vec_cross=T*r_cp;
    R_cr=[0 -vec_cross(3) vec_cross(2);vec_cross(3) 0 -vec_cross(1);-vec_cross(2) vec_cross(1) 0]; %��˾���
    Tau_vp=[[1 0 0;0 1 0; 0 0 1;] -R_cr];
    
    I_p0=T*I_p*T';%T*(I_p + m_p*[y_a^2+z_a^2 -x_a*y_a -x_a*z_a;-x_a*y_a x_a^2+z_a^2 -y_a*z_a;-x_a*z_a -y_a*z_a x_a^2+y_a^2])*T';
    M_pequ=Tau_vp'*Tau_vp*m_p+[zeros(3) zeros(3);zeros(3) I_p0];  %��Ч������������
    
    %�����˶��ܼ���1
    v_q=G'*[v_p;w_p];
    M_lequ=0;M_slequ=0;K_l=0;
    for i=1:6
        l_ui=LL(:,i)/lr;
        %����������ϵ��������ϵ�ı任  ע��������Ķ���   y��Ķ��壬I���i
        BO_l=cross([1;0;0],l_ui); %������ϵX�� �� �������� ��˵õ�һ����������ϵ������
        if norm(BO_l) == 0 
            BO_l = [0;1;0];
        else
            BO_l=BO_l/norm(BO_l);
        end
        T_li=[BO_l cross(l_ui,BO_l) l_ui];
                
        %����һ���Ի�����ϵΪ�ο���
        vec_cross=l_ui;
        X_ui=[0 -vec_cross(3) vec_cross(2);vec_cross(3) 0 -vec_cross(1);-vec_cross(2) vec_cross(1) 0]; %��˾���
        vec_cross=T*PAP(:,i);
        Tau_li=[1 0 0;0 1 0; 0 0 1;cross([1 0 0],vec_cross);cross([0 1 0],vec_cross);cross([0 0 1],vec_cross)]';
        Tau_wli=X_ui*(Tau_li-vec_e(:,i)*G(:,i)')/lr;
        
        I_l0=T_li*I_l*T_li';
        M_liequ1=Tau_wli'*I_l0*Tau_wli;
        Tau_vli=vec_e(:,i)*G(:,i)'-1/2*lr*X_ui*Tau_wli;
        M_liequ=Tau_vli'*Tau_vli*m_l(i)+M_liequ1;
        
        M_lequ=M_lequ+M_liequ;
        
        %���鶯��
        M_sliequ=G(:,i)*G(:,i)'*m_s;

        M_slequ=M_slequ+M_sliequ;
        %˿�ܶ���
    end
    
    M_equ=M_lequ+M_pequ;
    %�նȾ���
    Ks_sc=0.000201*210/lr*10^9;%������ѹ�ն� k=E*A/L
    %Ks_sc=0.43*10^9;
    Ks_ac=eye(6)*Ks_sc; %�������ָն�
    Ks_qc=GSPS*Ks_ac*GSPS'; 
    
    %�������ָն�
    K_qu=2/(1/(2*10^9) + 1/(1*10^9));
    Ks_qu=G*K_qu*G';
    
    Ks_q=inv(inv(Ks_qc)+inv(Ks_qu));%�ϲ������նȾ���
    
    %�����ն�
    %svd(Ks_q(:,[1,2,3]))
    %svd(Ks_q(:,[4,5,6]))
    MK=eig(Ks_q);
    
    %����Ƶ��
    In_Freq=eig(Ks_q,M_equ).^0.5/2/pi;
    %In_Freq=eig(M_equ\Ks_q).^0.5/2/pi;
    %eig(Ks_q/M_equ).^0.5/2/pi
end

