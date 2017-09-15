classdef PKM < handle
    properties
        %各关节初始位置
        S_init = zeros(3,6);
        U_init = zeros(3,6);
        P_dir = zeros(3,6);
        %各关节当前位置
        S_cur = zeros(3,6);
        U_cur = zeros(3,6);
        %最大最小行程
        q_min = zeros(1,6);
        q_max = zeros(1,6);
        %末端位姿
        pf_pe = zeros(1,6);
        %驱动量
        q = zeros(1,6);
        %速度雅可比
        jac = zeros(6);
    end
    
    methods
        function obj = PKM(r1,r2,h1,l,alpha,beta)
            %旋转矩阵
            RzB1  = [cosd(alpha) -sind(alpha) 0;
                     sind(alpha) cosd(alpha)  0;
                     0           0            1];
            Rz120 = [cosd(120) -sind(120) 0;
                     sind(120) cosd(120)  0;
                     0         0          1];
            RyP4 = [cosd(beta)  0 sind(beta);
                    0           1          0;
                    -sind(beta) 0 cosd(beta)];
            % 动平台S副位置
            S1_= [0 r1 h1]';
            S1 = RzB1 * S1_;
            S2 = Rz120 * S1;
            S3 = Rz120 * S2;
            S4 = [0 r2 0]';
            S5 = Rz120 * S4;
            S6 = Rz120 * S5;
            obj.S_init = [S1, S2, S3, S4, S5, S6];
            % P副方向
            P1_dir = [0 0 1]';
            P2_dir = P1_dir;
            P3_dir = P1_dir;
            P4_dir_ = [-1 0 0]';
            P4_dir = RyP4 * P4_dir_;
            P5_dir = Rz120 * P4_dir;
            P6_dir = Rz120 * P5_dir;
            obj.P_dir = [P1_dir, P2_dir, P3_dir, P4_dir, P5_dir, P6_dir];
            % U副初始位置
            obj.U_init = -l * P_dir + S_init;
        end
       %% 运动学反解
        function invKin( obj )
            % 运动后S副位置
            rotm = RotMat(obj.pf_pe(4),obj.pf_pe(5),obj.pf_pe(6));
            p = obj.pf_pe(1:3)';
            obj.S_cur = rotm * obj.S_init + p;
            UiSc = obj.S_cur - obj.U_init;
            for i = 1:6
                obj.q(i) = UiSc(:,i)'*obj.P_dir(:,i) - sqrt( limb^2 - (UiSc(:,i)'*UiSc(:,i) - (UiSc(:,i)'*obj.P_dir(:,i))^2));
                obj.U_cur(:,i) = obj.U_init(:,i) + obj.q(i) * obj.P_dir(:,i);
            end            
        end
       %% 速度雅可比
        function calVelJac( obj )
            Jinv = zeros(6);
            l = obj.S_cur - obj.U_cur; %连杆向量
            for i = 1:6
                liei = l(:,i)' * obj.P_dir(:, i);
                RSi = rotm * obj.S_init(:,i);
                Jinv(i, :) = [l(:,i)'/liei, (cross(RSi, l(:,i))/liei)'];
            end
            obj.jac = inv(Jinv);
        end
        
    end
end