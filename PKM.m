classdef PKM < handle
    properties
        %各关节初始位置
        S_init = zeros(3,6);
        U_init = zeros(3,6);
        P_dir = zeros(3,6);
        %连杆长度
        l = zeros(1,6);
        %驱动量
        q = zeros(6,1);
        %最大最小行程
        q_min = zeros(6,1);
        q_max = zeros(6,1);
        %末端位姿
        pose = zeros(6,1);
        rotm = eye(3);
        %各关节当前位置
        S_cur = zeros(3,6);
        U_cur = zeros(3,6);
        l_dir = zeros(3,6);
        %速度雅可比
        jac = zeros(6);
    end
    
    methods
        %% 初始化
        function obj = PKM( param_errors, ru, rl, thetau, thetal, l, qmin, qmax )
            %默认参数
            if nargin < 8 
                ru = 82.1 / 1000;
                rl = 240 / 1000;
                thetau = 86.0151 / 180 * pi;
                thetal = 26.9868 / 180 * pi;
                l = 260 / 1000;
                qmin = -0.005;
                qmax = 0.2;
                if nargin < 1
                    param_errors = zeros(60,1);
                end
            end
            %误差参数
            l_errors = param_errors(1:6);
            U_errors = zeros(3,6);
            S_errors = zeros(3,6);
            P_errors = zeros(3,6);
            for i = 1:6
                U_errors(:,i) = param_errors((9*i-2):(9*i));
                S_errors(:,i) = param_errors((9*i+1):(9*i+3));
                P_errors(:,i) = param_errors((9*i+4):(9*i+6));
            end
            %旋转矩阵
            Rz120 = [cosd(120) -sind(120) 0;
                     sind(120) cosd(120)  0;
                     0         0          1];
            % 动平台S副位置
            S1 = [ru * cos(pi/2 - thetau/2); ru * sin(pi/2 - thetau/2); 0];
            S2 = [ru * cos(pi/2 + thetau/2); ru * sin(pi/2 + thetau/2); 0];
            S3 = Rz120 * S1;
            S4 = Rz120 * S2;
            S5 = Rz120 * S3;
            S6 = Rz120 * S4;
            obj.S_init = [S1, S2, S3, S4, S5, S6] + S_errors;
            % U副初始位置
            U1 = [rl * cos(pi/2 - thetal/2); rl * sin(pi/2 - thetal/2); 0];
            U1(3) = -sqrt(l^2 - (U1(1) - S1(1))^2 - (U1(2) - S1(2))^2);
            U2 = [rl * cos(pi/2 + thetal/2); rl * sin(pi/2 + thetal/2); 0];
            U2(3) = -sqrt(l^2 - (U2(1) - S2(1))^2 - (U2(2) - S2(2))^2);
            U3 = Rz120 * U1;
            U4 = Rz120 * U2;
            U5 = Rz120 * U3;
            U6 = Rz120 * U4;
            obj.U_init = [U1, U2, U3, U4, U5, U6] + U_errors;
            % P副方向
            P1_dir = [0 0 1]';
            P2_dir = [0 0 1]';
            P3_dir = Rz120 * P1_dir;
            P4_dir = Rz120 * P2_dir;
            P5_dir = Rz120 * P3_dir;
            P6_dir = Rz120 * P4_dir;
            obj.P_dir = [P1_dir, P2_dir, P3_dir, P4_dir, P5_dir, P6_dir] + P_errors;
            % 连杆长度
            obj.l = l * ones(6,1) + l_errors;
            % 行程限制
            obj.q_min = qmin * ones(6,1);
            obj.q_max = qmax * ones(6,1);
        end
        %% 设置末端位姿
        function set.pose(obj,val)
            obj.pose = val;
            % disp('set function called');
            obj.invKin();
        end
        %% 运动学反解
        function invKin( obj )
            % 运动后S副位置
            % rotm = RotMat(obj.pose(4),obj.pose(5),obj.pose(6));
            alpha = obj.pose(4);
            beta = obj.pose(5);
            gamma = obj.pose(6);
            obj.rotm = [                                    cos(beta)*cos(gamma),                                     -cos(beta)*sin(gamma),              sin(beta);
                         sin(alpha)*sin(beta)*cos(gamma) + cos(alpha)*sin(gamma),  -sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma),  -sin(alpha)*cos(beta);
                        -cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma),   cos(alpha)*sin(beta)*sin(gamma) + sin(alpha)*cos(gamma),   cos(alpha)*cos(beta)];
            p = repmat(obj.pose(1:3),1,6);
            obj.S_cur = obj.rotm * obj.S_init + p;
            UiSc = obj.S_cur - obj.U_init;
            for i = 1:6
                obj.q(i) = UiSc(:,i)'*obj.P_dir(:,i) - sqrt(obj.l(i)^2 - (UiSc(:,i)'*UiSc(:,i) - (UiSc(:,i)'*obj.P_dir(:,i))^2));
                obj.U_cur(:,i) = obj.U_init(:,i) + obj.q(i) * obj.P_dir(:,i);
            end
            obj.l_dir = obj.S_cur - obj.U_cur;
            for i = 1:6
                obj.l_dir(:,i) = obj.l_dir(:,i) / norm(obj.l_dir(:,i));
            end
        end
        %% 速度雅可比
        function calVelJac( obj )
            Jinv = zeros(6);
            for i = 1:6
                liei = obj.l_dir(:,i)' * obj.P_dir(:, i);
                RSi = obj.rotm * obj.S_init(:,i);
                Jinv(i, :) = [obj.l_dir(:,i)'/liei, (cross(RSi, obj.l_dir(:,i))/liei)'];
            end
            % fprintf('det(J^-1) = %d\n',det(Jinv));
            obj.jac = inv(Jinv);
        end
        %% 判断是否在工作空间内
        function boolout = isInWorkspace( obj )
            boolout = true;
            obj.invKin();
            % 判断是否超行程
            if sum( obj.q >= obj.q_min & obj.q <= obj.q_max ) < 6
                boolout = false;
            end
        end
        %% 判断是否奇异
        function boolout = isSingular( obj )
            boolout = true;
            obj.calVelJac();
            if abs(det(obj.jac)) < 10e-3 || abs(det(obj.jac)) > 10e3
                boolout = false;
            end
        end
    end
end