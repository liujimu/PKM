classdef PKM < handle
    properties
        %设计参数
        rl = 0;
        thetal = 0;
        ru = 0;
        thetau = 0;
        ll = 0;
        %各关节初始位置
        S_init = zeros(3,6);
        U_init = zeros(3,6);
        P_dir = zeros(3,6);
        %连杆长度
        l = zeros(1,6);
        %驱动量
        q = zeros(6,1);
        home_errors = zeros(6,1);
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
        function obj = PKM( param_errors, rl, thetal, ru, thetau, l, qmin, qmax )
            %默认参数
            if nargin < 8
                qmin = -5;
                qmax = 260;
                if nargin < 6
                    l = 260;
                    if nargin < 5
                        ru = 82.1;
                        thetau = 86.0151 / 180 * pi;
                        if nargin < 3
                            %rl = 257;
                            %thetal = 25.1713 / 180 * pi;
                            rl = 240;
                            thetal = 26.9868 / 180 * pi;
                            if nargin < 1
                                param_errors = zeros(54,1);
                            end
                        end
                    end
                end
            end
            %设计参数
            obj.rl = rl;
            obj.thetal = thetal;
            obj.ru = ru;
            obj.thetau = thetau;
            obj.ll = l;
            %误差参数
            obj.home_errors = param_errors(1:6);%电机home位置误差
            l_errors = param_errors(7:12);  %杆长误差
            P_errors = zeros(3,6);          %导轨方向误差
            U_errors = zeros(3,6);          %U副初始位置误差
            S_errors = zeros(3,6);          %S副初始位置误差
            for i = 1:6
                P_errors(1:2,i) = param_errors((12+2*i-1):(12+2*i)); %只考虑x,y方向（直线方向只有两个独立变量）
                U_errors(1:2,i) = param_errors((24+2*i-1):(24+2*i)); %只考虑x,y方向
                S_errors(:,i) = param_errors((36+3*i-2):(36+3*i));
            end
            %旋转矩阵
            Rz120 = [cosd(120) -sind(120) 0;
                     sind(120) cosd(120)  0;
                     0         0          1];
            % 连杆长度
            obj.l = l * ones(6,1) + l_errors;
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
            %U1(3) = -sqrt(obj.l(1)^2 - (U1(1) - S1(1))^2 - (U1(2) - S1(2))^2);
            U2 = [rl * cos(pi/2 + thetal/2); rl * sin(pi/2 + thetal/2); 0];
            %U2(3) = -sqrt(obj.l(2)^2 - (U2(1) - S2(1))^2 - (U2(2) - S2(2))^2);
            U3 = Rz120 * U1;
            U4 = Rz120 * U2;
            U5 = Rz120 * U3;
            U6 = Rz120 * U4;
            obj.U_init = [U1, U2, U3, U4, U5, U6] + U_errors;
            for i = 1:6
                obj.U_init(3,i) = -sqrt(obj.l(i)^2 - (obj.U_init(1,i) - obj.S_init(1,i))^2 - (obj.U_init(2,i) - obj.S_init(2,i))^2);
            end
            % P副方向
            P1_dir = [0 0 1]';
            P2_dir = [0 0 1]';
            P3_dir = Rz120 * P1_dir;
            P4_dir = Rz120 * P2_dir;
            P5_dir = Rz120 * P3_dir;
            P6_dir = Rz120 * P4_dir;
            obj.P_dir = [P1_dir, P2_dir, P3_dir, P4_dir, P5_dir, P6_dir] + P_errors;
            for i = 1:6
                obj.P_dir(3,i) = sqrt(1 - obj.P_dir(1,i)^2 - obj.P_dir(2,i)^2);
            end
            % 行程限制
            obj.q_min = qmin * ones(6,1);
            obj.q_max = qmax * ones(6,1);
        end
        %% 设置末端位姿
        function setPose(obj,val)
            obj.pose = val;
            obj.invKin();
        end
        %% 运动学反解
        function invKin( obj )
            alpha = obj.pose(4);
            beta = obj.pose(5);
            gamma = obj.pose(6);
            % 123欧拉角
            obj.rotm = [                                    cos(beta)*cos(gamma),                                     -cos(beta)*sin(gamma),              sin(beta);
                         sin(alpha)*sin(beta)*cos(gamma) + cos(alpha)*sin(gamma),  -sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma),  -sin(alpha)*cos(beta);
                        -cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma),   cos(alpha)*sin(beta)*sin(gamma) + sin(alpha)*cos(gamma),   cos(alpha)*cos(beta)];
            p = repmat(obj.pose(1:3),1,6);
            % 运动后S副位置
            obj.S_cur = obj.rotm * obj.S_init + p;
            UiSc = obj.S_cur - obj.U_init;
            q = zeros(6,1);
            for i = 1:6
                radicand = obj.l(i)^2 - (UiSc(:,i)'*UiSc(:,i) - (UiSc(:,i)'*obj.P_dir(:,i))^2);
                % 判断被开方数是否小于0，避免虚数根
                if radicand < 0
                    % fprintf('Pose [%f, %f, %f, %f, %f, %f] is out of workspace\n',obj.pose(1),obj.pose(2),obj.pose(3),obj.pose(4),obj.pose(5),obj.pose(6));
                    q(i) = Inf;
                    obj.U_cur(:,i) = obj.U_init(:,i);
                else
                    q(i) = UiSc(:,i)'*obj.P_dir(:,i) - sqrt(radicand);
                    obj.U_cur(:,i) = obj.U_init(:,i) + q(i) * obj.P_dir(:,i);
                    UcSc_i = obj.S_cur(:,i) - obj.U_cur(:,i);
                    obj.l_dir(:,i) = UcSc_i / norm(UcSc_i);
                end
            end
            obj.q = q + obj.home_errors;
        end
        %% 速度雅可比
        function calVelJac( obj )
            % dX = J * dq
            % dX是笛卡尔空间速度，dq是关节空间速度
            % J = inv(G')
            Jinv = zeros(6);
            for i = 1:6
                liPi = obj.l_dir(:,i)' * obj.P_dir(:, i);
                RSi = obj.rotm * obj.S_init(:,i);
                Jinv(i, :) = [obj.l_dir(:,i)'/liPi, (cross(RSi, obj.l_dir(:,i))/liPi)'];
            end
            obj.jac = inv(Jinv);
        end
        %% 运动学正解
        function forKin( obj, qin, init_pose )
            if nargin == 2
                init_pose = zeros(6,1);
            end
            X = init_pose;
            dq = ones(6,1);
            eps = 10e-9;
            n = 0;
            while norm(dq) > eps
                X0 = X;
                obj.setPose(X0);
                q0 = obj.q;
                obj.calVelJac();
                X = X0 - obj.jac * (q0 - qin);
                dq = q0 - qin;
                n = n + 1;
            end
            %disp(n);
            obj.setPose(X);
        end
        %% 判断是否在工作空间内
        function boolout = isInWorkspace( obj, margin )
            if nargin == 1
                margin = 0;
            end
            boolout = true;
            obj.invKin();
            % 判断是否超行程
            if sum( obj.q >= (obj.q_min + margin.*ones(size(obj.q))) & obj.q <= (obj.q_max - margin.*ones(size(obj.q)))) < 6
                boolout = false;
            end
        end
        %% 判断是否奇异
        function boolout = isSingular( obj )
            boolout = false;
            obj.calVelJac();
            if abs(det(obj.jac)) < 10e-3 || abs(det(obj.jac)) > 10e3
                boolout = true;
            end
        end
        %% 计算连杆分别与动平台Z轴、P副的夹角
        function angle_out = checkSingularity( obj )
            angle = zeros(2,6);
            vz = obj.rotm * [0 0 1]';
            for i = 1:6
                angle(1,i) = acos(obj.l_dir(:,i)'*vz);
                angle(2,i) = acos(obj.l_dir(:,i)'*obj.P_dir(:,i));
            end
            angle_out = max(max(angle)) / pi * 180;
        end
        %% 计算平动工作空间的半径
        function r_out = getWorkspaceRadius( obj )
            theta = 1:60;
            n = length(theta);
            r_out = zeros(1,n);
            rmax = obj.ll + obj.ru - obj.rl;
            dr = 1;
            for i = 1:n
                inWorkspace = 0;
                r = rmax;
                while ~inWorkspace
                    r = r - dr;
                    pose_ = [r*sind(theta(i)); r*cosd(theta(i)); 0; 0; 0; 0];
                    obj.setPose(pose_);
                    inWorkspace = sum(isfinite(obj.q))==6;
                end
                r_out(i) = r;
            end
        end
        %% 根据末端载荷求关节力
        function joint_forces = getJointForces( obj, wrench )
            % τ = J' * F
            % τ是关节空间的力矢量，F是笛卡尔空间的力-力矩矢量
            obj.calVelJac();
            joint_forces = obj.jac' * wrench;
        end
        %% 计算刚度
        function K_p = getStiffness(obj)
            % K = inv(J') * k * inv(J)
            % k为连杆刚度，是一对角阵；K为刚度矩阵
            
            E_steel = 2.05e5; %钢的弹性模量，单位N/mm^2
            G_steel = 8e4; %钢的剪切模量，单位N/mm^2
            
            %Steward的雅可比
            G_l = zeros(6);
            for i = 1:6
                RSi = obj.rotm * obj.S_init(:,i);
                G_l(:,i) = [obj.l_dir(:,i); cross(RSi, obj.l_dir(:,i))];
            end
            %连杆的刚度
            d_li = 20; %连杆直径
            K_li = E_steel .* pi .* d_li^2 ./ obj.ll; %连杆拉压刚度
            K_l = K_li * eye(6);
            K_d = G_l * K_l * G_l';
            
            %丝杆的刚度
            d_s = 15; %丝杆小径，单位mm；
            l_s = 50 + obj.q;
            K_a = E_steel .* pi .* d_s^2 ./ ( 4*l_s ) .* eye(6);
            %K_a = diag([100 200 300 400 500 600]);
            obj.calVelJac();
            K_p = inv(obj.jac')*K_a*inv(obj.jac);
            
            disp('连杆刚度：')
            disp(K_d);
            disp('丝杆刚度：')
            disp(K_p);
            
            K = inv(inv(K_p) + inv(K_d));
        end
        %% 计算等效惯量
        function Mp_eq = getEquivalentInertia(obj)
            m_p = 5.78; %动平台质量，单位kg
            I_p = [2.405E+004 0          0;
                   0          1.424E+004 0;
                   0          0          1.424E+004]; %动平台转动惯量
            r_cmp = [0.0, 0.0, 18.28]';
            R_pc = eye(3);
            Rr_cmp = obj.rotm * r_cmp;
            Rr_cross = [0          -Rr_cmp(3) Rr_cmp(2);
                       Rr_cmp(3)  0          -Rr_cmp(1);
                       -Rr_cmp(2) Rr_cmp(1)  0         ];
            Ip_eq = obj.rotm * R_pc * I_p * R_pc' * obj.rotm';
            Mp_eq = m_p * [eye(3); -Rr_cross] * [eye(3), -Rr_cross] + [zeros(3), zeros(3); zeros(3), Ip_eq];
            
            m_li = 0.815; %连杆质量，单位kg
            I_li = [5498.42728 0         0;
                    0          5384.4886 0;
                    0          0         184.11219]; %连杆转动惯量
            m_sli = 0.614; %滑块质量，单位kg
                
            obj.calVelJac();
            G = inv(obj.jac)';
            M_eq = Mp_eq;
            for i = 1:6
                l_xi = cross([1 0 0]',obj.l_dir(:,i));
                R_li = [l_xi cross(obj.l_dir(:,i),l_xi) obj.l_dir(:,i)];
                lambda_i = G(:,i);
                rcmp_cross = [0         -r_cmp(3) r_cmp(2);
                              r_cmp(3)  0         -r_cmp(1);
                              -r_cmp(2) r_cmp(1)  0         ];
                li_cross = [0             -obj.l_dir(3) obj.l_dir(2);
                            obj.l_dir(3)  0             -obj.l_dir(1);
                            -obj.l_dir(2) obj.l_dir(1)  0           ];
                gamma_wi = li_cross / obj.l(i) * ([eye(3), -rcmp_cross] - obj.P_dir(:,i) * lambda_i');
                gamma_vi = obj.P_dir(:,i) * lambda_i' - 1/2 * obj.l(i) * li_cross * gamma_wi;
                Mli_eq = gamma_wi' * R_li * I_li * R_li' * gamma_wi + gamma_vi' * gamma_vi * m_li;
                
                Msli_eq = m_sli * (lambda_i * lambda_i');
                
                M_eq = M_eq + Mli_eq + Msli_eq;
            end
        end
    end
end