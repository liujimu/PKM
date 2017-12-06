%验证PKM中的雅可比公式是否正确
clear;
pkm = PKM();
pose = [0 0 0 0 0 0]';
pkm.pose = pose;
pkm.calVelJac();
%公式推导的雅可比
jac = pkm.jac;
h = 1e-8;
h12 = 1.2e-7;
dydx = zeros(6);
b = h * eye(6);
for i = 1:6
    br = b(:,i);
    pkm.pose = pose - 2 * br;
    y1 = pkm.q;
    pkm.pose = pose - br;
    y2 = pkm.q;
    pkm.pose = pose + br;
    y4 = pkm.q;
    pkm.pose = pose + 2 * br;
    y5 = pkm.q;
    dydx(:,i) = (8 * y4 - 8 * y2 + y1 - y5)' ./ h12;
end
%数值方法求得的雅可比的逆
inv_jac = dydx;
delta = jac*inv_jac - eye(6)