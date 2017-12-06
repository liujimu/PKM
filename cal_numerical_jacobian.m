function jac = cal_numerical_jacobian(fun, x)
% 雅可比矩阵的数值解法
    h = 1e-8;
    h12 = 1.2e-7;
    sz = length(x);
    f1 = fun(x);
    sz1 = length(f1);
    dydx = zeros(sz1, sz);
    b = h * eye(sz);
    for i = 1:sz
        br = b(:,i);
        y1 = fun(x - 2 * br);
        y2 = fun(x - br);
        y4 = fun(x + br);
        y5 = fun(x + 2 * br);
        dydx(:,i) = (8 * y4 - 8 * y2 + y1 - y5)' ./ h12;
    end
    jac = dydx;
end
