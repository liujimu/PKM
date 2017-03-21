function k = limb_stiffness( d, l )
%计算连杆刚度
%   k为连杆刚度
%   d为连杆直径
%   l为连杆长度
    E = 2.1e11; % 普通碳钢的弹性模量
    A = pi/4*d^2;
    k = E*A/l;
end

