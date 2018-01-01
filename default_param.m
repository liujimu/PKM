function param = default_param()
%生成6-PUS机构的尺寸参数
r1 = 100/1e3;       % 动平台竖直三角形分布半径
r2 = 140/1e3;       % 动平台倾斜三角形分布半径
h1 = -6/1e3;        % 动平台竖直三角形沿Z方向的偏置
limb = 154/1e3;     % 连杆长度
alpha = 30;         % 动平台竖直三角形与倾斜三角形的旋转角
beta = 0;          % 倾斜三角形P副方向与xy平面的夹角
d_l = 0.02;         % 连杆直径
d_s = 0.02;         % 丝杠小径
q_min = -0.2;       % 滑块的最小行程
q_max = 0.2;        % 滑块的最大行程
param = struct('r1',r1,'r2',r2,'h1',h1,'limb',limb,'alpha',alpha,'beta',beta,'d_l',d_l,'d_s',d_s,'q_min',q_min,'q_max',q_max);
end

