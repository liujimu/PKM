 clear;
close all
%%参数初始化
Z_MIN = -0.25;
Z_MAX = 0.25;
R_MIN = 0;
R_MAX = 0.5;
TH_MIN = 0;
TH_MAX = pi/180*120;

DZ = 0.01;
DR = 0.005;
DTH = pi/180*5;

z = Z_MIN:DZ:Z_MAX;
r = R_MIN:DR:R_MAX;
th = TH_MIN:DTH:TH_MAX;

rotm=RotMat(0,0,0);

for i = 1:length(z)
    for j = 1:length(th)
        for k = 1:length(r)
            p = [r(k).*cos(th(j)); r(k).*sin(th(j)); z(i)];
            isInWorkspace(i,j,k) = isinworkspace(p,rotm);
        end
    end
end

%% 查找工作空间边界
outerBoundary = [];
innerBoundary = [];
maxRadiusArray = zeros(size(z));
for i = 1:length(z)
    maxRadius = [];
    for j = 1:length(th)
        rIndex = find(isInWorkspace(i,j,:));
        if ~isempty(rIndex)
            maxIndex = max(rIndex);
            minIndex = min(rIndex);
            outerBoundary = [outerBoundary, [th(j); r(maxIndex); z(i)]];
            if(minIndex > 1)
                innerBoundary = [innerBoundary,[th(j); r(minIndex); z(i)]];
            end
            maxRadius = [maxRadius, r(maxIndex)];
        end
    end
    if ~isempty(maxRadius)
        maxRadiusArray(i) = min(maxRadius);
    end
end

%% 计算工作空间内切圆柱
d1 = abs(z - maxRadiusArray)./sqrt(2); % 与直线x-y=0的距离
d2 = abs(z + maxRadiusArray)./sqrt(2); % 与直线x+y=0的距离
[minD1, minIndex1] = min(d1);
[minD2, minIndex2] = min(d2);
rCylinder1 = sqrt(z(minIndex1)^2 + maxRadiusArray(minIndex1)^2)/sqrt(2);
rCylinder2 = sqrt(z(minIndex2)^2 + maxRadiusArray(minIndex2)^2)/sqrt(2);
rCylinder = min([rCylinder1 rCylinder2]);
fprintf('内切圆柱的最大半径为： %f\n',rCylinder);

%%画图
figure
hold on
% 内切圆柱
[X,Y,Z] = cylinder(rCylinder);
Z = 2*rCylinder.*Z - rCylinder;
surf(X,Y,Z)
theta = 0:DTH:2*pi;
rc = linspace(0,rCylinder,10);
[theta,rc] = meshgrid(theta,rc);
[X1,Y1] = pol2cart(theta,rc);
Z1 = rCylinder.*ones(size(X1));
Z2 = -rCylinder.*ones(size(X1));
surf(X1,Y1,Z1);
surf(X1,Y1,Z2);
% 点云
boundary = [outerBoundary, innerBoundary];
boundaryTotal=[boundary, ...
               boundary + repmat([pi/180*120 0 0]',1,size(boundary,2)), ...
               boundary + repmat([pi/180*240 0 0]',1,size(boundary,2))];
[X,Y,Z] = pol2cart(boundaryTotal(1,:),boundaryTotal(2,:),boundaryTotal(3,:));
scatter3(X,Y,Z,'.');
hold off
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

figure
plot(z,maxRadiusArray);
hold on
fill([-rCylinder -rCylinder rCylinder rCylinder], [0 rCylinder rCylinder 0], 'r');
hold off
xlabel('z');
ylabel('r_m');
axis equal