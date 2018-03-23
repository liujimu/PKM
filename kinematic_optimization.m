clear
ru = 104.216;
thetau = 51.163 / 180 * pi;
xu = ru * cos(pi/2 - thetau/2);
yu = ru * sin(pi/2 - thetau/2);
rl = 257.463;
%thetal = 20.132 / 180 * pi;
thetal = 2 * asin(xu / rl);
l = 332.5;
param_errors = zeros(54,1);
pkm = PKM( param_errors, rl, thetal, ru, thetau, l );
rmax = pkm.ll + pkm.ru - pkm.rl;
rw = pkm.getWorkspaceRadius();
plot(rw)
xlabel('\theta(deg)');
ylabel('r(mm)');
disp(min(rw));
pose = [160*sind(60); 160*cosd(60); 0; 0; 0; 0];
pkm.setPose(pose);
disp(pkm.q);
wrench = [0 0 1000 0 0 0]';
Fq = pkm.getJointForces(wrench);
disp(max(abs(Fq)))