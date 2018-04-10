clear
ru = 350;%130.2920;
thetau = 84.2426 / 180 * pi;
%xu = ru * cos(pi/2 - thetau/2);
%yu = ru * sin(pi/2 - thetau/2);
rl = 368.1854;
thetal = 24.9395 / 180 * pi;
%thetal = 2 * asin(xu / rl);
l = 388.25;
param_errors = zeros(54,1);
pkm = PKM( param_errors, rl, thetal, ru, thetau, l );
rmax = pkm.ll + pkm.ru - pkm.rl;
rw = pkm.getWorkspaceRadius();
plot(rw)
xlabel('\theta(deg)');
ylabel('r(mm)');
disp(min(rw));
%pose = [160*sind(60); 160*cosd(60); 0; 0; 0; 0];
pose = [0 -50 0 10/180*pi 0 0]';
pkm.setPose(pose);
disp(pkm.q);
wrench = [2000 0 -2000 0 670000 0]';
Fq = pkm.getJointForces(wrench);
disp(max(abs(Fq)))