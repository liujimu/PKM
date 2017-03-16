%%%%%%%%%%%%%%%%%%%%%%%%%% ROTATION MATRIX FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
% The three Euler angles that form the rotation matrix are defined by first    %
% rotating the mobile platform about the base z-axis by an angle -phi, then    %
% about the base z-axis by an angle theta, then about the base z-axis by an    %
% angle phi, and finally about the mobile z'-axis by an angle psi.             %
%                                                                              %
% Reference: (1) Bonev, I. A., and Ryu, J., "Orientation Workspace Analysis of %
%                6-DOF Parallel Manipulators", ASME Design Engineering         %
%                Technical Conferences (DETC'99), Las Vegas, NV, September     %
%                12-15, 1999.                                                  %
%                (http://wwwrobot.gmc.ulaval.ca/~bonev/DAC_8646.pdf)           %
%                                                                              %
% Created by: Ilian Bonev                e-mail: bonev@parallemic.org          %
% Last modified: July 20, 2000                                                 %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Rot_matrix = RotMat(phi, theta, psi)
 Rot_matrix = [ cos(phi)*cos(theta)*cos(psi-phi) - sin(phi)*sin(psi-phi), ...
                -cos(phi)*cos(theta)*sin(psi-phi) - sin(phi)*cos(psi-phi), ...
                cos(phi)*sin(theta); ...
                sin(phi)*cos(theta)*cos(psi-phi) + cos(phi)*sin(psi-phi), ...
                -sin(phi)*cos(theta)*sin(psi-phi) + cos(phi)*cos(psi-phi), ...
                sin(phi)*sin(theta);
                -sin(theta)*cos(psi-phi), sin(theta)*sin(psi-phi), cos(theta)];