R=eye(3);
P=[-0.2 0 0];
P_cross=[0 -P(3) P(2);
         P(3) 0 -P(1);
         -P(2) P(1) 0];
T_f=[R zeros(3);P_cross*R R];
F_end=[0 0 1000 0 0 0]';
F_mp=T_f*F_end;

pu=[0 0 0];
angu=[0 0 0];
[q, G_T]=InverseSolution_MOD(pu,angu);
J=inv(G_T);
F_act=J'*F_mp;