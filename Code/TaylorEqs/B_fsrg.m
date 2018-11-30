function B = B_fsrg(Flong,Fz_f,Iz,delta,lf,lr,m,mu,omegaB,vx,vy)
%B_FSRG
%    B = B_FSRG(FLONG,FZ_F,IZ,DELTA,LF,LR,M,MU,OMEGAB,VX,VY)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    29-Nov-2018 20:17:08

t2 = 1.0./m;
t3 = cos(delta);
t4 = lr.*omegaB;
t5 = t4+vy;
t6 = vx+1.0e-8;
t7 = 1.0./t6;
t8 = t5.*t7;
t9 = atan(t8);
t10 = -delta+t9;
t11 = sin(delta);
t12 = dirac(t10);
t13 = sign(t10);
t14 = 1.0./Iz;
B = reshape([0.0,0.0,0.0,t2.*(t3+1.0),t2.*t11,lf.*t11.*t14,0.0,0.0,0.0,-t2.*(Flong.*t11-Fz_f.*mu.*t3.*t13+Fz_f.*mu.*t11.*t12.*2.0),t2.*(Flong.*t3+Fz_f.*mu.*t3.*t12.*2.0+Fz_f.*mu.*t11.*t13),t14.*(Flong.*lf.*t3+Fz_f.*lf.*mu.*t3.*t12.*2.0+Fz_f.*lf.*mu.*t11.*t13)],[6,2]);
