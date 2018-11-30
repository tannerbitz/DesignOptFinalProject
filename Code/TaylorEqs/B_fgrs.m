function B = B_fgrs(C,Flong,Fz_f,Iz,R,delta,lf,lr,m,mu,omegaB,vx,vy)
%B_FGRS
%    B = B_FGRS(C,FLONG,FZ_F,IZ,R,DELTA,LF,LR,M,MU,OMEGAB,VX,VY)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    29-Nov-2018 20:39:16

t2 = 1.0./m;
t4 = lr.*omegaB;
t5 = t4+vy;
t6 = vx+1.0e-8;
t7 = 1.0./t6;
t8 = t5.*t7;
t9 = atan(t8);
t10 = delta-t9;
t3 = tan(t10);
t11 = t3.^2;
t12 = t11+1.0;
t13 = C.^2;
t14 = 1.0./Fz_f;
t15 = 1.0./mu;
t16 = t11.*t13;
t17 = R-2.0;
t18 = cos(delta);
t19 = 1.0./Fz_f.^2;
t20 = 1.0./mu.^2;
t21 = R.*(2.0./3.0);
t22 = t21-1.0;
t23 = sqrt(t16);
t24 = sin(delta);
t25 = C.*t3;
t26 = C.*t3.*t14.*t15.*t17.*t23.*(1.0./3.0);
t35 = C.*t3.*t11.*t13.*t19.*t20.*t22.*(1.0./9.0);
t27 = t25+t26-t35;
t28 = C.*t12;
t29 = C.*t12.*t14.*t15.*t17.*t23.*(1.0./3.0);
t30 = 1.0./sqrt(t16);
t31 = C.*t11.*t12.*t13.*t14.*t15.*t17.*t30.*(1.0./3.0);
t34 = C.*t11.*t12.*t13.*t19.*t20.*t22.*(1.0./3.0);
t32 = t28+t29+t31-t34;
t33 = 1.0./Iz;
B = reshape([0.0,0.0,0.0,t2.*(t18+1.0),t2.*t24,lf.*t24.*t33,0.0,0.0,0.0,-t2.*(Flong.*t24+t18.*t27+t24.*t32),t2.*(Flong.*t18+t18.*t32-t24.*t27),t33.*(Flong.*lf.*t18+lf.*t18.*t32-lf.*t24.*t27)],[6,2]);
