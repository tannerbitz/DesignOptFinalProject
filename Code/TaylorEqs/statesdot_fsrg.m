function statesdot = statesdot_fsrg(C,Flong,Fz_f,Fz_r,Iz,R,delta,lf,lr,m,mu,omegaB,varphi,vx,vy)
%STATESDOT_FSRG
%    STATESDOT = STATESDOT_FSRG(C,FLONG,FZ_F,FZ_R,IZ,R,DELTA,LF,LR,M,MU,OMEGAB,VARPHI,VX,VY)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    29-Nov-2018 20:39:14

t2 = cos(varphi);
t3 = sin(varphi);
t4 = 1.0./m;
t5 = sin(delta);
t6 = vx+1.0e-8;
t7 = 1.0./t6;
t8 = lr.*omegaB;
t9 = t8+vy;
t10 = t7.*t9;
t11 = atan(t10);
t12 = -delta+t11;
t13 = sign(t12);
t14 = cos(delta);
t15 = C.^2;
t18 = lf.*omegaB;
t16 = -t18+vy;
t17 = t16.^2;
t19 = C.*t7.*t16;
t20 = 1.0./Fz_r.^2;
t21 = 1.0./mu.^2;
t22 = R.*(2.0./3.0);
t23 = t22-1.0;
t24 = 1.0./t6.^3;
t25 = 1.0./Fz_r;
t26 = 1.0./mu;
t27 = R-2.0;
t28 = 1.0./t6.^2;
t29 = t15.*t17.*t28;
t30 = sqrt(t29);
t31 = C.*t7.*t16.*t25.*t26.*t27.*t30.*(1.0./3.0);
statesdot = [t2.*vx-t3.*vy;t3.*vx+t2.*vy;omegaB;t4.*(Flong+Flong.*t14+m.*omegaB.*vy+Fz_f.*R.*mu.*t5.*t13);-t4.*(t19+t31-Flong.*t5+m.*omegaB.*vx+Fz_f.*R.*mu.*t13.*t14-C.*t15.*t16.*t17.*t20.*t21.*t23.*t24.*(1.0./9.0));(lr.*(t19+t31-C.*t15.*t16.*t17.*t20.*t21.*t23.*t24.*(1.0./9.0))+Flong.*lf.*t5-Fz_f.*R.*lf.*mu.*t13.*t14)./Iz];
