function statesdot = statesdot_fgrg(C,Flong,Fz_f,Fz_r,Iz,R,delta,lf,lr,m,mu,omegaB,varphi,vx,vy)
%STATESDOT_FGRG
%    STATESDOT = STATESDOT_FGRG(C,FLONG,FZ_F,FZ_R,IZ,R,DELTA,LF,LR,M,MU,OMEGAB,VARPHI,VX,VY)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    29-Nov-2018 20:17:06

t2 = cos(varphi);
t3 = sin(varphi);
t4 = C.^2;
t5 = lr.*omegaB;
t6 = t5+vy;
t7 = vx+1.0e-8;
t8 = 1.0./t7;
t9 = t6.*t8;
t10 = atan(t9);
t11 = delta-t10;
t12 = tan(t11);
t13 = t12.^2;
t14 = 1.0./m;
t15 = cos(delta);
t16 = C.*t12;
t17 = 1.0./Fz_f.^2;
t18 = 1.0./mu.^2;
t19 = R.*(2.0./3.0);
t20 = t19-1.0;
t21 = 1.0./Fz_f;
t22 = 1.0./mu;
t23 = t4.*t13;
t24 = sqrt(t23);
t25 = R-2.0;
t26 = C.*t12.*t21.*t22.*t24.*t25.*(1.0./3.0);
t39 = C.*t4.*t12.*t13.*t17.*t18.*t20.*(1.0./9.0);
t27 = t16+t26-t39;
t28 = sin(delta);
t31 = lf.*omegaB;
t29 = -t31+vy;
t30 = t29.^2;
t32 = 1.0./Fz_r.^2;
t33 = 1.0./t7.^3;
t34 = C.*t4.*t18.*t20.*t29.*t30.*t32.*t33.*(1.0./9.0);
t35 = 1.0./Fz_r;
t36 = 1.0./t7.^2;
t37 = t4.*t30.*t36;
t38 = sqrt(t37);
statesdot = [t2.*vx-t3.*vy;t3.*vx+t2.*vy;omegaB;t14.*(Flong+Flong.*t15-t27.*t28+m.*omegaB.*vy);t14.*(t34+Flong.*t28+t15.*t27-C.*t8.*t29-m.*omegaB.*vx-C.*t8.*t22.*t25.*t29.*t35.*t38.*(1.0./3.0));(lr.*(-t34+C.*t8.*t29+C.*t8.*t22.*t25.*t29.*t35.*t38.*(1.0./3.0))+Flong.*lf.*t28+lf.*t15.*t27)./Iz];
