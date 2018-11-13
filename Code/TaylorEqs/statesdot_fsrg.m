function statesdot = statesdot_fsrg(C,Flong,Fz_f,Fz_r,Iz,R,delta,lf,lr,m,mu,omegaB,varphi,vx,vy)
%STATESDOT_FSRG
%    STATESDOT = STATESDOT_FSRG(C,FLONG,FZ_F,FZ_R,IZ,R,DELTA,LF,LR,M,MU,OMEGAB,VARPHI,VX,VY)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    13-Nov-2018 15:09:44

t2 = cos(varphi);
t3 = sin(varphi);
t4 = 1.0./m;
t5 = sin(delta);
t6 = vx+1.0e-8;
t7 = 1.0./t6;
t21 = lf.*omegaB;
t8 = -t21+vy;
t9 = t7.*t8;
t10 = atan(t9);
t11 = t10+1.0e-8;
t12 = tan(t11);
t13 = lr.*omegaB;
t14 = t13+vy;
t15 = t7.*t14;
t16 = atan(t15);
t17 = -delta+t16+1.0e-8;
t18 = sign(t17);
t19 = cos(delta);
t20 = C.^2;
t22 = t12.^2;
t23 = C.*t12;
t24 = 1.0./Fz_r;
t25 = 1.0./mu;
t26 = t20.*t22;
t27 = sqrt(t26);
t28 = R-2.0;
t29 = t24.*t25.*t27.*t28.*(1.0./3.0);
t30 = 1.0./Fz_r.^2;
t31 = 1.0./mu.^2;
t32 = R.*(2.0./3.0);
t33 = t32-1.0;
statesdot = [t2.*vx-t3.*vy;t3.*vx+t2.*vy;omegaB;t4.*(Flong+Flong.*t19+m.*omegaB.*vy+Fz_f.*mu.*t5.*t18);-t4.*(t23+t29-Flong.*t5+m.*omegaB.*vx+Fz_f.*mu.*t18.*t19-C.*t12.*t20.*t22.*t30.*t31.*t33.*(1.0./9.0));(lr.*(t23+t29-C.*t12.*t20.*t22.*t30.*t31.*t33.*(1.0./9.0))+Flong.*lf.*t5-Fz_f.*lf.*mu.*t18.*t19)./Iz];
