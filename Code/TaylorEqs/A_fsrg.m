function A = A_fsrg(C,Fz_f,Fz_r,Iz,R,delta,lf,lr,m,mu,omegaB,varphi,vx,vy)
%A_FSRG
%    A = A_FSRG(C,FZ_F,FZ_R,IZ,R,DELTA,LF,LR,M,MU,OMEGAB,VARPHI,VX,VY)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    13-Nov-2018 15:35:50

t2 = cos(varphi);
t3 = sin(varphi);
t4 = lr.*omegaB;
t5 = t4+vy;
t6 = vx+1.0e-8;
t7 = 1.0./t6.^2;
t8 = 1.0./m;
t9 = 1.0./t6;
t10 = t5.*t9;
t11 = atan(t10);
t12 = -delta+t11;
t13 = dirac(t12);
t14 = sin(delta);
t15 = t5.^2;
t16 = t7.*t15;
t17 = t16+1.0;
t18 = 1.0./t17;
t19 = m.*omegaB;
t20 = C.^2;
t23 = lf.*omegaB;
t21 = -t23+vy;
t22 = t21.^2;
t24 = 1.0./Fz_r.^2;
t25 = 1.0./mu.^2;
t26 = R.*(2.0./3.0);
t27 = t26-1.0;
t28 = 1.0./t6.^3;
t29 = cos(delta);
t30 = 1.0./Fz_r;
t31 = 1.0./mu;
t32 = R-2.0;
t33 = t7.*t20.*t22;
t34 = 1.0./sqrt(t33);
t35 = C.*t7.*t21;
t36 = 1.0./t6.^4;
t37 = t20.*t22.*t28.*t30.*t31.*t32.*t34.*(1.0./3.0);
t38 = 1.0./Iz;
t39 = C.*t9;
t40 = vy.*2.0;
t41 = t40-lf.*omegaB.*2.0;
t42 = t7.*t20.*t30.*t31.*t32.*t34.*t41.*(1.0./6.0);
t43 = C.*lf.*t20.*t22.*t24.*t25.*t27.*t28.*(1.0./3.0);
A = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t3.*vx-t2.*vy,t2.*vx-t3.*vy,0.0,0.0,0.0,0.0,t2,t3,0.0,Fz_f.*mu.*t5.*t7.*t8.*t13.*t14.*t18.*-2.0,t8.*(-t19+t35+t37+Fz_f.*mu.*t5.*t7.*t13.*t18.*t29.*2.0-C.*t20.*t21.*t22.*t24.*t25.*t27.*t36.*(1.0./3.0)),-t38.*(lr.*(t35+t37-C.*t20.*t21.*t22.*t24.*t25.*t27.*t36.*(1.0./3.0))-Fz_f.*lf.*mu.*t5.*t7.*t13.*t18.*t29.*2.0),-t3,t2,0.0,t8.*(t19+Fz_f.*mu.*t9.*t13.*t14.*t18.*2.0),-t8.*(t39+t42+Fz_f.*mu.*t9.*t13.*t18.*t29.*2.0-C.*t20.*t22.*t24.*t25.*t27.*t28.*(1.0./3.0)),t38.*(lr.*(t39+t42-C.*t20.*t22.*t24.*t25.*t27.*t28.*(1.0./3.0))-Fz_f.*lf.*mu.*t9.*t13.*t18.*t29.*2.0),0.0,0.0,1.0,t8.*(m.*vy+Fz_f.*lr.*mu.*t9.*t13.*t14.*t18.*2.0),-t8.*(t43+m.*vx-C.*lf.*t9+Fz_f.*lr.*mu.*t9.*t13.*t18.*t29.*2.0-lf.*t7.*t20.*t21.*t30.*t31.*t32.*t34.*(1.0./3.0)),-t38.*(lr.*(-t43+C.*lf.*t9+lf.*t7.*t20.*t21.*t30.*t31.*t32.*t34.*(1.0./3.0))+Fz_f.*lf.*lr.*mu.*t9.*t13.*t18.*t29.*2.0)],[6,6]);
