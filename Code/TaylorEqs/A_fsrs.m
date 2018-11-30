function A = A_fsrs(Fz_f,Fz_r,Iz,R,delta,lf,lr,m,mu,omegaB,varphi,vx,vy)
%A_FSRS
%    A = A_FSRS(FZ_F,FZ_R,IZ,R,DELTA,LF,LR,M,MU,OMEGAB,VARPHI,VX,VY)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    29-Nov-2018 20:39:17

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
t21 = lf.*omegaB;
t20 = -t21+vy;
t22 = t9.*t20;
t23 = atan(t22);
t24 = dirac(t23);
t25 = t20.^2;
t26 = t7.*t25;
t27 = t26+1.0;
t28 = 1.0./t27;
t29 = cos(delta);
t30 = 1.0./Iz;
A = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t3.*vx-t2.*vy,t2.*vx-t3.*vy,0.0,0.0,0.0,0.0,t2,t3,0.0,Fz_f.*R.*mu.*t5.*t7.*t8.*t13.*t14.*t18.*-2.0,t8.*(-t19+Fz_r.*R.*mu.*t7.*t20.*t24.*t28.*2.0+Fz_f.*R.*mu.*t5.*t7.*t13.*t18.*t29.*2.0),-t30.*(Fz_r.*R.*lr.*mu.*t7.*t20.*t24.*t28.*2.0-Fz_f.*R.*lf.*mu.*t5.*t7.*t13.*t18.*t29.*2.0),-t3,t2,0.0,t8.*(t19+Fz_f.*R.*mu.*t9.*t13.*t14.*t18.*2.0),-t8.*(Fz_r.*R.*mu.*t9.*t24.*t28.*2.0+Fz_f.*R.*mu.*t9.*t13.*t18.*t29.*2.0),t30.*(Fz_r.*R.*lr.*mu.*t9.*t24.*t28.*2.0-Fz_f.*R.*lf.*mu.*t9.*t13.*t18.*t29.*2.0),0.0,0.0,1.0,t8.*(m.*vy+Fz_f.*R.*lr.*mu.*t9.*t13.*t14.*t18.*2.0),-t8.*(m.*vx-Fz_r.*R.*lf.*mu.*t9.*t24.*t28.*2.0+Fz_f.*R.*lr.*mu.*t9.*t13.*t18.*t29.*2.0),-t30.*(Fz_r.*R.*lf.*lr.*mu.*t9.*t24.*t28.*2.0+Fz_f.*R.*lf.*lr.*mu.*t9.*t13.*t18.*t29.*2.0)],[6,6]);
