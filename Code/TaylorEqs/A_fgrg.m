function A = A_fgrg(C,Fz_f,Fz_r,Iz,R,delta,lf,lr,m,mu,omegaB,varphi,vx,vy)
%A_FGRG
%    A = A_FGRG(C,FZ_F,FZ_R,IZ,R,DELTA,LF,LR,M,MU,OMEGAB,VARPHI,VX,VY)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    30-Nov-2018 14:12:27

t2 = cos(varphi);
t3 = sin(varphi);
t4 = lf.*omegaB;
t5 = t4+vy;
t7 = vx+1.0e-8;
t10 = 1.0./t7;
t11 = t5.*t10;
t12 = atan(t11);
t13 = delta-t12;
t6 = tan(t13);
t8 = 1.0./t7.^2;
t9 = C.^2;
t14 = t6.^2;
t15 = t14+1.0;
t16 = t5.^2;
t17 = t8.*t16;
t18 = t17+1.0;
t19 = 1.0./t18;
t20 = 1.0./Fz_f;
t21 = 1.0./mu;
t22 = t9.*t14;
t23 = R-2.0;
t24 = 1.0./m;
t25 = sin(delta);
t26 = sqrt(t22);
t27 = 1.0./Fz_f.^2;
t28 = 1.0./mu.^2;
t29 = R.*(2.0./3.0);
t30 = t29-1.0;
t31 = 1.0./sqrt(t22);
t32 = C.*t5.*t8.*t15.*t19;
t33 = C.*t5.*t8.*t15.*t19.*t20.*t21.*t23.*t26.*(1.0./3.0);
t34 = C.*t5.*t8.*t9.*t14.*t15.*t19.*t20.*t21.*t23.*t31.*(1.0./3.0);
t58 = C.*t5.*t8.*t9.*t14.*t15.*t19.*t27.*t28.*t30.*(1.0./3.0);
t35 = t32+t33+t34-t58;
t38 = lr.*omegaB;
t36 = -t38+vy;
t37 = t36.^2;
t39 = 1.0./Fz_r;
t40 = 1.0./t7.^4;
t41 = t8.*t9.*t37;
t42 = cos(delta);
t43 = C.*t10.*t15.*t19;
t44 = C.*t10.*t15.*t19.*t20.*t21.*t23.*t26.*(1.0./3.0);
t45 = C.*t9.*t10.*t14.*t15.*t19.*t20.*t21.*t23.*t31.*(1.0./3.0);
t65 = C.*t9.*t10.*t14.*t15.*t19.*t27.*t28.*t30.*(1.0./3.0);
t46 = t43+t44+t45-t65;
t47 = sqrt(t41);
t48 = 1.0./Fz_r.^2;
t49 = 1.0./t7.^3;
t50 = 1.0./sqrt(t41);
t51 = C.*lf.*t10.*t15.*t19;
t52 = C.*lf.*t10.*t15.*t19.*t20.*t21.*t23.*t26.*(1.0./3.0);
t53 = C.*lf.*t9.*t10.*t14.*t15.*t19.*t20.*t21.*t23.*t31.*(1.0./3.0);
t69 = C.*lf.*t9.*t10.*t14.*t15.*t19.*t27.*t28.*t30.*(1.0./3.0);
t54 = t51+t52+t53-t69;
t55 = C.*t8.*t36;
t56 = C.*t8.*t21.*t23.*t36.*t39.*t47.*(1.0./3.0);
t57 = C.*t9.*t21.*t23.*t36.*t37.*t39.*t40.*t50.*(1.0./3.0);
t59 = 1.0./Iz;
t60 = C.*t10;
t61 = C.*t10.*t21.*t23.*t39.*t47.*(1.0./3.0);
t62 = vy.*2.0;
t63 = t62-lr.*omegaB.*2.0;
t64 = C.*t9.*t21.*t23.*t36.*t39.*t49.*t50.*t63.*(1.0./6.0);
t66 = C.*lr.*t10;
t67 = C.*lr.*t10.*t21.*t23.*t39.*t47.*(1.0./3.0);
t68 = C.*lr.*t9.*t21.*t23.*t37.*t39.*t49.*t50.*(1.0./3.0);
A = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t3.*vx-t2.*vy,t2.*vx-t3.*vy,0.0,0.0,0.0,0.0,t2,t3,0.0,-t24.*t25.*t35,t24.*(t55+t56+t57+t35.*t42-C.*t9.*t28.*t30.*t36.*t37.*t40.*t48.*(1.0./3.0)),-t59.*(lr.*(t55+t56+t57-C.*t9.*t28.*t30.*t36.*t37.*t40.*t48.*(1.0./3.0))-lf.*t35.*t42),-t3,t2,0.0,t24.*t25.*t46,-t24.*(t60+t61+t64+t42.*t46-C.*t9.*t28.*t30.*t37.*t48.*t49.*(1.0./3.0)),t59.*(lr.*(t60+t61+t64-C.*t9.*t28.*t30.*t37.*t48.*t49.*(1.0./3.0))-lf.*t42.*t46),0.0,0.0,1.0,t24.*t25.*t54,t24.*(t66+t67+t68-t42.*t54-C.*lr.*t9.*t28.*t30.*t37.*t48.*t49.*(1.0./3.0)),-t59.*(lr.*(t66+t67+t68-C.*lr.*t9.*t28.*t30.*t37.*t48.*t49.*(1.0./3.0))+lf.*t42.*t54)],[6,6]);
