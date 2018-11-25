function Gamma_ec = getGammaEc(X,Y,ax,ay,bx,by,cx,cy,dx,dy,theta)
%GETGAMMAEC
%    GAMMA_EC = GETGAMMAEC(X,Y,AX,AY,BX,BY,CX,CY,DX,DY,THETA)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    24-Nov-2018 16:42:57

t3 = theta.^2;
t10 = cy.*theta.*2.0;
t11 = dy.*t3.*3.0;
t2 = by+t10+t11;
t4 = cx.*theta.*2.0;
t5 = dx.*t3.*3.0;
t6 = bx+t4+t5+1.0e-8;
t7 = cy.*2.0;
t8 = dy.*theta.*6.0;
t9 = t7+t8;
t12 = 1.0./t6.^2;
t13 = t2.^2;
t14 = t12.*t13;
t15 = t14+1.0;
t16 = 1.0./t6;
t17 = cx.*2.0;
t18 = dx.*theta.*6.0;
t19 = t17+t18;
t20 = 1.0./sqrt(t15);
t21 = 1.0./t6.^3;
t25 = t2.*t9.*t12.*2.0;
t26 = t13.*t19.*t21.*2.0;
t22 = t25-t26;
t23 = 1.0./t15.^(3.0./2.0);
t24 = t9.*t16.*t20;
t27 = t24-t2.*t12.*t19.*t20-t2.*t16.*t22.*t23.*(1.0./2.0);
t28 = t22.*t23.*(1.0./2.0);
t29 = by.*theta;
t30 = cy.*t3;
t31 = dy.*t3.*theta;
t32 = -Y+ay+t29+t30+t31;
t33 = bx.*theta;
t34 = cx.*t3;
t35 = dx.*t3.*theta;
t36 = -X+ax+t33+t34+t35;
t37 = t22.^2;
t38 = 1.0./t15.^(5.0./2.0);
t39 = t19.^2;
t40 = bx+t4+t5;
t41 = t9.^2;
t42 = t12.*t41.*2.0;
t43 = 1.0./t6.^4;
t44 = t13.*t39.*t43.*6.0;
t45 = dy.*t2.*t12.*1.2e1;
t46 = t42+t44+t45-dx.*t13.*t21.*1.2e1-t2.*t9.*t19.*t21.*8.0;
Gamma_ec = reshape([0.0,0.0,t27,0.0,0.0,t28,t27,t28,t9.*t20-t2.*t22.*t23-t23.*t32.*t46.*(1.0./2.0)+t32.*t37.*t38.*(3.0./4.0)-dy.*t16.*t20.*t36.*6.0-t2.*t16.*t19.*t20-t9.*t16.*t20.*t40.*2.0+dx.*t2.*t12.*t20.*t36.*6.0+t2.*t12.*t19.*t20.*t40.*2.0+t9.*t12.*t19.*t20.*t36.*2.0+t2.*t16.*t22.*t23.*t40+t9.*t16.*t22.*t23.*t36-t2.*t20.*t21.*t36.*t39.*2.0+t2.*t16.*t23.*t36.*t46.*(1.0./2.0)-t2.*t16.*t36.*t37.*t38.*(3.0./4.0)-t2.*t12.*t19.*t22.*t23.*t36],[3,3]);
