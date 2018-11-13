function Gamma_ec = getGammaEc(X,Y,ax,ay,bx,by,cx,cy,dx,dy,theta)
%GETGAMMAEC
%    GAMMA_EC = GETGAMMAEC(X,Y,AX,AY,BX,BY,CX,CY,DX,DY,THETA)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    12-Nov-2018 19:11:04

t3 = theta.^2;
t4 = cy.*theta.*2.0;
t5 = dy.*t3.*3.0;
t2 = by+t4+t5;
t6 = t2.^2;
t7 = cx.*theta.*2.0;
t8 = dx.*t3.*3.0;
t9 = bx+t7+t8+1.0e-8;
t10 = 1.0./t9.^2;
t11 = t6.*t10;
t12 = t11+1.0;
t13 = 1.0./t12;
t14 = 1.0./t9;
t15 = 1.0./sqrt(t12);
t16 = cy.*2.0;
t17 = dy.*theta.*6.0;
t18 = t16+t17;
t19 = cx.*2.0;
t20 = dx.*theta.*6.0;
t21 = t19+t20;
t22 = 1.0./t9.^3;
t33 = t2.*t10.*t18.*2.0;
t34 = t6.*t21.*t22.*2.0;
t23 = t33-t34;
t24 = 1.0./t12.^(3.0./2.0);
t25 = bx.*theta;
t26 = cx.*t3;
t27 = dx.*t3.*theta;
t28 = -X+ax+t25+t26+t27;
t29 = by.*theta;
t30 = cy.*t3;
t31 = dy.*t3.*theta;
t32 = -Y+ay+t29+t30+t31;
t35 = t15.*t32;
t37 = t2.*t14.*t15.*t28;
t36 = t35-t37;
t38 = t2.*t15;
t39 = bx+t7+t8;
t40 = t2.*t14.*t23.*t24.*t28.*(1.0./2.0);
t41 = t2.*t10.*t15.*t21.*t28;
t43 = t23.*t24.*t32.*(1.0./2.0);
t44 = t2.*t14.*t15.*t39;
t45 = t14.*t15.*t18.*t28;
t42 = t38+t40+t41-t43-t44-t45;
t46 = t2.*t14.*t15.*t42.*2.0;
t47 = t23.*t24.*t36;
t48 = t47-t15.*t42.*2.0;
t49 = t23.^2;
t50 = 1.0./t12.^(5.0./2.0);
t51 = t21.^2;
t52 = t18.^2;
t53 = t10.*t52.*2.0;
t54 = 1.0./t9.^4;
t55 = t6.*t51.*t54.*6.0;
t56 = dy.*t2.*t10.*1.2e1;
t57 = t53+t55+t56-dx.*t6.*t22.*1.2e1-t2.*t18.*t21.*t22.*8.0;
Gamma_ec = reshape([t6.*t10.*t13.*2.0,t2.*t13.*t14.*-2.0,t46-t36.*(-t14.*t15.*t18+t2.*t10.*t15.*t21+t2.*t14.*t23.*t24.*(1.0./2.0)).*2.0,t2.*t13.*t14.*-2.0,t13.*2.0,t48,t46+t14.*t15.*t18.*t36.*2.0-t2.*t10.*t15.*t21.*t36.*2.0-t2.*t14.*t23.*t24.*t36,t48,t36.*(t15.*t18-t2.*t23.*t24-t24.*t32.*t57.*(1.0./2.0)+t32.*t49.*t50.*(3.0./4.0)-dy.*t14.*t15.*t28.*6.0-t2.*t14.*t15.*t21-t14.*t15.*t18.*t39.*2.0+dx.*t2.*t10.*t15.*t28.*6.0+t2.*t10.*t15.*t21.*t39.*2.0+t10.*t15.*t18.*t21.*t28.*2.0+t2.*t14.*t23.*t24.*t39+t14.*t18.*t23.*t24.*t28-t2.*t15.*t22.*t28.*t51.*2.0+t2.*t14.*t24.*t28.*t57.*(1.0./2.0)-t2.*t14.*t28.*t49.*t50.*(3.0./4.0)-t2.*t10.*t21.*t23.*t24.*t28).*2.0+t42.^2.*2.0],[3,3]);
