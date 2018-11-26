function el = getEl(X,Y,ax,ay,bx,by,cx,cy,dx,dy,theta)
%GETEL
%    EL = GETEL(X,Y,AX,AY,BX,BY,CX,CY,DX,DY,THETA)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    25-Nov-2018 21:36:21

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
t13 = 1.0./sqrt(t12);
el = t13.*(-X+ax+bx.*theta+cx.*t3+dx.*t3.*theta)+(t2.*t13.*(-Y+ay+by.*theta+cy.*t3+dy.*t3.*theta))./t9;
