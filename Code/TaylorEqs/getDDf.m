function DDf = getDDf(X0,Y0,ax,ay,bx,by,cx,cy,dx,dy,qc,ql,rF_long,rdelta,rv,theta0)
%GETDDF
%    DDF = GETDDF(X0,Y0,AX,AY,BX,BY,CX,CY,DX,DY,QC,QL,RF_LONG,RDELTA,RV,THETA0)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    25-Nov-2018 21:36:37

t3 = theta0.^2;
t5 = cy.*theta0.*2.0;
t6 = dy.*t3.*3.0;
t2 = by+t5+t6;
t4 = conj(theta0);
t7 = t2.^2;
t8 = cx.*theta0.*2.0;
t9 = dx.*t3.*3.0;
t10 = bx+t8+t9+1.0e-8;
t11 = 1.0./t10.^2;
t12 = t7.*t11;
t13 = t12+1.0;
t14 = conj(dy);
t15 = t4.^2;
t16 = conj(cy);
t17 = conj(bx);
t18 = conj(dx);
t19 = t15.*t18.*3.0;
t20 = conj(cx);
t21 = t4.*t20.*2.0;
t22 = conj(by);
t23 = t14.*t15.*3.0;
t24 = t4.*t16.*2.0;
t25 = t22+t23+t24;
t26 = t17+t19+t21+1.0e-8;
t27 = sqrt(t13);
t28 = conj(t27);
t29 = 1.0./t28;
t30 = t25.^2;
t31 = t16.*2.0;
t32 = t4.*t14.*6.0;
t33 = t31+t32;
t34 = 1.0./t26;
t35 = t20.*2.0;
t36 = t4.*t18.*6.0;
t37 = t35+t36;
t38 = 1.0./t26.^2;
t39 = conj(Y0);
t40 = conj(ay);
t41 = t15.*t16;
t42 = t4.*t14.*t15;
t43 = t4.*t22;
t44 = -t39+t40+t41+t42+t43;
t45 = t13.^(3.0./2.0);
t46 = conj(t45);
t47 = 1.0./t46;
t48 = 1.0./t26.^3;
t51 = t25.*t33.*t38.*2.0;
t52 = t30.*t37.*t48.*2.0;
t49 = t51-t52;
t53 = t17+t19+t21;
t54 = conj(X0);
t55 = conj(ax);
t56 = t15.*t20;
t57 = t4.*t15.*t18;
t58 = t4.*t17;
t59 = -t54+t55+t56+t57+t58;
t62 = t29.*t53;
t63 = t47.*t49.*t59.*(1.0./2.0);
t64 = t29.*t30.*t34;
t65 = t29.*t33.*t34.*t44;
t66 = t25.*t29.*t37.*t38.*t44;
t67 = t25.*t34.*t44.*t47.*t49.*(1.0./2.0);
t50 = t62-t63+t64+t65-t66-t67;
t69 = t25.*t29;
t70 = t44.*t47.*t49.*(1.0./2.0);
t71 = t25.*t29.*t34.*t53;
t72 = t29.*t33.*t34.*t59;
t73 = t25.*t29.*t37.*t38.*t59;
t74 = t25.*t34.*t47.*t49.*t59.*(1.0./2.0);
t60 = t69-t70-t71-t72+t73+t74;
t61 = conj(ql);
t68 = conj(qc);
t75 = t25.*t29.*t34.*t60.*t68.*2.0;
t76 = t75-t29.*t50.*t61.*2.0;
t77 = t30.*t38;
t78 = t77+1.0;
t79 = 1.0./t78;
t80 = t29.*t60.*t68.*-2.0-t25.*t29.*t34.*t50.*t61.*2.0;
t81 = t25.*t34.*t61.*t79.*2.0;
t82 = t81-t25.*t34.*t68.*t79.*2.0;
DDf = reshape([t50.^2.*t61.*2.0+t60.^2.*t68.*2.0,t76,t80,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t76,t61.*t79.*2.0+t30.*t38.*t68.*t79.*2.0,t82,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t80,t82,t68.*t79.*2.0+t30.*t38.*t61.*t79.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,rdelta,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,rF_long,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,rv],[13,13]);
