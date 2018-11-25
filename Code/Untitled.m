clear all;
syms x y theta  x0 y0 theta0 del_dx del_dy del_dtheta dec_dx dec_dy dec_dtheta ql qc el0 ec0

del = [del_dx, del_dy, del_dtheta];
dec = [dec_dx, dec_dy, dec_dtheta];

deltaX = [x-x0; y-y0; theta-theta0];

el2 = ql*(2*el0*del*deltaX + (del*deltaX)^2);
ec2  = qc*(2*ec0*dec*deltaX + (dec*deltaX)^2);

ff = el2 + ec2;
dff = gradient(ff,[x y theta]);
ddff = [gradient(dff(1),[x y theta]), gradient(dff(2), [x y theta]), gradient(dff(3), [x y theta])];


% plug in x,y,theta = 0 
dff = subs(dff, x, 0);
dff = subs(dff, y, 0);
dff = subs(dff, theta, 0);
