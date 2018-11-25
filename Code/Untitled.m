clear all;
syms theta X Y phi vx vy omegaB  delta F_long v x0 y0 theta0 del_dx del_dy del_dtheta dec_dx dec_dy dec_dtheta ql qc el0 ec0

optVar = [theta X Y phi vx vy omegaB delta F_long v]; 

del = [del_dx, del_dy, del_dtheta];
dec = [dec_dx, dec_dy, dec_dtheta];

deltaX = [X-x0; Y-y0; theta-theta0];

el2 = ql*(2*el0*del*deltaX + (del*deltaX)^2);
ec2  = qc*(2*ec0*dec*deltaX + (dec*deltaX)^2);

f = el2 + ec2 ;
df = gradient(f,optVar);
ddf = [gradient(df(1), optVar), gradient(df(2), optVar), gradient(df(3), optVar)];


% plug in x,y,theta = 0 
df = subs(df, X, 0);
df = subs(df, Y, 0);
df = subs(df, theta, 0);
