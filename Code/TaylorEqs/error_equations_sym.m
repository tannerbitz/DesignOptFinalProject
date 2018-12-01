clear all
syms theta ax bx cx dx ay by cy dy X Y 

eta = 1e-8;

fx = ax + bx*theta + cx*theta^2 + dx*theta^3;
fy = ay + by*theta + cy*theta^2 + dy*theta^3;
dfxdtheta = diff(fx, theta);
dfydtheta = diff(fy, theta);
Phi = atan(dfydtheta/ (dfxdtheta + eta));

ec = sin(Phi)*(X - fx) - cos(Phi)*(Y - fy);
el = -cos(Phi)*(X - fx) - sin(Phi)*(Y - fy);

matlabFunction(ec, 'File', 'getEc');
matlabFunction(el, 'File', 'getEl');

%% ec derivates
decdX = simplify(diff(ec, X));
decdY = simplify(diff(ec, Y));
decdtheta = simplify(diff(ec, theta));


%% el derivatives
deldX = simplify(diff(el, X));
deldY = simplify(diff(el, Y));
deldtheta = simplify(diff(el, theta));


%% Kenny's Addition
syms theta X Y phi vx vy omegaB  delta F_long v X0 Y0 theta0 ql qc el0 ec0 ddelta dF_long dv Ts rdelta rF_long rv gamma

optVar = [theta X Y phi vx vy omegaB delta F_long v ddelta dF_long dv]; 

%substite in decdX
decdX = subs(decdX, X, X0);
decdX = subs(decdX, Y, Y0);
decdX = subs(decdX, theta, theta0);

decdY = subs(decdY, X, X0);
decdY = subs(decdY, Y, Y0);
decdY = subs(decdY, theta, theta0);

decdtheta = subs(decdtheta, X, X0);
decdtheta = subs(decdtheta, Y, Y0);
decdtheta = subs(decdtheta, theta, theta0);

%substite in deldX
deldX = subs(deldX, X, X0);
deldX = subs(deldX, Y, Y0);
deldX = subs(deldX, theta, theta0);

deldY = subs(deldY, X, X0);
deldY = subs(deldY, Y, Y0);
deldY = subs(deldY, theta, theta0);

deldtheta = subs(deldtheta, X, X0);
deldtheta = subs(deldtheta, Y, Y0);
deldtheta = subs(deldtheta, theta, theta0);



% kennys code
del = [deldX, deldY, deldtheta];
dec = [decdX, decdY, decdtheta];

deltaX = [X-X0; Y-Y0; theta-theta0];

el2 = ql*(el0 + del*deltaX)^2;
ec2 =  qc*(ec0 + dec*deltaX)^2;

%ec2 = qc*ec^2;
%el2 = ql*el^2;

f = el2 + ec2 - gamma*v*Ts + rdelta*ddelta^2 + rF_long*dF_long^2 + rv*dv^2;
Df = gradient(f,optVar);

DDf = [gradient(Df(1), optVar), gradient(Df(2), optVar), gradient(Df(3), optVar),gradient(Df(4), optVar), gradient(Df(5), optVar), gradient(Df(6), optVar), ...
    gradient(Df(7), optVar), gradient(Df(8), optVar), gradient(Df(9), optVar),gradient(Df(10), optVar), gradient(Df(11), optVar), gradient(Df(12), optVar), ...
    gradient(Df(13), optVar)];

DDf = DDf';
%DDf = [DDf; zeros(length(optVar) - 3, length(optVar))];

%plug in x,y,theta = 0 
Df = subs(Df, X, 0);
Df = subs(Df, Y, 0);
Df = subs(Df, theta, 0);
Df = subs(Df, ddelta, 0);
Df = subs(Df, dF_long, 0);
Df = subs(Df, dv, 0);

matlabFunction(Df, 'File', 'getDf');
matlabFunction(DDf, 'File', 'getDDf');

