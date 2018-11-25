syms theta ax bx cx dx ay by cy dy X Y 

eta = 1e-8;

fx = ax + bx*theta + cx*theta^2 + dx*theta^3;
fy = ay + by*theta + cy*theta^2 + dy*theta^3;
dfxdtheta = diff(fx, theta);
dfydtheta = diff(fy, theta);
Phi = atan(dfydtheta/ (dfxdtheta + eta));

ec = sin(Phi)*(X - fx) - cos(Phi)*(Y - fy);
el = -cos(Phi)*(X - fx) - sin(Phi)*(Y - fy);

%% ec derivates
decdX = simplify(diff(ec, X));
decdY = simplify(diff(ec, Y));
decdtheta = simplify(diff(ec, theta));

% states = [X, Y, vx, vy, varphi, omegaB]

C_ec = [decdX;
        decdY;
        decdtheta];
matlabFunction(C_ec, 'File', 'getCec');

%% el derivatives
deldX = simplify(diff(el, X));
deldY = simplify(diff(el, Y));
deldtheta = simplify(diff(el, theta));

C_el = [deldX;
        deldY;
        deldtheta];
    
matlabFunction(C_el, 'File', 'getCel')