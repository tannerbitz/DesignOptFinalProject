syms theta ax bx cx dx ay by cy dy X Y 

eta = 1e-8;

fx = ax + bx*theta + cx*theta^2 + dx*theta^3;
fy = ay + by*theta + cy*theta^2 + dy*theta^3;
dfxdtheta = diff(fx, theta);
dfydtheta = diff(fy, theta);
Phi = atan(dfydtheta/ (dfxdtheta + eta));

ec = sin(Phi)*(X - fx) - cos(Phi)*(Y - fy);
el = -cos(Phi)*(X - fx) - sin(Phi)*(Y - fy);

ec2 = ec*ec;
el2 = el*el;

%% ec derivates
decdX = simplify(diff(ec2, X));
decdY = simplify(diff(ec2, Y));
decdtheta = simplify(diff(ec2, theta));

% states = [X, Y, vx, vy, varphi, omegaB]

Gamma_ec = [    diff(decdX, X),     diff(decdX, Y), diff(decdX, theta);
                diff(decdY, X),     diff(decdY, Y), diff(decdY, theta);
            diff(decdtheta, X), diff(decdtheta, Y), diff(decdtheta, theta)];


C_ec = [decdX;
        decdY;
        decdtheta];

matlabFunction(Gamma_ec, 'File', 'getGammaEc')
matlabFunction(C_ec, 'File', 'getCec');

%% el derivatives
deldX = simplify(diff(el2, X));
deldY = simplify(diff(el2, Y));
deldtheta = simplify(diff(el2, theta));


Gamma_el = [    diff(deldX, X),     diff(deldX, Y), diff(deldX, theta);
                diff(deldY, X),     diff(deldY, Y), diff(deldY, theta);
            diff(deldtheta, X), diff(deldtheta, Y), diff(deldtheta, theta)];

C_el = [deldX;
        deldY;
        deldtheta];

matlabFunction(Gamma_el, 'File', 'getGammaEl')
matlabFunction(C_el, 'File', 'getCel')