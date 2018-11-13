syms X Y vx vy varphi omegaB theta delta lf lr R mu Fz_f Fz_r Flong m Iz C

eps = 10^-8;
alphaf = atan((vy + lr*omegaB)/(vx+eps)) - delta;
alphar = atan((vy - lf*omegaB)/(vx+eps));

falphaf = C*tan(alphaf);
falphar = C*tan(alphar);

%% no slip
Flat_f = -falphaf + (2 - R)/(3*mu*Fz_f)*sqrt(falphaf.^2) - (1 - 2/3*R)/(3*mu*Fz_f)^2 * falphaf^3;
Flat_r = -falphar + (2 - R)/(3*mu*Fz_r)*sqrt(falphar.^2) - (1 - 2/3*R)/(3*mu*Fz_r)^2 * falphar^3;


Xdot = vx*cos(varphi) - vy*sin(varphi);
Ydot = vx*sin(varphi) + vy*cos(varphi);
varphidot = omegaB;
vxdot = 1/m*(Flong - Flat_f*sin(delta) + Flong*cos(delta) + m*vy*omegaB);
vydot = 1/m*(Flat_r + Flat_f*cos(delta) + Flong*sin(delta) - m*vx*omegaB);
omegaBdot = 1/Iz*(Flat_f*lf*cos(delta) + Flong*sin(delta)*lf - Flat_r*lr);


statesdot = [Xdot;
            Ydot;
            varphidot;
            vxdot;
            vydot;
            omegaBdot];

A = [diff(statesdot, X), ...
     diff(statesdot, Y), ...
     diff(statesdot, varphi), ...
     diff(statesdot, vx), ...
     diff(statesdot, vy), ...
     diff(statesdot, omegaB)];

B = [diff(statesdot, Flong), ...
     diff(statesdot, delta)];


matlabFunction(A, 'File', 'A_fgrg');
matlabFunction(B, 'File', 'B_fgrg');
matlabFunction(statesdot, 'File', 'statesdot_fgrg');


%% front slip, rear no slip
Flat_f = -sign(alphaf)*mu*Fz_f;
Flat_r = -falphar + (2 - R)/(3*mu*Fz_r)*sqrt(falphar.^2) - (1 - 2/3*R)/(3*mu*Fz_r)^2 * falphar^3;


Xdot = vx*cos(varphi) - vy*sin(varphi);
Ydot = vx*sin(varphi) + vy*cos(varphi);
varphidot = omegaB;
vxdot = 1/m*(Flong - Flat_f*sin(delta) + Flong*cos(delta) + m*vy*omegaB);
vydot = 1/m*(Flat_r + Flat_f*cos(delta) + Flong*sin(delta) - m*vx*omegaB);
omegaBdot = 1/Iz*(Flat_f*lf*cos(delta) + Flong*sin(delta)*lf - Flat_r*lr);


statesdot = [Xdot;
            Ydot;
            varphidot;
            vxdot;
            vydot;
            omegaBdot];

A = [diff(statesdot, X), ...
     diff(statesdot, Y), ...
     diff(statesdot, varphi), ...
     diff(statesdot, vx), ...
     diff(statesdot, vy), ...
     diff(statesdot, omegaB)];

B = [diff(statesdot, Flong), ...
     diff(statesdot, delta)];


matlabFunction(A, 'File', 'A_fsrg');
matlabFunction(B, 'File', 'B_fsrg');
matlabFunction(statesdot, 'File', 'statesdot_fsrg');

%% front no slip, rear slip
Flat_f = -falphaf + (2 - R)/(3*mu*Fz_f)*sqrt(falphaf.^2) - (1 - 2/3*R)/(3*mu*Fz_f)^2 * falphaf^3;
Flat_r = -sign(alphar)*mu*Fz_r;

Xdot = vx*cos(varphi) - vy*sin(varphi);
Ydot = vx*sin(varphi) + vy*cos(varphi);
varphidot = omegaB;
vxdot = 1/m*(Flong - Flat_f*sin(delta) + Flong*cos(delta) + m*vy*omegaB);
vydot = 1/m*(Flat_r + Flat_f*cos(delta) + Flong*sin(delta) - m*vx*omegaB);
omegaBdot = 1/Iz*(Flat_f*lf*cos(delta) + Flong*sin(delta)*lf - Flat_r*lr);


statesdot = [Xdot;
            Ydot;
            varphidot;
            vxdot;
            vydot;
            omegaBdot];

A = [diff(statesdot, X), ...
     diff(statesdot, Y), ...
     diff(statesdot, varphi), ...
     diff(statesdot, vx), ...
     diff(statesdot, vy), ...
     diff(statesdot, omegaB)];

B = [diff(statesdot, Flong), ...
     diff(statesdot, delta)];


matlabFunction(A, 'File', 'A_fgrs');
matlabFunction(B, 'File', 'B_fgrs');
matlabFunction(statesdot, 'File', 'statesdot_fgrs');


%% front slip, rear slip
Flat_f = -sign(alphaf)*mu*Fz_f;
Flat_r = -sign(alphar)*mu*Fz_r;


Xdot = vx*cos(varphi) - vy*sin(varphi);
Ydot = vx*sin(varphi) + vy*cos(varphi);
varphidot = omegaB;
vxdot = 1/m*(Flong - Flat_f*sin(delta) + Flong*cos(delta) + m*vy*omegaB);
vydot = 1/m*(Flat_r + Flat_f*cos(delta) + Flong*sin(delta) - m*vx*omegaB);
omegaBdot = 1/Iz*(Flat_f*lf*cos(delta) + Flong*sin(delta)*lf - Flat_r*lr);


statesdot = [Xdot;
            Ydot;
            varphidot;
            vxdot;
            vydot;
            omegaBdot];

A = [diff(statesdot, X), ...
     diff(statesdot, Y), ...
     diff(statesdot, varphi), ...
     diff(statesdot, vx), ...
     diff(statesdot, vy), ...
     diff(statesdot, omegaB)];

B = [diff(statesdot, Flong), ...
     diff(statesdot, delta)];


matlabFunction(A, 'File', 'A_fsrs');
matlabFunction(B, 'File', 'B_fsrs');
matlabFunction(statesdot, 'File', 'statesdot_fsrs');





