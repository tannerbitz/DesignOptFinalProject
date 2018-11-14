addpath ./TaylorEqs/

theta = 0.01;
rc = RaceTrack('track1.txt', 12);
rc.computeRaceTrack();
[theta_out,X_out,Y_out,lapCount] = rc.getXY(theta, 0)

[ec, el] = rc.getErrors(X_out, Y_out, 0.05);

ec2 = ec^2;
el2 = el^2;

theta_vec = 0:0.5:5;
[theta1_out,X1_out,Y1_out,lapCount] = rc.getXY(theta_vec, 0)
[ax,bx,cx,dx,ay,by,cy,dy] = rc.getCubicPolynomial(theta1_out, X1_out, Y1_out);

C_ec = Cec(X_out, Y_out, ax, ay, bx, by, cx, cy, dx, dy, theta_out);
C_el = Cel(X_out, Y_out, ax, ay, bx, by, cx, cy, dx, dy, theta_out);
Gamma_ec = GammaEc(X_out, Y_out, ax, ay, bx, by, cx, cy, dx, dy, theta_out);
Gamma_el = GammaEl(X_out, Y_out, ax, ay, bx, by, cx, cy, dx, dy, theta_out);

vec = [X_out; Y_out; theta_out];
ec2_hat = C_ec'*vec + vec'*Gamma_ec*vec
el2_hat = C_el'*vec + vec'*Gamma_el*vec



