function A = A_fgrg(C,Fz_f,Fz_r,Iz,R,delta,lf,lr,m,mu,omegaB,varphi,vx,vy)

A = [ 0, 0, - vy*cos(varphi) - vx*sin(varphi),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            cos(varphi),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            -sin(varphi),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0;
      0, 0,   vx*cos(varphi) - vy*sin(varphi),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            sin(varphi),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             cos(varphi),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0;
      0, 0,                                 0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        1;
      0, 0,                                 0,                                                                                                                                                                                                                                                                 (sin(delta)*((C^3*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2*((2*R)/3 - 1)*(vy + lr*omegaB)*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1))/(3*Fz_f^2*mu^2*((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)^2) - (C*(vy + lr*omegaB)*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1))/(((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)^2) + (C^2*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))*(vy + lr*omegaB)*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1)*(R - 2))/(3*Fz_f*mu*(C^2*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2)^(1/2)*((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)^2)))/m,                                                                                                                                                                                                                               -(sin(delta)*((C^3*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2*((2*R)/3 - 1)*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1))/(3*Fz_f^2*mu^2*((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)) - (C*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1))/(((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)) + (C^2*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1)*(R - 2))/(3*Fz_f*mu*(C^2*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2)^(1/2)*((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000))) - m*omegaB)/m,                                                                                                                                                                                                                                            (m*vy - sin(delta)*((C^3*lr*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2*((2*R)/3 - 1)*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1))/(3*Fz_f^2*mu^2*((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)) - (C*lr*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1))/(((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)) + (C^2*lr*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1)*(R - 2))/(3*Fz_f*mu*(C^2*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2)^(1/2)*((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000))))/m;
      0, 0,                                 0, -(m*omegaB + cos(delta)*((C^3*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2*((2*R)/3 - 1)*(vy + lr*omegaB)*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1))/(3*Fz_f^2*mu^2*((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)^2) - (C*(vy + lr*omegaB)*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1))/(((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)^2) + (C^2*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))*(vy + lr*omegaB)*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1)*(R - 2))/(3*Fz_f*mu*(C^2*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2)^(1/2)*((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)^2)) - (C*(vy - lf*omegaB))/(vx + 1/100000000)^2 + (C^3*((2*R)/3 - 1)*(vy - lf*omegaB)^3)/(3*Fz_r^2*mu^2*(vx + 1/100000000)^4) - (C^2*(vy - lf*omegaB)^2*(R - 2))/(3*Fz_r*mu*(vx + 1/100000000)^3*((C^2*(vy - lf*omegaB)^2)/(vx + 1/100000000)^2)^(1/2)))/m,          (cos(delta)*((C^3*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2*((2*R)/3 - 1)*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1))/(3*Fz_f^2*mu^2*((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)) - (C*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1))/(((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)) + (C^2*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1)*(R - 2))/(3*Fz_f*mu*(C^2*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2)^(1/2)*((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000))) - C/(vx + 1/100000000) + (C^3*((2*R)/3 - 1)*(vy - lf*omegaB)^2)/(3*Fz_r^2*mu^2*(vx + 1/100000000)^3) - (C^2*(2*vy - 2*lf*omegaB)*(R - 2))/(6*Fz_r*mu*(vx + 1/100000000)^2*((C^2*(vy - lf*omegaB)^2)/(vx + 1/100000000)^2)^(1/2)))/m,    (cos(delta)*((C^3*lr*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2*((2*R)/3 - 1)*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1))/(3*Fz_f^2*mu^2*((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)) - (C*lr*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1))/(((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)) + (C^2*lr*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1)*(R - 2))/(3*Fz_f*mu*(C^2*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2)^(1/2)*((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000))) - m*vx + (C*lf)/(vx + 1/100000000) - (C^3*lf*((2*R)/3 - 1)*(vy - lf*omegaB)^2)/(3*Fz_r^2*mu^2*(vx + 1/100000000)^3) + (C^2*lf*(vy - lf*omegaB)*(R - 2))/(3*Fz_r*mu*(vx + 1/100000000)^2*((C^2*(vy - lf*omegaB)^2)/(vx + 1/100000000)^2)^(1/2)))/m;
      0, 0,                                 0,   -(lr*((C*(vy - lf*omegaB))/(vx + 1/100000000)^2 - (C^3*((2*R)/3 - 1)*(vy - lf*omegaB)^3)/(3*Fz_r^2*mu^2*(vx + 1/100000000)^4) + (C^2*(vy - lf*omegaB)^2*(R - 2))/(3*Fz_r*mu*(vx + 1/100000000)^3*((C^2*(vy - lf*omegaB)^2)/(vx + 1/100000000)^2)^(1/2))) + lf*cos(delta)*((C^3*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2*((2*R)/3 - 1)*(vy + lr*omegaB)*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1))/(3*Fz_f^2*mu^2*((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)^2) - (C*(vy + lr*omegaB)*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1))/(((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)^2) + (C^2*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))*(vy + lr*omegaB)*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1)*(R - 2))/(3*Fz_f*mu*(C^2*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2)^(1/2)*((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)^2)))/Iz, (lr*(C/(vx + 1/100000000) - (C^3*((2*R)/3 - 1)*(vy - lf*omegaB)^2)/(3*Fz_r^2*mu^2*(vx + 1/100000000)^3) + (C^2*(2*vy - 2*lf*omegaB)*(R - 2))/(6*Fz_r*mu*(vx + 1/100000000)^2*((C^2*(vy - lf*omegaB)^2)/(vx + 1/100000000)^2)^(1/2))) + lf*cos(delta)*((C^3*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2*((2*R)/3 - 1)*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1))/(3*Fz_f^2*mu^2*((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)) - (C*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1))/(((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)) + (C^2*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1)*(R - 2))/(3*Fz_f*mu*(C^2*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2)^(1/2)*((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000))))/Iz, -(lr*((C*lf)/(vx + 1/100000000) - (C^3*lf*((2*R)/3 - 1)*(vy - lf*omegaB)^2)/(3*Fz_r^2*mu^2*(vx + 1/100000000)^3) + (C^2*lf*(vy - lf*omegaB)*(R - 2))/(3*Fz_r*mu*(vx + 1/100000000)^2*((C^2*(vy - lf*omegaB)^2)/(vx + 1/100000000)^2)^(1/2))) - lf*cos(delta)*((C^3*lr*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2*((2*R)/3 - 1)*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1))/(3*Fz_f^2*mu^2*((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)) - (C*lr*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1))/(((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000)) + (C^2*lr*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))*(tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2 + 1)*(R - 2))/(3*Fz_f*mu*(C^2*tan(delta - atan((vy + lr*omegaB)/(vx + 1/100000000)))^2)^(1/2)*((vy + lr*omegaB)^2/(vx + 1/100000000)^2 + 1)*(vx + 1/100000000))))/Iz];
 
end