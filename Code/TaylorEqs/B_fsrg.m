function B = B_fsrg(Flong,Fz_f,Iz,delta,lf,lr,m,mu,omegaB,vx,vy)

B = [                  0,                                                                                                                                                                                           0;
                       0,                                                                                                                                                                                           0;
                       0,                                                                                                                                                                                           0;
      (cos(delta) + 1)/m,          -(Flong*sin(delta) - Fz_f*mu*sign(atan((vy + lr*omegaB)/(vx + 1/100000000)) - delta)*cos(delta) + 2*Fz_f*mu*dirac(atan((vy + lr*omegaB)/(vx + 1/100000000)) - delta)*sin(delta))/m;
            sin(delta)/m,           (Flong*cos(delta) + 2*Fz_f*mu*dirac(atan((vy + lr*omegaB)/(vx + 1/100000000)) - delta)*cos(delta) + Fz_f*mu*sign(atan((vy + lr*omegaB)/(vx + 1/100000000)) - delta)*sin(delta))/m;
      (lf*sin(delta))/Iz, (Flong*lf*cos(delta) + 2*Fz_f*lf*mu*dirac(atan((vy + lr*omegaB)/(vx + 1/100000000)) - delta)*cos(delta) + Fz_f*lf*mu*sign(atan((vy + lr*omegaB)/(vx + 1/100000000)) - delta)*sin(delta))/Iz];
 
end