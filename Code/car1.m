function car = car1()
car = SimModel();

% vehicle physical properties - Suburu BRZ
car.mass   = 1270;                                            % vehicle mass (kg)
car.length = 2.57;                                            % vehicle length (m)
car.width  = 1.54;                                            % vehicle width    (m)
car.dm     = .50;                                             % front/back load balance
car.Iz     = 2*1/12*car.mass*(car.length^2 + car.width^2);    % body moment of inertia about z axis (kg*m^2), 2 is a fundge factor
car.hcg    = .5;                                              % center of gravity height      (m)
car.wr     = .3048;                                           % wheel radius  (m)
car.Iw     = 1;                                               % tire/wheel/drive train rotational intertia  (kg*m^2)
car.mu_s   = 1;                                               % kinetic road/tire friction coefficient
car.mu_k   = 0.8;                                             % static road/tire friction coefficient
car.C      = 120000;                                          % tire stiffness   (N/rad)
car.F_long_accel   = 1000;                                    % max (positive) acceleration force on each wheel (N)
car.F_long_brake   = -50000;                                   % max (negative) braking force on each wheel (N)
car.delta_max = 45*pi/180;
car.delta_min = -45*pi/180;
car.v_max = 100;
car.v_min = 0;

car.refA   = 50.6*77.8*.0254^2;                               % aero reference area (m^2)
car.Cd     = .27;                                             % drag coefficient
car.rhoAir = 1.225;                                           % density of air (kg/m^3)

% vehicle initial conditions
car.X      = 0;     % global X coordinate
car.Y      = 0;     % global Y coordinate
car.phi    = 0;  % vehicle heading relative to global X axis  [rad]
car.vx     = 0;   % body x velocity, initialize at 1 m/s
car.vy     = 0;   % body y velocity
car.omegaB = 0;   % body rotation rate


car.Fx = 0 ;     % total x-force on vehicle from previous interation
car.Fy = 0;       % total y-force on vehicle from previous interation
end
