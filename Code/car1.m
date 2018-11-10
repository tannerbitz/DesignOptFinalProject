function car = car1()
car = SimModel();

% vehicle physical properties - Suburu BRZ
car.mass   = 1270;                                            % vehicle mass (kg)
car.length = 2.57;                                            % vehicle length (m)
car.width  = 1.54;                                            % vehicle width    (m)
car.dm     = .55;                                             % front/back load balance
car.Iz     = 2*1/12*car.mass*(car.length^2 + car.width^2);    % body moment of inertia about z axis (kg*m^2), 2 is a fundge factor
car.hcg    = .5;                                              % center of gravity height      (m)
car.wr     = .3048;                                           % wheel radius  (m)
car.Iw     = .8;                                              % tire/wheel/drive train rotational intertia
car.maxT   = 205;                                             % maximum drivetrain torque    (N*m)
car.mu     = 1;                                               % road/tire friction coefficient (kg*m^2)
car.C      = 120000;                                          % tire stiffness   (N/rad)

car.refA   = 50.6*77.8*.0254^2;                               % aero reference area (m^2)
car.Cd     = .27;                                             % drag coefficient

% vehicle initial conditions
car.X      = 0;     % global X coordinate
car.Y      = 0;     % global Y coordinate
car.phi    = 0;  % vehicle heading relative to global X axis  [rad]
car.vx     = 0;   % body x velocity
car.vy     = 0;   % body y velocity
car.omegaB = 0;   % body rotation rate
car.omegaW_fl = 0; % front left wheel rotation rate
car.omegaW_fr = 0; % front right wheel rotation rate
car.omegaW_rl = 0;  % front left wheel rotation rate
car.omegaW_rr = 0; % front right wheel rotation rate

car.Fx = 0 ;     % total x-force on vehicle from previous interation
car.Fy = 0;       % total y-force on vehicle from previous interation


end