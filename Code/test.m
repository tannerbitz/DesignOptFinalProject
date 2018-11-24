close all; clear all;

endTime = 10;
Ts = .1;
N = 40;      % horizon length
inputFile = 'track1.txt';
width = 12;

time1 = 0:Ts:endTime;

% compute racetrack using knots from inputsFile
rt = RaceTrack(inputFile,width);
rt.computeRaceTrack();
X0 = rt.X(1);
Y0 = rt.Y(1);

% intialize sim model physical properties and states
car = car1();  
car.X = X0;
car.Y = Y0;
car.vx = 10;

% initalize Model Pridictive Controller
mpc = MPC(Ts,N,X0,Y0,rt, car);

% initialize vehicle state histories vectors
X_hist = [];
Y_hist = [];
phi_hist = [];
vx_hist = [];
vy_hist = [];
omegaB_hist = [];
F_long_hist = [];
delta_hist  = [];
thetaA_hist = [];

lapCount = 0;

%[A, B] = mpc.getIneqCons(car);
[lb,ub] = mpc.getBounds(car,rt);

for n = 1:length(time1)
    t = time1(n);
    
    % calculate cubic approximation to track over the the prediction
    % horizon
    theta1 = mpc.thetaA;
    [theta1,X1,Y1,lapCount] = rt.getXY(theta1,lapCount);
    [ax,bx,cx,dx,ay,by,cy,dy] = rt.getCubicPolynomial(theta1,X1,Y1);
    

    %linearize prediction and error models 
    mpc.linearize(car, ax, ay, bx, by, cx, cy, dx, dy)
   
    % solve optimization problem for control inputs
    opt0 = zeros(1+mpc.N*3,1);
    opt0(1) = mpc.thetaA(1);
    opt0(2:3:end) = mpc.F_long;
    opt0(3:3:end) = mpc.delta;
    opt0(4:3:end) = mpc.v;
    
    %opt = fmincon(@mpc.cost2, opt0, A, B);
    ub(1) = max(theta1);
    opt = fmincon(@mpc.cost2, opt0, [], [],[],[],lb,ub);
    mpc.setStates(opt);
    
    
    % update simulation model states using ode45
    car.update(Ts,mpc.delta(1),mpc.F_long(1)/2);
    car.plotCar();

    % save history of states
    X_hist = [X_hist car.X];
    Y_hist = [Y_hist car.Y];
    phi_hist = [phi_hist car.phi];
    vx_hist = [vx_hist car.vx];
    vy_hist = [vy_hist car.vy];
    omegaB_hist = [omegaB_hist car.omegaB];
    F_long_hist = [F_long_hist mpc.F_long(1)];
    delta_hist  = [delta_hist mpc.delta(1)];
    thetaA_hist = [thetaA_hist, mpc.thetaA(1)];
    
    % set the new linearization states for the MPC controller
    States = [car.X,car.Y, car.phi,car.vx,car.vy,car.omegaB];
    mpc.setLinPoints(States);
    
end
%%
close all 
figure
plot(time1, vx_hist);
hold on
plot(time1, vy_hist);
hold off
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Velocity Components in Body Frame')
legend('v_x', 'v_y')

figure
plot(time1, X_hist)
hold on
plot(time1, Y_hist)
hold off
xlabel('Time (s)')
ylabel('Position (m)')
title('X and Y Parametric Movements')
legend('X', 'Y')

figure
plot(time1, F_long_hist);
xlabel('Time (s)')
ylabel('Force (N)')
title('Longitudinal Force on Each Wheel')


figure
plot(time1, delta_hist*180/pi);
xlabel('Time (s)')
ylabel('Steering Angle (deg)')
title('Steering Angle')

figure
plot(time1, phi_hist)
hold on
plot(time1, omegaB_hist)
hold off
xlabel('Time (s)')
ylabel('Angle/Anglular Rate')
title('Global Angle and Angular Rate')
legend({'$\varphi$', '$\omega_B$'}, 'Interpreter', 'latex')
