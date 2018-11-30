close all; clear all;

endTime = 1;
Ts = .02;
N = 40;      % horizon length
inputFile = 'track1.txt';
width = 12;
nVarsPerIter = 13;

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
mpc = MPC(Ts,N,car);

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
    

    %linearize dynamics model 
    mpc.linearizeModel(car)
   
    % make optVar vector
    optVar0 = zeros(nVarsPerIter*N, 1);
    optVar0(1:nVarsPerIter:end) = mpc.thetaA;
    for i = 1:N
        optVar0(2+(i-1)*nVarsPerIter:+7+(i-1)*nVarsPerIter) = mpc.States(:,i);
    end
    optVar0(8:nVarsPerIter:end)  = mpc.delta;
    optVar0(9:nVarsPerIter:end)  = mpc.F_long;
    optVar0(10:nVarsPerIter:end) = mpc.v;
    optVar0(11:nVarsPerIter:end) = mpc.ddelta;
    optVar0(12:nVarsPerIter:end) = mpc.dF_long;
    optVar0(13:nVarsPerIter:end) = mpc.dv;
    
    % solve optimization problem to yield inputs [delta, Flong]
    [Aeq, Beq] = mpc.getEqualityCons();
    [G,H] = mpc.getCostMatrices(ax, ay, bx, by, cx, cy, dx, dy);
    
    optVar = quadprog(H,G,[],[],Aeq,Beq,lb,ub);
    
    % update variable vectors
    mpc.thetaA(1,:) = optVar(1:nVarsPerIter:end) ;
    for i = 1:N
         mpc.States(:,i) = optVar(2+(i-1)*nVarsPerIter:+7+(i-1)*nVarsPerIter);
    end
    mpc.delta(1,:) = optVar(8:nVarsPerIter:end);
    mpc.F_long(1,:) = optVar(9:nVarsPerIter:end);
    mpc.v(1,:) = optVar(10:nVarsPerIter:end) ;
    mpc.ddelta(1,:) = optVar(11:nVarsPerIter:end);
    mpc.dF_long(1,:) = optVar(12:nVarsPerIter:end);
    mpc.dv(1,:) = optVar(13:nVarsPerIter:end);
    
    plot(mpc.States(1,:), mpc.States(2,:), 'b.-')
    % update simulation model states using ode45
    car.update(Ts,mpc.delta(1),mpc.F_long(1)/2);
    car.plotCar();
    drawnow();

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
%close all 
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
