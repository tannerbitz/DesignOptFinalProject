close all; clear all;

endTime = 20;
Ts = .1;
N = 40;      % horizon length
inputFile = 'track1.txt';
width = 12;


rt = RaceTrack(inputFile,width);
rt.computeRaceTrack();
X0 = rt.X(1);
Y0 = rt.Y(1);

car = car1();  % intialize sim model physical properties and states
car.X = X0;
car.Y = Y0;

mpc = MPC(Ts,N);

time1 = 0:Ts:endTime;

% initialize vehicle state histories vectors
X = [];
Y = [];
phi = [];
vx = [];
vy = [];
omegaB = [];


F_long = 1000;
delta = 10;

lapCount = 0;
for n = 1:length(time1)
    t = time1(n);
    
    
    % calculate cubic approximation to track over the the prediction
    % horizon
    theta1 = mpc.thetaA;
    [theta1,X1,Y1,lapCount] = rt.getXY(theta1,lapCount);
    [ax,bx,cx,dx,ay,by,cy,dy] = rt.getCubicPolynomial(theta1,X1,Y1);
    
    
    
    
    
    % update simulation model states using ode45
    car.update(Ts,delta*pi/180,F_long);
    car.plotCar();

    
    % save history of states
    X = [X car.X];
    Y = [Y car.Y];
    phi = [phi car.phi];
    vx = [vx car.vx];
    vy = [vy car.vy];
    omegaB = [omegaB car.omegaB];
    %f_long = [f_long ];
    %delta  = [delta]

end

figure
plot(time1,vx,time1,vy);
hold on

