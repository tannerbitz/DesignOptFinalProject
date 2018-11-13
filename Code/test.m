close all; clear all;

endTime = 20;
Ts = .1;
N = 40;      % horizon length
inputFile = 'track1.txt';
width = 12;


rt = RaceTrack(inputFile,width);
rt.computeRaceTrack();
car = car1();
mpc = MPC(Ts,N);

time1 = 0:Ts:endTime;

% initialize vehicle state histories vectors
X = [];
Y = [];
phi = [];
vx = [];
vy = [];
omegaB = [];

%initialize prediction horizon vectors


F_long = 1000;
delta = 10;

lapCount = 0;
for n = 1:length(time1)
    t = time1(n);
    
    theta1 = MPC.theta;
    [theta1,X1,Y1,lapCount] = rt.getXY(theta1,lapCount);
    [ax,bx,cx,dx,ay,by,cy,dy] = rt.getCubicPolynomial(theta1,X1,Y1);
    
    
    % update vehicle states using ode45
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
plot(time1,vx);
hold on


start = 50;
ind = start:4:start+81;
aL = rt.arcLen(ind);
x = rt.X(ind);
y = rt.Y(ind);


aLL = aL(1):.1:aL(end); 
xx = spline(aL,x,aLL);
yy = spline(aL,y,aLL);
figure(1)
plot(x,y,'o',xx,yy,'linewidth',5)

for i = 1:length(aL)
   V(i,:) = [1 aL(i) aL(i)^2 aL(i)^3];
end

coef_x = inv(V'*V)*V'*x';
coef_y = inv(V'*V)*V'*y';

aLL =aL(1):.01:aL(end);

xx = coef_x(1) + coef_x(2).*aLL + coef_x(3).*aLL.^2 + coef_x(4).*aLL.^3;
yy = coef_y(1) + coef_y(2).*aLL + coef_y(3).*aLL.^2 + coef_y(4).*aLL.^3;

figure(1)
plot(x,y,'o',xx,yy,'linewidth',5)
