close all; clear all;

N = 40;      % horizon length
inputFile = 'track1.txt';
width = 12;
rt = RaceTrack(inputFile,width);
rt.computeRaceTrack();

car = car1();
car.plotCar()

dt = .1;
<<<<<<< HEAD
time1 = 0:dt:10;
=======
time1 = 0:dt:20;
>>>>>>> 20c9ccf719cc3a794831c4c74486f6c7c590f2f5

% initialize vehicle state histories vectors
X = [];
Y = [];
phi = [];
vx = [];
vy = [];
omegaB = [];

<<<<<<< HEAD
torque = 2000;
for n = 1:length(time1)
    t = time1(n);
    if  t > 5
        torque = .000;
    end
    car.update(dt,5*pi/180,torque);
=======

%initialize prediction horizon vectors

force = 1000;
delta = 10;

T = [];
for n = 1:length(time1)
    t = time1(n);
    
    
    
    
    
    
    % update vehicle states using ode45
    car.update(dt,delta*pi/180,force);
>>>>>>> 20c9ccf719cc3a794831c4c74486f6c7c590f2f5
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

<<<<<<< HEAD
start = 100;
ind = start:20:start+121;
=======
start = 50;
ind = start:4:start+81;
>>>>>>> 20c9ccf719cc3a794831c4c74486f6c7c590f2f5
aL = rt.arcLen(ind);
x = rt.X(ind);
y = rt.Y(ind);

<<<<<<< HEAD
aLL = aL(1):.1:aL(end); 
xx = spline(aL,x,aLL);
yy = spline(aL,y,aLL);
figure(1)
plot(x,y,'o',xx,yy,'linewidth',5)

=======
for i = 1:length(aL)
   V(i,:) = [1 aL(i) aL(i)^2 aL(i)^3];
end

coef_x = inv(V'*V)*V'*x';
coef_y = inv(V'*V)*V'*y';
>>>>>>> 20c9ccf719cc3a794831c4c74486f6c7c590f2f5

aLL =aL(1):.01:aL(end);

xx = coef_x(1) + coef_x(2).*aLL + coef_x(3).*aLL.^2 + coef_x(4).*aLL.^3;
yy = coef_y(1) + coef_y(2).*aLL + coef_y(3).*aLL.^2 + coef_y(4).*aLL.^3;

figure(1)
plot(x,y,'o',xx,yy,'linewidth',5)
