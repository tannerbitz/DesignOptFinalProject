close all; clear all;
inputFile = 'track1.txt';
width = 15;
rt = RaceTrack(inputFile,width);
rt.computeRaceTrack();

car = car1();
car.X = rt.X(1);
car.Y = rt.Y(1);
car.plotCar()

dt = .2;
time1 = 0:dt:20;

X = [];
Y = [];
phi = [];
vx = [];
vy = [];
omegaB = [];
omegaW_fl = [];
omegaW_fr = [];
omegaW_rl = [];
omegaW_rr = [];

force = 500/2;
delta = -1;

T = [];
for n = 1:length(time1)
    t = time1(n);
    T = [T ];
    if mod(100,n)
        t
    end
    t = time1(n);
        if  t > 7
            delta = 1;
            force = 1000/2;
        elseif t > 5
            delta = -2;
            force = 500/2;
        end
    tic
    car.update(dt,delta*pi/180,force);
    T = [T toc];
    car.plotCar();

    X = [X car.X];
    Y = [Y car.Y];
    phi = [phi car.phi];
    vx = [vx car.vx];
    vy = [vy car.vy];
    omegaB = [omegaB car.omegaB];

end

figure
plot(time1,vx);
hold on

start = 350;
ind = start:10:start+101;
start = 100;
ind = start:40:start+121;
aL = rt.arcLen(ind);
x = rt.X(ind);
y = rt.Y(ind);

for i = 1:length(aL)
   V(i,:) = [1 aL(i) aL(i)^2 aL(i)^3];
end

coef_x = inv(V'*V)*V'*x';
coef_y = inv(V'*V)*V'*y';

aLL =aL(1):.01:aL(end);

xx = coef_x(1) + coef_x(2).*aLL + coef_x(3).*aLL.^2 + coef_x(4).*aLL.^3;
yy = coef_y(1) + coef_y(2).*aLL + coef_y(3).*aLL.^2 + coef_y(4).*aLL.^3;

figure(1)
% plot(x,y,'o',xx,yy,'linewidth',5)
