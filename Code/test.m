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
plot(time1,vy);
plot(time1,T);

legend('vx', 'vy')


