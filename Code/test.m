close all; clear all;
inputFile = 'track1.txt';
width = 12;
rt = RaceTrack(inputFile,width);
rt.computeRaceTrack();

car = car1();
car.plotCar()

dt = .1;
time1 = 0:dt:10;

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

torque = 2000;
for n = 1:length(time1)
    t = time1(n)
    if  t > 5
        torque = .000;
    end
    car.update(dt,0*pi/180,torque);
    car.plotCar();
    
    X = [X car.X];
    Y = [Y car.Y];
    phi = [phi car.phi];
    vx = [vx car.vx];
    vy = [vy car.vy];
    omegaB = [omegaB car.omegaB];
    omegaW_fl = [omegaW_fl car.omegaW_fl];
    omegaW_fr = [omegaW_fr car.omegaW_fr];
    omegaW_rl = [omegaW_rl car.omegaW_rl];
    omegaW_rr = [omegaW_rr car.omegaW_rr];
    
end

figure
plot(time1,vx);
hold on
plot(time1, omegaW_rl*car.wr, '.')
plot(time1,vy);

legend('vx', 'omegaW*wr', 'vy')


