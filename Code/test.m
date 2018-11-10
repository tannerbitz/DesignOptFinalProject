close all; clear all;
inputFile = 'track1.txt';
width = 12;
rt = RaceTrack(inputFile,width);
rt.computeRaceTrack();

car = car1();
car.plotCar()

dt = .1;
time1 = 0:dt:10;
vx = [];
torque = 1;
for n = 1:length(time1)
    t = time1(n);
    if  t > .1
        torque = .001;
    end
    car.update(dt,20*pi/180,torque);
    car.plotCar();
    
    vx = [vx car.vx];
end

figure
plot(time1,vx);

