Ts = 0.05;

rt = RaceTrack('track1.txt', 12);
rt.computeRaceTrack();
X0 = rt.X(1);
Y0 = rt.Y(1);

% Inputs
delta = 5*pi/180;
Flong = 1000; 

% intialize sim model physical properties and states
car1 = car1();  
car1.X = X0;
car1.Y = Y0;
car1.updateOneTrack(Ts,delta, Flong);



printCar = @(car) fprintf(['X: %6.2d\n' ...
                             'Y: %6.2d\n' ...
                             'varphi: %6.2d\n' ...
                             'vx: %6.2d\n' ...
                             'vy: %6.2d\n' ...
                             'omegaB: %6.2d\n'], ...
                             car.X, car.Y, car.phi, car.vx, car.vy, car.omegaB);
                         
printCar(car1)

car2 = car1();
car2.X = X0;
car2.Y = Y0;
car2.updateOneTrackLinear(Ts, delta, Flong)

printCar(car2)
