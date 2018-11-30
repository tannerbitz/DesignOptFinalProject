close all;clear all;
R = .7;
mu = 1.0;

Fz_f = 300*9.8;

alphaf = [-10:.1:10]*pi/180;
C = 120000;


for i = 1:length(alphaf)
    falphaf = C*tan(alphaf(i));
    if abs(alphaf(i)) < atan2(3*mu*Fz_f,C)
        f_lat(i) = -falphaf + (2 - R)/(3*mu*Fz_f)*sqrt(falphaf^2)*falphaf - (1 - 2/3*R)/(3*mu*Fz_f)^2 * falphaf^3;
    else
        f_lat(i) = -sign(alphaf(i))*mu*R*Fz_f;
    end
end

plot(alphaf*180/pi,f_lat)