close all; clear all;
inputFile = 'track1.txt';
width = .1;
rt = RaceTrack(inputFile,width);

rt.computeRaceTrack();

plot(rt.X(1), rt.Y(1), 'o')
