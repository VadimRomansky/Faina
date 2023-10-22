clear;

%data = importdata('../output4.dat');
%data = importdata('../distribution.dat');
%data = importdata('../anisotropicCompton.dat');
%data = importdata('../differentialFlux.dat');
%data = importdata('../bremsstrahlung.dat');
data = importdata('../outputCompton.dat');

N = size(data,1);

%approx(1:N) = 0;
%for i=1:N,
%    approx(i) = 0.5*data(1,2)*(1 + (cos(data(i,1))^2));
%end;



set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{E}');
xlabel ('E эрг');
ylabel ('F_{E} см^{-2} с^{-1}');

mc2 = (9.1*10^-28) * (3*10^10)^2;

loglog(data(1:N,1)/(1.6*10^-9),data(1:N,2),'red','LineWidth',2,'Marker','+');
%plot(data(1:N,1),data(1:N,3),'green','LineWidth',2,'Marker','+');
%plot(data(1:N,1),data(1:N,4),'magenta','LineWidth',2,'Marker','+');
%plot(data(1:N,1),approx(1:N),'blue','LineWidth',2,'Marker','+');
%plot(data(1:N,1)/(mc2),data(1:N,2),'red','LineWidth',2,'Marker','+');
