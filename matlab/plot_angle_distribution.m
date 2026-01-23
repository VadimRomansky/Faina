clear;

data = importdata('../angledistribution.dat');

N = size(data,1);



set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
set(gca, 'YScale', 'log');
%set(gca, 'XScale', 'log');
title ('F_{\mu}');
xlabel ('\mu');
ylabel ('F_{\mu}');

mc2 = (9.1*10^-28) * (3*10^10)^2;
hplank = 6.626E-27;


loglog(data(1:N,1),data(1:N,2),'red','LineWidth',2,'Marker','+');
loglog(data(1:N,1),data(1:N,3),'green','LineWidth',2,'Marker','+');
loglog(data(1:N,1),data(1:N,4),'blue','LineWidth',2,'Marker','+');
loglog(data(1:N,1),data(1:N,5),'black','LineWidth',2,'Marker','+');
loglog(data(1:N,1),data(1:N,6),'magenta','LineWidth',2,'Marker','+');
loglog(data(1:N,1),data(1:N,7),'yellow','LineWidth',2,'Marker','+');

