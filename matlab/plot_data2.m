clear;

%data = importdata('../output4.dat');
data = importdata('../output/pdf_AD.dat');
data1 = importdata('../output/pdf_MC.dat');
%data = importdata('../anisotropicCompton.dat');
%data = importdata('../differentialFlux.dat');
%data = importdata('../outputSynch3.dat');
p = importdata('../output/p_grid.dat');

N = size(data,1);
Np = size(p,1);
Nx = N/Np;

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

plot(p(1:Np),data(1:Np),'Color','red');
plot(p(1:Np),data1(1:Np),'--','Color','red');

plot(p(1:Np),data(130*Np + 1:130*Np + Np),'Color','green');
plot(p(1:Np),data1(130*Np + 1:130*Np + Np),'--','Color','green');

plot(p(1:Np),data(152*Np + 1:152*Np + Np),'Color','blue');
plot(p(1:Np),data1(152*Np + 1:152*Np + Np),'--','Color','blue');

%loglog(data(1:N,1),data(1:N,2),'red','LineWidth',2,'Marker','+');
%loglog(data1(1:N,1),data1(1:N,2),'blue','LineWidth',2,'Marker','+');
%plot(data(1:N,1),data(1:N,3),'green','LineWidth',2,'Marker','+');
%plot(data(1:N,1),data(1:N,4),'magenta','LineWidth',2,'Marker','+');
%plot(data(1:N,1),approx(1:N),'blue','LineWidth',2,'Marker','+');
%plot(data(1:N,1)/(mc2),data(1:N,2),'red','LineWidth',2,'Marker','+');
