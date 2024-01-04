clear;

%data = importdata('../output4.dat');
%data = importdata('../distribution.dat');
%data = importdata('../anisotropicCompton.dat');
%data = importdata('../differentialFlux.dat');
%data = importdata('../bremsstrahlung.dat');
data = importdata('../output.dat');
data7 = importdata('../output7.dat');

N = size(data,1);
N7 = size(data7,1);

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
title ('E F_{E}');
xlabel ('E keV');
ylabel ('E F_{E} erg cm^{-2} s^{-1}');

data1(1:N) = 0;
for i = 1:N,
    data1(i) = data(i,1)*data(i,2);
end;

for i = 1:N7,
    data7(i,2) = data7(i,1)*data7(i,2)/(1.6*10^-9);
end;

mc2 = (9.1*10^-28) * (3*10^10)^2;
hplank = 6.626E-27;

%startPower = 200;
%endPower = 300;
%Fpa(1:N) = 0;

%Fpa(startPower) = data1(startPower);
%Fpa(endPower) = data1(endPower);
%polyfitx(1:endPower-startPower + 1) = 0;
%polyfity(1:endPower-startPower + 1) = 0;

%for i = 1:endPower-startPower + 1,
%    polyfitx(i) = log(data(i+startPower-1,1));
%    %polyfitx(i) = log((me*energy(i+startPower - 1)));
%    polyfity(i) = log(data1(i+startPower-1));
%end;
%p = polyfit(polyfitx, polyfity, 1);

%for i = startPower-5:endPower+5,
%    Fpa(i) = exp(polyval(p, log(data(i,1))));
%end;

loglog(data(1:N,1),data1(1:N),'red','LineWidth',2,'Marker','+');
loglog(data7(1:N7,1)/(1.6*10^-9), data7(1:N7,2),'blue','LineWidth',2,'Marker','+');
%loglog(data(1:N,1)/(1.6*10^-9),Fpa(1:N),'green','LineWidth',2);
%plot(data(1:N,1),data(1:N,3),'green','LineWidth',2,'Marker','+');
%plot(data(1:N,1),data(1:N,4),'magenta','LineWidth',2,'Marker','+');
%plot(data(1:N,1),approx(1:N),'blue','LineWidth',2,'Marker','+');
%plot(data(1:N,1)/(mc2),data(1:N,2),'red','LineWidth',2,'Marker','+');
