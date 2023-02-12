clear;

radiation = importdata('../outputPionE.dat');

N = size(radiation,1);

r = 5.17E16;
d = 40*3*10^24;
B = 0.346;
n = 2912;
me = 0.91*10^-27;
c = 3*10^10;
c1 = 6.27*10^18;
c5 = 7.25*10^-24;
c6 = 7.97*10^-41;
g = 3;
E0 = me*c*c;
N0 = n*(g-1)*(E0^(g-1));
f = 0.5;
s = 4*f*r/3;

h = 6.626*10^-27;


set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(2);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
set(gca, 'xlim',[0.01,30000]);
set(gca, 'ylim',[1E-14,1E-10]);
title ('F_{E}');
xlabel ('E GeV');
%ylabel ('cm^{-2} erg^{-1} s^{-1}');
ylabel ('E^2 dN/dE TeV cm^{-2} s^{-1}');
fe2(1:N) = 0;
for i = 1:N,
    fe2(i) = radiation(i,1)*radiation(i,1)*radiation(i,2)*1.6E-12/1E24;
end;

startPower = 120;
endPower = 150;

polyfitx(1:endPower-startPower + 1) = 0;
polyfity(1:endPower-startPower + 1) = 0;

for i = 1:endPower-startPower + 1,
    polyfitx(i) = log(radiation(i+startPower - 1,1));
    %polyfitx(i) = log((me*energy(i+startPower - 1)));
    polyfity(i) = log(fe2(i+startPower - 1));
end;
p = polyfit(polyfitx, polyfity, 1);

%ap = exp(log(Fpa(startPower)) - gammap*log((me*energy(startPower)+m)));

Fpa(1:N) = 0;
for i = startPower-5:endPower+5,
    Fpa(i) = exp(polyval(p, log(radiation(i,1))));
    %Fpa(i) = exp(polyval(p, log(me*energy(i))));
end;

plot(radiation(1:N,1)/1E9,fe2(1:N),'red','LineWidth',2);
plot(radiation(1:N,1)/1E9,Fpa(1:N),'blue','LineWidth',2);