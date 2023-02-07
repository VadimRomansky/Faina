clear;

sigma = importdata('outputSigma.dat');

N = size(sigma,1);

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
title ('F_{E}');
xlabel ('E GeV');
%ylabel ('cm^{-2} erg^{-1} s^{-1}');
ylabel ('\sigma mb/GeV');

plot(sigma(1:N,1)/1E9,sigma(1:N,2)*1E27/(1.6E-3),'red','LineWidth',2);
