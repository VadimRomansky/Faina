clear;

data = importdata('../output1.dat');

N = size(data,1);



set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{E}');
xlabel ('E эрг');
ylabel ('F_{E} см^{-2} с^{-1}');

m_e = 9.1*10^-28;
c=3*10^10;
plot(data(1:N,1)/(m_e*c*c),data(1:N,2),'red','LineWidth',2);
