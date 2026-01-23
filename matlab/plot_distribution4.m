clear;

data = importdata('../distribution.dat');
data1 = importdata('../distribution1.dat');
data2 = importdata('../distribution2.dat');
data3 = importdata('../distribution3.dat');
data4 = importdata('../distribution4.dat');


N = size(data1,1);




set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{E}');
xlabel ('E erg');
ylabel ('F_{E}');


hplank = 6.626E-27;

me = 9.1*10^-28;
mp = 1.6*10^-24;
c =3*10^10;
m=me;



loglog(data(1:N,1)/(m*c*c),data(1:N,2),'yellow','LineWidth',2);
loglog(data1(1:N,1)/(m*c*c),data1(1:N,2),'red','LineWidth',2);
loglog(data2(1:N,1)/(m*c*c),data2(1:N,2),'green','LineWidth',2);
loglog(data3(1:N,1)/(m*c*c),data3(1:N,2),'blue','LineWidth',2);
loglog(data4(1:N,1)/(m*c*c),data4(1:N,2),'magenta','LineWidth',2);
