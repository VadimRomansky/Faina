clear;

data1 = importdata('../output2.dat');
data2 = importdata('../output3.dat');

data3 = importdata('../output1.dat');

N1 = size(data1,1);
N2 = size(data2,1);
N3 = size(data3,1);



set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');


plot(data1(1:N1,1),data1(1:N1,2),'red','LineWidth',2,'Marker','+');
plot(data2(1:N2,1),data2(1:N2,2),'blue','LineWidth',2,'Marker','+');
plot(data3(1:N3,1),data3(1:N3,2),'green','LineWidth',2,'Marker','+');
