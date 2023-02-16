clear;

data = importdata('../output.dat');

N = size(data,1);



set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
set(gca, 'YScale', 'log');
%set(gca, 'XScale', 'log');


plot(data(1:N,1),data(1:N,2),'red','LineWidth',2);
