clear;

radiation = importdata('../css161010.dat');

N = size(radiation,1);
Nr = size(radiation,2);

css0 = importdata('../examples_data/css_data/coppejans69.txt');
css1 = importdata('../examples_data/css_data/coppejans99.txt');
css2 = importdata('../examples_data/css_data/coppejans357.txt');

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

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{\nu}');
xlabel ('{\nu} GHz');
ylabel ('mJy');

loglog(radiation(1:N,1),radiation(1:N,2),'red','LineWidth',2);
loglog(radiation(1:N,1),radiation(1:N,3),'green','LineWidth',2);
%loglog(radiation(1:N,1),radiation(1:N,4),'blue','LineWidth',2);


%errorbar(css0(:,5),css0(:,7),css0(:,8),'red','LineWidth',2,'LineStyle','none');
errorbar(css1(:,5),css1(:,7),css1(:,8),'green','LineWidth',2,'LineStyle','none');
errorbar(css2(:,5),css2(:,7),css2(:,8),'blue','LineWidth',2,'LineStyle','none');

legend('69 days', '99 days', '357 days');

xlim([0.1 200]);
ylim([0.01 20]);