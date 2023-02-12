clear;

radiation = importdata('../outputBremNu.dat');

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

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{\nu}');
xlabel ('{\nu} GHz');
ylabel ('mJy');

loglog(radiation(1:N,1),radiation(1:N,2),'red','LineWidth',2);

%loglog(radiation(1:N,1),radiation(1:N,5),'magenta','LineWidth',2);
%plot(radiation(1:N,1),radiation(1:N,6),'Color',[1.0,0.6,0],'LineWidth',2);
%plot(radiation(1:N,1),radiation(1:N,7),'black','LineWidth',2);

%loglog(aprx(1:4),apry(1:4),'--o','Color','red','LineWidth',2);
%loglog(mayx(1:3),mayy(1:3),'--o','Color','green','LineWidth',2);
%loglog(junx(1:4),juny(1:4),'--o','Color','blue','LineWidth',2);
%loglog(augx(1:5),augy(1:5),'--o','Color','magenta','LineWidth',2);
%plot(octx(1:3),octy(1:3),'--o','Color',[1.0,0.6,0],'LineWidth',2);
%plot(decx(1:3),decy(1:3),'--o','Color','black','LineWidth',2);

%errorbar(cssx1,cssy1,cssError1,'red','LineWidth',2);
%errorbar(cssx2,cssy2,cssError2,'green','LineWidth',2);
%errorbar(cssx3,cssy3,cssError3,'blue','LineWidth',2);

%errorbar(aprx,apry,aprerr,'red','LineWidth',2);
%errorbar(mayx,mayy,mayerr,'green','LineWidth',2);
%errorbar(junx,juny,junerr,'blue','LineWidth',2);
%errorbar(augx,augy,augerr,'magenta','LineWidth',2);

%legend('theory', 'shevalier', 'observation');

%xlim([0.1 100]);
%ylim([0.05 50]);

%legend('april','may','june','august','Location','northwest');
%legend('99 days','162 days ','357 days','Location','northwest');

grid ;