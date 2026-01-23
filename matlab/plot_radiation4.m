clear;

radiation0 = importdata('../outputSynch0.dat');
radiation1 = importdata('../outputSynch1.dat');
radiation2 = importdata('../outputSynch2.dat');
radiation3 = importdata('../outputSynch3.dat');

css0 = importdata('../examples_data/css_data/coppejans69.txt');
css1 = importdata('../examples_data/css_data/coppejans99.txt');
css2 = importdata('../examples_data/css_data/coppejans162.txt');
css3 = importdata('../examples_data/css_data/coppejans357.txt');

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{\nu}');
xlabel ('{\nu} GHz');
ylabel ('mJy');

loglog(radiation0(:,1),radiation0(:,2),'red','LineWidth',2);
loglog(radiation1(:,1),radiation1(:,2),'green','LineWidth',2);
loglog(radiation2(:,1),radiation2(:,2),'blue','LineWidth',2);
loglog(radiation3(:,1),radiation3(:,2),'magenta','LineWidth',2);
%loglog(compton(1:N,1),compton(1:N,2),'blue','LineWidth',2);

%loglog(radiation(1:N,1),radiation(1:N,5),'magenta','LineWidth',2);
%plot(radiation(1:N,1),radiation(1:N,6),'Color',[1.0,0.6,0],'LineWidth',2);
%plot(radiation(1:N,1),radiation(1:N,7),'black','LineWidth',2);

%loglog(aprx(1:4),apry(1:4),'--o','Color','red','LineWidth',2);
%loglog(mayx(1:3),mayy(1:3),'--o','Color','green','LineWidth',2);
%loglog(junx(1:4),juny(1:4),'--o','Color','blue','LineWidth',2);
%loglog(augx(1:5),augy(1:5),'--o','Color','magenta','LineWidth',2);
%plot(octx(1:3),octy(1:3),'--o','Color',[1.0,0.6,0],'LineWidth',2);
%plot(decx(1:3),decy(1:3),'--o','Color','black','LineWidth',2);

errorbar(css0(:,5),css0(:,7),css0(:,8),'red','LineWidth',2,'LineStyle','none');
errorbar(css1(:,5),css1(:,7),css1(:,8),'green','LineWidth',2,'LineStyle','none');
errorbar(css2(:,5),css2(:,7),css2(:,8),'blue','LineWidth',2,'LineStyle','none');
errorbar(css3(:,5),css3(:,7),css3(:,8),'magenta','LineWidth',2,'LineStyle','none');
%loglog(radiation4(1:Nn,1),radiation4(1:Nn,2),'blue','LineWidth',2);
%errorbar(cssx2,cssy2,cssError2,'green','LineWidth',2);
%errorbar(cssx3,cssy3,cssError3,'blue','LineWidth',2);

%errorbar(aprx,apry,aprerr,'red','LineWidth',2);
%errorbar(mayx,mayy,mayerr,'green','LineWidth',2);
%errorbar(junx,juny,junerr,'blue','LineWidth',2);
%errorbar(augx,augy,augerr,'magenta','LineWidth',2);

%legend('theory', 'shevalier', 'observation');

xlim([0.1 100]);
%ylim([0.05 50]);

%legend('april','may','june','august','Location','northwest');
%legend('99 days','162 days ','357 days','Location','northwest');

grid ;

