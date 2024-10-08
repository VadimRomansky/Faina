clear;

radiation = importdata('../outputSynch1.dat');
%compton = importdata('../outputCompt.dat');
%radiation2 = importdata('../outputSynch2.dat');
%radiation3 = importdata('../outputSynch3.dat');

radiationObserved = importdata('../examples_data/AT2020xnd_data/bright26.dat');

N = size(radiation,1);
Nr = size(radiation,2);

Nn = size(radiationObserved,1);

%CSS161010
cssx1(1:4) = [1.5, 3.0, 6.1, 9.87 ];
cssy1(1:4) = [ 1.5, 4.3, 6.1, 4.2 ];
cssError1(1:4) = [ 0.1, 0.2, 0.3, 0.2 ];

cssx2(1:5) = [ 1.5, 2.94, 6.1, 9.74, 22.0 ];
cssy2(1:5) = [ 4.7, 2.9, 2.3, 1.74, 0.56 ];
cssError2(1:5) = [ 0.6, 0.2, 0.1, 0.09, 0.03 ];

cssx3(1:6) = [ 0.33, 0.61, 1.5, 3.0, 6.05, 10.0 ];
%todo 0? or what?
cssy3(1:6) = [ 0.375, 0.79, 0.27, 0.17, 0.07, 0.032 ];
cssError3(1:6) = [ 0.375, 0.09, 0.07, 0.03, 0.01, 0.008 ];
%

%SN2009bb
jul12x(1:2) = 0;
jul12x(1) = 0.325;
jul12x(2) = 0.61;
jul12y(1:2) = 0;
jul12y(1) = 4.4;
jul12y(2) = 2.0;

jan12x(1:3) = 0;
jan12x(1) = 0.325;
jan12x(2) = 0.61;
jan12x(3) = 1.28;
jan12y(1:3) = 0;
jan12y(1) = 5.9;
jan12y(2) = 1.4;
jan12y(3) = 0.7;

apr11x(1:3) = 0;
apr11x(1) = 0.325;
apr11x(2) = 0.61;
apr11x(3) = 1.28;
apr11y(1:3) = 0;
apr11y(1) = 4.6;
apr11y(2) = 1.4;
ap11y(3) = 0.7;

decx(1:3) = 0;
decx(1) = 0.325;
decx(2) = 0.61;
decx(3) = 1.28;
decy(1:3) = 0;
decy(1) = 6.1;
decy(2) = 1.9;
decy(3) = 0.9;

octx(1:3) = 0;
octx(1) = 0.325;
octx(2) = 0.61;
octx(3) = 1.28;
octy(1:3) = 0;
octy(1) = 12;
octy(2) = 6.3;
octy(3) = 4.4;

augx(1:5) = 0;
augx(1) = 0.332;
augx(2) = 0.617;
augx(3) = 1.43;
augx(4) = 4.86;
augx(5) = 8.46;
augy(1:5) = 0;
augy(1) = 3.3;
augy(2) = 7.9;
augy(3) = 8.68;
augy(4) = 2.47;
augy(5) = 1.084;
augerr(1:5) = 0;
augerr(1) = 0.7;
augerr(2) = 0.3;
augerr(3) = 0.341;
augerr(4) = 0.155;
augerr(5) = 0.091;

augmax = 0.886;
augmaxy = 11.2;

junx(1:4) = 0;
junx(1) = 0.617;
junx(2) = 1.43;
junx(3) = 4.86;
junx(4) = 8.46;
juny(1:4) = 0;
juny(1) = 3.0;
juny(2) = 9.465;
juny(3) = 5.877;
juny(4) = 3.203;
junerr(1:4) = 0;
junerr(1) = 0.3;
junerr(2) = 0.186;
junerr(3) = 0.11;
junerr(4) = 0.074;

junmaxx = 1.65;
junmaxy = 13.2;

mayx(1:3) = 0;
mayx(1) = 1.43;
mayx(2) = 4.86;
mayx(3) = 8.46;
mayy(1:3) = 0;
mayy(1) = 4.93;
mayy(2) = 12.2;
mayy(3) = 6.82;
mayerr(1:3) = 0;
mayerr(1) = 0.122;
mayerr(2) = 0.165;
mayerr(3) = 0.098;

maymaxx = 2.96;
maymaxy = 15.2;

aprx(1:4) = 0;
aprx(1) = 1.43;
aprx(2) = 4.86;
aprx(3) = 8.46;
aprx(4) = 22.5;
apry(1:4) = 0;
apry(1) = 1.3;
apry(2) = 12.86;
apry(3) = 17.57;
apry(4) = 5.2;
aprerr(1:4) = 0;
aprerr(1) = 0.1;
aprerr(2) = 0.161;
aprerr(3) = 0.088;
aprerr(4) = 0.19;

aprmaxx = 6.50;
aprmaxy = 19.3;
%

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

radiationLinear(Nn) = 0;
power = 2.0;
for i = 1:Nn,
    radiationLinear(i) = radiationObserved(1,7)*(radiationObserved(i,5)/radiationObserved(1,5))^power;
end;


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
%loglog(radiation2(1:N,1),radiation2(1:N,2),'green','LineWidth',2);
%loglog(radiation3(1:N,1),radiation3(1:N,2),'blue','LineWidth',2);
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

%errorbar(cssx1,cssy1,cssError1,'red','LineWidth',2);
%errorbar(cssx2,cssy2,cssError2,'green','LineWidth',2);
%errorbar(cssx3,cssy3,cssError3,'blue','LineWidth',2);

errorbar(radiationObserved(1:Nn,5), radiationObserved(1:Nn,7), radiationObserved(1:Nn,8),'blue','LineWidth',2);
loglog(radiationObserved(1:Nn,5), radiationLinear(1:Nn),'green');
%loglog(radiation4(1:Nn,1),radiation4(1:Nn,2),'blue','LineWidth',2);
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

%nufnuSynch(1:N)=0;
%nufnuCompt(1:N)=0;
%for i=1:N,
%    nufnuSynch(i) = radiation(i,1)*radiation(i,2)*10^-26*10^9;
%    nufnuCompt(i) = compton(i,1)*compton(i,2)*10^-26*10^9;
%end;

%figure(2);
%hold on;
%set(gca, 'YScale', 'log');
%set(gca, 'XScale', 'log');
%title ('\nu F_{\nu}');
%xlabel ('{\E} eV');
%ylabel ('\nu F_{\nu} erg cm^{-2} s^{-1}');

%loglog(radiation(1:N,1)*10^9*h/(1.6*10^-12),nufnuSynch(1:N),'red','LineWidth',2);
%loglog(compton(1:N,1)*10^9*h/(1.6*10^-12),nufnuCompt(1:N),'blue','LineWidth',2);

%dlmwrite('radiation.dat',radiation,'delimiter',' ');