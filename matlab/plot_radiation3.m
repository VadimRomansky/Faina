clear;

radiation = importdata('../outputSynch3.dat');

N = size(radiation,1);

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{\nu}');xlabel ('{\nu} GHz');
ylabel ('mJy');

loglog(radiation(1:N,1),radiation(1:N,2),'red','LineWidth',2);
loglog(radiation(1:N,1),radiation(1:N,3),'green','LineWidth',2);
loglog(radiation(1:N,1),radiation(1:N,4),'blue','LineWidth',2);

ratioleft = radiation(1,3)/radiation(1,2);
ratioright = radiation(N,4)/radiation(N,2);


grid ;
