clear;

data = importdata('../distribution.dat');

Ee=importdata('../examples_data/gamma0.2_theta0-90/Ee3.dat');
Fs=importdata('../examples_data/gamma0.2_theta0-90/Fs3.dat');

Ee2=importdata('../examples_data/gamma0.3_theta0-90/Ee3.dat');
Fs2=importdata('../examples_data/gamma0.3_theta0-90/Fs3.dat');

Ee3=importdata('../examples_data/gamma0.5_theta0-90/Ee3.dat');
Fs3=importdata('../examples_data/gamma0.5_theta0-90/Fs3.dat');

Ee4=importdata('../examples_data/gamma1.5_theta0-90/Ee3.dat');
Fs4=importdata('../examples_data/gamma1.5_theta0-90/Fs3.dat');

N = size(data,1);

N2 = size(Ee,2);



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



loglog(data(1:N,1)/(m*c*c),data(1:N,2),'red','LineWidth',2,'Marker','+');

figure(2);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{E}');
xlabel ('E erg');
ylabel ('F_{E}');

loglog((Ee(1:N2)+1),Fs(1:N2),'red');
loglog((Ee2(1:N2)+1),Fs2(1:N2),'green');
loglog((Ee3(1:N2)+1),Fs3(1:N2),'blue');
loglog((Ee4(1:N2)+1),Fs4(1:N2),'magenta');
