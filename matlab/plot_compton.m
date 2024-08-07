clear;


radiation1 = importdata('../outputCompt.dat');
%radiation2 = importdata('../output1.dat');
%radiation3 = importdata('../output2.dat');


%radiation = importdata('../outputNu.dat');


Nnu = size(radiation1,1);


figure(2);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('E F_{E}');
xlabel ('E eV');
ylabel ('E F_{E} erg cm^{-2} s^{-1}');

xlim([0.0001,10^17])
ylim([10^-40,10^-32])

plot(radiation1(1:Nnu,1),radiation1(1:Nnu,2),'red','LineWidth',2);
plot(radiation1(1:Nnu,1),radiation1(1:Nnu,3),'green','LineWidth',2);
plot(radiation1(1:Nnu,1),radiation1(1:Nnu,4),'blue','LineWidth',2);
plot(radiation1(1:Nnu,1),radiation1(1:Nnu,5),'magenta','LineWidth',2);
plot(radiation1(1:Nnu,1),radiation1(1:Nnu,6),'cyan','LineWidth',2);
%plot(radiation2(1:Nnu,1),radiation2(1:Nnu,2),'blue','LineWidth',2);
%plot(radiation3(1:Nnu,1),radiation3(1:Nnu,2),'green','LineWidth',2);
%plot(radiation1(1:Nnu,1),Fa(1:Nnu),'blue','LineWidth',2);
legend('Jones','Isotropic KN','Anisotropic KN1','Anisotropic KN2','Anisotropic KN3');
%grid ;

