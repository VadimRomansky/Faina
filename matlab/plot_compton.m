clear;


radiationE = importdata('../output.dat');


%radiation = importdata('../outputNu.dat');


Nnu = size(radiationE,1);


figure(2);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('E F_{E}');
xlabel ('E eV');
ylabel ('E F_{E} erg cm^{-2} s^{-1}');

plot(radiationE(1:Nnu,1),radiationE(1:Nnu,2),'red','LineWidth',2);
%plot(radiation1(1:Nnu,1),Fa(1:Nnu),'blue','LineWidth',2);
%legend('Dubus','Uvarov/mc^2','\nu^{-1.5}');
%grid ;

