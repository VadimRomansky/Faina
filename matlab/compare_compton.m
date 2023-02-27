clear;


radiationEu = importdata('../outputEu.dat');
radiationEd = importdata('../outputEd.dat');
radiationEt = importdata('../outputEt.dat');
%radiationEa = importdata('../outputEa.dat');

radiationUvarov = importdata('../compt_test1_KN.dat');
radiationMine = importdata('../output.dat');


Nnu = size(radiationEu,1);
N1 = size(radiationUvarov,1);
N2 = size(radiationMine,1);



startPower = 90;
endPower = 110;

Fa(1:Nnu) = 0;

Fa(startPower) = radiationEd(startPower,2);
Fa(endPower) = radiationEd(endPower,2);

%gammap = log(Fpa(startPower)/Fpa(endPower))/log((me*energy(startPower)+m)/(me*energy(endPower)+m));

polyfitx(1:endPower-startPower + 1) = 0;
polyfity(1:endPower-startPower + 1) = 0;

for i = 1:endPower-startPower + 1,
    polyfitx(i) = log(radiationEd(i+startPower - 1,1));
    polyfity(i) = log(radiationEd(i+startPower - 1,2));
end;
p = polyfit(polyfitx, polyfity, 1);

%ap = exp(log(Fpa(startPower)) - gammap*log((me*energy(startPower)+m)));

for i = startPower-2:endPower+2,
    %Fpa(i) = ap*((me*energy(i)+m)^gammap);
    Fa(i) = exp(polyval(p, log(radiationEd(i))));
end;

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('EF(E)');
xlabel ('E eV');
ylabel ('E F_{E} erg/{s cm^{2}}');


relation = radiationEd(startPower,2)/radiationEu(startPower,2);

plot(radiationEd(1:Nnu,1),radiationEd(1:Nnu,2),'red','LineWidth',2);
plot(radiationEu(1:Nnu,1),radiationEu(1:Nnu,2),'green','LineWidth',2);
plot(radiationEt(1:Nnu,1),radiationEt(1:Nnu,2),'blue','LineWidth',2);
%plot(radiationEa(1:Nnu,1),radiationEa(1:Nnu,2)/(4*3.14),'magenta','LineWidth',2);
plot(radiationEd(1:Nnu,1),Fa(1:Nnu),'black','LineWidth',2);
legend('klein nisina','kang jones','thompson','anisotropic K-N','Power Law');
grid ;

figure(2)
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
plot(radiationUvarov(1:N1,1), radiationUvarov(1:N1,2),'red','LineWidth',2);
plot(radiationMine(1:N2,1), radiationMine(1:N2,2),'blue','LineWidth',2);
