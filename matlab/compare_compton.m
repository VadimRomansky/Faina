clear;


%radiationEu = importdata('../outputEu.dat');
%radiationEd = importdata('../outputEd.dat');
%radiationEt = importdata('../outputEt.dat');
%radiationEa = importdata('../outputEa.dat');

radiationUvarov = importdata('../compt_test1_KN.dat');
radiationMineJones = importdata('../output1.dat');
radiationMine1 = importdata('../output2.dat');
radiationMine2 = importdata('../output3.dat');
radiationMine3 = importdata('../output4.dat');
radiationMine4 = importdata('../output5.dat');


%Nnu = size(radiationEu,1);
N1 = size(radiationUvarov,1);
N2 = size(radiationMine1,1);



% startPower = 90;
% endPower = 110;
% 
% Fa(1:Nnu) = 0;
% 
% Fa(startPower) = radiationEd(startPower,2);
% Fa(endPower) = radiationEd(endPower,2);
% 
% gammap = log(Fpa(startPower)/Fpa(endPower))/log((me*energy(startPower)+m)/(me*energy(endPower)+m));
% 
% polyfitx(1:endPower-startPower + 1) = 0;
% polyfity(1:endPower-startPower + 1) = 0;
% 
% for i = 1:endPower-startPower + 1,
%     polyfitx(i) = log(radiationEd(i+startPower - 1,1));
%     polyfity(i) = log(radiationEd(i+startPower - 1,2));
% end;
% p = polyfit(polyfitx, polyfity, 1);
% 
% ap = exp(log(Fpa(startPower)) - gammap*log((me*energy(startPower)+m)));
% 
% for i = startPower-2:endPower+2,
%     Fpa(i) = ap*((me*energy(i)+m)^gammap);
%     Fa(i) = exp(polyval(p, log(radiationEd(i))));
% end;

% figure(1);
% hold on;
% set(gca, 'YScale', 'log');
% set(gca, 'XScale', 'log');
% title ('EF(E)');
% xlabel ('E eV');
% ylabel ('E F_{E} erg/{s cm^{2}}');
% 
% 
% relation = radiationEd(startPower,2)/radiationEu(startPower,2);
% 
% plot(radiationEd(1:Nnu,1),radiationEd(1:Nnu,2),'red','LineWidth',2);
% plot(radiationEu(1:Nnu,1),radiationEu(1:Nnu,2),'green','LineWidth',2);
% plot(radiationEt(1:Nnu,1),radiationEt(1:Nnu,2),'blue','LineWidth',2);
% plot(radiationEa(1:Nnu,1),radiationEa(1:Nnu,2)/(4*3.14),'magenta','LineWidth',2);
% plot(radiationEd(1:Nnu,1),Fa(1:Nnu),'black','LineWidth',2);
% legend('klein nisina','kang jones','thompson','anisotropic K-N','Power Law');
% grid ;

startPower = 25;
endPower = 35;

polyfitxJones(1:endPower-startPower + 1) = 0;
polyfityJones(1:endPower-startPower + 1) = 0;
for i = 1:endPower-startPower + 1,
    polyfitxJones(i) = log(radiationMineJones(i+startPower - 1,1));
    polyfityJones(i) = log(radiationMineJones(i+startPower - 1,2));
end;
pJones = polyfit(polyfitxJones, polyfityJones, 1);
appJones(1:N2) = 0;
for i = startPower-5:endPower+5,
    appJones(i) = exp(polyval(pJones, log(radiationMineJones(i,1))));
end;

polyfitx1(1:endPower-startPower + 1) = 0;
polyfity1(1:endPower-startPower + 1) = 0;
for i = 1:endPower-startPower + 1,
    polyfitx1(i) = log(radiationMine1(i+startPower - 1,1));
    polyfity1(i) = log(radiationMine1(i+startPower - 1,2));
end;
p1 = polyfit(polyfitx1, polyfity1, 1);
app1(1:N2) = 0;
for i = startPower-5:endPower+5,
    app1(i) = exp(polyval(p1, log(radiationMine1(i,1))));
end;

polyfitx3(1:endPower-startPower + 1) = 0;
polyfity3(1:endPower-startPower + 1) = 0;
for i = 1:endPower-startPower + 1,
    polyfitx3(i) = log(radiationMine3(i+startPower - 1,1));
    polyfity3(i) = log(radiationMine3(i+startPower - 1,2));
end;
p3 = polyfit(polyfitx3, polyfity3, 1);
app3(1:N2) = 0;
for i = startPower-5:endPower+5,
    app3(i) = exp(polyval(p3, log(radiationMine3(i,1))));
end;

figure(2)
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
for i = 1:N2,
    %radiationMine1(i,2) = radiationMine1(i,2)*4*3.14*radiationMine1(i,1);
end;
plot(radiationUvarov(1:N1,1), radiationUvarov(1:N1,2)/(4*3.14),'cyan','LineWidth',2);
plot(radiationMineJones(1:N2,1), radiationMineJones(1:N2,2),'blue','LineWidth',2);
plot(radiationMine1(1:N2,1), radiationMine1(1:N2,2),'green','LineWidth',2);
plot(radiationMine2(1:N2,1), radiationMine2(1:N2,2),'red','LineWidth',2);
plot(radiationMine3(1:N2,1), radiationMine3(1:N2,2),'magenta','LineWidth',2);
plot(radiationMine4(1:N2,1), radiationMine4(1:N2,2),'yellow','LineWidth',2);

plot(radiationMineJones(1:N2,1), appJones(1:N2),'--','Color','blue','LineWidth',2);
plot(radiationMine1(1:N2,1), app1(1:N2),'--','Color','green','LineWidth',2);
plot(radiationMine3(1:N2,1), app3(1:N2),'--','Color','magenta','LineWidth',2);
legend('uvarov','kang jones','anisotropic KN','isotropic KN','anisotropic photons KN \theta = \pi/18','anisotropic photons KN \theta = 0.9\pi', strcat('fit jones\gamma = ',num2str(pJones(1))), strcat('fit KN\gamma = ',num2str(p1(1))), strcat('fit KN \gamma = ',num2str(p3(1))));
