clear;

radiationUvarov = importdata('../compt_test1_KN.dat');
radiationMineJones = importdata('../output1.dat');
radiationMine1 = importdata('../output2.dat');
radiationMine2 = importdata('../output3.dat');
radiationMine3 = importdata('../output4.dat');
radiationMine4 = importdata('../output5.dat');

N2 = size(radiationMine1,1);

for i = 1:N2,
    radiationMineJones(i,2) = radiationMineJones(i,2)*radiationMineJones(i,1)*1.6*10^-12;
    radiationMine1(i,2) = radiationMine1(i,2)*radiationMine1(i,1)*1.6*10^-12;
    radiationMine2(i,2) = radiationMine2(i,2)*radiationMine2(i,1)*1.6*10^-12;
    radiationMine3(i,2) = radiationMine3(i,2)*radiationMine3(i,1)*1.6*10^-12;
    radiationMine4(i,2) = radiationMine4(i,2)*radiationMine4(i,1)*1.6*10^-12;
end;

figure(2)
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
for i = 1:N2,
    %radiationMine1(i,2) = radiationMine1(i,2)*4*3.14*radiationMine1(i,1);
end;
plot(radiationMineJones(1:N2,1), radiationMineJones(1:N2,2),'blue','LineWidth',2);
plot(radiationMine1(1:N2,1), radiationMine1(1:N2,2),'green','LineWidth',2);
plot(radiationMine2(1:N2,1), radiationMine2(1:N2,2),'red','LineWidth',2);
%plot(radiationMine3(1:N2,1), radiationMine3(1:N2,2),'magenta','LineWidth',2);
%plot(radiationMine4(1:N2,1), radiationMine4(1:N2,2),'yellow','LineWidth',2);

%legend('uvarov','kang jones','anisotropic KN','isotropic KN','anisotropic photons KN \theta = \pi/100','anisotropic photons KN \theta = 0.9\pi', strcat('fit jones\gamma = ',num2str(pJones(1))), strcat('fit KN\gamma = ',num2str(p1(1))), strcat('fit KN \gamma = ',num2str(p3(1))));
