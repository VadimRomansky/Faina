clear;

me = 9.1*10^-28;
mp = 1.6*10^-24;
kB = 1.38*10^-16;
c =3*10^10;
m=mp;

E1 = importdata('../examples_data/gamma0.66_PICMC/PICnew/Ee3.dat');
F1 = importdata('../examples_data/gamma0.66_PICMC/PICnew/Fs3.dat');
N1 = size(E1,2);
dE1(1:N1)=0;
FEE1(1:N1)=0;

for i = 1:N1,
    %E1(i) = m*c*c*(1.0 + 1.2*E1(i));
    E1(i) = m*c*c*(1.0 + E1(i));
end;
dE1(1)=0;
for i = 2:N1,
    dE1(i) = E1(i)-E1(i-1);
end;
norm1 = 0;
for i = 1:N1,
    norm1 = norm1 + F1(i)*dE1(i);
end;

for i=1:N1,
    F1(i) = F1(i)/norm1;
end;

for i = 1:N1,
    FEE1(i)=F1(i)*E1(i)*E1(i);
end;

T1 = 1.2E12;
Fjuttner(1:N1) = 0;
theta = kB*T1/(m*c*c);
bes = besselk(2, 1/theta);
for i = 1:N1,
        gam = E1(i)/(m*c*c);
        beta = sqrt(1.0 - 1.0/(gam*gam));
        exp1 = exp(-gam/theta);       
        Fjuttner(i) = ((1.0/(theta*bes))*exp1*gam*gam*beta)/(m*c*c);
end;

normt = 0;
for i = 1:N1,
    normt = normt + Fjuttner(i)*dE1(i);
    %FEE1(i)=F1(i)*E1(i)*E1(i);
end;

%MC_F = importdata('../examples_data/Grafik_u0_03_B0_003/GLE_pdf_sf8.dat');
%MC_F2 = importdata('../examples_data/Grafik_u0_03_B0_000003/GLE_pdf_sf8.dat');
%MC_F3 = importdata('../examples_data/Grafik_u0_03_B0_000003/GLE_pdf_sf8.dat');
MC_F = importdata('../examples_data/gamma0.66_PICMC/MC0.75/GLE_pdf_pf_306.dat');
MC_F2 = importdata('../examples_data/gamma0.66_PICMC/MC0.75/GLE_pdf_pf_306.dat');
MC_F3 = importdata('../examples_data/gamma0.66_PICMC/MC0.75/GLE_pdf_pf_306.dat');
N2 = size(MC_F,1);


E2(1:N2,1:3)=0;
dE2(1:N2,1:3)=0;
F2(1:N2,1:3)=0;
P2(1:N2,1:3)=0;
dE2(1:N2,1:3)=0;
FEE2(1:N2,1:3)=0;
norm2(1:3) = 0;
dE2(1,1)=0;
dE2(1,2)=0;
dE2(1,3)=0;
for i = 1:N2,
    P2(i,1)=(10^MC_F(i,1))*mp*c;
    P2(i,2)=(10^MC_F2(i,1))*mp*c;
    P2(i,3)=(10^MC_F3(i,1))*mp*c;
    %P2(i,1)=MC_F(i,1)*mp*c;
    %P2(i,2)=MC_F2(i,1)*mp*c;
    %P2(i,3)=MC_F3(i,1)*mp*c;
    F2(i,1) = MC_F(i,2);
    F2(i,2) = MC_F2(i,2);
    F2(i,3) = MC_F3(i,2);
    for j = 1:3,
        E2(i,j)=sqrt(P2(i,j)*P2(i,j)*c*c + mp*mp*c*c*c*c);
        F2(i,j)=F2(i,j)*E2(i,j)/(P2(i,j)*P2(i,j)*P2(i,j)*c*c);
        if(i > 1)
            dE2(i,j) = E2(i,j)-E2(i-1,j);
        end;
        norm2(j) = norm2(j) + (F2(i,j))*dE2(i,j);
    end;
end;


for i=1:N2,
    for j = 1:3,
        F2(i,j) = F2(i,j)/norm2(j);
        FEE2(i,j)=F2(i,j)*E2(i,j)*E2(i,j);
    end;
end;

Tmin = 2E11;
Tmax = 1E13;

Tleft = Tmin;
Tright = Tmax;

index1 = 50;
index2 = 100;

Fjuttner2(1:N2)=0;
for j = 1:20,
    T1 = Tleft + (Tright - Tleft)/3;
    T2 = Tleft + (Tright - Tleft)*2/3;
    s1 = 0;
    s2 = 0;
    theta = kB*T1/(m*c*c);
    bes = besselk(2, 1/theta);
    for i = index1:index2,
        gam = E1(i)/(m*c*c);
        beta = sqrt(1.0 - 1.0/(gam*gam));
        exp1 = exp(-gam/theta);       
        Fjuttner(i) = ((1.0/(theta*bes))*exp1*gam*gam*beta)/(m*c*c);
        s1 = s1 + ((Fjuttner(i) - F1(i))^2)*dE1(i);
    end;
    theta = kB*T2/(m*c*c);
    bes = besselk(2, 1/theta);
    for i = index1:index2,
        gam = E1(i)/(m*c*c);
        beta = sqrt(1.0 - 1.0/(gam*gam));
        exp1 = exp(-gam/theta);       
        Fjuttner(i) = ((1.0/(theta*bes))*exp1*gam*gam*beta)/(m*c*c);
        s2 = s2 + ((Fjuttner(i) - F1(i))^2)*dE1(i);
    end;
    if(s1 < s2)
        Tright = T2;
    else 
        Tleft = T1;
    end;
end;
Tpic = (Tleft + Tright)/2;
%T = 10^12;
theta = kB*Tpic/(m*c*c);
bes = besselk(2, 1/theta);
Fshifted(1:N1) = 0;
for i = 1:N1,   
    gam = E1(i)/(m*c*c);
    beta = sqrt(1.0 - 1.0/(gam*gam));
    exp1 = exp(-gam/theta);       
    Fjuttner(i) = ((1.0/(theta*bes))*exp1*gam*gam*beta)/(m*c*c);
    Fshifted(i) = juttner_shifted_integrated(gam, 0.127, sqrt(1.0 - 1.0/2.25));
end;


Tmin = 2E11;
Tmax = 1E13;

Tleft = Tmin;
Tright = Tmax;

index1 = 35;
index2 = 46;

for j = 1:20,
    T1 = Tleft + (Tright - Tleft)/3;
    T2 = Tleft + (Tright - Tleft)*2/3;
    s1 = 0;
    s2 = 0;
    theta = kB*T1/(m*c*c);
    bes = besselk(2, 1/theta);
    for i = index1:index2,
        gam = E2(i,3)/(m*c*c);
        beta = sqrt(1.0 - 1.0/(gam*gam));
        exp1 = exp(-gam/theta);       
        Fjuttner2(i) = ((1.0/(theta*bes))*exp1*gam*gam*beta)/(m*c*c);
        s1 = s1 + ((Fjuttner2(i) - F2(i,3))^2)*dE2(i);
    end;
    theta = kB*T2/(m*c*c);
    bes = besselk(2, 1/theta);
    for i = index1:index2,
        gam = E2(i,3)/(m*c*c);
        beta = sqrt(1.0 - 1.0/(gam*gam));
        exp1 = exp(-gam/theta);       
        Fjuttner2(i) = ((1.0/(theta*bes))*exp1*gam*gam*beta)/(m*c*c);
        s2 = s2 + ((Fjuttner2(i) - F2(i,3))^2)*dE2(i);
    end;
    if(s1 < s2)
        Tright = T2;
    else 
        Tleft = T1;
    end;
end;
Tmc = (Tleft + Tright)/2;
%Tmc = 10*10^11;
theta = kB*Tmc/(m*c*c);
bes = besselk(2, 1/theta);
for i = 1:N2,   
    gam = E2(i,3)/(m*c*c);
    beta = sqrt(1.0 - 1.0/(gam*gam));
    exp1 = exp(-gam/theta);       
    Fjuttner2(i) = ((1.0/(theta*bes))*exp1*gam*gam*beta)/(m*c*c);
    Fshifted(i) = juttner_shifted_integrated(gam, 0.127, sqrt(1.0 - 1.0/2.25));
end;

endPower1 = 120;
startPower1 = 100;
polyfitx1(1:endPower1-startPower1 + 1) = 0;
polyfity1(1:endPower1-startPower1 + 1) = 0;
Fpa1(1:N1) = 0;

for i = 1:endPower1-startPower1 + 1,
    polyfitx1(i) = log(E1(i+startPower1 - 1)-m*c*c);
    polyfity1(i) = log(F1(i+startPower1 - 1));
end
p1 = polyfit(polyfitx1, polyfity1, 1);

for i = startPower1-20:endPower1+20,
    Fpa1(i) = exp(polyval(p1, log(E1(i)-m*c*c)));
end;

endPower2 = 160;
startPower2 = 135;
polyfitx2(1:endPower2-startPower2 + 1) = 0;
polyfity2(1:endPower2-startPower2 + 1) = 0;
Fpa2(1:N1) = 0;

for i = 1:endPower2-startPower2 + 1,
    polyfitx2(i) = log(E1(i+startPower2 - 1)-m*c*c);
    polyfity2(i) = log(F1(i+startPower2 - 1));
end
p2 = polyfit(polyfitx2, polyfity2, 1);

for i = startPower2-20:endPower2+20,
    Fpa2(i) = exp(polyval(p2, log(E1(i)-m*c*c)));
end;

endPower3 = 60;
startPower3 = 40;
polyfitx3(1:endPower3-startPower3 + 1) = 0;
polyfity3(1:endPower3-startPower3 + 1) = 0;
Fpa3(1:N2) = 0;

for i = 1:endPower3-startPower3 + 1,
    polyfitx3(i) = log(E2(i+startPower3 - 1)-m*c*c);
    polyfity3(i) = log(F2(i+startPower3 - 1));
end
p3 = polyfit(polyfitx3, polyfity3, 1);

for i = startPower3-10:endPower3+10,
    Fpa3(i) = exp(polyval(p3, log(E2(i)-m*c*c)));
end;

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{E}');
xlabel ('E');
ylabel ('F_{E}');

plot(E1(1:N1)/(m*c*c)-1, F1(1:N1),'red','LineWidth',2);
plot(E1(1:N1)/(m*c*c)-1, Fjuttner(1:N1),'--','Color','magenta','LineWidth',2);
%plot(E1(1:N1)/(m*c*c)-1, Fpa1(1:N1),':','Color', 'green', 'LineWidth',4);
%plot(E1(1:N1)/(m*c*c)-1, Fpa2(1:N1),'-.','Color', '#FFA500', 'LineWidth',4);
plot(E2(1:N2,3)/(m*c*c)-1, F2(1:N2,3),'blue','LineWidth',2);
plot(E2(1:N2)/(m*c*c)-1, Fjuttner2(1:N2),'--','Color','cyan','LineWidth',2);
%plot(E2(1:N2)/(m*c*c)-1, Fpa3(1:N2),':','Color','#000080','LineWidth',4);
%legend('PIC','T=6.6E11', 'p=-2.09','p=-2.37', 'Monte-Carlo','T=1.1E12','p=-2.12','Location','southwest');
grid;