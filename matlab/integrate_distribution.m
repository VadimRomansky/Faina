clear;

me = 9.1*10^-28;
mp = 1.6*10^-24;
c =3*10^10;
m=mp;

pev = 1.6*10^-12*10^15;

energy1 = 10*pev;
energy2 = 100*pev;

minEcr = 1.6*10^-24*c*c;
maxEcr = 10^10*minEcr;
Ncr = 1000;
Ecr(1:Ncr) = 0;
Fcr(1:Ncr) = 0;
factor = power(maxEcr/minEcr, 1/(Ncr - 1));
Ecr(1)=minEcr;
for i = 2:Ncr,
    Ecr(i) = Ecr(i-1)*factor;
end;
Fcr(1)=1;
for i = 2:Ncr,
    p=1.0;
    if(Ecr(i) < 3*1.6*1000) 
        p = 2.7;
    elseif (Ecr(i) < 3*1.6*1000000)
        p=3.0;
    else
        p=2.6;
    end;
    Fcr(i) = Fcr(i-1)*power(factor,-p);
end;

Ucr = 1.6*10^-12;
U = 0;
for i = 2:Ncr,
    U = U + Fcr(i)*Ecr(i)*(Ecr(i) - Ecr(i-1));
end;
norm = U/Ucr;
for i = 1:Ncr,
    Fcr(i) = Fcr(i)/norm;
end;

pevindex1 = 1;
for i = 1:Ncr,
    if(Ecr(i) > energy1)
        pevindex1 = i;
        break;
    end;
end;
pevindex2 = 1;
for i = 1:Ncr,
    if(Ecr(i) > energy2)
        pevindex2 = i;
        break;
    end;
end;

realPevenergy = 0;
for i = pevindex1:pevindex2,
    realPevenergy = realPevenergy + Fcr(i)*(Ecr(i)-m*c*c)*(Ecr(i)-Ecr(i-1));
end;


E1 = importdata('../examples_data/gamma1.5_combined_protons/Ee3.dat');
F1 = importdata('../examples_data/gamma1.5_combined_protons/Fs3.dat');
N1 = size(E1,2);
dE1(1:N1)=0;
FEE1(1:N1)=0;

E2 = importdata('./Ee3.dat');
F2 = importdata('./Fs3.dat');
N2 = size(E2,2);
dE2(1:N2)=0;
FEE2(1:N2)=0;


for i = 1:N1,
    E1(i) = m*c*c*(1.0 + E1(i));
end;
dE1(1)=0;
for i = 2:N1,
    dE1(i) = E1(i)-E1(i-1);
end;
for i = 1:N1,
    FEE1(i)=F1(i)*E1(i)*E1(i);
end;

norm = 0;
for i = 1:N1,
    norm = norm + F1(i)*dE1(i);
end;
for i = 1:N1,
    F1(i) = F1(i)/norm;
end;

for i = 1:N2,
    E2(i) = m*c*c*(1.0 + E2(i));
end;
dE2(1)=0;
for i = 2:N2,
    dE2(i) = E2(i)-E2(i-1);
end;
for i = 1:N2,
    FEE2(i)=F2(i)*E2(i)*E2(i);
end;

norm = 0;
for i = 1:N2,
    norm = norm + F2(i)*dE2(i);
end;
for i = 1:N2,
    F2(i) = F2(i)/norm;
end;

pevindex = 1;
for i = 1:N1,
    if(E1(i) > energy1)
        pevindex = i;
        break;
    end;
end;
pevupindex = 1;
for i = 1:N1,
    if(E1(i) > energy2)
        pevupindex = i;
        break;
    end;
end;

concentration = 80;
R = 1.9E17;
volume = 3.14*R*R*R*0.1;
pevenergy = 0;
for i = pevindex:pevupindex,
    pevenergy = pevenergy + F1(i)*(E1(i)-m*c*c)*dE1(i);
end;
pevenergy = 10*pevenergy*concentration*volume;

totalenergy = 0;
for i = 2:N1,
    totalenergy = totalenergy + F1(i)*(E1(i)-m*c*c)*dE1(i);
end;
totalenergy = totalenergy*concentration*volume;


q=4.84*10^-10;
lgalaxy = 3000*3*10^18;
Rgalaxy = 30000*3*10^18;
B=3*10^-6;

Dpev= 10^30;
time = lgalaxy^2/(4*Dpev);
Vgalaxy = 3.14*Rgalaxy*Rgalaxy*lgalaxy;
densityGalaxy = pevenergy/Vgalaxy;

timeYears = time/3E7;

timeExp = (time/(realPevenergy/densityGalaxy))/3E7;

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
plot(Ecr(1:Ncr)/(1.6*10^-12), Fcr(1:Ncr),'red','LineWidth',2);

figure(2);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
plot(E1(1:N1)/(1.6*10^-12), F1(1:N1),'red','LineWidth',2);
plot(E2(1:N1)/(1.6*10^-12), F2(1:N1),'blue','LineWidth',2);
