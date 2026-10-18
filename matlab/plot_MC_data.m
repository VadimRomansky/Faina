clear;

me = 9.1*10^-28;
mp = 1.6*10^-24;
kB = 1.38*10^-16;
c =3*10^10;
m=mp;


frontDistribution = importdata('../examples_data/gamma0.7_MC1/GLE_pdf_sf.dat');
fullDistribution = importdata('../examples_data/gamma0.7_MC1/pdf_sf.dat');
Beff = importdata('../examples_data/gamma0.7_MC1/Beff.dat');
n = importdata('../examples_data/gamma0.7_MC1/number_dens.dat');

Np = size(frontDistribution,1);
Nx = size(n,1);


E2(1:Np)=0;
dE2(1:Np)=0;
F2(1:Np)=0;
P2(1:Np)=0;
dE2(1:Np)=0;
FEE2(1:Np)=0;
norm2 = 0;
dE2(1)=0;
for i = 1:Np,
    P2(i)=frontDistribution(i,1)*mp*c;
    F2(i) = frontDistribution(i,2);
        E2(i)=sqrt(P2(i)*P2(i)*c*c + mp*mp*c*c*c*c);
        F2(i)=F2(i)*E2(i)/(P2(i)*P2(i)*P2(i)*c*c);
        if(i > 1)
            dE2(i) = E2(i)-E2(i-1);
        end;
        norm2 = norm2 + (F2(i))*dE2(i);
end;


for i=1:Np,
        F2(i) = F2(i)/norm2;
        FEE2(i)=F2(i)*E2(i)*E2(i);
end;

concentration(1:Nx)=0;

for k=1:Nx,
for i = 1:Np,
    p=fullDistribution((k-1)*Np + i,1)*mp*c;
    f = fullDistribution((k-1)*Np + i,2);
        E=sqrt(P2(i)*P2(i)*c*c + mp*mp*c*c*c*c);
        f=f*E2(i)/(P2(i)*P2(i)*P2(i)*c*c);
        concentration(k) = concentration(k) + f*dE2(i);
end;
end;


figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{E}');
xlabel ('E');
ylabel ('F_{E}');

%plot(E1(1:N1)/(m*c*c)-1, F1(1:N1),'red','LineWidth',2);
%plot(E1(1:N1)/(m*c*c)-1, Fjuttner(1:N1),'cyan','LineWidth',2);
%loglog(E2(1:N2,1)/(m*c*c)-1, F2(1:N2,1),'magenta','LineWidth',2);
%loglog(E2(1:N2,2)/(m*c*c)-1, F2(1:N2,2),'green','LineWidth',2);
loglog(E2(1:Np)/(m*c*c)-1, F2(1:Np),'blue','LineWidth',2);
grid;


figure(2);
hold on;
set(gca, 'YScale', 'log');
%set(gca, 'XScale', 'log');
title ('n');
xlabel ('x');
ylabel ('n');

%plot(E1(1:N1)/(m*c*c)-1, F1(1:N1),'red','LineWidth',2);
%plot(E1(1:N1)/(m*c*c)-1, Fjuttner(1:N1),'cyan','LineWidth',2);
%loglog(E2(1:N2,1)/(m*c*c)-1, F2(1:N2,1),'magenta','LineWidth',2);
%loglog(E2(1:N2,2)/(m*c*c)-1, F2(1:N2,2),'green','LineWidth',2);
plot(n(1:Nx,2), concentration(1:Nx),'blue','LineWidth',2);
plot(n(1:Nx,2), n(1:Nx,3),'red','LineWidth',2);
grid;

ratio = concentration(20)/n(20,3);