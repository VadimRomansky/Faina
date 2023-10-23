clear;

me = 9.1*10^-28;
mp = 1.6*10^-24;
c =3*10^10;
m=mp;

E1 = importdata('../examples_data/gamma1.5_theta0-90_protons/Ee3.dat');
F1 = importdata('../examples_data/gamma1.5_theta0-90_protons/Fs3.dat');
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
for i = 1:N1,
    FEE1(i)=F1(i)*E1(i)*E1(i);
end;

%MC_F = importdata('../examples_data/pdf_sf_gamma1.5/GLE_pdf_sf_B0_00003.dat');
MC_F = importdata('../examples_data/pdf_sf_gamma1.5/GLE_pdf_pf_306.dat');
N2 = size(MC_F,1);


E2(1:N2)=0;
F2(1:N2)=0;
P2(1:N2)=0;
dE2(1:N2)=0;
FEE2(1:N2)=0;
for i = 1:N2,
    P2(i)=(10^MC_F(i,1))*mp*c;
    E2(i)=sqrt(P2(i)*P2(i)*c*c + mp*mp*c*c*c*c);
    F2(i)=MC_F(i,2)*E2(i)/(P2(i)*P2(i)*P2(i)*c*c);
end;

dE2(1)=0;
for i=2:N2,
    dE2(i)=E2(i)-E2(i-1);
end;
for i=1:N2,
    FEE2(i)=F2(i)*E2(i)*E2(i);
end;

Nconcat1 = 185;
Econcat1 = E1(Nconcat1);
Nconcat2=1;
for i=1:N2,
    if(E2(i) > Econcat1)
        Nconcat2=i;
        break;
    end;
end;

N3 = Nconcat1 + (N2 - Nconcat2 + 1);
E3(1:N3)=0;
F3(1:N3)=0;
dE3(1:N3)=0;
FEE3(1:N3)=0;
for i=1:Nconcat1,
    E3(i)=E1(i);
    F3(i)=F1(i);
end;
tempF = exp(log(F1(Nconcat1)) + (log(F1(Nconcat1)/F1(Nconcat1-1))/log(E1(Nconcat1)/E1(Nconcat1-1)))*log(E2(Nconcat2)/E1(Nconcat1)));
%tempF = F1(Nconcat1) + ((F1(Nconcat1) - F1(Nconcat1-1))/(E1(Nconcat1)-E1(Nconcat1-1))*(E2(Nconcat2)-E1(Nconcat1)));
for i=Nconcat1+1:N3,
    E3(i)=E2(i - Nconcat1-1 + Nconcat2);
    F3(i)=F2(i - Nconcat1-1 + Nconcat2)*tempF/F2(Nconcat2);
    %F3(i)=tempF;
end;
dE3(1)=0;
for i=2:N3,
    dE3(i)=E3(i)-E3(i-1);
end;
norm = 0;
for i = 1:N3,
    norm = norm + F3(i)*dE3(i);
end;
for i = 1:N3,
    F3(i) = F3(i)/norm;
end;
for i = 1:N3,
    FEE3(i)=F3(i)*E3(i)*E3(i);
end;



set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_{E}');
xlabel ('E/m_p c^2');
ylabel ('F_{E} E^2');


plot(E1(1:N1)/(mp*c*c)-1, FEE1(1:N1),'red','LineWidth',2);
plot(E2(1:N2)/(mp*c*c)-1, FEE2(1:N2),'green','LineWidth',2);
plot(E3(1:N3)/(mp*c*c)-1, FEE3(1:N3),'blue','LineWidth',2);

E3kin(1:N3)=0;
F3kin(1:N3)=0;
for i=1:N3,
    E3kin(i) = E3(i)/(m*c*c) - 1.0;
    F3kin(i) = F3(i)*m*c*c;
end;


dlmwrite('Ee3.dat',E3kin,'delimiter',' ');
dlmwrite('Fs3.dat',F3kin,'delimiter',' ');
