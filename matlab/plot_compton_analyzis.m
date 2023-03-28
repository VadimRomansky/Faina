clear;

N = 1000;

E0min = 0.01;
E0max = 10;
E0(1:N) = 0.1;
factor = (E0max/E0min)^(1.0/(N-1));
E0(1) = E0min;
for i = 2:N,
    E0(i) = E0(i-1)*factor;
end;

E1min = 1E7;
E1max = 1E10;
E1(1:N) = 1E10;
factor = (E1max/E1min)^(1.0/(N-1));
E1(1) = E1min;
for i = 2:N,
    E1(i) = E1(i-1)*factor;
end;

gmin = 1E8;
gmax = 2E10;
gam(1:N) = 2E10;
factor = (gmax/gmin)^(1.0/(N-1));
gam(1) = gmin;
for i = 2:N,
    gam(i) = gam(i-1)*factor;
end;

index0 = 10;
index1 = 10;
indexg = 10;

s1(1:N)=0;
s2(1:N)=0;
s3(1:N)=0;
q1(1:N)=0;
q2(1:N)=0;
q3(1:N)=0;



set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');

for i = 1:N,
    index0 = i;
    e1 = 1E9;
    g0 = 1E10;
    G = 4*E0(index0)*g0;
    q = e1/((g0 - e1)*G);
    q1(i) = q;
    s1(i) = 2*q*log(q) + 1 + q - 2*q*q +0.5*q*q*(1-q)*G*G/(1+q*G);
end;
plot(E0(1:N),s1(1:N),'red','LineWidth',2);
plot(E0(1:N),q1(1:N),'blue','LineWidth',2);

figure(2);
hold on;
%set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');

index0 = 10;
for i = 1:N,
    index1 = i;
    e0 = 1.0;
    g0 = 1E10;
    G = 4*e0*g0;
    q = E1(index1)/((g0 - E1(index1))*G);
    q2(i) = q;
    s2(i) = 2*q*log(q) + 1 + q - 2*q*q +0.5*q*q*(1-q)*G*G/(1+q*G);
end;
plot(E1(1:N),s2(1:N),'red','LineWidth',2);
plot(E1(1:N),q2(1:N),'blue','LineWidth',2);

figure(3);
hold on;
%set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
index1 = 10;
for i = 1:N,
    indexg = i;
    e0 = 1.0;
    e1 = 1E9;
    G = 4*e0*gam(indexg);
    q = e1/((gam(indexg) - e1)*G);
    q3(i) = q;
    s3(i) = 2*q*log(q) + 1 + q - 2*q*q +0.5*q*q*(1-q)*G*G/(1+q*G);
end;
plot(gam(1:N),s3(1:N),'red','LineWidth',2);
plot(gam(1:N),q3(1:N),'blue','LineWidth',2);
