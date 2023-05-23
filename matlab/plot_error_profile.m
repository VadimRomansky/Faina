clear;

Np1 = 1;
Np2 = 2;

param1 = importdata(strcat('../parameter',num2str(Np1),'.dat'));
param2 = importdata(strcat('../parameter',num2str(Np2),'.dat'));
error = importdata(strcat('../error_optimized_',num2str(Np1),'_',num2str(Np2),'.dat'));

N = size(param1,1);

sqrterror(1:N,1:N)=0;
for i = 1:N,
    for j = 1:N,
        sqrterror(i,j) = sqrt(error(i,j));
    end;
end;


figure(1);
colormap Jet;
[X, Y] = meshgrid(param2, param1);
surf(X, Y, sqrterror);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')
shading interp;
title ('error');
xlabel ('n');
ylabel ('\sigma');
zlabel ('error');
grid ;