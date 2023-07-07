clear;
data = importdata('../image.dat');
profile = importdata('../outputRadial.dat');
R = 259.5;
R0 = R*(1 - 2*0.05);
R0=0;
Nr = size(data,1);
Nphi = size(data,2);
Nprofile = size(profile,1);

Ndata = 20;
observedR(1:Ndata) = [250, 250.5, 251, 251.5, 252, 252.5, 253, 253.5, 254, 254.5, 255, 255.5, 256, 256.5, 257, 257.5, 258, 258.5, 259, 259.5];
observedFlux(1:Ndata) = [0.35, 0.32, 0.35, 0.37, 0.4, 0.45, 0.5, 0.6, 0.75, 0.8, 0.9, 1.0, 0.95, 0.8, 0.55, 0.45, 0.25, 0.2, 0.15, 0.15 ];
observedError(1:Ndata) = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1];

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

image(1:Nphi+1,1:Nr) = 0;
for i = 1:Nr,
    for j = 1:Nphi,
        image(j,i) = data(i,j);
    end;
    image(Nphi+1,i) = data(i,j);
end;

radius(1:Nr)=0;
for i = 1:Nr,
    radius(i) = R0 + (R-R0)*(i-0.5)/Nr;
end;

figure(1)
[r,t] = meshgrid(radius,pi/Nphi:2*pi/Nphi:2*pi+pi/Nphi);
x = r.*cos(t);
y = r.*sin(t);
%contourf(x,y,image);
surf(x,y,image);
shading interp;
view(0,90)
hold on
%plot([zeros(1,13); 90*cosd(0:30:360)], [zeros(1,13); 90*sind(0:30:360)],'k')
%plot(90*((0:0.33:1)'*cosd(0:10:360))', 90*((0:0.33:1)'*sind(0:10:360))','k')
colorbar
set(colorbar,'FontSize',16)
axis equal
%set(gca, 'Box','off', 'XColor','none', 'YColor','none',  'Color','none')
hold off

data1(1:Nr)=0;
phipoint = 1;
for i = 1:Nr,
    data1(i) = data(i, phipoint)*10^30;
end;
figure(2)
hold on;
plot(radius(1:Nr),data1(1:Nr));
plot(profile(1:Nprofile,1), profile(1:Nprofile,2)*10^30);
errorbar(observedR(1:Ndata),observedFlux(1:Ndata),observedError(1:Ndata),'red','LineWidth',2);
