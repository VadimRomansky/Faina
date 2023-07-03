clear;
data = importdata('../image.dat');
R = 258;
R0 = R*(1 - 2*0.05);
%R0=0;
Nr = size(data,1);
Nphi = size(data,2);

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
    data1(i) = data(i, phipoint);
end;
figure(2)
plot(radius(1:Nr),data1(1:Nr));
hold on;
