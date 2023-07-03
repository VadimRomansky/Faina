clear;
Bdata = importdata('../B.dat');

R = 258;
R0 = R*(1 - 2*0.25);
%R0=0;

Nr = 200;
Nz = 400;
Nphi = 8;

B(1:Nr, 1:Nz, 1:Nphi)=0;

for i = 1:Nr,
    for j = 1:Nz,
        for k = 1:Nphi,
            B(i,j,k) = Bdata((i-1)*Nz*Nphi + (j-1)*Nphi + k);
        end;
    end;
end;

zpoint = fix(Nz/2);
phipoint = 1;

radius(1:Nr)=0;
for i = 1:Nr,
    radius(i) = R0 + (R-R0)*(i-0.5)/Nr;
end;

B2d(1:Nphi+1,1:Nr)=0;
B1d(1:Nr) = 0;
for i = 1:Nr,
        for k = 1:Nphi,
            B2d(k,i) = B(i,zpoint,k);
        end;
        B2d(Nphi+1,i) = B(i, zpoint,Nphi);
        B1d(i) = B(i,zpoint, phipoint);
end;

figure(1)
[r,t] = meshgrid(radius,pi/Nphi:2*pi/Nphi:2*pi+pi/Nphi);
x = r.*cos(t);
y = r.*sin(t);
%contourf(x,y,image);
surf(x,y,B2d);
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


figure(2);
hold on;
plot(radius(1:Nr),B1d(1:Nr));
