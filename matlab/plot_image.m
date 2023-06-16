clear;
data = importdata('../image.dat');
R = 1.4*10^17;
Nr = size(data,1);
Nphi = size(data,2);

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

Az(1:Nphi*Nr) = 0;
Rd(1:Nphi*Nr) = 0;
data2(1:Nphi*Nr) = 0;
sz(1:Nphi*Nr) = 1;
for i = 1:Nr,
    for j = 1:Nphi,
        Az((i-1)*Nphi + j) = (j+0.5)*2*pi/Nphi;
        Rd((i-1)*Nphi + j) = (i+0.5)*R/Nr;
        data2((i-1)*Nphi + j) = data(i,j);
    end;
end;

cn = 100;                                             % Number Of Colors
cm = colormap(jet(cn));
figure(1)
%polarscatter(Az, Rd, sz, cm(fix(data2),:), 'filled')
polarscatter(Az, Rd, sz, 'red', 'filled')
grid on