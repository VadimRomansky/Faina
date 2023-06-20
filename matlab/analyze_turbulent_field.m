clear;
Bxdata = importdata("../Bx.dat");
Bydata = importdata("../By.dat");
Bzdata = importdata("../Bz.dat");

Nr = 200;
Nz = 200;
Nphi = 4;

Bx(1:Nr, 1:Nz, 1:Nphi)=0;
By(1:Nr, 1:Nz, 1:Nphi)=0;
Bz(1:Nr, 1:Nz, 1:Nphi)=0;

Bxmean = 0;
Bymean = 0;
Bzmean = 0;

volume = 0;

for i = 1:Nr,
    for j = 1:Nz,
        for k = 1:Nphi,
            Bx(i,j,k) = Bxdata((i-1)*Nz*Nphi + (j-1)*Nphi + k);
            By(i,j,k) = Bydata((i-1)*Nz*Nphi + (j-1)*Nphi + k);
            Bz(i,j,k) = Bzdata((i-1)*Nz*Nphi + (j-1)*Nphi + k);

            Bxmean = Bxmean + Bx(i,j,k)*Bx(i,j,k)*i;
            Bymean = Bymean + By(i,j,k)*By(i,j,k)*i;
            Bzmean = Bzmean + Bz(i,j,k)*Bz(i,j,k)*i;

            volume = volume + i;
        end;
    end;
end;

Bxmean = Bxmean/volume;
Bymean = Bymean/volume;
Bzmean = Bzmean/volume;

relation = 2*Bzmean/(Bxmean + Bymean);


