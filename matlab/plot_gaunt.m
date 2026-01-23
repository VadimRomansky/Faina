clear;

etamin = 0.0001;
etamax = 100;
thetamin = 0.001;
thetamax = 1000;

logetamin = -4;
logetamax = 2;
logthetamin = -3;
logthetamax = 3;

N = 100;

dlogeta = (logetamax - logetamin)/(N+1);
dlogtheta = (logthetamax - logthetamin)/(N+1);

gaunt(N,N)=0;
for i = 1:N,
    logeta = logetamin + (i-0.5)*dlogeta;
    eta = 10^logeta;
    for j = 1:N,
        logtheta = logthetamin + (j-0.5)*dlogtheta;
        theta = 10^logtheta;
        if(logeta > 0)
            if(logtheta > 0)
                gaunt(i,j) = sqrt(3/(3.14*theta));
            else
                gaunt(i,j) = (sqrt(3)/3.14)*log(4/(1.78*theta));
            end;
        else
            if(logtheta > 0)
                if(logtheta > -logeta)
                    gaunt(i,j) = sqrt(12/(eta*theta));
                else
                    gaunt(i,j) = 1.0;
                end;
            else
                if(-logtheta > -0.5*logeta)
                    gaunt(i,j) = (sqrt(3)/3.14)*log(4*sqrt(eta)/((1.78^2.5)*theta));
                else
                    gaunt(i,j) = 1.0;
                end;
            end
        end;
    end;
end;

figure(1);
colormap Jet;
[X, Y] = meshgrid(logthetamin+(1:N)*dlogtheta, logetamin+(1:N)*dlogeta);
surf(X, Y, gaunt);
shading interp;
title ('g');
xlabel ('theta');
ylabel ('eta');
zlabel ('g');
grid ;
