function [result] = juttner_shifted(G, mu,theta,beta)
beta2 = beta*beta;
G2 = G*G;
mu2 = mu*mu;
Gsqrt = sqrt(G2 - 1);
gamma = 1.0/sqrt(1.0 - beta2);
gamma2 = gamma*gamma;
bes = besselk(2, 1/theta);

A = (G2-1)*(gamma2*mu2+beta2*gamma2+1-mu2)+beta2*gamma2-2*beta*gamma2*mu*G*Gsqrt;
J = (G-Gsqrt*beta*mu)*gamma*Gsqrt/sqrt(A*(A+1));
newG = sqrt(A+1);
result = (2.0*pi*newG*sqrt(newG*newG-1)/(4*pi*theta*bes))*exp(-newG/theta)*J;
end