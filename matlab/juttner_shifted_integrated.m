function [result] = juttner_shifted_integrated(G,theta,beta)
result = 0;
N = 1000;
dmu = 2/N;
for i = 1:N,
    mu = -1 + dmu*(i-0.5);
    result = result + juttner_shifted(G, mu, theta, beta)*dmu;
end;
end