clc
clear
close all
format long e

x = 0: 0.01 : 50;

% conical

r_cr1 = 1e-3;
alpha1 = 0.117*pi;
r1 = r_cr1 + x.*r_cr1.*tan(alpha1);
S1 = pi*r1.^2./2./r_cr1^2;

% hyperbolic

r_cr2 = 3e-3;
alpha2 = 1/18*pi;
r2 = r_cr2.*sqrt(1 + (x.*r_cr2).^2.*(tan(alpha2)/r_cr2)^2);
S2 = pi*r2.^2./2./r_cr2^2;

% F4

a = 0.3599;
bb = 0.2277;
cc = 0.1884;
d = 0.0184;
e = 0.1447;
r = a - bb - d/e;
r_cr3 = a - bb - d/e;
xx = x.*r_cr3;
r3 = a - bb.*exp(-cc.*xx.^2) - d./(xx.^2 + e);
S3 = pi.*r3.^2./2./r_cr3^2;

figure(1)
semilogy(x, [S1; S2; S3])
legend('conical','hyperbolic','F4')
xlabel('x/r^*')
ylabel('S(x)/r^*, m^2')

figure(2)
plot(x, [r1; r2; r3])
legend('conical','hyperbolic','F4')
xlabel('x/r^*')
ylabel('r(x), m')