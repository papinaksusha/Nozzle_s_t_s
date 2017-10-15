function [VT, VV_N2, VV_O2, VV_NO] = k_ssh(sp,T)

%% The rate coefficients of the VV, VV' and VT transitions,SSH theory

% 1. wave_vibr.pdf
% 2. G.G. Chernyi, S.A. Losev, S.O. Macheret, B.V. Potapkin. Physical and 
% Chemical Processes in Gas Dynamics: Cross Sections and Rate Constants, Vol.1

% INPUT:
%
% sp -- chemical specie sp = 1, 2, 3 (1 = N2, 2 = O2, 3 = NO)
% T  -- temperature

% OUTPUT:
%
% VT   -- k_{c, i+1->i}^M, , m^3/s, 
%         size(VT) = I(SW_O, c), 5
% VV_c -- k_{sp, i+1->i}^{c, m->m+1}, m^3/s,
%         size(VV_c) = I(SW_O, sp), I(SW_O, c)

global M C K H W WX SW_O I

d = [3.621 3.458 3.47 3.298 2.75].*1e-10;                                  % m, molecular diameter 
r0 = [3.621 3.458 3.47 3.298 2.75].*1e-10;                                 % m, the distance at which the interparticle Lennard-Jones potential is zero
r0 = diag(r0);
eps = [97.5 107.4 119 71.4 80];                                            % eps/k, K, eps -- the minimum value of the potential (the depth of the potential well)
eps = diag(eps);
%re = [1.097 1.207 1.151].*1e-10;% для Z0 зависящего от диаметра

for y = 1 : 5
    for j = 1 : 5
        eps(y,j) = sqrt(eps(y,y)*eps(j,j));
        r0(y,j) = 0.5*(r0(y,y) + r0(j,j));
    end
end

% VT

Z0 = 3;                                                                    % the steric factor

mu = M(sp).*M./(M(sp) + M);                                                % kg, the reduced oscillator mass
omega = 2*pi*C*(W(sp) - 2*WX(sp));                                         % c^-1, the angular frequency of the oscillator  
alpha = 17.5./r0(sp,:);                                                    % m-1
hi = ((0.5*pi^2*omega^2/K/T).*mu./(alpha.^2)).^(1/3);
R = (0.5.*sqrt(1 + hi.*T./eps(sp,:)) + 0.5).^(-1/3);                       % (r/r0)^2

%r = r0(sp,:).*(0.5.*sqrt(1 + hi.*T./eps(sp,:) + 0.5)).^(-1/6);
%Z0 = (alpha.*re(sp)).^2.*exp(-3/8.*alpha.*re(sp).^2./r); % так приблизительно сходится P10 c wavevibr

P10 = 1.294.*R./Z0./(1 + 1.1.*eps(sp,:)./T).*8.*pi^3.*mu.*omega./...       % the averaged probability for the VT transition A(1) + B -> A(0) + B
      alpha.^2./H.*sqrt(4*pi/3.*hi).*exp(-3.*hi + H*omega/4/pi/K/T + eps(sp,:)./T);
Zn = (d(sp) + d).^2.*sqrt(pi*K*T/2./mu);                                   % collision frequency per unit number density
k10 = P10.*Zn;

i = {0 : I(SW_O,1) - 1, 0 : I(SW_O,2) - 1, 0 : I(SW_O,3) - 1};
dE = H*C.*WX(sp);
E1 = H*C.*(W(sp) - 2*WX(sp));

if SW_O == 1
        VT = (cell2mat(i(sp)) + 1)'*k10;
else
        gamma0 = 2*pi^2.*E1./alpha./H.*sqrt(0.5.*mu./K./T);
        g_i = (2*pi^2.*(E1 - 2.*cell2mat(i(sp)).*dE))'*(alpha.*H./sqrt(0.5.*mu./K./T)).^(-1);
      %  g_i = pi*(W(sp) - 2.*(cell2mat(i(sp)) + 1).*WX(sp))'*(alpha.*sqrt(0.5.*mu./K./T)).^(-1); % из справочника
        deltaVT = zeros(I(SW_O,sp),5);
        l = find(g_i >= 20);
        ll = find(g_i < 20);
        deltaVT(l) = 4.*gamma0(fix((l - 1)./I(SW_O,sp)) + 1).^(2/3).*dE./E1;
        deltaVT(ll) = 4/3.*gamma0(fix((ll - 1)./I(SW_O,sp)) + 1).*dE./E1;
        VT = (cell2mat(i(sp)) + 1)'*k10.*exp((ones(1,5)'*cell2mat(i(sp)))'.*deltaVT).*exp(-(ones(1,5)'*cell2mat(i(sp)))'.*H*C*WX(sp)/K/T);

end

% VV & VV'

m_r(1) = M(4)*M(4)/(M(4) + M(4));
m_r(2) = M(5)*M(5)/(M(5) + M(5));
m_r(3) = M(4)*M(5)/(M(4) + M(5));
lambda1 = 0.5;
lambda2 = 0.5;

Q1001 = lambda1^2*lambda2^2*4.*alpha(1:3).^2.*K.*T./omega^2./m_r(sp);      % the averaged probability for the VV transition A(1) + B(0) -> A(0) + B(1)
k1001 = Q1001.*Zn(1:3);

VV = cell(1,3);                                                            

if SW_O == 1
    for l = 1:3
       VV(l) = {(cell2mat(i(sp)) + 1)'*(cell2mat(i(l)) + 1).*k1001(l)};
    end
else
    deltaVV = 8/3*pi^2*dE/H/alpha(sp)*sqrt(mu(sp)/2/K/T);
    %deltaVV = 4/3*pi*WX(sp)./alpha(sp)*sqrt(mu(sp)/2/K/T);
    for l = 1:3
        if (l == sp)
            A = cell2mat(i(sp))'*ones(1,I(SW_O,sp)) - ...
                ones(I(SW_O,sp),1)*cell2mat(i(sp));
            VV(l) = {(cell2mat(i(sp)) + 1)'*(cell2mat(i(l)) + 1).*k1001(l).* ...
                exp(-deltaVV.*abs(A)).*(1.5 - 0.5.*exp(-deltaVV.*abs(A))).* ...
                exp(A'.*H.*C.*WX(sp)./K./T)};
        else
            deltaVVs = 8/3*pi^2*C*WX(l)/alpha(l)*sqrt(mu(l)/2/K/T);
            %deltaVVs = 4/3*pi*WX(l)./alpha(l)*sqrt(mu(l)/2/K/T);
            p = (W(sp) - W(l) - 2*(WX(l) - WX(sp)))/2/WX(sp);
            A = deltaVVs.*ones(I(SW_O,sp) , 1)*cell2mat(i(l)) - ...
            deltaVV.*cell2mat(i(sp))'*ones(1 , I(SW_O,l)) + ...
            deltaVV*p.*ones(I(SW_O,sp) , I(SW_O,l));
            B = ones(I(SW_O,sp) , 1)*cell2mat(i(l)).*H.*C*WX(sp) - ...
                cell2mat(i(sp))'*ones(1, I(SW_O,l)).*H.*C.*WX(l);
            VV(l) = {(cell2mat(i(sp)) + 1)'*(cell2mat(i(l)) + 1).*k1001(l).* ...
                    exp(-abs(A)).*exp(deltaVV.*p).*exp(B./K./T)};
        end
     end
end

VV_N2 = cell2mat(VV(1));
VV_O2 = cell2mat(VV(2));
VV_NO = cell2mat(VV(3));

end