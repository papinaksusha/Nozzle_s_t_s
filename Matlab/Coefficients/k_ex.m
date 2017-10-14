function f = k_ex(sp,T)

% Exchange rate coefficients, Warnatz model k = A*(i+1)*T^b*exp[-(Ea-eps_i/k/T)*Th(E-eps_i)]
%
% INPUT:
%
% sp - chemical specie (c = N2, O2)
% T - temperature
%
% OUTPUT:
%
% f(:, model) = k_{ci, NO}^{at1->at2}, m^3/s,
%               size(f) = I(SW_O, sp), 5

global K NA SW_O I M H C W WX

i = 0 : I(SW_O, sp);
f = zeros(I(SW_O,sp) + 1,5); % 4 модели

e = H*C.*(W(sp).*i - WX(sp).*i(sp) - WX(sp).*i.^2);   % Vibrational energy of moleculese = e_i_c;

%Ea = [319 37.37]./NA.*1e3;  % Activation energy, J из Варнатца
Ea = [3.2 0.33].*1.6e-19; % Лосев, справочник, Т 2 стр. 97

%% 1. Русанов, Фридман

r0 = [(3.621 + 2.75)/2, (3.458 + 3.298)/2 ].*1e-10;
mu = [M(1)*M(5)/(M(1) + M(5)) , M(2)*M(4)/(M(2) + M(4))];
A = pi*r0(sp)^2.*sqrt(8*K*T/pi/mu(sp));
alpha = [0.51, 0.24];
f(:,1) = A.*exp(-(Ea(sp)-alpha(sp).*e)./K./T.*(Ea(sp) - alpha(sp).*e >= 0));

%% 2. Полак, Гольденберг

betta = [0.9, 0.46];
gamma = [0.52, 0.12];
f(:,2) = A.*exp(-(Ea(sp) - gamma(sp).*e)./betta(sp)./K./T.*(Ea(sp) - gamma(sp).*e >= 0));

%% 3. Warnatz model
        % Warnatz J., Riedel U., Schmidt R. Different levels of air dissociation
        % and Physico-chimitchesKie protsessi v gazovoy dynamiKe, T.2, p. 96 of 370
        A = [4.168e12 1.151e9].*1e-6./NA;                                          % Model constant, M^3/s
                                                % Activation energy, J               
        betta = [0,1];
        f(:,3) = A(sp).*(i + 1).*T.^betta(sp).*exp(-(Ea(sp) - e)./...
            K./T.*(Ea(sp) - e >= 0));
%% 4.  Капителли как Синицын
       
    if sp == 2
        f(:,4) = 8.4e12*exp(-19450/T)/NA*1e-6;
    else
        v = 0 : I(SW_O,sp);
        Ev = 3395.*v.*(1 - 6.217e-2.*(v + 1));
        f(1:9,4) = exp(-23.044680).*(Ev(1:9) + 3000).^(-0.419312)/T^(-0.37836).*...
                   exp(38370*0.99243/T + Ev(1:9)*0.989385/T).*1e-6;
        f(10:13,4) = exp(1.423118).*(Ev(10:13) + 3000).^(-3.42306)/T^(-1.4234).*...
                     exp(38370*(-0.919692)/T + Ev(10:13)*0.917323/T).*1e-6;
        f(14:end,4) = exp(-96.75885).*real((Ev(14:end) + 3000).^(6.480504))/T^(-0.27937).*...
                      exp(38370*(-0.037869)/T + Ev(14:end)*0.019647/T.*1e-6);
      %  f(1:13,4) = 1.3e-10*exp(-38000/T)*1e-6;
       % f(14:end,4) = 1e-11*1e-6;

    end
%% 5. Без обменных реакци

    f(:,5) = 0;

end
