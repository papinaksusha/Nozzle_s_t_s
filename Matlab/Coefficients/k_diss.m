function f = k_diss(sp,T)
% Dissociation rate coefficients, Treanor-Marrone model
%
% INPUT:
%
% sp - chemical specie
% T - temperature
% 
% OUTPUT:
%
% f(:, 1:3) -- k_i_diss_^mol, , m^3/s,  
%              size(f(:, 1:3)) = I(SW_O, sp), 3
% f(:, 4:5) -- k_i_diss_^at, , m^3/s,  
%              size(f(:, 1:2)) = I(SW_O, sp), 2

global D K I SW_O NA H C W WX

% Constants for Arrenius law from
% International Workshop on Radiation of High Temperature Gases in Atmospheric
% Entry. Part II // 30 Sep.- 10 Oct.2004. Porquerolles, France.

% как брал Синицын
A = [1.6e16 3.7e15; 
     9.99e15 1.99e15;
     0.3e12 0.41e12]./NA; % A{(sp,:) = [Aat Amol]
N = [-1.6 -1.6;
     -1.5 -1.5;
     0.5 -1];

%Оля

% A = [1.6e16 3.7e15; 
%      9.99e15 1.99e15;
%      0.41e13 0.41e13]./NA; % A{(sp,:) = [Aat Amol]
% N = [-1.6 -1.6;
%      -1.5 -1.5;
%      -1 -1];

 
i = 0 : I(SW_O, sp);
e = H*C.*(W(sp).*i - WX(sp).*i - WX(sp).*i.^2); 

Z_vibr = sum(exp(-e./K./T));                                               % Equilibrium vibrational partition function                                               
k_c_diss_eq = A(sp,:).*T.^N(sp,:).*exp(-D(sp)/K/T);                              % Thermal equilibrium dissociation rate coefficient, m^3/s
U = D(sp)/6/K;                                                          % Parameter (U = inf, U = D/6k, U = 3T)
%U = Inf;
Z_ci = Z_vibr/sum(exp(e./K./U)).*exp(e./K*(1/T + 1/U));             % Nonequilibrium factor Z_ci
f = zeros(I(SW_O,sp) + 1, 5);
f(:,1:3) = Z_ci'*k_c_diss_eq(2)*ones(1,3);
f(:,4:5) = Z_ci'*k_c_diss_eq(1)*ones(1,2);

end

