clc
clear 
close all
format long e
addpath(('./MAT/'))
addpath(('./ODE_systems/')) % переписать кросс-платформенно
addpath(('./Coefficients/'))

load('.\MAT\1T.mat')

global K H C W WX I SW_O SW_N

p_cr = 101325.*[1 1 10 10 100 100];
T_cr = [5000 7000 5000 7000 5000 7000];
N = 10000;
hx = 50/N;
X = (0:hx:50)';
xr = [0 1 2 3 5 10 15 20 25 30 40 50];

i_N2 = 0 : I(SW_O, 1);
i_O2 = 0 : I(SW_O, 2);
i_NO = 0 : I(SW_O, 3);

e_i_N2 = H*C.*(W(1).*i_N2 - WX(1).*i_N2 - WX(1).*i_N2.^2);
e_i_O2 = H*C.*(W(2).*i_O2 - WX(2).*i_O2 - WX(2).*i_O2.^2);
e_i_NO = H*C.*(W(3).*i_NO - WX(3).*i_NO - WX(3).*i_NO.^2);

for i = 1: 6
    n_N2 = nn2(i, :);
    n_O2 = no2(i, :);
    n_NO = nno(i, :);
    n_N = nn(i, :);
    n_O = no(i, :);
    TT = (T(i, :))';
    vv = (v(i, :))';
    
    n_i_N2_1t = zeros(length(TT), I(SW_O, 1) + 1);
    n_i_O2_1t = zeros(length(TT), I(SW_O, 2) + 1);
    u_N2_1t = zeros(length(xr), I(SW_O, 1) + 1);
    u_O2_1t = zeros(length(xr), I(SW_O, 2) + 1);
    
    for j = 1 : length(TT)
        n_i_N2_1t(j,:) = n_N2(j)/sum(exp(-e_i_N2/K/TT(j))).*exp(-e_i_N2/K/TT(j)); 
        n_i_O2_1t(j,:) = n_O2(j)/sum(exp(-e_i_O2/K/TT(j))).*exp(-e_i_O2/K/TT(j)); 
    end
    
    for k = 1 : length(xr)
        ind = find(X == xr(k));
        u_N2_1t(k,:) = n_i_N2_1t(ind, :);
        u_O2_1t(k,:) = n_i_O2_1t(ind, :);
    end
    
    filename = strcat('1T_NOZ', num2str(SW_N), '_', num2str(p_cr(i)/101325), '_' , num2str(T_cr(i)), '_', 'OSC', ...
            num2str(SW_O), '_', '.mat');
    save(filename, 'X', 'n_N2', 'n_O2', 'n_NO', 'n_N', 'n_O', 'vv', 'TT', 'n_i_N2_1t', 'n_i_O2_1t', 'u_N2_1t', 'u_O2_1t')
end
 