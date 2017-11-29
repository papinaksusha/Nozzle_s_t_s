%initial values

x_N = 0.005;
x = 0 : x_N : 50;
%x = [0,50];
T_cr = 5000;
p_cr = 1*101325;
n_cr = p_cr/K/T_cr; 
n_N2_cr = 0.79;
n_O2_cr = 0.21;

i_N2 = 0 : I(SW_O, 1);
i_O2 = 0 : I(SW_O, 2);
i_NO = 0 : I(SW_O, 3);

e_i_N2 = H*C.*(W(1).*i_N2 - WX(1).*i_N2 - WX(1).*i_N2.^2);
e_i_O2 = H*C.*(W(2).*i_O2 - WX(2).*i_O2 - WX(2).*i_O2.^2);
e_i_NO = H*C.*(W(3).*i_NO - WX(3).*i_NO - WX(3).*i_NO.^2);

Z_vibr_N2 = sum(exp(-e_i_N2/K/T_cr));
Z_vibr_O2 = sum(exp(-e_i_O2/K/T_cr));

init = [n_N2_cr n_O2_cr 0 0 0 1 1];

v_cr = v_critical([n_N2_cr n_O2_cr 0 0 0] , T_cr);
v_cr = v_cr + v_cr*0.1;

v_cr = 1360.6;

%options = odeset('AbsTol', 1e-52, 'RelTol', 2.3e-14, 'OutputFcn', @odeplot, 'OutputSel', l_T);
options=odeset('AbsTol', 1e-54, 'RelTol', 2.3e-14,'Stats', 'on', 'OutputFcn', @odeplot, 'BDF', 'on', 'OutputSel', 7);
%%
[X,Y] = Nozzle_5_1T(x , init , options , T_cr , p_cr , v_cr);

toc

n_N2 = Y(: , 1)./sum(Y(: , 1 : 5) , 2);
n_O2 = Y(: , 2)./sum(Y(: , 1 : 5) , 2);
n_NO = Y(: , 3)./sum(Y(: , 1 : 5) , 2);
n_N = Y(: , 4)./sum(Y(: , 1 : 5) , 2);
n_O = Y(: , 5)./sum(Y(: , 1 : 5) , 2);
v = Y(: , 6).*v_cr;
T = Y(: , 7).*T_cr;

figure(1)
plot(X , T);


figure(4)
semilogy(X, [n_N2 n_O2 n_NO n_N n_O])
legend('N_2','O_2','NO','N','O');
ylabel('n_c/n');
xlabel('x/_r*');
xlim([0,5]);
%ylim([1e-5,1]);

for i = 1 : length(T)
n_i_N2_1t(i,:) = n_N2(i)/sum(exp(-e_i_N2/K/T(i))).*exp(-e_i_N2/K/T(i)); 
n_i_O2_1t(i,:) = n_O2(i)/sum(exp(-e_i_O2/K/T(i))).*exp(-e_i_O2/K/T(i)); 
end

xr = [0 1 2 3 5 10 15 20 25 30 40 50];

for i = 1 : length(xr)
    ind = find(X == xr(i));
    u_N2_1t(i,:) = n_i_N2_1t(ind, :);
    u_O2_1t(i,:) = n_i_O2_1t(ind, :);
end

filename = strcat('./MAT/1T_NOZ', num2str(SW_N), '_', num2str(p_cr/101325), '_' , num2str(T_cr), '_', 'OSC', ...
           num2str(SW_O), '_', '.mat');
save(filename, 'X', 'n_N2', 'n_O2', 'n_NO', 'n_N', 'n_O', 'v', 'T', 'n_i_N2_1t', 'n_i_O2_1t', 'u_N2_1t', 'u_O2_1t')