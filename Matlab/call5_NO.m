%initial values
l_N2 = I(SW_O , 1) + 1;
l_O2 = I(SW_O , 2) + 1;
l_NO = I(SW_O , 3) + 1;
l_mol =  l_N2 + l_O2 + l_NO;
l_c = l_mol + 2;
l_v = l_c + 1;
l_T = l_v + 1;

x_N = 0.005;
x = 0 : x_N : 50;
%x = [0,50];
T_cr = 7000;
p_cr = 100*101325;
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
Z_vibr_NO = sum(exp(-e_i_NO/K/T_cr));

init = zeros(l_T , 1);

init(1 : l_N2) = n_N2_cr/Z_vibr_N2.*exp(-e_i_N2./K./T_cr);
init(l_N2 + 1 : l_N2 + l_O2) = n_O2_cr/Z_vibr_O2.*exp(-e_i_O2./K./T_cr);
init(l_v : l_T) = 1;
v_cr = v_critical([sum(init(1 : l_N2)) sum(l_N2 + 1 : l_N2 + l_O2) ...
                   sum(init(l_N2 + l_O2 + 1 : l_mol)) init(l_mol + 1) init(l_c)] , T_cr);
v_cr = v_cr + v_cr*0.2;

%v_cr = 1607.84;
%v_cr = [1360.6 1607.84 1360.6 1607.84 1360.6 1607.84];

%options = odeset('AbsTol', 1e-52, 'RelTol', 2.3e-14, 'OutputFcn', @odeplot, 'OutputSel', l_T);
options=odeset('AbsTol', 1e-54, 'RelTol', 2.3e-14,'Stats', 'on', 'OutputFcn', @odeplot, 'BDF', 'on', 'OutputSel', l_T);

[X,Y] = Nozzle_5_NO_full(x , init , options , T_cr , p_cr , v_cr);

toc

n_i_N2 = Y(: , 1 : l_N2)./(sum(Y(: , 1 : l_c) , 2)*ones(1 , l_N2));
n_i_O2 = Y(: , l_N2 + 1 : l_N2 + l_O2)./(sum(Y(: , 1 : l_c) , 2)*ones(1 , l_O2));
n_i_NO = Y(: , l_N2 + l_O2 + 1 : l_mol)./(sum(Y(: , 1 : l_c) , 2)*ones(1 , l_NO));
n_N2 = sum(n_i_N2 , 2);
n_O2 = sum(n_i_O2 , 2);
n_NO = sum(n_i_NO , 2);
n_N = Y(: , l_mol + 1)./sum(Y(: , 1 : l_c) , 2);
n_O = Y(: , l_c)./sum(Y(: , 1 : l_c) , 2);
v = Y(: , l_v).*v_cr;
T = Y(: , l_T).*T_cr;

figure(1)
plot(X , T);
%save('MAT/1_5000_war_anhar_full.mat','X','Y');

%%

xr = [0 1 2 3 5 10 15 20 25 30 40 50];

colours = colormap(jet(length(xr)));

u_N2 = zeros(length(xr),length(i_N2));
u_O2 = zeros(length(xr),length(i_O2));
u_NO = zeros(length(xr),length(i_O2));

figure(2)

for i = 1 : length(xr)
    ind = find(X == xr(i));
    u_N2(i,:) = n_i_N2(ind, :);
    semilogy(i_N2, u_N2(i,:), 'color', colours(i,:)), hold on
end

ylabel('n_i/n');
xlim([0, I(SW_O,1)]);
%ylim([1e-15,1])
xlabel('i');
title('N_2');
legend('x/r = 0', 'x/r = 1', 'x/r = 2', 'x/r = 3', 'x/r = 5', 'x/r = 10', ...
       'x/r = 15', 'x/r = 20', 'x/r = 25', 'x/r = 30', 'x/r = 40', 'x/r =50');
hold off

figure(3)

for i = 1 : length(xr)
    ind = find(X == xr(i));
    u_O2(i,:) = n_i_O2(ind, :);
    semilogy(i_O2, u_O2(i,:), 'color', colours(i,:)), hold on
end

ylabel('n_i/n');
xlim([0,I(SW_O,2)]);
%ylim([1e-15,1])
xlabel('i');
title('O_2');
legend('x/r = 0', 'x/r = 1', 'x/r = 2', 'x/r = 3', 'x/r = 5', 'x/r = 10', ...
       'x/r = 15', 'x/r = 20', 'x/r = 25', 'x/r = 30', 'x/r = 40', 'x/r =50');
hold off

figure(4)

for i = 1 : length(xr)
    ind = find(X == xr(i));
    u_NO(i,:) = n_i_NO(ind, :);
    semilogy(i_O2, u_NO(i,:), 'color', colours(i,:)), hold on
end

ylabel('n_i/n');
xlim([0,I(SW_O,3)]);
%ylim([1e-15,1])
xlabel('i');
title('NO');
legend('x/r = 0', 'x/r = 1', 'x/r = 2', 'x/r = 3', 'x/r = 5', 'x/r = 10', ...
       'x/r = 15', 'x/r = 20', 'x/r = 25', 'x/r = 30', 'x/r = 40', 'x/r =50');
hold off

figure(5)
semilogy(X, [n_N2 n_O2 n_NO n_N n_O])
legend('N_2','O_2','NO','N','O');
ylabel('n_c/n');
xlabel('x/_r*');
xlim([0,5]);
%ylim([1e-5,1]);


filename = strcat('./MAT/NO_NOZ', num2str(SW_N), '_', num2str(p_cr/101325), '_' , num2str(T_cr), '_', 'OSC', ...
           num2str(SW_O), '_', 'EX', num2str(EX_MODEL), '_', 'REC', num2str(REC), '_', '.mat');
save(filename, 'X', 'n_i_N2', 'n_i_O2', 'n_N2', 'n_O2', 'n_NO', 'n_N', 'n_O', 'v', 'T','u_N2','u_O2','e_i_N2', 'e_i_O2', 'i_N2', 'i_O2')