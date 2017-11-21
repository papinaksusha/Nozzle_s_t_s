%initial values
l_N2 = I(SW_O , 1) + 1;
l_O2 = I(SW_O , 2) + 1;
l_mol =  l_N2 + l_O2 + 1;
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

init = zeros(l_T , 1);

init(1 : l_N2) = n_N2_cr/Z_vibr_N2.*exp(-e_i_N2./K./T_cr);
init(l_N2 + 1 : l_N2 + l_O2) = n_O2_cr/Z_vibr_O2.*exp(-e_i_O2./K./T_cr);
init(l_v : l_T) = 1;
v_cr = v_critical([sum(init(1 : l_N2)) sum(l_N2 + 1 : l_N2 + l_O2) ...
                   init(l_mol) init(l_mol + 1) init(l_c)] , T_cr);
v_cr = v_cr + v_cr*0.5;

v_cr = 5000;

%options = odeset('AbsTol', 1e-52, 'RelTol', 2.3e-14, 'OutputFcn', @odeplot, 'OutputSel', l_T);
options=odeset('AbsTol', 1e-54, 'RelTol', 2.3e-14,'Stats', 'on', 'OutputFcn', @odeplot, 'BDF', 'on', 'OutputSel', l_T);

[X,Y] = Nozzle_5_full(x , init , options , T_cr , p_cr , v_cr);

toc

n_i_N2 = Y(: , 1 : l_N2)./(sum(Y(: , 1 : l_c) , 2)*ones(1 , l_N2));
n_i_O2 = Y(: , l_N2 + 1 : l_N2 + l_O2)./(sum(Y(: , 1 : l_c) , 2)*ones(1 , l_O2));
n_N2 = sum(n_i_N2 , 2);
n_O2 = sum(n_i_O2 , 2);
n_NO = Y(: , l_mol)./sum(Y(: , 1 : l_c) , 2);
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

figure(2)

for i = 1 : length(xr)
    for g = 1 : length(i_N2)
        u_N2(i,g) = interp1q(X,n_i_N2(:,g),xr(i));
    end
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
    for g = 1 : length(i_O2)
        u_O2(i,g) = interp1q(X,n_i_O2(:,g),xr(i));
    end
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
semilogy(X, [n_N2 n_O2 n_NO n_N n_O])
legend('N_2','O_2','NO','N','O');
ylabel('n_c/n');
xlabel('x/_r*');
xlim([0,5]);
%ylim([1e-5,1]);


filename = strcat('./MAT/NOZ', num2str(SW_N), '_', num2str(p_cr/101325), '_' , num2str(T_cr), '_', 'OSC', ...
           num2str(SW_O), '_', 'EX', num2str(EX_MODEL), '_', 'REC', num2str(REC), '_', '.mat');
save(filename, 'X', 'n_i_N2', 'n_i_O2', 'n_N2', 'n_O2', 'n_NO', 'n_N', 'n_O', 'v', 'T','u_N2','u_O2','e_i_N2', 'e_i_O2', 'i_N2', 'i_O2')
%% conservation laws

N_init = [sum(init(1 : l_N2)) sum(init(l_N2 + 1 : l_N2 + l_O2)) ...
          init(l_mol : l_c)'].*n_cr;
n_i_N2_d = Y(: , 1 : l_N2).*n_cr;
n_i_O2_d = Y(: , l_N2 + 1 : l_N2 + l_O2).*n_cr;
n_N2_d = sum(n_i_N2_d , 2);
n_O2_d = sum(n_i_O2_d , 2);
n_NO_d = Y(: , l_mol).*n_cr;
n_N_d = Y(: , l_mol + 1).*n_cr;
n_O_d = Y(: , l_c).*n_cr;
N = [n_N2_d n_O2_d n_NO_d n_N_d n_O_d];
rho = sum(ones(length(X),1)*M.*N , 2);

%  continuity equation rho*v*S = const (S - conical)

r_cr = 1e-3;
r = r_cr.*(1 + X.*tan(0.117*pi)); 

Q_cr = sum(M.*N_init)*v_cr*r_cr^2;
Q = rho.*r.^2.*v;
eps1 = max(abs(Q - Q_cr)./Q_cr);

% energy equation

R_c = K*NA./MOLAR;
rho_c = (ones(length(X),1)*M.*N);
Y_c =  rho_c./(rho*ones(1,5));

h_c(:,1) = 3.5*R_c(1).*T + 1./rho_c(:,1).*sum(ones(length(X),1)*(e_i_N2 + H*C*(W(1)/2 - WX(1)/4)).*n_i_N2_d , 2);
h_c(:,2) = 3.5*R_c(2).*T + 1./rho_c(:,2).*sum(ones(length(X),1)*(e_i_O2 + H*C*(W(2)/2 - WX(2)/4)).*n_i_O2_d , 2);
h_c(:,3) = 3.5*R_c(3).*T + H*C*(W(3)/2 - WX(3)/4)/M(3) + (D(1)/2 + D(2)/2 - D(3))/M(3);
h_c(:,4) = 2.5*R_c(4).*T + D(1)/2;
h_c(:,5) = 2.5*R_c(5).*T + D(2)/2;

H = v.^2./2 + sum(Y_c.*h_c , 2);
eps2 = max(abs(H - H(1))./H(1));