function [f1,f2] = Nozzle_5_1T (x,init,options,T_cr,p_cr,v_cr)


global K H C W WX M D I SW_O SW_N NA THETA_R

n_cr = p_cr/K/T_cr;

[f1,f2] = ode15s(@nozzle_5,x,init,options);

    function dy = nozzle_5(x,y)
    % y(1) = n_N2, y(2) = n_O2, y(3) = n_NO, y(4) = n_N, y(5) = n_O, y(6) = v, y(7) = T;   
    
    dy = zeros(7,1);
    b = zeros(7,1);
    A = zeros(7,7);
    
    % nozzle geometry
    
    switch SW_N
        case 1
              r_cr = 1e-3;
              alpha = 0.117*pi;
              S(1) = (1 + x*tan(alpha))^2;
              S(2) = 2*tan(alpha)*(1 + x*tan(alpha));
        case 2
              r_cr = 3e-3;
              alpha = 1/18*pi;
              S(1) = 1+x^2*(tan(alpha))^2;
              S(2) = 2*x*(tan(alpha))^2;
        case 3
              a = 0.3599;
              bb = 0.2277;
              cc = 0.1884;
              d = 0.0184;
              e = 0.1447;
              r_cr = a - bb - d/e;
              S(1) = 1/r_cr^2*(a-bb*exp(-cc*x^2*r_cr^2)-d/(x^2*r_cr^2+e))^2; 
              S(2) = 2/r_cr^2*(a-bb*exp(-cc*x^2*r_cr^2)-d/(x^2*r_cr^2+e))*(bb*cc*2*...
                     r_cr^2*x*exp(-cc*x^2*r_cr^2)+2*d*x*r_cr^2/(x^2*r_cr^2+e)^2);
    end
    
    % dimensionless variables
    
    v = y(6);
    T = y(7);
    M_cr = sum(M); 
    MM = M./M_cr;
    
   % n = [n_N2 n_O2 n_NO n_N n_O];
   % rho = sum([MM(1).*n_N2' MM(2).*n_O2' MM(3)*n_NO MM(4)*n_N MM(5)*n_O]);
    
    i_N2 = 0 : I(SW_O, 1);
    i_O2 = 0 : I(SW_O, 2);
    i_NO = 0 : I(SW_O, 3);

    e_i_N2 = H*C/K/T_cr.*(W(1).*i_N2 - WX(1).*i_N2 - WX(1).*i_N2.^2);    % 1*48
    e_i_O2 = H*C/K/T_cr.*(W(2).*i_O2 - WX(2).*i_O2 - WX(2).*i_O2.^2);
    e_i_NO = H*C/K/T_cr.*(W(3).*i_NO - WX(3).*i_NO - WX(3).*i_NO.^2);
    
    
    % vibrational energy of 0th levels
     
    e_0_N2 = H*C/K/T_cr*(W(1)/2 - WX(1)/4);
    e_0_O2 = H*C/K/T_cr*(W(2)/2 - WX(2)/4);
    e_0_NO = H*C/K/T_cr*(W(3)/2 - WX(3)/4);
    
    e_NO = (D(1)/2 + D(2)/2 - D(3))/K/T_cr;
    e_N = D(1)/2/K/T_cr;
    e_O = D(2)/2/K/T_cr;
    
    E_vibr_N2 = sum((e_i_N2 + e_0_N2).*exp(-e_i_N2./T))/sum(exp(-e_i_N2./T));
    E_vibr_O2 = sum((e_i_O2 + e_0_O2).*exp(-e_i_O2./T))/sum(exp(-e_i_O2./T));
    E_vibr_NO = sum((e_i_NO + e_0_NO).*exp(-e_i_NO./T))/sum(exp(-e_i_NO./T));
    
    C_vibr_N2 = 1/T^2*( sum((e_i_N2 + e_0_N2).*e_i_N2.*exp(-e_i_N2./T))/sum(exp(-e_i_N2./T))- ...
                sum((e_i_N2 + e_0_N2).*exp(-e_i_N2./T))*sum(e_i_N2.*exp(-e_i_N2./T))/(sum(exp(-e_i_N2./T)))^2 );
    C_vibr_O2 = 1/T^2*( sum((e_i_O2 + e_0_O2).*e_i_O2.*exp(-e_i_O2./T))/sum(exp(-e_i_O2./T))- ...
                sum((e_i_O2 + e_0_O2).*exp(-e_i_O2./T))*sum(e_i_O2.*exp(-e_i_O2./T))/(sum(exp(-e_i_O2./T)))^2 );
    C_vibr_NO = 1/T^2*( sum((e_i_NO + e_0_NO).*e_i_NO.*exp(-e_i_NO./T))/sum(exp(-e_i_NO./T))- ...
                sum((e_i_NO + e_0_NO).*exp(-e_i_NO./T))*sum(e_i_NO.*exp(-e_i_NO./T))/(sum(exp(-e_i_NO./T)))^2 );     
 
%% Матрица M ( M(x,y)*y' = f(x,y))

    A(1, :) = [y(6) 0 0 0 0 y(1) 0];
    A(2, :) = [0 y(6) 0 0 0 y(2) 0];
    A(3, :) = [0 0 y(6) 0 0 y(3) 0];
    A(4, :) = [0 0 0 y(6) 0 y(4) 0];
    A(5, :) = [0 0 0 0 y(6) y(5) 0];
    A(6, 1:5) = y(7);
    A(6,6) = M_cr*v_cr^2/K/T_cr*y(6)*sum(MM'.*y(1:5));
    A(6,7) = sum(y(1:5));
    A(7,1) = 5/2*y(7) + E_vibr_N2;
    A(7,2) = 5/2*y(7) + E_vibr_O2;
    A(7,3) = 5/2*y(7) + E_vibr_NO + e_NO;
    A(7,4) = 3/2*y(7) + e_N;
    A(7,5) = 3/2*y(7) + e_O;
    A(7,6) = 1/y(6)*( y(1)*(7/2*y(7) + E_vibr_N2) + ...
             y(2)*(7/2*y(7) + E_vibr_O2) + ...
             y(3)*(7/2*y(7) + E_vibr_NO + e_NO) + ... 
             y(4)*(5/2*y(7) + e_N) +...
             y(5)*(5/2*y(7) + e_O));
    A(7,7) = 5/2*(y(1) + y(2) + y(3)) + y(1)*C_vibr_N2 + y(2)*C_vibr_O2 + y(3)*C_vibr_NO + ...
             3/2*(y(4)+y(5));
    
%% Правые части

%    k_ar = Arrenius_coeff(y(7)*T_cr, 1, e_i_N2*T_cr, e_i_O2*T_cr, e_i_NO*T_cr);
%     
%     R(1) = n_cr^2*( y(3)*y(4)*k_ar(2) - y(1)*y(5)*k_ar(1) );
%     R(2) = n_cr^2*( y(3)*y(5)*k_ar(4) - y(2)*y(4)*k_ar(3) );
%     R(3) = -R(1) - R(2);
%     R(4) = n_cr^2*( n_cr*y(4)^2*( y(1)*k_ar(10) + y(2)*k_ar(11) + y(3)*k_ar(12) + y(4)*k_ar(13) + y(5)*k_ar(14) ) - ...
%                            y(1)*( y(1)*k_ar(5) + y(2)*k_ar(6) + y(3)*k_ar(7) + y(4)*k_ar(8) + y(5)*k_ar(9) ) );
%     
%     R(5) = n_cr^2*( n_cr*y(5)^2*( y(1)*k_ar(20) + y(2)*k_ar(21) + y(3)*k_ar(22) + y(4)*k_ar(23) + y(5)*k_ar(24) ) - ...
%                            y(2)*( y(1)*k_ar(15) + y(2)*k_ar(16) + y(3)*k_ar(17) + y(4)*k_ar(18) + y(5)*k_ar(19) ) );
%     
%     R(6) = n_cr^2*( n_cr*y(4)*y(5)*( y(1)*k_ar(30) + y(2)*k_ar(31) + y(3)*k_ar(32) + y(4)*k_ar(33) + y(5)*k_ar(34) ) - ...
%                               y(1)*( y(1)*k_ar(25) + y(2)*k_ar(26) + y(3)*k_ar(27) + y(4)*k_ar(28) + y(5)*k_ar(29) ) );
%     R(7) = -R(1) + R(2);
%     R(8) = -2*R(4) - R(6);
%     R(9) =  -R(2) + R(1);
%     R(10) = -2*R(5) - R(6);
   R = React1T(y(1)*n_cr,y(2)*n_cr,y(3)*n_cr,y(4)*n_cr,y(5)*n_cr,y(7)*T_cr,e_i_N2*K*T_cr,e_i_O2*T_cr,e_i_NO*T_cr,2);
sum(R);
%% International Workship
% [A_c_diss_n2, A_c_diss_o2, A_c_diss_no, A_c_diss_n, A_c_diss_o]

         An2_diss(1:5) = 2.5e13/NA;
        Ao2_diss(1:5) = 9.1e12/NA; 
        Ano_diss(1:5) = 4.1e12/NA;
        An2o_non = 7.4e5/NA;
        Anoo_o2n = 3e5/NA;
        
% [N_c_diss_n2, N_c_diss_o2, N_c_diss_no, N_c_diss_n, N_c_diss_o]    

        Nn2_diss = [-1 -1 -1 -1 -1];
        No2_diss = [-1 -1 -1 -1 -1];
        Nno_diss = [-1 -1 -1 -1 -1];
        Nn2o_non = 0.5;
        Nnoo_o2n = 0.5;
        
        En2_diss = [D(1) D(1) D(1) D(1) D(1)];
        Eo2_diss = [D(2) D(2) D(2) D(2) D(2)];
        Eno_diss = [D(3) D(3) D(3) D(3) D(3)];
        En2o_non = 37940*K;
        Enoo_o2n = 19460*K;
        
% Forward

        kn2_diss = An2_diss.*(T.*T_cr).^Nn2_diss.*exp(-En2_diss./K./T./T_cr);
        ko2_diss = Ao2_diss.*(T.*T_cr).^No2_diss.*exp(-Eo2_diss./K./T./T_cr);
        kno_diss = Ano_diss.*(T.*T_cr).^Nno_diss.*exp(-Eno_diss./K./T./T_cr);
        kn2o_non = An2o_non*(T.*T_cr)^Nn2o_non*exp(-En2o_non/K/T/T_cr);
        knoo_o2n = Anoo_o2n*(T.*T_cr)^Nnoo_o2n*exp(-Enoo_o2n/K/T/T_cr);
        
% Reverse

        Zint_n2 = sum(exp(-e_i_N2./T))*T*T_cr/2/THETA_R(1); 
        Zint_o2 = sum(exp(-e_i_O2./T))*T*T_cr/2/THETA_R(2);    
        Zint_no = sum(exp(-e_i_NO./T))*T*T_cr/THETA_R(3);

        kn2_rec = kn2_diss.*(2/M(4))^(3/2)*H^3*(2*pi*K*T*T_cr)^(-3/2)*Zint_n2*exp(D(1)/K/T/T_cr);
        ko2_rec = ko2_diss.*(2/M(5))^(3/2)*H^3*(2*pi*K*T*T_cr)^(-3/2)*Zint_o2*exp(D(2)/K/T/T_cr);
        kno_rec = kno_diss.*(M(3)/M(4)/M(5))^(3/2)*H^3*(2*pi*K*T*T_cr)^(-3/2)*Zint_no*exp(D(3)/K/T/T_cr);
        knon_n2o = kn2o_non*(2*M(5)/M(3))^(3/2)*Zint_n2/Zint_no*exp((D(1)-D(3))/K/T/T_cr);
        ko2n_noo = knoo_o2n*(2*M(4)/M(3))^(-3/2)*Zint_no/Zint_o2*exp((D(3)-D(2))/K/T/T_cr);

    R = zeros(1, 10);
    
    n_N2 = y(1)*n_cr;
    n_O2 = y(2)*n_cr;
    n_NO = y(3)*n_cr;
    n_N = y(4)*n_cr;
    n_O = y(5)*n_cr;
    
    R(1) = n_NO*n_N*knon_n2o - n_N2*n_O*kn2o_non; % RN2 2-2
    R(2) = n_NO*n_O*knoo_o2n - n_O2*n_N*ko2n_noo; % RO2 2-2
    R(3) = -R(1) - R(2); % RNO 2-2
    
    n = [n_N2 n_O2 n_NO n_N n_O];
    
    R(4) = sum(n.*(n_N^2.*kn2_rec - n_N2*kn2_diss)); % RN2 2-3
    R(5) = sum(n.*(n_O^2.*ko2_rec - n_O2*ko2_diss)); % RO2 2-3
    R(6) = sum(n.*(n_N*n_O.*kno_rec - n_NO*kno_diss)); % RNO 2-3
    
    R(7) = -R(1) + R(2);% RN 2-2
    R(8) = -R(2) + R(1);% RO 2-2
    R(9) = -2*R(4) - R(6);% RN 2-3
    R(10) = -2*R(5) - R(6);% RO 2-3
%sum(R)
    b(1) = r_cr/n_cr/v_cr*( R(1) + R(4) ) - y(1)*y(6)/S(1)*S(2);
    b(2) = r_cr/n_cr/v_cr*( R(2) + R(5) ) - y(2)*y(6)/S(1)*S(2);
    b(3) = r_cr/n_cr/v_cr*( R(3) + R(6) ) - y(3)*y(6)/S(1)*S(2);
    b(4) = r_cr/n_cr/v_cr*( R(7) + R(9) ) - y(4)*y(6)/S(1)*S(2);
    b(5) = r_cr/n_cr/v_cr*( R(8) + R(10) ) - y(5)*y(6)/S(1)*S(2);
    b(6) = 0;
    b(7) = -( y(1)*(7/2*y(7) + E_vibr_N2) + ...
              y(2)*(7/2*y(7) + E_vibr_O2) + ...
              y(3)*(7/2*y(7) + E_vibr_NO + e_NO) + ... 
              y(4)*(5/2*y(7) + e_N) +...
              y(5)*(5/2*y(7) + e_O))/S(1)*S(2);  
    dy = A\b;

    end
end

