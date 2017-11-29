function coefficient = Arrenius_coeff(T,sw,e_n2,e_o2,e_no)

global D M THETA_R K NA H 

Dn2 = D(1);
Do2 = D(2);
Dno = D(3);

mn2 = M(1);
mo2 = M(2);
mno = M(3);
mn = M(4);
mo = M(5);

th_n2_rot = THETA_R(1);
th_o2_rot = THETA_R(2);
th_no_rot = THETA_R(3);

switch sw
    case 1
%% DSMC, Boyd
% [A_c_diss_n2, A_c_diss_o2, A_c_diss_no, A_c_diss_n, A_c_diss_o]

        An2_diss = [4.1e-12 1.5e-11 1.5e-11 1e-11 4e-12];
        Ao2_diss = [1.3e-10 5.33e-11 1.1e-10 1.1e-10 1.5e-10]; 
        Ano_diss = [2.1e-10 2e-10 1e-10 4e-10 4e-10];
        An2o_non = 0.8e-16;
        Ao2n_noo = 4e-15;
        
% [N_c_diss_n2, N_c_diss_o2, N_c_diss_no, N_c_diss_n, N_c_diss_o]    

        Nn2_diss = [-0.62 -0.68 -0.68 -0.68 -0.54];
        No2_diss = [-1 -1 -1 -1 -1.05];
        Nno_diss = [-1 -1 -1 -1.1 -1.1];
        Nn2o_non = 0;
        No2n_noo = -0.39;
        
        En2_diss = [15.67e-19 15.67e-19 15.67e-19 15.67e-19 15.67e-19];
        Eo2_diss = [8.197e-19 8.197e-19 8.197e-19 8.197e-19 8.197e-19];
        Eno_diss = [10.43e-19 10.43e-19 10.43e-19 10.43e-19 10.43e-19];
        En2o_non = 5.175e-19;
        Eo2n_noo = 0.2e-19;
        
% Forward

        kn2_diss = An2_diss.*T.^Nn2_diss.*exp(-En2_diss./K./T);
        ko2_diss = Ao2_diss.*T.^No2_diss.*exp(-Eo2_diss./K./T);
        kno_diss = Ano_diss.*T.^Nno_diss.*exp(-Eno_diss./K./T);
        kn2o_non = An2o_non.*T.^Nn2o_non.*exp(-En2o_non./K./T);
        ko2n_noo = Ao2n_noo.*T.^No2n_noo.*exp(-Eo2n_noo./K./T);
        
% Reverse

        Zint_n2 = sum(exp(-e_n2./T))*T/2/th_n2_rot; 
        Zint_o2 = sum(exp(-e_o2./T))*T/2/th_o2_rot;    
        Zint_no = sum(exp(-e_no./T))*T/th_no_rot;

        kn2_rec = kn2_diss.*(mn2/mn/mn)^(3/2)*H^3*(2*pi*K*T)^(-3/2)*Zint_n2*exp(Dn2/K/T);
        ko2_rec = ko2_diss.*(mo2/mo/mo)^(3/2)*H^3*(2*pi*K*T)^(-3/2)*Zint_o2*exp(Do2/K/T);
        kno_rec = kno_diss.*(mno/mn/mo)^(3/2)*H^3*(2*pi*K*T)^(-3/2)*Zint_no*exp(Dno/K/T);
        knon_n2o = kn2o_non*(mn2*mo/mno/mn)^(3/2)*Zint_n2/Zint_no*exp((Dn2-Dno)/K/T);
        knoo_o2n = ko2n_noo*(mo2*mn/mno/mo)^(3/2)*Zint_o2/Zint_no*exp((Do2-Dno)/K/T);
    case 2
%% International Workship
% [A_c_diss_n2, A_c_diss_o2, A_c_diss_no, A_c_diss_n, A_c_diss_o]

        An2_diss = [2.5e13/NA 2.5e13/NA 2.5e13/NA 2.5e13/NA 2.5e13/NA];
        Ao2_diss = [9.1e12/NA 9.1e12/NA 9.1e12/NA 9.1e12/NA 9.1e12/NA]; 
        Ano_diss = [4.1e12/NA 4.1e12/NA 4.1e12/NA 4.1e12/NA 4.1e12/NA];
        An2o_non = 7.4e5/NA;
        Anoo_o2n = 3e5/NA;
        
% [N_c_diss_n2, N_c_diss_o2, N_c_diss_no, N_c_diss_n, N_c_diss_o]    

        Nn2_diss = [-1 -1 -1 -1 -1];
        No2_diss = [-1 -1 -1 -1 -1];
        Nno_diss = [-1 -1 -1 -1 -1];
        Nn2o_non = 0.5;
        Nnoo_o2n = 0.5;
        
        En2_diss = [Dn2 Dn2 Dn2 Dn2 Dn2];
        Eo2_diss = [Do2 Do2 Do2 Do2 Do2];
        Eno_diss = [Dno Dno Dno Dno Dno];
        En2o_non = 37940*K;
        Enoo_o2n = 19460*K;
        
% Forward

        kn2_diss = An2_diss.*T.^Nn2_diss.*exp(-En2_diss./K./T);
        ko2_diss = Ao2_diss.*T.^No2_diss.*exp(-Eo2_diss./K./T);
        kno_diss = Ano_diss.*T.^Nno_diss.*exp(-Eno_diss./K./T);
        kn2o_non = An2o_non*T^Nn2o_non*exp(-En2o_non/K/T);
        knoo_o2n = Anoo_o2n*T^Nnoo_o2n*exp(-Enoo_o2n/K/T);
        
% Reverse

        Zint_n2 = sum(exp(-e_n2./T))*T/2/th_n2_rot; 
        Zint_o2 = sum(exp(-e_o2./T))*T/2/th_o2_rot;    
        Zint_no = sum(exp(-e_no./T))*T/th_no_rot;

        kn2_rec = kn2_diss.*(mn2/mn/mn)^(3/2)*H^3*(2*pi*K*T)^(-3/2)*Zint_n2*exp(Dn2/K/T);
        ko2_rec = ko2_diss.*(mo2/mo/mo)^(3/2)*H^3*(2*pi*K*T)^(-3/2)*Zint_o2*exp(Do2/K/T);
        kno_rec = kno_diss.*(mno/mn/mo)^(3/2)*H^3*(2*pi*K*T)^(-3/2)*Zint_no*exp(Dno/K/T);
        knon_n2o = kn2o_non*(mn2*mo/mno/mn)^(3/2)*Zint_n2/Zint_no*exp((Dn2-Dno)/K/T);
        ko2n_noo = knoo_o2n*(mo2*mn/mno/mo)^(-3/2)*Zint_no/Zint_o2*exp((Dno-Do2)/K/T);

    case 3
%% Капителли
        ko2n_noo = 3.2e-18.*(T./300).*exp(-3150./T);
        kn2o_non = 3e-16.*exp(-38370./T);
        knon_n2o = 1.8e-17*exp(T./300)^(1/2);
        knoo_o2n = 7.5e-18*exp(T./300)*exp(-19500./T);
        kn2diss = 5.4e-14.*(1-exp(-3354./T)).*exp(-113200./T);
        kn2_diss = [kn2diss kn2diss kn2diss 6.6*kn2diss 6.6*kn2diss];
        ko2diss = 6.1e-15.*(1-exp(-2240./T)).*exp(-59380./T);
        ko2_diss = [ko2diss 5.9*ko2diss ko2diss ko2diss 21*ko2diss];
        knodiss = 8.7e-15.*exp(-75994./T);
        kno_diss = [knodiss knodiss 20*knodiss 20*knodiss 20*knodiss];
        kn2rec = 1.8e-45*exp(435/T);
        kn2_rec = [8.3e-44*exp(500/T) kn2rec kn2rec 3*kn2rec 3*kn2rec];
        ko2rec = 4e-45.*(300./T).^0.41;
        ko2_rec = [1e-45*(300./T).^0.41 ko2rec 0.17*ko2rec 0.8*ko2rec 3.6*ko2rec];
        kno_rec = [1e-44.*(300./T).^(1/2) 1e-44.*(300./T).^(1/2) 1.8e-43.*(300./T) 1.8e-43.*(300./T) 1.8e-43.*(300./T)];
end

coefficient = [kn2o_non knon_n2o ko2n_noo knoo_o2n kn2_diss kn2_rec ko2_diss ko2_rec kno_diss kno_rec];

end
