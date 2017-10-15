function v_cr = v_critical(n,T)

global M K NA W WX C H SW_O MOLAR I

i_N2 = 0 : I(SW_O, 1);
i_O2 = 0 : I(SW_O, 2);
i_NO = 0 : I(SW_O, 3);

e_i_N2 = H*C.*(W(1).*i_N2 - WX(1).*i_N2 - WX(1).*i_N2.^2);
e_i_O2 = H*C.*(W(2).*i_O2 - WX(2).*i_O2 - WX(2).*i_O2.^2);
e_i_NO = H*C.*(W(3).*i_NO - WX(3).*i_NO - WX(3).*i_NO.^2);

M_inverse = sum(M.*n./MOLAR);
R_bar = K*NA*M_inverse;

Ctr = 1.5*R_bar;
Crot = K*sum(n(1:3));


e_0_N2 = H*C*(W(1)/2 - WX(1)/4);
e_0_O2 = H*C*(W(2)/2 - WX(2)/4);
e_0_NO = H*C*(W(3)/2 - WX(3)/4);

Z_vibr_N2 = sum(exp(-e_i_N2./K./T));
Z_vibr_O2 = sum(exp(-e_i_O2./K./T));
Z_vibr_NO = sum(exp(-e_i_NO./K./T));
E_vibr_N2 = sum((e_i_N2 + e_0_N2).*exp(-e_i_N2./K./T))/Z_vibr_N2;
E_vibr_O2 = sum((e_i_O2 + e_0_O2).*exp(-e_i_O2./K./T))/Z_vibr_O2;
E_vibr_NO = sum((e_i_NO + e_0_NO).*exp(-e_i_NO./K./T))/Z_vibr_NO;
C_vibr_N2 = T^(-2)/K*(sum((e_i_N2 + e_0_N2).*e_i_N2.*exp(-e_i_N2./K./T))/Z_vibr_N2 - ...
                 E_vibr_N2*sum(e_i_N2.*exp(-e_i_N2./K./T))/Z_vibr_N2);
C_vibr_O2 = T^(-2)/K*(sum((e_i_O2 + e_0_O2).*e_i_O2.*exp(-e_i_O2./K./T))/Z_vibr_O2 - ...
                 E_vibr_O2*sum(e_i_O2.*exp(-e_i_O2./K./T))/Z_vibr_O2);
C_vibr_NO = T^(-2)/K*(sum((e_i_NO + e_0_NO).*e_i_NO.*exp(-e_i_NO./K./T))/Z_vibr_NO - ...
                 E_vibr_NO*sum(e_i_NO.*exp(-e_i_NO./K./T))/Z_vibr_NO);
     
Cvibr = sum([C_vibr_N2 C_vibr_O2 C_vibr_NO].*n(1:3));
Cv = Ctr + Crot + Cvibr;

Cp = R_bar + Cv;

kappa = Cp/Cv;

v_cr = sqrt(kappa*R_bar/sum(M.*n)*T);

end
