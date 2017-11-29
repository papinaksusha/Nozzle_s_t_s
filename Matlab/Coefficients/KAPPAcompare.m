global SW_O WX W H C K

TT = [2500, 5000, 10000, 20000];

k_VT_N2  = zeros(I(SW_O,1), 5, length(TT));
k_VT_O2  = zeros(I(SW_O,2), 5, length(TT));

for  i = 1 : length(TT)
   [k_VT_N2(:,:,i), ~, ~, ~] = k_ssh(1,TT(i));
   [k_VT_O2(:,:,i), ~, ~, ~] = k_ssh(2,TT(i));
end

l_N2 = I(SW_O , 1) + 1;
l_O2 = I(SW_O , 2) + 1;
i_N2 = 0 : l_N2 - 2;
i_O2 = 0 : l_O2 - 2;

for j = 1 : length(TT)
    
    k_VT_N2_r(:,:,j) = k_VT_N2(:,:,j).*exp(- H*C/K.*(W(1) - 2*WX(1) - 2*WX(1).*i_N2)'*ones(1 , 5)./TT(j)); 
    k_VT_O2_r(:,:,j) = k_VT_O2(:,:,j).*exp(- H*C/K.*(W(2) - 2*WX(2) - 2*WX(2).*i_O2)'*ones(1 , 5)./TT(j));
    
end

% N2

k_VT_N2_N_1_0 = squeeze(k_VT_N2(1,4,:))

k_VT_N2_N_5_4 = squeeze(k_VT_N2(5,4,:))

k_VT_N2_N_r_0_1 = squeeze(k_VT_N2_r(1,4,:))

k_VT_N2_N_r_4_5 = squeeze(k_VT_N2_r(5,4,:))

% O2

k_VT_O2_O2_1_0 = squeeze(k_VT_O2(1,2,:))

k_VT_O2_O2_5_4 = squeeze(k_VT_O2(5,2,:))

k_VT_O2_O2_r_0_1 = squeeze(k_VT_O2_r(1,2,:))

k_VT_O2_O2_r_4_5 = squeeze(k_VT_O2_r(5,2,:))
