%%% VV probability

global SW_O WX H C K  W

SW = SW_O;

SW_O = 2;
wx = WX;
WX = [1432 1198 1407.5];

TT = [500, 1500, 2000, 7000];

K_ssh_VV_N2_N2 = zeros(I(SW_O,1), I(SW_O,1), length(TT));
K_ssh_VV_N2_O2 = zeros(I(SW_O,1), I(SW_O,2), length(TT));

K_ssh_VV_O2_N2 = zeros(I(SW_O,2), I(SW_O,1), length(TT));
K_ssh_VV_O2_O2 = zeros(I(SW_O,2), I(SW_O,2), length(TT));


for  i = 1 : length(TT)
   [~, K_ssh_VV_N2_N2(:,:,i), K_ssh_VV_N2_O2(:,:,i), ~] = k_ssh(1,TT(i));
   [~, K_ssh_VV_O2_N2(:,:,i), K_ssh_VV_O2_O2(:,:,i), ~] = k_ssh(2,TT(i));
end

for j = 1 : length(TT)
    figure(fig)
    image(K_ssh_VV_N2_N2(:, :, j),'CDataMapping','scaled')
    title(strcat('N_2 + N_2, T = ', num2str(TT(j)), ' K'))
    colorbar
    fig = fig + 1;
    
end

for j = 1 : length(TT)
    figure(fig)
    image(K_ssh_VV_N2_O2(:, :, j),'CDataMapping','scaled')
    title(strcat('N_2 + O_2, T = ', num2str(TT(j)), ' K'))
    colorbar
    fig = fig + 1;
    
end

for j = 1 : length(TT)
    figure(fig)
    image(K_ssh_VV_O2_N2(:, :, j),'CDataMapping','scaled')
    title(strcat('O_2 + N_2, T = ', num2str(TT(j)), ' K'))
    colorbar
    fig = fig + 1;
    
end

for j = 1 : length(TT)
    figure(fig)
    image(K_ssh_VV_O2_O2(:, :, j),'CDataMapping','scaled')
    title(strcat('O_2 + O_2, T = ', num2str(TT(j)), ' K'))
    colorbar
    fig = fig + 1;
    
end

% reverse coefficients

l_N2 = I(SW_O , 1) + 1;
l_O2 = I(SW_O , 2) + 1;
i_N2 = 0 : l_N2 - 2;
i_O2 = 0 : l_O2 - 2;

for j = 1 : length(TT)
    
    K_ssh_VV_N2_N2_r(:,:,j) = K_ssh_VV_N2_N2(:,:,j).*exp(-2*H*C*WX(1)/K.*(i_N2'*ones(1 , l_N2 - 1) - ... 
                                   ones(l_N2 - 1, 1)*i_N2)./TT(j));
    K_ssh_VV_O2_O2_r(:,:,j) = K_ssh_VV_O2_O2(:,:,j).*exp(-2*H*C*WX(2)/K.*(i_O2'*ones(1 , l_O2 - 1) - ...
                                   ones(l_O2 - 1, 1)*i_O2)./TT(j));
    K_ssh_VV_N2_O2_r(:,:,j) =  K_ssh_VV_N2_O2(:,:,j).*exp(-H*C/K.*((W(1) - WX(1).*(i_N2 + 2))'*ones(1 , l_O2 - 1) - ...
                                    ones(l_N2 - 1 , 1)*(W(2) - WX(2).*(i_O2 + 2)))./TT(j));
    K_ssh_VV_O2_N2_r(:,:,j) =  K_ssh_VV_O2_N2(:,:,j).*exp(-H*C/K.*((W(2) - WX(2).*(i_O2 + 2))'*ones(1 , l_N2 - 1) + ...
                                    ones(l_O2 - 1 , 1)*(W(1) - WX(1).*(i_N2 + 2)))./TT(j));   
                               
end


for j = 1 : length(TT)
    figure(fig)
    image(K_ssh_VV_N2_N2_r(:, :, j),'CDataMapping','scaled')
    title(strcat('N_2 + N_2, reverse, T = ', num2str(TT(j)), ' K'))
    colorbar
    fig = fig + 1;
    
end

for j = 1 : length(TT)
    figure(fig)
    image(K_ssh_VV_N2_O2_r(:, :, j),'CDataMapping','scaled')
    title(strcat('N_2 + O_2, reverse, T = ', num2str(TT(j)), ' K'))
    colorbar
    fig = fig + 1;
    
end

for j = 1 : length(TT)
    figure(fig)
    image(K_ssh_VV_O2_N2_r(:, :, j),'CDataMapping','scaled')
    title(strcat('O_2 + N_2, reverse , T = ', num2str(TT(j)), ' K'))
    colorbar
    fig = fig + 1;
    
end

for j = 1 : length(TT)
    figure(fig)
    image(K_ssh_VV_O2_O2_r(:, :, j),'CDataMapping','scaled')
    title(strcat('O_2 + O_2, reverse T = ', num2str(TT(j)), ' K'))
    colorbar
    fig = fig + 1;
    
end
 


SW_O = SW;
WX = wx;


