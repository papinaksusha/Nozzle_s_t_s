%% TEST SSH 

global SW_O WX

SW = SW_O;

SW_O = 2;
wx = WX;
WX = [1432 1198 1407.5];

N_T = 100;
TT = 2000 : N_T : 14000;

K_ssh_VT_N2  = zeros(I(SW_O,1), 5, length(TT));
K_ssh_VV_N2_N2 = zeros(I(SW_O,1), I(SW_O,1), length(TT));
K_ssh_VV_N2_O2 = zeros(I(SW_O,1), I(SW_O,2), length(TT));
K_ssh_VV_N2_NO = zeros(I(SW_O,1), I(SW_O,3), length(TT));

K_ssh_VT_O2  = zeros(I(SW_O,2), 5, length(TT));
K_ssh_VV_O2_N2 = zeros(I(SW_O,2), I(SW_O,1), length(TT));
K_ssh_VV_O2_O2 = zeros(I(SW_O,2), I(SW_O,2), length(TT));
K_ssh_VV_O2_NO = zeros(I(SW_O,2), I(SW_O,3), length(TT));


for  i = 1 : length(TT)
   [K_ssh_VT_N2(:,:,i), K_ssh_VV_N2_N2(:,:,i), K_ssh_VV_N2_O2(:,:,i), K_ssh_VV_N2_NO(:,:,i)] = k_ssh(1,TT(i));
   [K_ssh_VT_O2(:,:,i), K_ssh_VV_O2_N2(:,:,i), K_ssh_VV_O2_O2(:,:,i), K_ssh_VV_O2_NO(:,:,i)] = k_ssh(2,TT(i));
end

% comparing with Olya's dissertation
%% VT

f = squeeze(K_ssh_VT_N2(1,:,:));

% 1
figure(fig)
plot(TT, log10(f));
title('VT N2 1->0');
legend('N2','O2','NO','N','O')
xlabel('T');
ylabel('lg(k_{N2,1->0}^{M})')
fig = fig + 1;

% 2
figure(fig)
plot(TT, log10(squeeze(K_ssh_VT_O2(1,:,:))));
title('VT O2 1->0');
legend('N2','O2','NO','N','O')
xlabel('T');
ylabel('lg(k_{O2,1->0}^{M})')
fig = fig + 1;

% 3
figure(fig)
plot(TT, log10(squeeze(K_ssh_VT_N2(21,:,:))));
title('VT N2 21->20');
legend('N2','O2','NO','N','O')
xlabel('T');
ylabel('lg(k_{N2,21->20}^{M})')
fig = fig + 1;

% 4
figure(fig)
plot(TT, log10(squeeze(K_ssh_VT_O2(21,:,:))));
title('VT O2 21->20');
legend('N2','O2','NO','N','O')
xlabel('T');
ylabel('lg(k_{O2,21->20}^{M})')
fig = fig + 1;

% 5
figure(fig)
plot(TT, log10(squeeze(K_ssh_VT_N2(47,:,:))));
title('VT N2 47->46');
legend('N2','O2','NO','N','O')
xlabel('T');
ylabel('lg(k_{N2,47->46}^{M})')
fig = fig + 1;

% 6
figure(fig)
plot(TT, log10(squeeze(K_ssh_VT_O2(36,:,:))));
title('VT O2 36->35');
legend('N2','O2','NO','N','O')
xlabel('T');
ylabel('lg(k_{O2,36->35}^{M})')
fig = fig + 1;

% 7
figure(fig)
hold on
subplot(2,2,1) 
plot(0 : I(SW_O,1) - 1, log10(squeeze(K_ssh_VT_N2(:,:,(5000-TT(1))/N_T+1))))
title('VT N2 T = 5000 K');
xlabel('i');
ylabel('lg(k_{N_2,i->i-1}^{M})');
xlim([0,I(SW_O,1) - 1]);
subplot(2,2,2) 
plot(0 : I(SW_O,2) - 1, log10(squeeze(K_ssh_VT_O2(:,:,(5000-TT(1))/N_T+1))))
title('VT O2 T = 5000 K');
xlabel('i');
ylabel('lg(k_{O_2,i->i-1}^{M})');
xlim([0,I(SW_O,2) - 1]);
subplot(2,2,3) 
plot(0 : I(SW_O,1) - 1, log10(squeeze(K_ssh_VT_N2(:,:,(8000-TT(1))/N_T+1))))
title('VT N2 T = 8000 K');
xlabel('i');
xlim([0,I(SW_O,1) - 1]);
ylabel('lg(k_{N_2,i->i-1}^{M})');
subplot(2,2,4) 
plot(0 : I(SW_O,2) - 1, log10(squeeze(K_ssh_VT_O2(:,:,(8000-TT(1))/N_T+1))))
title('VT O2 T = 8000 K');
xlim([0,I(SW_O,2) - 1]);
xlabel('i');
ylabel('lg(k_{O_2,i->i-1}^{M})');
hold off
fig = fig + 1;

%% VV

% 8
figure(fig)
hold on
plot(TT, log10(squeeze(K_ssh_VV_N2_N2(1,1,:))))
plot(TT, log10(squeeze(K_ssh_VV_N2_N2(11,1,:))))
plot(TT, log10(squeeze(K_ssh_VV_N2_N2(21,1,:))))
plot(TT, log10(squeeze(K_ssh_VV_N2_N2(31,1,:))))
title('VV N2-N2');
legend('i = 1','i = 11','i = 21','i = 31')
xlabel('T');
ylabel('lg(k_{N_2,i->i-1}^{N_2,0->1})');
hold off
fig = fig +1;

% 9
figure(fig)
hold on
plot(TT, log10(squeeze(K_ssh_VV_O2_O2(1,36,:))))
plot(TT, log10(squeeze(K_ssh_VV_O2_O2(11,36,:))))
plot(TT, log10(squeeze(K_ssh_VV_O2_O2(21,36,:))))
plot(TT, log10(squeeze(K_ssh_VV_O2_O2(31,36,:))))
title('VV O2-O2');
legend('i = 1','i = 11','i = 21','i = 31')
xlabel('T');
ylabel('lg(k_{O_2,i->i-1}^{O_2,35->36})');
hold off
fig = fig + 1;

% 10
figure(fig)
hold on
plot(0 : I(SW_O,1) - 1, log10(squeeze(K_ssh_VV_N2_N2(:,1,(5000-TT(1))/N_T+1))))
plot(0 : I(SW_O,1) - 1, log10(squeeze(K_ssh_VV_N2_N2(:,1,(8000-TT(1))/N_T+1))))
plot(0 : I(SW_O,1) - 1, log10(squeeze(K_ssh_VV_N2_N2(:,1,(14000-TT(1))/N_T+1))))
title('VV N2-N2');
legend('T = 5000','T = 8000','T = 14000')
xlabel('i');
xlim([0, I(SW_O,1)]);
ylabel('lg(k_{N_2,i->i-1}^{N_2,0->1})');
hold off
fig = fig + 1;

% 11
figure(fig)
hold on
plot(0 : I(SW_O,2) - 1, log10(squeeze(K_ssh_VV_O2_O2(:,36,(5000-TT(1))/N_T+1))))
plot(0 : I(SW_O,2) - 1, log10(squeeze(K_ssh_VV_O2_O2(:,36,(8000-TT(1))/N_T+1))))
plot(0 : I(SW_O,2) - 1, log10(squeeze(K_ssh_VV_O2_O2(:,36,(14000-TT(1))/N_T+1))))
title('VV O2-O2');
legend('T = 5000','T = 8000','T = 14000')
xlabel('i');
xlim([0, I(SW_O,2)]);
ylabel('lg(k_{O_2,i->i-1}^{O_2,35->36})');
hold off
fig = fig + 1;

% VV'

% 12
figure(fig)
hold on
plot(TT, log10(squeeze(K_ssh_VV_N2_O2(1,1,:))))
plot(TT, log10(squeeze(K_ssh_VV_N2_O2(11,1,:))))
plot(TT, log10(squeeze(K_ssh_VV_N2_O2(21,1,:))))
plot(TT, log10(squeeze(K_ssh_VV_N2_O2(31,1,:))))
title('VV" N2-O2');
legend('i = 1','i = 11','i = 21','i = 31')
xlabel('T');
ylabel('lg(k_{N_2,i->i-1}^{O_2,0->1})');
hold off
fig = fig + 1;

% 13
figure(fig)
hold on
plot(TT, log10(squeeze(K_ssh_VV_O2_N2(1,1,:))))
plot(TT, log10(squeeze(K_ssh_VV_O2_N2(11,1,:))))
plot(TT, log10(squeeze(K_ssh_VV_O2_N2(21,1,:))))
plot(TT, log10(squeeze(K_ssh_VV_O2_N2(31,1,:))))
title('VV" O2-N2');
legend('i = 1','i = 11','i = 21','i = 31')
xlabel('T');
ylabel('lg(k_{O_2,i->i-1}^{N_2,35->36})');
hold off
fig = fig + 1;

% 14
figure(fig)
hold on
plot(0 : I(SW_O,1) - 1, log10(squeeze(K_ssh_VV_N2_O2(:,1,(5000-TT(1))/N_T+1))))
plot(0 : I(SW_O,1) - 1, log10(squeeze(K_ssh_VV_N2_O2(:,1,(8000-TT(1))/N_T+1))))
plot(0 : I(SW_O,1) - 1, log10(squeeze(K_ssh_VV_N2_O2(:,1,(14000-TT(1))/N_T+1))))
title('VV N2-O2');
legend('T = 5000','T = 8000','T = 14000')
xlabel('i');
xlim([0, I(SW_O,1)]);
ylabel('lg(k_{N_2,i->i-1}^{O_2,0->1})');
hold off
fig = fig + 1;

% 15
figure(fig)
hold on
plot(0 : I(SW_O,2) - 1, log10(squeeze(K_ssh_VV_O2_N2(:,1,(5000-TT(1))/N_T+1))))
plot(0 : I(SW_O,2) - 1, log10(squeeze(K_ssh_VV_O2_N2(:,1,(8000-TT(1))/N_T+1))))
plot(0 : I(SW_O,2) - 1, log10(squeeze(K_ssh_VV_O2_N2(:,1,(14000-TT(1))/N_T+1))))
title('VV" O2-N2');
legend('T = 5000','T = 8000','T = 14000')
xlabel('i');
xlim([0, I(SW_O,2)]);
ylabel('lg(k_{O_2,i->i-1}^{N_2,35->36})');
hold off
fig = fig + 1;
%%
% 16
figure(fig)
semilogy(TT, squeeze(K_ssh_VT_N2(9,1,:))), hold on
semilogy(TT, squeeze(K_ssh_VT_N2(9,4,:)))
semilogy(TT, squeeze(K_ssh_VV_N2_N2(9,3,:)))
legend('VT whith molecule','VT with atom','VV');
hold off
fig = fig + 1;

%%
% comparing with Ildar's dissertation


% 17
figure(fig)
semilogy(0 : I(SW_O,1) - 1, squeeze(K_ssh_VT_N2(:,4,(2500 - TT(1))/N_T+1)))
xlabel('i');
xlim([0, I(SW_O,1)]);
title('VT N_2 + N, T = 2500 K')
grid on
fig = fig + 1;

% 18
figure(fig)
semilogy(0 : I(SW_O,2) - 1, squeeze(K_ssh_VT_O2(:,5,(2500 - TT(1))/N_T+1)))
xlabel('i');
xlim([0, I(SW_O,2)]);
title('VT O_2 + O, T = 2500 K')
grid on
fig = fig + 1;

% 19
f = figure(fig);
semilogy(0 : I(SW_O,1) - 1, squeeze(K_ssh_VT_N2(:,4,(5000 - TT(1))/N_T+1)))
xlabel('i');
xlim([0, I(SW_O,1)]);
title('VT N_2 + N, T = 5000 K')
grid on;
fig = fig + 1;

% 20
figure(fig)
semilogy(0 : I(SW_O,2) - 1, squeeze(K_ssh_VT_O2(:,5,(5000 - TT(1))/N_T+1)))
xlabel('i');
xlim([0, I(SW_O,2)]);
title('VT O_2 + O, T = 5000 K')
grid on
fig = fig + 1;

% 21
f = figure(fig);
semilogy(0 : I(SW_O,1) - 1, squeeze(K_ssh_VT_N2(:,4,(10000 - TT(1))/N_T+1)))
xlabel('i');
xlim([0, I(SW_O,1)]);
title('VT N_2 + N, T = 10000 K')
grid on;
fig = fig + 1;

% 22
figure(fig)
semilogy(0 : I(SW_O,2) - 1, squeeze(K_ssh_VT_O2(:,5,(10000 - TT(1))/N_T+1)))
xlabel('i');
xlim([0, I(SW_O,2)]);
title('VT O_2 + O, T = 10000 K')
grid on
fig = fig + 1;

% 23
figure(fig)
hold on
plot(0 : I(SW_O,2) - 1, squeeze(K_ssh_VV_O2_N2(:,10,(2500-TT(1))/N_T+1)))
title('VV" O2-N2, T = 2500 K');
xlabel('i');
xlim([0, I(SW_O,2)]);
ylabel('k_{O_2,i->i-1}^{N_2,9->10}');
hold off
fig = fig + 1;

% 24
figure(fig)
hold on
plot(0 : I(SW_O,2) - 1, squeeze(K_ssh_VV_O2_N2(:,10,(5000-TT(1))/N_T+1)))
title('VV" O2-N2, T = 5000 K');
xlabel('i');
xlim([0, I(SW_O,2)]);
ylabel('k_{O_2,i->i-1}^{N_2,9->10}');
hold off
fig = fig + 1;

% 25
figure(fig)
hold on
plot(0 : I(SW_O,2) - 1, squeeze(K_ssh_VV_O2_N2(:,10,(10000-TT(1))/N_T+1)))
title('VV" O2-N2, T = 10000 K');
xlabel('i');
xlim([0, I(SW_O,2)]);
ylabel('k_{O_2,i->i-1}^{N_2,9->10}');
hold off

SW_O = SW;
WX = wx;
