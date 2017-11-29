function Relax = React1T(nn2,no2,nno,nn,no,T,e_n2,e_o2,e_no,sw)
% Релаксационные члены Rn2, Ro2, Rno, Rn, Ro


k = Arrenius_coeff(T, sw, e_n2, e_o2, e_no);

Relax(1) = nno*nn*k(2) - nn2*no*k(1);  % Rn2 2-2
Relax(2) = nno*no*k(4) - no2*nn*k(3); % Ro2 2-2
Relax(3) = nn2*no*k(1) - nno*nn*k(2) + no2*nn*k(3) - nno*no*k(4); % Rno 2-2
Relax(4) = nn2*(nn^2*k(10)-nn2*k(5))+no2*(nn^2*k(11)-nn2*k(6))+nno*(nn^2*k(12)-nn2*k(7))+nn*(nn^2*k(13)-nn2*k(8))+no*(nn^2*k(14)-nn2*k(9)); %Rn2 2-3
Relax(5) = nn2*(no^2*k(20)-no2*k(15))+no2*(no^2*k(21)-no2*k(16))+nno*(no^2*k(22)-no2*k(17))+nn*(no^2*k(23)-no2*k(18))+no*(no^2*k(24)-no2*k(19)); %Ro2 2-3
Relax(6) = nn2*(nn*no*k(30)-nno*k(25))+no2*(nn*no*k(31)-nno*k(26))+nno*(nn*no*k(32)-nno*k(27))+nn*(nn*no*k(33)-nno*k(28))+no*(nn*no*k(34)-nno*k(29)); %Rno 2-3
Relax(7) = -Relax(1) + Relax(2); % Rn 2-2
Relax(8) = -2*Relax(4)- Relax(6); % Rn 2-3
Relax(9) = -Relax(2) + Relax(1); % Ro 2-2
Relax(10) = -2*Relax(5) - Relax(6); %Ro 2-3


end