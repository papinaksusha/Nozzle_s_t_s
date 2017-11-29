function [f1 , f2] = Nozzle_5_NO_full(x , init , options , T_cr , p_cr , v_cr)

global K H C M THETA_R D SW_N SW_O

n_cr = p_cr/K/T_cr;

[f1,f2] = ode15s(@nozzle_5_NO , x , init , options);

    function dy = nozzle_5_NO(x , y)

	l_N2 = I(SW_O , 1) + 1;
	l_O2 = I(SW_O , 2) + 1;
	l_NO = I(SW_O , 3) + 1;
	l_mol =  l_N2 + l_O2 + l_NO;
	l_c = l_mol + 2;
	l_v = l_c + 1;
	l_T = l_v + 1;

    
    A = zeros(l_T , l_T); 
    b = zeros(l_T , 1);
    S = zeros(2);                                                          % S(1) - nozzle cross-section, S(2) - its derivative 
    
    % nozzle geometry
    
    switch SW_N
        case 1
             r_cr = 1e-3;
             alpha = 0.117*pi;
             S(1) = (1 + x*tan(alpha))^2;
             S(2) = 2*tan(alpha)*(1 + x*tan(alpha));
    end
    
    % dimensionless variables
    
    n_N2 = y(1 : l_N2);
    n_O2 = y(l_N2 + 1 : l_N2 + l_O2);
    n_NO = y(l_N2 + l_O2 + 1 : l_mol);
    n_N = y(l_mol + 1);
    n_O = y(l_c);
    v = y(l_v);
    T = y(l_T);
    M_cr = sum(M); 
    M = M./M_cr;
    
    n = [n_N2' n_O2' n_NO' n_N n_O];
    rho = sum([M(1).*n_N2' M(2).*n_O2' M(3)*n_NO' M(4)*n_N M(5)*n_O]);
    
    % vibrational energy of molecules

    i_N2 = 0 : l_N2 - 1;
	i_O2 = 0 : l_O2 - 1;
	i_NO = 0 : l_NO - 1;
    
	e_i_N2 = H*C./K./T_cr.*(W(1).*i_N2 - WX(1).*i_N2 - WX(1).*i_N2.^2);
	e_i_O2 = H*C./K./T_cr.*(W(2).*i_O2 - WX(2).*i_O2 - WX(2).*i_O2.^2);
	e_i_NO = H*C./K./T_cr.*(W(3).*i_NO - WX(3).*i_NO - WX(3).*i_NO.^2);

	 e_i = [e_i_N2 e_i_O2 e_i_NO 0 0];
    
     % vibrational energy of 0th levels
     
    e_0_N2 = H*C/K/T_cr*(W(1)/2 - WX(1)/4);
    e_0_O2 = H*C/K/T_cr*(W(2)/2 - WX(2)/4);
    e_0_NO = H*C/K/T_cr*(W(3)/2 - WX(3)/4);

    e_0 = [e_0_N2.*ones(1 , l_N2) e_0_O2.*ones(1 , l_O2) e_0_NO.*ones(1 , l_NO) 0 0];
     
     % formation energy
     
    e_NO = (D(1)/2 + D(2)/2 - D(3))/K/T_cr;
    e_N = D(1)/2/K/T_cr;
    e_O = D(2)/2/K/T_cr;
     
    e_f = [zeros(1 , l_N2 + l_O2) e_NO.*ones(1 , l_NO) e_N e_O];
     
    %% Left part, matrix A (A(y)*y = b)
   
    % kinetic equations
    
    A(1 : l_c , 1 : l_c) = diag(v.*ones(1 , l_c)); 
    A(1 : l_c , l_v) = n; 
    
    % momentum equation
    
    A(l_v, 1 : l_c) = T; 
    A(l_v, l_v) = M_cr*v_cr^2/K/T_cr * v *rho;
    A(l_v, l_T) = sum(n);
    
    % energy equation
    
    T_energy1 = [2.5*T.*ones(1 , l_mol) 1.5*T 1.5*T];
    T_energy2 = [3.5*T.*ones(1 , l_mol) 2.5*T 2.5*T];
    
    A(l_T, 1 : l_c) = T_energy1 + e_i + e_0 + e_f; 
    A(l_T, l_v) = 1/v*sum(n.*(T_energy2  + e_i + e_0 + e_f));
    A(l_T, l_T) = 2.5*(sum(n_N2) + sum(n_O2) + sum(n_NO)) + 1.5*(n_N +n_O);
   
    AA = sparse(A);
    
    %% Right part, vector b
    
    %%% reaction rate coefficients and source terms
        
    % dimensional variables
    
    n_N2_d = [0 n_N2'.*n_cr 0]'*ones(1 , 5); %50*5
    n_O2_d = [0 n_O2'.*n_cr 0]'*ones(1 , 5);
    n_NO_d = [0 n_NO'.*n_cr 0]'*ones(1 , 5);
    n_N_d = n_N*n_cr;
    n_O_d = n_O*n_cr;
    T_d = T*T_cr;
    
    %%% vibrational energy transitions (SSH)
    
    [k_N2_VT, k_N2_N2_VV, k_N2_O2_VV, k_N2_NO_VV] = k_ssh(1 , T_d);
    [k_O2_VT, k_O2_N2_VV, k_O2_O2_VV, k_O2_NO_VV] = k_ssh(2 , T_d);
    [k_NO_VT, k_NO_N2_VV, k_NO_O2_VV, k_NO_NO_VV] = k_ssh(3 , T_d);
        
    % VT

    i_N2 = 0 : l_N2 - 2;
    i_O2 = 0 : l_O2 - 2;
    i_NO = 0 : l_NO - 2;

    k_N2_VT_r = k_N2_VT.*exp(- H*C/K/T_cr.*(W(1) - 2*WX(1) - 2*WX(1).*i_N2)'*ones(1 , 5)./T); 
    k_O2_VT_r = k_O2_VT.*exp(- H*C/K/T_cr.*(W(2) - 2*WX(2) - 2*WX(2).*i_O2)'*ones(1 , 5)./T);
    k_NO_VT_r = k_NO_VT.*exp(- H*C/K/T_cr.*(W(3) - 2*WX(3) - 2*WX(3).*i_NO)'*ones(1 , 5)./T);
    
       
    k_N2_VT = [zeros(1 , 5); k_N2_VT; zeros(1 , 5)]; %49*5
    k_O2_VT = [zeros(1 , 5); k_O2_VT; zeros(1 , 5)];
    k_NO_VT = [zeros(1 , 5); k_NO_VT; zeros(1 , 5)];
    k_N2_VT_r = [zeros(1 , 5); k_N2_VT_r; zeros(1 , 5)]; %49*5
    k_O2_VT_r = [zeros(1 , 5); k_O2_VT_r; zeros(1 , 5)];
    k_NO_VT_r = [zeros(1 , 5); k_NO_VT_r; zeros(1 , 5)];
    
    n_c_N2 = ones(l_N2 , 1)*[sum(n_N2) sum(n_O2) sum(n_NO) n_N n_O].*n_cr; %  %*  48*5
    
    R_N2_VT = sum(n_c_N2.*(n_N2_d(1 : end - 2 , :).*k_N2_VT_r(1 : end - 1 , :) - ...
                           n_N2_d(2 : end - 1 , :).*k_N2_VT(1 : end - 1 , :) + ...
                           n_N2_d(3 : end , :).*k_N2_VT(2 : end , :) - ...
                           n_N2_d(2 : end - 1 , :).*k_N2_VT_r(2 : end , :)) , 2); %48*1
                     
    n_c_O2 = ones(l_O2 , 1)*[sum(n_N2) sum(n_O2) sum(n_NO) n_N n_O].*n_cr; 
    
    R_O2_VT = sum(n_c_O2.*(n_O2_d(1 : end - 2 , :).*k_O2_VT_r(1 : end - 1 , :) - ...
                           n_O2_d(2 : end - 1 , :).*k_O2_VT(1 : end - 1 , :) + ...
                           n_O2_d(3 : end , :).*k_O2_VT(2 : end , :) - ...
                           n_O2_d(2 : end - 1 , :).*k_O2_VT_r(2 : end , :)) , 2);
    
    n_c_NO = ones(l_NO , 1)*[sum(n_N2) sum(n_O2) sum(n_NO) n_N n_O].*n_cr; 
                       
    R_NO_VT = sum(n_c_NO.*(n_NO_d(1 : end - 2 , :).*k_NO_VT_r(1 : end - 1 , :) - ...
                           n_NO_d(2 : end - 1 , :).*k_NO_VT(1 : end - 1 , :) + ...
                           n_NO_d(3 : end , :).*k_NO_VT(2 : end , :) - ...
                           n_NO_d(2 : end - 1 , :).*k_NO_VT_r(2 : end , :)) , 2);                 
    % VV 
    
    n_N2_d = [0 n_N2'.*n_cr 0]; %1*50
    n_O2_d = [0 n_O2'.*n_cr 0];
    n_NO_d = [0 n_NO'.*n_cr 0];
    
    k_N2_N2_VV_r = k_N2_N2_VV.*exp(-((e_i_N2(2 : end) - e_i_N2(1 : end - 1))'*ones(1 , l_N2 - 1) + ... 
                                   ones(l_N2 - 1 , 1)*(e_i_N2(1 : end - 1) - e_i_N2(2 : end)))./T);
    k_O2_O2_VV_r = k_O2_O2_VV.*exp(-((e_i_O2(2 : end) - e_i_O2(1 : end - 1))'*ones(1 , l_O2 - 1) + ... 
                                   ones(l_O2 - 1 , 1)*(e_i_O2(1 : end - 1) - e_i_O2(2 : end)))./T);
    k_NO_NO_VV_r = k_NO_NO_VV.*exp(-((e_i_NO(2 : end) - e_i_NO(1 : end - 1))'*ones(1 , l_NO - 1) + ... 
                                   ones(l_NO - 1 , 1)*(e_i_NO(1 : end - 1) - e_i_NO(2 : end)))./T);

    k_N2_N2_VV = [zeros(1 , l_N2 - 1) ; k_N2_N2_VV; zeros(1 , l_N2 - 1)]; %49*47
    k_O2_O2_VV = [zeros(1 , l_O2 - 1) ; k_O2_O2_VV; zeros(1 , l_O2 - 1)];
    k_NO_NO_VV = [zeros(1 , l_NO - 1) ; k_NO_NO_VV; zeros(1 , l_NO - 1)];
    k_N2_N2_VV_r = [zeros(1 , l_N2 - 1) ; k_N2_N2_VV_r; zeros(1 , l_N2 - 1)]; %49*47
    k_O2_O2_VV_r = [zeros(1 , l_O2 - 1) ; k_O2_O2_VV_r; zeros(1 , l_O2 - 1)];
    k_NO_NO_VV_r = [zeros(1 , l_NO - 1) ; k_NO_NO_VV_r; zeros(1 , l_NO - 1)];
    
    R_N2_VV = sum(n_N2_d(1 : end - 2)'*n_N2_d(3 : end - 1).*k_N2_N2_VV_r(1 : end - 1, :) - ...
                  n_N2_d(2 : end - 1)'*n_N2_d(2 : end - 2).*k_N2_N2_VV(1 : end - 1, :) + ...
                  n_N2_d(3 : end)'*n_N2_d(2 : end - 2).*k_N2_N2_VV(2 : end, :) - ...
                  n_N2_d(2 : end - 1)'*n_N2_d(3 : end - 1).*k_N2_N2_VV_r(2 : end, :) , 2); % 48*1
    
    R_O2_VV = sum(n_O2_d(1 : end - 2)'*n_O2_d(3 : end - 1).*k_O2_O2_VV_r(1 : end - 1, :) - ...
                  n_O2_d(2 : end - 1)'*n_O2_d(2 : end - 2).*k_O2_O2_VV(1 : end - 1, :) + ...
                  n_O2_d(3 : end)'*n_O2_d(2 : end - 2).*k_O2_O2_VV(2 : end, :) - ...
                  n_O2_d(2 : end - 1)'*n_O2_d(3 : end - 1).*k_O2_O2_VV_r(2 : end, :) , 2);
       
    R_NO_VV = sum(n_NO_d(1 : end - 2)'*n_NO_d(3 : end - 1).*k_NO_NO_VV_r(1 : end - 1, :) - ...
                  n_NO_d(2 : end - 1)'*n_NO_d(2 : end - 2).*k_NO_NO_VV(1 : end - 1, :) + ...
                  n_NO_d(3 : end)'*n_NO_d(2 : end - 2).*k_NO_NO_VV(2 : end, :) - ...
                  n_NO_d(2 : end - 1)'*n_NO_d(3 : end - 1).*k_NO_NO_VV_r(2 : end, :) , 2);
    
    % VV'
              
    k_N2_O2_VV_r =  k_N2_O2_VV.*exp(-((e_i_N2(2 : end) - e_i_N2(1 : end - 1))'*ones(1 , l_O2 - 1) + ...
                                    ones(l_N2 - 1 , 1)*(e_i_O2(1 : end - 1) - e_i_O2(2 : end)))./T);
    k_N2_NO_VV_r =  k_N2_NO_VV.*exp(-((e_i_N2(2 : end) - e_i_N2(1 : end - 1))'*ones(1 , l_NO - 1) + ...
                                    ones(l_N2 - 1, 1)*(e_i_NO(1 : end - 1) - e_i_NO(2 : end)))./T);
    k_O2_N2_VV_r =  k_O2_N2_VV.*exp(-((e_i_O2(2 : end) - e_i_O2(1 : end - 1))'*ones(1 , l_N2 - 1) + ...
                                    ones(l_O2 - 1, 1)*(e_i_N2(1 : end - 1) - e_i_N2(2 : end)))./T); 
    k_O2_NO_VV_r =  k_O2_NO_VV.*exp(-((e_i_O2(2 : end) - e_i_O2(1 : end - 1))'*ones(1 , l_NO - 1) + ...
                                    ones(l_O2 - 1, 1)*(e_i_NO(1 : end - 1) - e_i_NO(2 : end)))./T);                            
    k_NO_N2_VV_r =  k_NO_N2_VV.*exp(-((e_i_NO(2 : end) - e_i_NO(1 : end - 1))'*ones(1 , l_N2 - 1) + ...
                                    ones(l_NO - 1, 1)*(e_i_N2(1 : end - 1) - e_i_N2(2 : end)))./T); 
    k_NO_O2_VV_r =  k_NO_O2_VV.*exp(-((e_i_NO(2 : end) - e_i_NO(1 : end - 1))'*ones(1 , l_O2 - 1) + ...
                                    ones(l_NO - 1, 1)*(e_i_O2(1 : end - 1) - e_i_O2(2 : end)))./T); 
    
    k_N2_O2_VV = [zeros(1 , l_O2 - 1);  k_N2_O2_VV; zeros(1 , l_O2 - 1)];
    k_N2_O2_VV_r = [zeros(1 , l_O2 - 1);  k_N2_O2_VV_r; zeros(1 , l_O2 - 1)];
    k_N2_NO_VV = [zeros(1 , l_NO - 1);  k_N2_NO_VV; zeros(1 , l_NO - 1)];
    k_N2_NO_VV_r = [zeros(1 , l_NO - 1);  k_N2_NO_VV_r; zeros(1 , l_NO - 1)];
    k_O2_N2_VV = [zeros(1 , l_N2 - 1);  k_O2_N2_VV ; zeros(1 , l_N2 - 1)];
    k_O2_N2_VV_r = [zeros(1 , l_N2 - 1);  k_O2_N2_VV_r; zeros(1 , l_N2 - 1)];
    k_O2_NO_VV = [zeros(1 , l_NO - 1);  k_O2_NO_VV ; zeros(1 , l_NO - 1)];
    k_O2_NO_VV_r = [zeros(1 , l_NO - 1);  k_O2_NO_VV_r; zeros(1 , l_NO - 1)];
    k_NO_N2_VV = [zeros(1 , l_N2 - 1);  k_NO_N2_VV; zeros(1 , l_N2 - 1)];
    k_NO_N2_VV_r = [zeros(1 , l_N2 - 1);  k_NO_N2_VV_r; zeros(1 , l_N2 - 1)];
    k_NO_O2_VV = [zeros(1 , l_O2 - 1);  k_NO_O2_VV; zeros(1 , l_O2 - 1)];
    k_NO_O2_VV_r = [zeros(1 , l_O2 - 1);  k_NO_O2_VV_r; zeros(1 , l_O2 - 1)];
    
    R_N2_VV_s = sum(n_N2_d(1 : end - 2)'*n_O2_d(3 : end - 1).*k_N2_O2_VV_r(1 : end - 1, :) - ...
                  n_N2_d(2 : end - 1)'*n_O2_d(2 : end - 2).*k_N2_O2_VV(1 : end - 1, :) + ...
                  n_N2_d(3 : end)'*n_O2_d(2 : end-2).*k_N2_O2_VV(2 : end, :) - ...
                  n_N2_d(2 : end - 1)'*n_O2_d(3 : end - 1).*k_N2_O2_VV_r(2 : end, :) , 2) + ...
                  sum(n_N2_d(1 : end - 2)'*n_NO_d(3 : end - 1).*k_N2_NO_VV_r(1 : end - 1, :) - ...
                  n_N2_d(2 : end - 1)'*n_NO_d(2 : end - 2).*k_N2_NO_VV(1 : end - 1, :) + ...
                  n_N2_d(3 : end)'*n_NO_d(2 : end-2).*k_N2_NO_VV(2 : end, :) - ...
                  n_N2_d(2 : end - 1)'*n_NO_d(3 : end - 1).*k_N2_NO_VV_r(2 : end, :) , 2);
    
    R_O2_VV_s = sum(n_O2_d(1 : end - 2)'*n_N2_d(3 : end - 1).*k_O2_N2_VV_r(1 : end - 1, :) - ...
                  n_O2_d(2 : end - 1)'*n_N2_d(2 : end - 2).*k_O2_N2_VV(1 : end - 1, :) + ...
                  n_O2_d(3 : end)'*n_N2_d(2 : end-2).*k_O2_N2_VV(2 : end, :) - ...
                  n_O2_d(2 : end - 1)'*n_N2_d(3 : end - 1).*k_O2_N2_VV_r(2 : end, :) , 2) + ...
                  sum(n_O2_d(1 : end - 2)'*n_NO_d(3 : end - 1).*k_O2_NO_VV_r(1 : end - 1, :) - ...
                  n_O2_d(2 : end - 1)'*n_NO_d(2 : end - 2).*k_O2_NO_VV(1 : end - 1, :) + ...
                  n_O2_d(3 : end)'*n_NO_d(2 : end-2).*k_O2_NO_VV(2 : end, :) - ...
                  n_O2_d(2 : end - 1)'*n_NO_d(3 : end - 1).*k_O2_NO_VV_r(2 : end, :) , 2);
       
   % ���������!!!!!
    R_NO_VV_s = sum(n_NO_d(1 : end - 2)'*n_N2_d(3 : end - 1).*k_NO_N2_VV_r(1 : end - 1, :) - ...
                  n_NO_d(2 : end - 1)'*n_N2_d(2 : end - 2).*k_NO_N2_VV(1 : end - 1, :) + ...
                  n_NO_d(3 : end)'*n_N2_d(2 : end-2).*k_NO_N2_VV(2 : end, :) - ...
                  n_NO_d(2 : end - 1)'*n_N2_d(3 : end - 1).*k_NO_N2_VV_r(2 : end, :) , 2) + ...
                  sum(n_NO_d(1 : end - 2)'*n_O2_d(3 : end - 1).*k_NO_O2_VV_r(1 : end - 1, :) - ...
                  n_NO_d(2 : end - 1)'*n_O2_d(2 : end - 2).*k_NO_O2_VV(1 : end - 1, :) + ...
                  n_NO_d(3 : end)'*n_O2_d(2 : end-2).*k_NO_O2_VV(2 : end, :) - ...
                  n_NO_d(2 : end - 1)'*n_O2_d(3 : end - 1).*k_NO_O2_VV_r(2 : end, :) , 2);
    
    R_N2_vibr = R_N2_VT + R_N2_VV + R_N2_VV_s; 
    R_O2_vibr = R_O2_VT + R_O2_VV + R_O2_VV_s;
    R_NO_vibr = R_NO_VT + R_NO_VV + R_NO_VV_s; 
    
    % Chemical exchange
    
    n_N2_d = n_N2.*n_cr;
    n_O2_d = n_O2.*n_cr;
    n_NO_d = n_NO.*n_cr;
    
    TT = 100 : 100 : 100000;
    
    k_ex_N2 = interp3(0 : l_NO - 1 , 0 : l_N2 - 1 , TT , k_ex_N2_STELLAR , 0 : l_NO - 1 , 0 : l_N2 - 1, T_d, 'spline');
    k_ex_O2 = interp3(0 : l_NO - 1 , 0 : l_O2 - 1 , TT , k_ex_O2_STELLAR , 0 : l_NO - 1 , 0 : l_O2 - 1, T_d, 'spline');
    
%     k_ex_N2loop = zeros(l_N2 , l_NO);
%     k_ex_O2loop = zeros(l_O2 , l_NO);
%     for q = 1 : l_NO
%         for v1 = 1 : l_N2
%             k_vector1 = reshape(k_ex_N2_STELLAR(v1,q,:),[1 1e3]);
%             k_ex_N2loop(v1,q) = interp1(TT,k_vector1,T_d,'spline');
%         end
%         for v2 = 1 : l_O2
%             k_vector2 = reshape(k_ex_O2_STELLAR(v2,q,:),[1 1e3]);
%             k_ex_O2loop(v2,q) = interp1(TT,k_vector2,T_d,'spline');
%         end
%     end

    k_ex_N2_r = k_ex_N2.*(M(1)*M(5)/M(3)/M(4))^(1.5).*...
                THETA_R(3)/THETA_R(1)*0.5.*exp((ones(l_N2 , 1)*e_i_NO - ...
                e_i_N2'*ones(1 , l_NO))./T + D(1)/K/T_d - D(3)/K/T_d);
    
    k_ex_O2_r = k_ex_O2.*(M(2)*M(4)/M(3)/M(5))^(1.5).*...
                THETA_R(3)/THETA_R(2)*0.5.*exp((ones(l_O2 , 1)*e_i_NO - ...
                e_i_O2'*ones(1 , l_NO))./T + D(2)/K/T_d - D(3)/K/T_d);
    
    R_N2_ex = sum(ones(l_N2 , 1)*n_NO_d'*n_N_d.*k_ex_N2_r - n_N2_d*ones(1, l_NO).*n_O_d.*k_ex_N2 , 2); 
    R_O2_ex = sum(ones(l_O2 , 1)*n_NO_d'*n_O_d.*k_ex_O2_r - n_O2_d*ones(1, l_NO).*n_N_d.*k_ex_O2 , 2);
    R_NO_ex = sum( - ones(l_N2 , 1)*n_NO_d'*n_N_d.*k_ex_N2_r + n_N2_d*ones(1, l_NO).*n_O_d.*k_ex_N2 , 1)' - ...
              sum( - ones(l_O2 , 1)*n_NO_d'*n_O_d.*k_ex_O2_r + n_O2_d*ones(1, l_NO).*n_N_d.*k_ex_O2 , 1)';
    R_N_ex = - sum(R_N2_ex) + sum(R_O2_ex);
    R_O_ex = sum(R_N2_ex) - sum(R_O2_ex);
    
    % Dissociation
    
    n_N2_d = n_N2.*n_cr*ones(1,5); 
    n_O2_d = n_O2.*n_cr*ones(1,5);
    n_NO_d = n_NO.*n_cr*ones(1,5);
    
    k_diss_N2 = k_diss(1,T_d);
    k_diss_O2 = k_diss(2,T_d);
    k_diss_NO = k_diss(3,T_d);
    
    k_rec_N2 = k_diss_N2.*(M(1)/M(4)^2)^(1.5).*H^3.*(2*pi*K*T_d)^(-1.5).*...
                T_d./THETA_R(1).*0.5.*exp(-e_i_N2'*ones(1,5)./T + D(1)/K/T_d); 
    k_rec_O2 = k_diss_O2.*(M(2)/M(5)^2)^(1.5).*H^3.*(2*pi*K*T_d)^(-1.5).*...
               T_d./THETA_R(2).*0.5.*exp(-e_i_O2'*ones(1,5)./T + D(2)/K/T_d);
    k_rec_NO = k_diss_NO.*(M(3)/M(1)/M(2))^(1.5).*H^3.*(2*pi*K*T_d)^(-1.5).*...
                T_d./THETA_R(3).*exp(-e_i_NO'*ones(1,5)./T + D(3)/K/T_d); 
    
    R_N2_diss = sum(n_c_N2.*(n_N_d^2.*k_rec_N2 - n_N2_d.*k_diss_N2) , 2);
    R_O2_diss = sum(n_c_O2.*(n_O_d^2.*k_rec_O2 - n_O2_d.*k_diss_O2) , 2);
    R_NO_diss = sum(n_c_NO.*(n_N_d*n_O_d.*k_rec_NO - n_NO_d.*k_diss_NO) , 2);  
    R_N_diss = -2*sum(R_N2_diss) - sum(R_NO_diss);
    R_O_diss = -2*sum(R_O2_diss) - sum(R_NO_diss);
    
    R_N2 = r_cr/n_cr/v_cr.*(R_N2_vibr + R_N2_ex + R_N2_diss);
    R_O2 = r_cr/n_cr/v_cr.*(R_O2_vibr + R_O2_ex + R_O2_diss);
    R_NO = r_cr/n_cr/v_cr.*(R_NO_vibr + R_NO_ex + R_NO_diss);
    R_N = r_cr/n_cr/v_cr.*(R_N_ex + R_N_diss);
    R_O = r_cr/n_cr/v_cr.*(R_O_ex + R_O_diss);
   
    %disp(sum(R_N2) + sum(R_O2) + sum(R_NO) + R_N + R_O);
    
    b(1 : l_N2) = R_N2 - v.*n_N2./S(1).*S(2); 
    b(l_N2 + 1 : l_N2 + l_O2) = R_O2 - v.*n_O2./S(1).*S(2);
    b(l_N2 + l_O2 + 1 : l_mol) = R_NO - v.*n_NO./S(1).*S(2);
    b(l_mol + 1) = R_N - v*n_N/S(1)*S(2);
    b(l_c) = R_O - v*n_O/S(1)*S(2);
    b(l_v) = 0;
    b(l_T) = -S(2)/S(1)*sum(n.*(T_energy2  + e_i + e_0 + e_f));

    dy = AA\b;
    
    end

end
