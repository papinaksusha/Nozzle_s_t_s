%% 5-components air mixture N2/O2/NO/N/O, state-to-state kinetics in nozzle  

clc
clear 
close all
format long e
addpath(('./MAT/'))
addpath(('./ODE_systems/')) % переписать кросс-платформенно
addpath(('./Coefficients/'))

global NA K H C W WX M THETA_R D I SW_O SW_N EX_MODEL MOLAR REC VV VT

%% CONSTs

% J = kg*m^2/s^2, Pa = kg/m/s^2
% [N2 O2 NO N O]

tic

NA = 6.022141e23;                                                          % mol^-1, Avogadro constant 
K = 1.380648e-23;                                                          % J/K, Boltzmann constant
H = 6.626070e-34;                                                          % J*s, Planck constant
C = 299792458;                                                             % m/s, speed of light 
W = [235857 158019 190420];                                                % m^-1, spectroscopic constant
MOLAR = 1e-3.*[28.0134 31.99880 30.00610 14.0067 15.9994];                 % kg/mol, molar mass
M = 1.6605402e-24.*MOLAR;                                                  % kg, molecular mass
THETA_R = [2.86 2.07 2.42];                                                % K, characteristic rotational emperature
D = H*C.*[7.871e6 4.126e6 5.24e6];                                         % J, dissociation energy
VV = 1;
VT = 1;

SW_O = 2;                                                                  % switch on oscillator, 1 -  harmonic oscillator; 2 -  anharmonic oscillator
SW_N = 1;                                                                  % switcn on nozzle
EX_MODEL = 2;                                                              % exchange rate coefficients model
REC = 1;                                                                   % rec = 0 - without recombination   

switch SW_O 
    case 1
        WX = [0 0 0];   
    case 2
        WX = [1432 1198 1407.5];                                           % m^-1, spectroscopic constant
end
             
I = [33 26 27; 47 36 39];                                                  % Max vibrational levels                                                 
                                                                           % I(1,:) - harmonic oscillator, I(2,:) - anharmonic oscillator 
%%

fig = 1;

% TEST_SSH;
% TEST_SSH_ALEX
% TEST_DISS_TM;
% TEST_EXCHANGE
 VVprobability_compare
 
% call5;
% call2;