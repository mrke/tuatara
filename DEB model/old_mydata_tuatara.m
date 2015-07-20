%% mydata_tuatara
% Bas Kooijman at 2010/05/30; modified 2011/04/25
% Egernia cunninghami estimation by Michael Kearney, 18th July 2012

% the parameter estimation method for the standard DEB model is described in
%   http://www.bio.vu.nl/thb/deb/deblab/add_tuatara/add_tuatara.pdf
%   http://www.bio.vu.nl/thb/research/bib/LikaFrei2011.html
% copy this template; replace 'tuatara' by the name of your species and insert data
%   required at several places in this script
%   edit real data (not pseudo-data)
% this template is at completeness level 2.5
%   using one zero- and one-univariate data matrix
%   you can add more zero- and uni-variate data
%   and adapt the predict_tuatara-file for the computation of the expected values
% all rates and times must have an actual body temperature; 
%   T_ref (see below) is the temperature at which the parameter values are given
% if possible, only use units K, d, cm, g, mol, J
% provide the source of the data that you enter
% if you have less or more data: adapt predict_tuatara
% if you are satisfied with the result: 
%   copy parameter estimates, COMPLETENESS and FIT to the file pars_tuatara
%   run pars_tuatara and judge implied properties critically
% if all seems fine, submit your mydata-, predict- and pars-files 
%   to the add_tuatara collection: bas@bio.vu.nl

COMPLETE = 2.5;   % judge the level using LikaFrei2011; adjust if you have less or more data
FIT = 9.01;           % compute after having obtained the estimates
                      % this script shows the value on the screen
                      
  global dwm tT % pass d_O, w_O, mu_O directly to predict_tuatara

%% set data

% zero-variate data
% real data
% typically depend on scaled functional response f, see initial estimates below
%   here assumed to be equal for all real data
%   if they differ, the real values can be inserted in column 1 or 4 of  Data (see below)
ab = 183;      %  1 d, age at birth at f (age 0 is at onset of embryo development) 
  T_ab = 273 + 21; % K, temperature for ab 
  % observed age at birth is frequently larger than ab, because of diapauzes during incubation
ap = 12*365;     %  2 d, age at puberty at f 
  T_ap = 273 + 18.5; % K, temperature for ap 
  % observed age at puberty is fruently larger tha ap, because allocation
  % to reproduction starts below the first eggs appear
Lb = 5.3;    %  3 cm, physical length at birth at f 
Lp = 17;    %  4 cm, physical length at puberty at f 
Li = 22.0;    %  5 cm, ultimate physical length at f 
Wb = 4.7*0.3;  %  6 g, dry weight at birth at f 
Wp = 163*0.3;    %  7 g, dry weight at puberty at f 
Wi = 357*0.3;    %  8 g, ultimate dry weight at f 
Ri = 10/(4*365.);     %  9 #/d, maximum reproduction rate at f (for individual of max length) 
  T_Ri = 273 + 18.5; % K, temperature for Ri temperature for ap 
tm = 100*365;     % 10 y, life span at f (accounting for aging only) 
  T_tm = 273 + 18.5; % K, temperature for tm 
  
% pseudo-data from pars_tuatara at T_ref; don't change these data
%   they supplement the real data, if necessary, to provide enough information
%   for all parameters to be estimated
v = 0.02;     % 11 cm/d, energy conductance
kap = 0.5;    % 12 -, allocation fraction to soma = growth + somatic maintenance
kap_R = 0.95; % 13 -, reproduction efficiency
p_M = 18;     % 14 J/d.cm^3, [p_M] vol-specific somatic maintenance
p_T =  0;     % 15 J/d.cm^2, {p_T} surface-specific som maintenance
k_J = 0.002;  % 16 1/d, < k_M = p_M/E_G, maturity maint rate coefficient
kap_G = .8;   % 17 -, growth efficiency

% pack data
data = [ab; ap; Lb; Lp; Li; Wb; Wp; Wi; Ri; tm; % 01:10 real data
        v; kap; kap_R; p_M; p_T; k_J; kap_G];   % 11:17 pseudo data

% nmregr and nmcvregr want to have 3 columns: independent var, dependent var, weight coeff
%   prepend independent variable (generally not used for zero-variate data)
%   append weight coefficients inverse to squared data, but modify according to insight
% weight coefficients for WLS and ML criteria differ
%Data = data(:,[1 1 1]); Data(:,3) = 1; % nmvcregr, nrvcregr (ML criterion)
Data = [data(:,[1 1]), min(100,1 ./ max(1e-6, data) .^ 2)]; % nmregr, nrregr (WLS criterion)
Data(1:10,3) = 10 * Data(1:10,3); % give real data more weight
%Data(6:8,3) = 10 * Data(6:8,3);   % give weight data more weight
%Data([1,9],3) = 10 * Data([1,9],3);  % more weight to repro and age at birth
Data(17,3) = 20 * Data(17,3);  % more weight to kap_G
%Data(2,3) = 20 * Data(2,3);  % more weight to age at puberty
%Data(12,3) = 10 * Data(12,3);  % more weight to kappa
% insert temperature data for rates and times in the first column
Data([1 2 9 10], 1) = [T_ab; T_ap; T_Ri; T_tm];

txt_data = {... % for presentation of predictions
    '1 ab, d, age at birth ';
    '2 ap, d, age at puberty ';
    '3 Lb, cm, physical length at birth ';
    '4 Lp, cm, physical length at puberty ';
    '5 Li, cm, ultimate physical length ';
    '6 Wb, g, dry weight at birth ';
    '7 Wp, g, dry weight at puberty ';
    '8 Wi, g, ultimate dry weight ';
    '9 Ri, #/d, maximum reprod rate ';
   '10 am, d, life span ';
   '11 v, cm/d, energy conductance ';
   '12 kap, -, allocation fraction to soma ';
   '13 kap_R, -, reproduction efficiency ';
   '14 [p_M], J/d.cm^3, vol-spec som maint ';
   '15 {p_T}, J/d.cm^2, sur-spec som maint ';
   '16 k_J, 1/d, maturity maint rate coefficient ';
   '17 kap_G, -, growth efficiency'};

% uni-variate data at f = 0.8 and T = 273 + 25
%  This info is copied directly in predict_my_pat
%  tL = [0     50  100 200 300 400 500 600;    % d, time since birth
%        0.45  1.1 1.7 2.7 3.4 4.0 4.5 4.9]'; % cm, physical length at f and T
%  tL = tL(:,[1 2 2]); tL(:,3) = 1; % append weight coefficients for ML criterion
  %tL = [tL, 10./tL(:,2).^2]; % append weight coefficients for WLS criterion

%% conversion coefficients (selected copy-paste from pars_tuatara)
%  don't change these values, unless you have a good reason

% chemical indices
%       X     V     E     P
n_O = [1.00, 1.00, 1.00, 1.00;  % C/C, equals 1 by definition
       1.80, 1.80, 1.80, 1.80;  % H/C  these values show that we consider dry-mass
       0.50, 0.50, 0.50, 0.50;  % O/C
       0.15, 0.15, 0.15, 0.15]; % N/C
%       C     H     O     N
n_M = [ 1     0     0     0;    % C/C, equals 0 or 1
        0     2     0     3;    % H/C
        2     1     2     0;    % O/C
        0     0     0     1];   % N/C
    
% specific densities
%       X     V     E     P
d_O = [0.3;  0.3;  0.3;  0.3];    % g/cm^3, specific densities for organics
% For a specific density of wet mass of 1 g/cm^3,
%   a specific density of d_E = d_V = 0.1 g/cm^3 means
%   a dry-over-wet weight ratio of 0.1
%   modify this for your pet

% chemical potentials
%        X     V     E     P
mu_O = [525; 500;  550;  480] * 1000; % J/mol, chemical potentials for organics

% molecular weights
w_O = n_O' * [12; 1; 16; 14];  % g/mol, mol-weights for organics

% pack coefficients
dwm = [d_O, w_O, mu_O]; % g/cm^3, g/mol, kJ/mol spec density, mol weight, chem pot

% temp profile for puberty
tT = [0:12; 273 + [25 23 21 18 16 14 12 14 16 18 20 22 24]]';

%% parameters: initial values at T_ref
% edit these values such that predictions are not too far off

T_ref  = 293;      % 1 K, temp for which rate pars are given; don't change this vulue
T_A  = 11689;       % 2 K, Arrhenius temp Wilson cited in Shine 1971 (unpub thesis)
f = 1.0;           % 3 scaled functional response
z = 5.892;             % 4 -, zoom factor; for z = 1: L_m = 1 cm
del_M = 0.2688;      % 5 -, shape coefficient
F_m = 6.5;         % 6 l/d.cm^2, {F_m} max spec searching rate
kap_X = 0.8;       % 7 -, digestion efficiency of food to reserve
v = 0.02805;          % 8 cm/d, energy conductance
kap = 0.5522;         % 9 -, alloaction fraction to soma = growth + somatic maintenance
kap_R = 0.95;      %10 -, reproduction efficiency
p_M = 12.41;          %11 J/d.cm^3, [p_M] vol-specific somatic maintenance
p_T =  0;          %12 J/d.cm^2, {p_T} surface-specific som maintenance
E_G = 7932;        %14 J/cm^3, [E_G], spec cost for structure
k_M = p_M/ E_G;
k_J = k_M;       %13 1/d, < k_M = p_M/E_G, maturity maint rate coefficient
E_Hb = 1.808e+004;       %15 J, E_H^b maturity threshold at birth
E_Hp = 6.385e+005;         %16 J, E_H^p maturity threshold at puberty
h_a = 3.535e-011;        %17 1/d^2, Weibull aging acceleration
s_G = 0;           %18 -, Gompertz stress coefficient
T_L  =   285; % K, lower boundary tolerance range
T_H  =   304; % K, upper boundary tolerance range
T_AL = 50000; % K, Arrhenius temp for lower boundary
T_AH = 65000; % K, Arrhenius temp for upper boundary

% pack parameters and fix T_ref and f and possibly other at well
%   in second column: 0 = fix; 1 = release
pars = [T_ref 0; T_A 0; f    0; z     1; del_M 1; F_m  0;  
        kap_X 0; v   1; kap  1; kap_R 0; p_M   1; p_T  0;  
        k_J   0; E_G 1; E_Hb 1; E_Hp  1; h_a   1; s_G  0; T_L 0; T_H 0; T_AL 0; T_AH 0];

txt_pars = { ...    % for presentation of parameter estimates
  'T_ref, K'; 'T_A, K';        'f, -'; 
  'z, -';     'del_M, -';      '{F_m}, l/d.cm^2'; 
  'kap_X, -'; 'v, cm/d';       'kap, -'; 
  'kap_R, -'; '[p_M], J/d.cm^3'; '{p_T}, J/d.cm^2'; 
  'k_J, 1/d'; '[E_G], J/cm^3'; 'E_Hb, J'; 
  'E_Hp, J';  'h_a, 1/d^2';    's_G, -'; 'T_L'; 'T_H'; 'T_AL'; 'T_AH'};

%% estimate parameters

nmregr_options('default'); % set options for parameter estimation
nmregr_options('max_step_number',5e3); % set options for parameter estimation
nmregr_options('max_fun_evals',2e4);   % set options for parameter estimation

% first out-comment parameter estimation, and run this script
%   manually tune initial estimates till the predictions are not too far off
%   then activate the parameter estimation keeping the more certain parameters fixed
%   copy the result in the initial estimates, release the more certain parameters
%   repeat the estimation, possibly after fine-tuning the weight coefficients
%   (if a prediction is too far off, increase the weight coefficient for that data point)

%pars = nmvcregr('predict_tuatara', pars, Data, tL);  % ML estimate parameters using overwrite
%pars = nrvcregr('predict_tuatara', pars, Data, tL);  % ML estimate parameters using overwrite
%pars = nmregr('predict_tuatara', pars, Data);   % WLS estimate parameters using overwrite
%pars = nrregr('predict_tuatara', pars, Data);   % WLS estimate parameters using overwrite
sd = 0 * pars(:,1);                                 % initiate standard deviations

%if 0 % start with 0, but choose 1 if initial parameter values are the resulting ones
  %[cov cor sd] = pvcregr('predict_tuatara', pars, Data); % get standard deviation for ML
  [cov cor sd] = pregr('predict_tuatara', pars, Data); % get standard deviation for WLS
%end

%% get FIT

Data(:,3) = 0; Data(1:10,3) = 1; % give unit weight to real data, zero to pseudo-data
%[mean_relative_error relative_error] = mrevc('predict_tuatara', pars, Data, tL); % ML-method
 [mean_relative_error relative_error] = mre('predict_tuatara', pars, Data); % WLS-method
FIT = 10 * (1 - mean_relative_error) % get mark for goodness of fit 
% copy FIT to header and to pars_tuatara; mrevc and mre gives mean relative error

%% get predictions

t = linspace(0,650,100)'; % times for plotting length data
[Edata] = predict_tuatara(pars(:,1), Data); % notice use of first column of pars only

%% get initial reserve
mu_E = dwm(3,3);
p_Am = z * p_M/ kap;
E_m = p_Am/ v;
g = E_G/ kap/ E_m;
J_E_Am = p_Am/ mu_E;
M_Hb = E_Hb/ mu_E;
U_Hb = M_Hb/ J_E_Am;
V_Hb = U_Hb/ (1 - kap);
pars_E0 = [V_Hb; g; k_J; k_M; v]; % pars for initial_scaled_reserve
[U_E0 L_b info] = initial_scaled_reserve(f, pars_E0); % d cm^2, initial scaled reserve
E_0 = p_Am * U_E0;    % J, initial reserve (of embryo)
%% present results

printpar(txt_data, data, Edata, 'data and expectations'); % for zero-variate data
fprintf('\n') % insert a blank line in screen output
printpar(txt_pars, pars, sd);

close all % remove existing figures, else you get more and more if you retry

% figure % one figure to show results of uni-variate data
% plot(tL(:,1), tL(:,2), '.r', t, EL, 'g')
% xlabel('time since birth, d')
% ylabel('length, cm')

% if you are satisfied by the results in the MATLAB window, 
%   copy the resulting estimates in the initial values and in pars_tuatara