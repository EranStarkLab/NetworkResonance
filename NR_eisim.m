% NR_eisim              an E-cell and an I-cell that receives input to one of the cells
%
% call                  [ t, Ve_1, Vi_1, finput_e, finput_i ] = NR_eisim( MDL )
%
% gets                  MDL         model.  4:  I-to-E, chirp to I. 
%                                               w/o E-E, I-I, or E-to-I connections
%                                           5:  isolated E, chirp to E.
%                                           6:  I-to-E, chirp to I; gamma-resonant I
%                                               w/o E-E, I-I, or E-to-I connections
%                                           7:  isolated gamma I, chirp to I.
%
% model can be overloaded with parameters. In that case, MDL is a
% two-element cell array: { MDL, MDL_params }
%
%                       MDL_params  5-element vector: 
%                                   [ Ain Iapp G_I2E D_i D_e ]
%                                   where   Ain     chirp amplitude [muA/cm^2]
%                                           Iapp    bias current [muA/cm^2]
%                                           G_I2E   conductance of the inhibitory synapse [mS/cm^2]
%                                           D_i     SD of the I-cell noise [mV] 
%                                           D_e     SD of the E-cell noise [mV] 
%
% returns               t           time [ms]
%                       Ve_1        E-cell membrane potential [mV]
%                       Vi_1        I-cell membrane potential [mV]
%                       finput_e    chirp input to E-cell [muA/cm^2]
%                       finput_i    chirp input to I-cell [muA/cm^2]
%
% calls                 nothing

% references:
% Basic model:          [BTLBK-2012]
% I-cells:              [WB-1996]
% E-cells:              Modified Traub-Miles model + M-current [OWCK-2003]
%                       E-cell M-currents: [BTLBK-2012]
%                       E-cell h-current: [PMJ-2002], [ZKPFH-2012]
%                       E-cell AHP-current: [J-2001]
% Synaptic connections: [BTLBK-2012]
%
% Numerical method: Modified-Euler (Runge-Kutta II)

% 06-aug-12 HGR

% last update
% 01-jul-22

function [ t, Ve_1, Vi_1, finput_e, finput_i ] = NR_eisim( MDL )

% [C] = muF/cm^2
% [I] = muA/cm^2
% [G] = mS/cm^2
% [t] = msec
% [V] = mV
% [freq] = Hz

%------------------------------------------------------------------------
% arguments
%------------------------------------------------------------------------
nargs                           = nargin;
if nargs < 1 || isempty( MDL )
    MDL                         = 4;
end
if isa( MDL, 'cell' )
    MDL_params                  = MDL{ 2 };
    MDL                         = MDL{ 1 };
else
    MDL_params                  = [];
end
if ~ismember( MDL, 4 : 7 ) 
    error( 'MDL should be 4, 5, 6, or 7' )
end

%------------------------------------------------------------------------
% Activation/Inactivation functions
%------------------------------------------------------------------------

% E-cells
minfek                          = @(v) (0.32*(54+v)./(1-exp(-(v+54)/4)))/((0.32*(54+v)./(1-exp(-(v+54)/4)))+(0.28*(v+27)./(exp((v+27)/5)-1)));
hinfek                          = @(v) (0.128*exp(-(50+v)/18)) ./((0.128*exp(-(50+v)/18))+(4./(1+exp(-(v+27)/5))));
htauek                          = @(v) 1.0 ./((0.128*exp(-(50+v)/18))+(4./(1+exp(-(v+27)/5))));
ninfek                          = @(v) (0.032*(v+52)./(1-exp(-(v+52)/5))) ./ ((0.032*(v+52)./(1-exp(-(v+52)/5)))+(0.5*exp(-(57+v)/40)));
ntauek                          = @(v) 1.0 ./ ((0.032*(v+52)./(1-exp(-(v+52)/5)))+(0.5*exp(-(57+v)/40)));
qinfek                          = @(v) 1./(1+exp(-(v+35)/10));
qtauek                          = @(v) 400./(3.3*exp((v+35)/20)+exp(-(v+35)/20));
pinfe                           = @(v) 1./(1+exp(-(v+38)/6.67));
rinfe                           = @(v) 1./(1 + exp((v+82.9)/12.4));
rtaue                           = @(v) 1.5*(exp(0.033*(v+75))./(0.011*(1+exp(0.083*(v+75)))));

% I-cells
minfwb                          = @(v) (-0.2*(v+35)./(exp(-0.1*(v+35))-1)) ./((-0.2*(v+35)./(exp(-0.1*(v+35))-1))+(4*exp(-(v+60)/18)));
hinfwb                          = @(v) (0.07*(exp(-(v+58)/20))) ./((0.07*(exp(-(v+58)/20)))+(1.0./(exp(-0.1*(v+28))+1)));
htauwb                          = @(v) 1.0 ./((0.07*(exp(-(v+58)/20)))+(1.0./(exp(-0.1*(v+28))+1)));
ninfwb                          = @(v) (-0.01*(v+34)./(exp(-0.1*(v+34))-1)) ./((-0.01*(v+34)./(exp(-0.1*(v+34))-1))+(0.125*exp(-(v+44)/80)));
ntauwb                          = @(v) 1.0 ./((-0.01*(v+34)./(exp(-0.1*(v+34))-1))+(0.125*exp(-(v+44)/80)));
qinfwb                          = @(v) 1./(1+exp(-(v+35)/10));
qtauwb                          = @(v) 200./(3.3*exp((v+35)/20)+exp(-(v+35)/20));

% Sigmoid function
H                               = @(v) 0.5*(1 + tanh(v/4));

%------------------------------------------------------------------------
% parameters
%------------------------------------------------------------------------
% E-cell parameters
C_e                             = 1.0;
Gna_e                           = 100;
Gk_e                            = 80;
Gl_e                            = 0.1;
Gm_e                            = 0.0;
Gp_e                            = 0.0;
Gh_e                            = 0.485;
Ena_e                           = 50;
Ek_e                            = -100;
El_e                            = -67;
Eh_e                            = -33;
Iapp_e                          = -2.7;

% I-cell parameters
C_i                             = 1.0;
Gna_i                           = 35;
Gk_i                            = 9;
Gl_i                            = 0.1;
Gm_i                            = 4;
Ena_i                           = 55;
Ek_i                            = -90;
El_i                            = -65;
phi_i                           = 5;
Iapp_i                          = -0.5;

% Synaptic parameters
taurse_e                        = 0.1;     % AMPA
taudec_e                        = 3.0;     % AMPA
taurse_i                        = 0.3;     % GABA_A
taudec_i                        = 9.0;     % GABA_A
Eex                             = 0.0;
Ein                             = -80;

% Connectivity parameters
% Gxy = Connection from X-cell to Y-cell
Gee                             = 0.0;
Gei                             = 0.15;
Gii                             = 0.05;
Gie                             = 0.27;

% noise
D_e                             = 0.1;
D_i                             = 0.1;

%----------------------------------------------------------------
% Time and Sinusoidal Input definitions
%----------------------------------------------------------------

% initialize output
Ve_1                            = [];
Vi_1                            = [];

% build chirp input
Ain_e                           = 0.0;
Ain_i                           = 1.35;

Fs                              = 40000;     % 40 kHz -> 0.025 ms
duration                        = 10; % [s]
f0                              = 0;
f1                              = 40;
if MDL == 7
    f1                          = 80;
end

t                               = (1/Fs:1/Fs:duration)';
Fzap                            = sin(2*pi*((f1-f0)/(2*t(end)).*t.^2+f0.*t+270/360));
dt                              = 1/Fs*1000;
t                               = dt:dt:duration*1000;
finput_i                        = Ain_i * Fzap;
finput_e                        = Ain_e * Fzap;

% generate white noise
wne_1                           = randn(1,length(t));
wni_1                           = randn(1,length(t));

%----------------------------------------------------------------
% simultation
%----------------------------------------------------------------
if ( MDL == 4 ) || ( MDL == 5 )
    
    if MDL == 4 % remove all synapses except the one from I to E
        
        % only I-to-E (no synapses E-to-I, I-to-I)
        Gei                     = 0;                                        % no feedback from PYR to INT
        Gii                     = 0;                                        % no recurrent inhibition
        Ain_i0                  = Ain_i;
        
        if isempty( MDL_params ) || length( MDL_params ) < 3
            Ain_i               = 0.5;
            Iapp_i              = -0.5;
            Gie                 = 0.27;
        else
            Ain_i               = MDL_params( 1 );
            Iapp_i              = MDL_params( 2 );
            Gie                 = MDL_params( 3 );
        end
        
        % rescale the ZAP to the I-cell by the requested Ain_i
        finput_i                = finput_i / Ain_i0 * Ain_i;
        
    else % MDL 5: give the input directly to the PYR
        
        % no connectivity at all
        Gei                     = 0;
        Gii                     = 0;
        Gie                     = 0;
        Gee                     = 0;
        
        if isempty( MDL_params ) || length( MDL_params ) < 2
            Ain_e               = 0.2;
            Iapp_e              = -2.7;
        else
            Ain_e               = MDL_params( 1 );
            Iapp_e              = MDL_params( 2 );
        end
        
        % swap the zap to the E and rescale by the requested Ain_e
        finput_e                = finput_i / Ain_i * Ain_e;
        
        % to prevent I spikes, no noise to the I, hyperpolarizing current
        wni_1( : )              = 0;
        Ain_i                   = 0;
        Iapp_i                  = -1;
        finput_i( : )           = Ain_i;
        
    end
    
    % set noise level for E- and I-cells manually
    if length( MDL_params ) > 3
        D_i                     = MDL_params( 4 );
    end
    if length( MDL_params ) > 4
        D_e                     = MDL_params( 5 );
    end
    
    % Initial conditions
    Ve_1                        = zeros(1,length(t));
    he_1                        = zeros(1,length(t));
    ne_1                        = zeros(1,length(t));
    qe_1                        = zeros(1,length(t));
    re_1                        = zeros(1,length(t));
    Se_1                        = zeros(1,length(t));
    
    Vi_1                        = zeros(1,length(t));
    hi_1                        = zeros(1,length(t));
    ni_1                        = zeros(1,length(t));
    Si_1                        = zeros(1,length(t));
    
    Ve_1(1)                     = -65.157;
    he_1(1)                     = 0.99355;
    ne_1(1)                     = 0.050507;
    qe_1(1)                     = 0.032955;
    re_1(1)                     = 0.19177;
    Se_1(1)                     = 0.0;
    
    Vi_1(1)                     = El_i;
    hi_1(1)                     = 0.94645;
    ni_1(1)                     = 0.036746;
    Si_1(1)                     = 0.0;
    
    % Computation of the solution
    for j                       = 1:length(t)-1
        
        k1_ve1 = Iapp_e - Gl_e*(Ve_1(j)-El_e);
        k1_ve1 = k1_ve1 - Gna_e*he_1(j)*minfek(Ve_1(j))^3*(Ve_1(j)-Ena_e);
        k1_ve1 = k1_ve1 - Gk_e*ne_1(j)^4*(Ve_1(j)-Ek_e);
        k1_ve1 = k1_ve1 - Gm_e*qe_1(j)*(Ve_1(j)-Ek_e);
        k1_ve1 = k1_ve1 - Gp_e*pinfe(Ve_1(j))*(Ve_1(j)-Ena_e);
        k1_ve1 = k1_ve1 - Gh_e*re_1(j)*(Ve_1(j)-Eh_e);
        k1_ve1 = k1_ve1 - Gee*Se_1(j)*(Ve_1(j)-Eex);
        k1_ve1 = k1_ve1 - Gie*Si_1(j)*(Ve_1(j)-Ein);
        k1_ve1 = k1_ve1 + finput_e(j);
        k1_ve1 = k1_ve1 / C_e;
        k1_he1 = (hinfek(Ve_1(j))-he_1(j))/htauek(Ve_1(j));
        k1_ne1 = (ninfek(Ve_1(j))-ne_1(j))/ntauek(Ve_1(j));
        k1_qe1 = (qinfek(Ve_1(j))-qe_1(j))/qtauek(Ve_1(j));
        k1_re1 = (rinfe(Ve_1(j))-re_1(j))/rtaue(Ve_1(j));
        k1_se1 = H(Ve_1(j))*(1-Se_1(j))/taurse_e-Se_1(j)/taudec_e;
        %
        k1_vi1 = Iapp_i - Gl_i*(Vi_1(j)-El_i);
        k1_vi1 = k1_vi1 - Gna_i*hi_1(j)*minfwb(Vi_1(j))^3*(Vi_1(j)-Ena_i);
        k1_vi1 = k1_vi1 - Gk_i*ni_1(j)^4*(Vi_1(j)-Ek_i);
        k1_vi1 = k1_vi1 - Gei*Se_1(j)*(Vi_1(j)-Eex);
        k1_vi1 = k1_vi1 - Gii*Si_1(j)*(Vi_1(j)-Ein);
        k1_vi1 = k1_vi1 + finput_i(j);
        k1_vi1 = k1_vi1 / C_i;
        k1_hi1 = phi_i*(hinfwb(Vi_1(j))-hi_1(j))/htauwb(Vi_1(j));
        k1_ni1 = phi_i*(ninfwb(Vi_1(j))-ni_1(j))/ntauwb(Vi_1(j));
        k1_si1 = H(Vi_1(j))*(1-Si_1(j))/taurse_i-Si_1(j)/taudec_i;
        %
        ave_1 = Ve_1(j) + k1_ve1*dt;
        ave_1 = ave_1 + D_e*dt^0.5*wne_1(j);
        ahe_1 = he_1(j) + k1_he1*dt;
        ane_1 = ne_1(j) + k1_ne1*dt;
        aqe_1 = qe_1(j) + k1_qe1*dt;
        are_1 = re_1(j) + k1_re1*dt;
        ase_1 = Se_1(j) + k1_se1*dt;
        %
        avi_1 = Vi_1(j) + k1_vi1*dt;
        avi_1 = avi_1 + D_i*dt^0.5*wni_1(j);
        ahi_1 = hi_1(j) + k1_hi1*dt;
        ani_1 = ni_1(j) + k1_ni1*dt;
        asi_1 = Si_1(j) + k1_si1*dt;
        %
        k2_ve1 = Iapp_e - Gl_e*(ave_1-El_e);
        k2_ve1 = k2_ve1 - Gna_e*ahe_1*minfek(ave_1)^3*(ave_1-Ena_e);
        k2_ve1 = k2_ve1 - Gk_e*ane_1^4*(ave_1-Ek_e);
        k2_ve1 = k2_ve1 - Gm_e*aqe_1*(ave_1-Ek_e);
        k2_ve1 = k2_ve1 - Gp_e*pinfe(ave_1)*(ave_1-Ena_e);
        k2_ve1 = k2_ve1 - Gh_e*are_1*(ave_1-Eh_e);
        k2_ve1 = k2_ve1 - Gee*ase_1*(ave_1-Eex);
        k2_ve1 = k2_ve1 - Gie*asi_1*(ave_1-Ein);
        k2_ve1 = k2_ve1 + finput_e(j+1);
        k2_ve1 = k2_ve1 / C_e;
        k2_he1 = (hinfek(ave_1)-ahe_1)/htauek(ave_1);
        k2_ne1 = (ninfek(ave_1)-ane_1)/ntauek(ave_1);
        k2_qe1 = (qinfek(ave_1)-aqe_1)/qtauek(ave_1);
        k2_re1 = (rinfe(ave_1)-are_1)/rtaue(ave_1);
        k2_se1 = H(ave_1)*(1-ase_1)/taurse_e-ase_1/taudec_e;
        %
        k2_vi1 = Iapp_i - Gl_i*(avi_1-El_i);
        k2_vi1 = k2_vi1 - Gna_i*ahi_1*minfwb(avi_1)^3*(avi_1-Ena_i);
        k2_vi1 = k2_vi1 - Gk_i*ani_1^4*(avi_1-Ek_i);
        k2_vi1 = k2_vi1 - Gei*ase_1*(avi_1-Eex);
        k2_vi1 = k2_vi1 - Gii*asi_1*(avi_1-Ein);
        k2_vi1 = k2_vi1 + finput_i(j+1);
        k2_vi1 = k2_vi1 / C_i;
        k2_hi1 = phi_i*(hinfwb(avi_1)-ahi_1)/htauwb(avi_1);
        k2_ni1 = phi_i*(ninfwb(avi_1)-ani_1)/ntauwb(avi_1);
        k2_si1 = H(avi_1)*(1-asi_1)/taurse_i-asi_1/taudec_i;
        %
        Ve_1(j+1) = Ve_1(j) + (k1_ve1 + k2_ve1)*dt/2;
        Ve_1(j+1) = Ve_1(j+1) + D_e*dt^0.5*wne_1(j);
        he_1(j+1) = he_1(j) + (k1_he1 + k2_he1)*dt/2;
        ne_1(j+1) = ne_1(j) + (k1_ne1 + k2_ne1)*dt/2;
        qe_1(j+1) = qe_1(j) + (k1_qe1 + k2_qe1)*dt/2;
        re_1(j+1) = re_1(j) + (k1_re1 + k2_re1)*dt/2;
        Se_1(j+1) = Se_1(j) + (k1_se1 + k2_se1)*dt/2;
        %
        Vi_1(j+1) = Vi_1(j) + (k1_vi1 + k2_vi1)*dt/2;
        Vi_1(j+1) = Vi_1(j+1) + D_i*dt^0.5*wni_1(j);
        hi_1(j+1) = hi_1(j) + (k1_hi1 + k2_hi1)*dt/2;
        ni_1(j+1) = ni_1(j) + (k1_ni1 + k2_ni1)*dt/2;
        Si_1(j+1) = Si_1(j) + (k1_si1 + k2_si1)*dt/2;
        
    end
    
elseif ( MDL == 6 || MDL == 7 )
    
    % reduce Ain_i
    Ain_i0                      = Ain_i;
    Ain_i                       = MDL_params( 1 );
    Iapp_i                      = MDL_params( 2 );
    if length( MDL_params ) > 2
        Gie                     = MDL_params( 3 );
    end
    
    % remove E-to-I, I-to-I
    Gei                         = 0;                                        % no feedback from PYR to INT
    Gii                         = 0;                                        % no recurrent inhibition
    if MDL == 7                                                             % by default, noiseless, isolated gamma I-cell
        Gie                     = 0;
        D_i                     = 0;
    end
    
    % set noise level for I-cells manually
    if length( MDL_params ) > 3
        D_i                     = MDL_params( 4 );
    end
    
    % rescale the ZAP to the I-cell by the requested Ain_i
    finput_i                    = finput_i / Ain_i0 * Ain_i;
    
    % Initial conditions
    Ve_1                        = zeros(1,length(t));
    he_1                        = zeros(1,length(t));
    ne_1                        = zeros(1,length(t));
    qe_1                        = zeros(1,length(t));
    re_1                        = zeros(1,length(t));
    Se_1                        = zeros(1,length(t));
    
    Vi_1                        = zeros(1,length(t));
    hi_1                        = zeros(1,length(t));
    ni_1                        = zeros(1,length(t));
    qi_1                        = zeros(1,length(t));
    Si_1                        = zeros(1,length(t));
    
    Ve_1(1)                     = -64.692;
    he_1(1)                     = 0.9927;
    ne_1(1)                     = 0.054578;
    qe_1(1)                     = 0.034574;
    re_1(1)                     = 0.18454;
    Se_1(1)                     = 0.0;
    
    Vi_1(1)                     = -67.254; 
    hi_1(1)                     = 0.85187;
    ni_1(1)                     = 0.069048;
    qi_1(1)                     = 0.03812;   
    Si_1(1)                     = 0.0;
    
    % Computation of the solution
    for j=1:length(t)-1
        
        k1_ve1 = Iapp_e - Gl_e*(Ve_1(j)-El_e);
        k1_ve1 = k1_ve1 - Gna_e*he_1(j)*minfek(Ve_1(j))^3*(Ve_1(j)-Ena_e);
        k1_ve1 = k1_ve1 - Gk_e*ne_1(j)^4*(Ve_1(j)-Ek_e);
        k1_ve1 = k1_ve1 - Gm_e*qe_1(j)*(Ve_1(j)-Ek_e);
        k1_ve1 = k1_ve1 - Gp_e*pinfe(Ve_1(j))*(Ve_1(j)-Ena_e);
        k1_ve1 = k1_ve1 - Gh_e*re_1(j)*(Ve_1(j)-Eh_e);
        k1_ve1 = k1_ve1 - Gee*Se_1(j)*(Ve_1(j)-Eex);
        k1_ve1 = k1_ve1 - Gie*Si_1(j)*(Ve_1(j)-Ein);
        k1_ve1 = k1_ve1 + finput_e(j);
        k1_ve1 = k1_ve1 / C_e;
        k1_he1 = (hinfek(Ve_1(j))-he_1(j))/htauek(Ve_1(j));
        k1_ne1 = (ninfek(Ve_1(j))-ne_1(j))/ntauek(Ve_1(j));
        k1_qe1 = (qinfek(Ve_1(j))-qe_1(j))/qtauek(Ve_1(j));
        k1_re1 = (rinfe(Ve_1(j))-re_1(j))/rtaue(Ve_1(j));
        k1_se1 = H(Ve_1(j))*(1-Se_1(j))/taurse_e-Se_1(j)/taudec_e;
        %
        k1_vi1 = Iapp_i - Gl_i*(Vi_1(j)-El_i);
        k1_vi1 = k1_vi1 - Gna_i*hi_1(j)*minfwb(Vi_1(j))^3*(Vi_1(j)-Ena_i);
        k1_vi1 = k1_vi1 - Gk_i*ni_1(j)^4*(Vi_1(j)-Ek_i);
        k1_vi1 = k1_vi1 - Gm_i*qi_1(j)*(Vi_1(j)-Ek_i);
        k1_vi1 = k1_vi1 - Gei*Se_1(j)*(Vi_1(j)-Eex);
        k1_vi1 = k1_vi1 - Gii*Si_1(j)*(Vi_1(j)-Ein);
        k1_vi1 = k1_vi1 + finput_i(j);
        k1_vi1 = k1_vi1 / C_i;
        k1_hi1 = phi_i*(hinfwb(Vi_1(j))-hi_1(j))/htauwb(Vi_1(j));
        k1_ni1 = phi_i*(ninfwb(Vi_1(j))-ni_1(j))/ntauwb(Vi_1(j));
        k1_qi1 = phi_i*(qinfwb(Vi_1(j))-qi_1(j))/qtauwb(Vi_1(j));
        k1_si1 = H(Vi_1(j))*(1-Si_1(j))/taurse_i-Si_1(j)/taudec_i;
        %
        ave_1 = Ve_1(j) + k1_ve1*dt;
        ave_1 = ave_1 + D_e*dt^0.5*wne_1(j);
        ahe_1 = he_1(j) + k1_he1*dt;
        ane_1 = ne_1(j) + k1_ne1*dt;
        aqe_1 = qe_1(j) + k1_qe1*dt;
        are_1 = re_1(j) + k1_re1*dt;
        ase_1 = Se_1(j) + k1_se1*dt;
        avi_1 = Vi_1(j) + k1_vi1*dt;
        avi_1 = avi_1 + D_i*dt^0.5*wni_1(j);
        ahi_1 = hi_1(j) + k1_hi1*dt;
        ani_1 = ni_1(j) + k1_ni1*dt;
        aqi_1 = qi_1(j) + k1_qi1*dt;
        asi_1 = Si_1(j) + k1_si1*dt;
        %
        k2_ve1 = Iapp_e - Gl_e*(ave_1-El_e);
        k2_ve1 = k2_ve1 - Gna_e*ahe_1*minfek(ave_1)^3*(ave_1-Ena_e);
        k2_ve1 = k2_ve1 - Gk_e*ane_1^4*(ave_1-Ek_e);
        k2_ve1 = k2_ve1 - Gm_e*aqe_1*(ave_1-Ek_e);
        k2_ve1 = k2_ve1 - Gp_e*pinfe(ave_1)*(ave_1-Ena_e);
        k2_ve1 = k2_ve1 - Gh_e*are_1*(ave_1-Eh_e);
        k2_ve1 = k2_ve1 - Gee*ase_1*(ave_1-Eex);
        k2_ve1 = k2_ve1 - Gie*asi_1*(ave_1-Ein);
        k2_ve1 = k2_ve1 + finput_e(j+1);
        k2_ve1 = k2_ve1 / C_e;
        k2_he1 = (hinfek(ave_1)-ahe_1)/htauek(ave_1);
        k2_ne1 = (ninfek(ave_1)-ane_1)/ntauek(ave_1);
        k2_qe1 = (qinfek(ave_1)-aqe_1)/qtauek(ave_1);
        k2_re1 = (rinfe(ave_1)-are_1)/rtaue(ave_1);
        k2_se1 = H(ave_1)*(1-ase_1)/taurse_e-ase_1/taudec_e;
        %
        k2_vi1 = Iapp_i - Gl_i*(avi_1-El_i);
        k2_vi1 = k2_vi1 - Gna_i*ahi_1*minfwb(avi_1)^3*(avi_1-Ena_i);
        k2_vi1 = k2_vi1 - Gk_i*ani_1^4*(avi_1-Ek_i);
        k2_vi1 = k2_vi1 - Gm_i*aqi_1*(avi_1-Ek_i);
        k2_vi1 = k2_vi1 - Gei*ase_1*(avi_1-Eex);
        k2_vi1 = k2_vi1 - Gii*asi_1*(avi_1-Ein);
        k2_vi1 = k2_vi1 + finput_i(j+1);
        k2_vi1 = k2_vi1 / C_i;
        k2_hi1 = phi_i*(hinfwb(avi_1)-ahi_1)/htauwb(avi_1);
        k2_ni1 = phi_i*(ninfwb(avi_1)-ani_1)/ntauwb(avi_1);
        k2_qi1 = phi_i*(qinfwb(avi_1)-aqi_1)/qtauwb(avi_1);
        k2_si1 = H(avi_1)*(1-asi_1)/taurse_i-asi_1/taudec_i;
        %
        Ve_1(j+1) = Ve_1(j) + (k1_ve1 + k2_ve1)*dt/2;
        Ve_1(j+1) = Ve_1(j+1) + D_e*dt^0.5*wne_1(j);
        he_1(j+1) = he_1(j) + (k1_he1 + k2_he1)*dt/2;
        ne_1(j+1) = ne_1(j) + (k1_ne1 + k2_ne1)*dt/2;
        qe_1(j+1) = qe_1(j) + (k1_qe1 + k2_qe1)*dt/2;
        re_1(j+1) = re_1(j) + (k1_re1 + k2_re1)*dt/2;
        Se_1(j+1) = Se_1(j) + (k1_se1 + k2_se1)*dt/2;
        %
        Vi_1(j+1) = Vi_1(j) + (k1_vi1 + k2_vi1)*dt/2;
        Vi_1(j+1) = Vi_1(j+1) + D_i*dt^0.5*wni_1(j);
        hi_1(j+1) = hi_1(j) + (k1_hi1 + k2_hi1)*dt/2;
        ni_1(j+1) = ni_1(j) + (k1_ni1 + k2_ni1)*dt/2;
        qi_1(j+1) = qi_1(j) + (k1_qi1 + k2_qi1)*dt/2;
        Si_1(j+1) = Si_1(j) + (k1_si1 + k2_si1)*dt/2;
        
    end
    
end

return

% EOF
