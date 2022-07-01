% NR_syntransmit        synaptic transmission simulation (LIF/HH model)
%
% call          [ st1, V, t, S, avars ] = NR_syntransmit( st )
% 
% gets          st              spike times; [ms]
%                               three options supported:
%                               (1) N-element cell array: each cell considered a different pre-synaptic train
%                               (2) vector: considered as a single spike train
%                               (3) two-column matrix: each row is considered as clu/res pair
%
% optional arguments (name/value pairs):
%
%               model           {'LIF'}, 'Inap', 'I-cell', 'E-cell', 'LIF_C'
%               dt              {0.1};  [ms]    0.1   	for LIF, Inap, LIC_C, or noiseless I-cell
%                                               0.02    for noisy I-cell or any E-cell 
%               C               {1};    [muF/cm^2]
%               synmodel        {'null'}, or 'DF', 'D', 'F'
%
%               Dn              {2};    [SD];       of noise
%               Tmax            {1000}; [ms]; should match Iext (if given)
%               Ibias           {};     [muA/cm^2]
%                                   Inap:       -1.8
%                                   LIF:        0       (@ Dn=0, Ibias >= 5: repetitive spont. spikes)
%                                   I-cell:     -2.0    (@ Dn=0, Ibias >= -0.53: onset spike; >= -0.5: repetitive spont. spikes)
%                                   E-cell:     -2.7    (@ Dn=0, Ibias >= -2.51: onset spike; >= -2.3: repetitive spont. spikes)
%               Iext            {[]};   [muA/cm^2]
%                                   if empty:       Iapp = Ibias
%                                   if specified:   Iapp = Ibias + Iext
%                                               must match Tmax (can be specified) and dt (hard coded internally)
%               calcIC          {0}
%               graphics        {1} 
%               verbose         {1}
%
%           LIF parameters:
%
%               Gl_LIF          {0.5}   [mS/cm^2]   leak conductance
%               El_LIF          {-60}   [mV]        leak reversal potential
%               Vth_LIF      	{-50}   [mV]        spike threshold
%               Vreset_LIF      {-70}   [mV]        post-spike reset
%               Vpeak_LIF       {50}    [mV]        spike peak (for plotting the spike manually)
%
%           I_Na,h parameters:
%
%               Gl_Inap         {0.1}   [mS/cm^2]   leak conductance
%               El_Inap         {-65}   [mV]        leak reversal potential
%               Vth_Inap     	{-50}   [mV]        spike threshold
%               Vreset_Inap     {-70}   [mV]        Vm post-spike reset
%               Vpeak_Inap      {50}    [mV]        spike peak (for plotting the spike manually)
%
%               Gp_Inap         {0.1}   [mS/cm^2]   persistent sodium conductance
%               Gh_Inap         {1}     [mS/cm^2]   h-current conductance
%               rreset_Inap     {0}     []          h-current post-spike reset
%
%           LIF_C parameters
%
%               Gc_LIF_C        {0.08}  [mS/cm^2]   Calcium conductance
%
%           spike dynamics parameters:
%
%               SD_In           {0.1}   [ms]        duration of every input spike
%               SD_Out          {0.1}   [ms]        duration of every output spike (relevant to LIF model only)
%               ARP             {1}     [ms]        duration of refractory period (relevant to LIF model only)
%
%           connectivity parameters - chemical synapses:
%
%               tau_rise_S      {0.1}   [ms]        rise time for synaptic activation 
%               tau_fall_S      {3.0}   [ms]        fall time for synaptic activation 
%               Es              {0}     [mV]        synaptic reversal potential 
%               Gs              {0.12}  [mS/cm^2]   synaptic conductance
%
%           connectivity parameters - plasticity of chemical synapses:
%
%               tau_rise_D      {0.1}   [ms]        rise time for synaptic depression 
%               tau_fall_D      {100}   [ms]        fall time for synaptic depression 
%               tau_rise_F      {0.2}   [ms]        rise time for synaptic facilitation 
%               tau_fall_F      {300}   [ms]        fall time for synaptic facilitation
%
%               all connectivity parameters can be set at the level of
%               individual interaction (synapses). if that is desired, the size of the
%               array should correspond to N. 
%
% returns       st1             output spike train; [ms]
%               V               membrane potential of the post-synaptic cell; [mV]
%               t               time vector; [ms]
%               S               synaptic activation parameters (for each synapse); [prob]
%               avars           activation variables; [prob] 
%                                   LIF:        none
%                                   Inap:       r
%                                   I-cell:     h, n
%                                   E-cell:     h, n, r
%                                   LIF_C:      c
%
% note          E-cell runs are slower. dt is: 
%               	Inap, LIF:          dt = 0.1 ms
%                   I-cell, Dn <= 5:    dt = 0.1 ms 
%                   I-cell, Dn > 5:     dt = 0.02 ms 
%                   E-cell:             dt = 0.02 ms
%
% calls         ParseArgPairs
%               plot_raster, plotTraces
%
% see also      NR_einet

%-----------------------------------------------------------------------
% Post-synaptic cellular models:
%
% (0) 'LIF'
% with noise and refractory period (1D + reset)
% 
% C dV/dt = I_ext - gL * (V-E_L) - gS * S * (V-Es) + gN * eta
%
% if V > Vth, then V = Vreset (for spike duration, then ARP)
%
%-----------------------------------------------------------------------
% (1) 'Inap'
% with noise and spiking mechanism (2D + reset)
%
% C dV/dt = I_ext - gL * (V-E_L) - gS * S * (V-Es) + gN * eta
%           - g_P * pinf( V ) * ( V - E_Na ) 
%           - g_h * r * ( V - E_h )
% dr/dt = ( rinf( V )- r ) / rtau
%
% if V > Vth, then V = Vreset; r = rreset (for spike duration, then ARP)
%
%-----------------------------------------------------------------------
% (2) 'I-cell'
% with noise and dynamics on V, h, and n (3D)
% Wang and Buzsaki, 1996, JNS
%
% C dV/dt = I_ext - gL * (V-E_L) - gS * S * (V-Es) + gN * eta
%           - g_Na * h * minf^3 * ( V - E_Na )
%           - g_K * n^4 * ( V - E_K )
% dh/dt = ( hinf( V )- h ) / htau( V )
% dn/dt = ( ninf( V )- n ) / ntau( V )
%
%-----------------------------------------------------------------------
% (3) 'E-cell'
% with noise and dynamics on V, h, n, and h-current (4D)
% basic model:
% Olufsen et al., 2003, J. Comp. Neurosci.
% h-current:
% Poolos et al., 2002, Nat. Neurosci. 
% Zemankovics et al., 2010, J. Physiol
%
% C dV/dt = I_ext - gL * (V-E_L) - gS * S * (V-Es) + gN * eta
%           - g_Na * h * minf^3 * ( V - E_Na )
%           - g_K * n^4 * ( V - E_K )
%           - g_h * r * ( V - E_h )
% dh/dt = ( hinf( V )- h ) / htau( V )
% dn/dt = ( ninf( V )- n ) / ntau( V )
% dr/dt = ( rinf( V )- r ) / rtau( V )
%
%-----------------------------------------------------------------------
% (4) 'LIF_C'
% with noise and refractory period, and voltage-dependent calcium dynamics (3D + reset)
% Stark et al, 2022, PLOS Comput. Biol.
% 
% C dV/dt = I_ext - gL * (V-E_L) - gCa * c * (V-E_Ca) - gS * S * (V-Es) + gN * eta
%   dc/dt =  Nc * ( 1 - c ) / tau_act_c - c / tau_inact_c
%  dNc/dt = -Nc / tau_dec_n
%
% if V > Vth, then V = Vreset; N = N_reset (for spike duration, then ARP)
%
%-----------------------------------------------------------------------
% Chemical synaptic model:
% Borgers et al., 2012, PLOS Comput. Biol.
%
% dS/dt = H(Vpre) * (1-S)/tau_rise_S - S/tau_fall_S
% H(V) = (1+tanh(V - Vth)/4)/2;
% 
% possibly different parameters (tau_rise_S, tau_fall_S, Gs, Es) for different synapses
%
%-----------------------------------------------------------------------
% Synaptic plasticity model:
% Stark et al, 2022, PLOS Comput. Biol.
%
% dS/dt = H(Vpre) * (1-S)/tau_rise_S - S/tau_fall_S
% dD/dt = -H(Vpre) * D/tau_rise_D - (1-D)/tau_fall_D
% dF/dt = H(Vpre) * (1-F)/tau_rise_F - F/tau_fall_F
% H(V) = (1+tanh(V - Vth)/4)/2;
% 
% C dV/dt = I_ext - gL*(V-E_L) - gS*S*D*F*(V-Es) + eta
% 
% possibly different parameters (tau_rise_D/F, tau_fall_D/F) for different synapses
%
%-----------------------------------------------------------------------

% 24-aug-20 ES

% last update:
% 30-jun-22

function [ st1, V, t, S, avar ] = NR_syntransmit( st, varargin )

% [C] = muF/cm^2
% [I] = muA/cm^2
% [G] = mS/cm^2
% [t] = ms
% [V] = mV

%-------------------------------------------------------------------------
% constants
%-------------------------------------------------------------------------
% time steps required for stability
MAX_DT_INAP                     = 0.1;                                      % [ms]
MAX_DT_LIF                      = 0.1;                                      % [ms]
MAX_DT_I_CELL_WITH_NOISE        = 0.02;                                     % [ms]; high noise - above 5 mV
MAX_DT_E_CELL                   = 0.02;                                     % [ms]

% membrane potential parameters - pre synaptic cell
% used only for converting input spike trains to presynaptic potentials
El_pre                          = -60;                                      % [mV], to "analogize" presynaptic spikes
Vpeak_pre                       = 50;                                       % [mV], to "analogize" presynaptic spikes
Vth_pre                         = 0;                                        % [mV], sets inflection point of presynaptic sigmoid

%-------------------------------------------------------------------------
% default parameters
%-------------------------------------------------------------------------
% timing parameters
dt_DFLT                         = 0.1;                                      % [ms]; E-cell model requires 0.02 or lower to maintain stability
TMAX_DFLT                       = 1000;                                     % [ms]

% default parameters
C_DFLT                          = 1.0;                                      % [muF/cm^2]

% input parameters 
SD_In_DFLT                      = 0.1;                                      % [ms], duration of each input spike

% default noise level
Dn_DFLT                         = 2;                                        % [mV]

% default synaptic parameters
tau_rise_S_DFLT                 = 0.1;                                      % [ms], AMPA
tau_fall_S_DFLT                 = 3.0;                                      % [ms], AMPA
Es_DFLT                         = 0.0;                                      % [mV], AMPA

% default connectivity parameters
Gs_DFLT                         = 0.12;                                     % [mS/cm^2]; from pre to post

% default synaptic plasticity parameters
tau_rise_D_DFLT                 = 0.1;                                      % [ms]; depression rise time
tau_fall_D_DFLT                 = 100;                                      % [ms]; depression decay time
tau_rise_F_DFLT                 = 0.2;                                      % [ms]; facilitation rise time
tau_fall_F_DFLT                 = 300;                                      % [ms]; facilitation decay time

% RP (for LIF model only)
ARP_DFLT                        = 1;                                        % [ms]; ARP
SD_Out_DFLT                    	= 0.1;                                      % [ms]; duration of each output spike

% LIF parameters
Gl_LIF_DFLT                     = 0.5;                                      % [mS/cm^2]; or 0.1; leak conductance
El_LIF_DFLT                     = -60;                                      % [mV]; leak reversal potential
Vth_LIF_DFLT                    = -50;                                      % [mV]; spike threshold
Vreset_LIF_DFLT                 = -70;                                      % [mV]; or -60; post-spike reset
Vpeak_LIF_DFLT                  = 50;                                       % [mV]; spike peak (for plotting the spike manually)

% Inap parameters        
Gl_Inap_DFLT                    = 0.1;                                      % [mS/cm^2]; or 0.1; leak conductance
El_Inap_DFLT                    = -65;                                      % [mV]; leak reversal potential
Vth_Inap_DFLT                   = -50;                                      % [mV]; spike threshold
Vreset_Inap_DFLT                = -70;                                      % [mV]; or -60; post-spike reset
Vpeak_Inap_DFLT                 = 50;                                       % [mV]; spike peak (for plotting the spike manually)
Gp_Inap_DFLT                    = 0.1;                                      % [mS/cm^2]; persistent sodium conductance
Gh_Inap_DFLT                    = 1;                                        % [mS/cm^2]; h-current conductance
rreset_Inap_DFLT                = 0;                                        % h-current gating variable post-spike reset

% LIF_C parameters
Gc_LIF_C_DFLT                   = 0.08;                                     % [mS/cm^2]; Calcium conductance

%-------------------------------------------------------------------------
% functions
%-------------------------------------------------------------------------
% presynaptic sigmoid
H                               = @(v) ( 1 + tanh( ( v - Vth_pre ) / 4 ) ) / 2;

% Inap 
pinfnap                         = @(v) 1.0 ./ ( 1 + exp( -(v+38) / 6.5 ) );
rinfnap                         = @(v) 1.0 ./ ( 1 + exp( (v+79.2) / 9.78 ) ); 

% I-cell
minfwb                          = @(v) (-0.2*(v+35)./(exp(-0.1*(v+35))-1)) ./((-0.2*(v+35)./(exp(-0.1*(v+35))-1))+(4*exp(-(v+60)/18)));
hinfwb                          = @(v) (0.07*(exp(-(v+58)/20))) ./((0.07*(exp(-(v+58)/20)))+(1.0./(exp(-0.1*(v+28))+1)));
htauwb                          = @(v) 1.0 ./((0.07*(exp(-(v+58)/20)))+(1.0./(exp(-0.1*(v+28))+1)));
ninfwb                          = @(v) (-0.01*(v+34)./(exp(-0.1*(v+34))-1)) ./((-0.01*(v+34)./(exp(-0.1*(v+34))-1))+(0.125*exp(-(v+44)/80)));
ntauwb                          = @(v) 1.0 ./((-0.01*(v+34)./(exp(-0.1*(v+34))-1))+(0.125*exp(-(v+44)/80)));

% E-cell
minfek                          = @(v) (0.32*(54+v)./(1-exp(-(v+54)/4)))/((0.32*(54+v)./(1-exp(-(v+54)/4)))+(0.28*(v+27)./(exp((v+27)/5)-1)));
hinfek                          = @(v) (0.128*exp(-(50+v)/18)) ./((0.128*exp(-(50+v)/18))+(4./(1+exp(-(v+27)/5))));
htauek                          = @(v) 1.0 ./((0.128*exp(-(50+v)/18))+(4./(1+exp(-(v+27)/5))));
ninfek                          = @(v) (0.032*(v+52)./(1-exp(-(v+52)/5))) ./ ((0.032*(v+52)./(1-exp(-(v+52)/5)))+(0.5*exp(-(57+v)/40)));
ntauek                          = @(v) 1.0 ./ ((0.032*(v+52)./(1-exp(-(v+52)/5)))+(0.5*exp(-(57+v)/40)));
rinfek                          = @(v) 1./(1 + exp((v+82.9)/12.4));
rtauek                          = @(v) 1.5*(exp(0.033*(v+75))./(0.011*(1+exp(0.083*(v+75)))));

%-------------------------------------------------------------------------
% arguments
%-------------------------------------------------------------------------
nargs                           = nargin;

% spike train(s)
if nargs < 1 || isempty( st )                                               % no input
    N                           = 1;
    st                          = { [] };
elseif isa( st, 'cell' )                                                    % cell array: each cell considered a different pre-synaptic train
    N                        = length( st );
elseif isa( st, 'double' ) || isa( st, 'single' )
    if min( size( st ) ) == 1
        st                      = st( : );
    end
    ncols                       = size( st, 2 );
    if ncols == 1                                                           % vector: considered as a single spike train
        N                       = 1;
        st                      = { st };
    elseif ncols == 2                                                       % two-column matrix: considered as clu/res pair
        clu                     = st( :, 1 );
        res                     = st( :, 2 );
        uclu                    = unique( clu );
        N                       = length( uclu );
        st                      = cell( 1, N );
        for i                   = 1 : N
            st{ i }             = res( clu == uclu( i ) );
        end
    end
end

% other arguments
[ model, dt, C, synmodel ...
    , Gl_LIF, El_LIF, Vth_LIF, Vreset_LIF, Vpeak_LIF ...
    , Gl_Inap, El_Inap, Vth_Inap, Vreset_Inap, Vpeak_Inap ...
    , Gp_Inap, Gh_Inap, rreset_Inap ...
    , Gc_LIF_C ...
    , tau_rise_S, tau_fall_S, Es, Gs ...
    , tau_rise_D, tau_fall_D, tau_rise_F, tau_fall_F ...
    , Dn ...
    , Iext, Ibias ...
    , SD_Out, ARP, calcIC, Tmax, SD_In ... 
    , maxGB, graphics, verbose ]         = ParseArgPairs(...
    { 'model', 'dt', 'C', 'synmodel' ...
    , 'Gl_LIF', 'El_LIF', 'Vth_LIF', 'Vreset_LIF', 'Vpeak_LIF' ...
    , 'Gl_Inap', 'El_Inap', 'Vth_Inap', 'Vreset_Inap', 'Vpeak_Inap' ...
    , 'Gp_Inap', 'Gh_Inap', 'rreset_Inap' ...
    , 'Gc_LIF_C' ...
    , 'tau_rise_S', 'tau_fall_S', 'Es', 'Gs' ... 
    , 'tau_rise_D', 'tau_fall_D', 'tau_rise_F', 'tau_fall_F' ...
    , 'Dn' ...
    , 'Iext', 'Ibias' ...
    , 'SD_Out', 'ARP', 'calcIC', 'Tmax', 'SD_In' ...
    , 'maxGB', 'graphics', 'verbose' }...
    , { 'LIF', dt_DFLT, C_DFLT, 'null' ...
    , Gl_LIF_DFLT, El_LIF_DFLT, Vth_LIF_DFLT, Vreset_LIF_DFLT, Vpeak_LIF_DFLT ...
    , Gl_Inap_DFLT, El_Inap_DFLT, Vth_Inap_DFLT, Vreset_Inap_DFLT, Vpeak_Inap_DFLT ...
    , Gp_Inap_DFLT, Gh_Inap_DFLT, rreset_Inap_DFLT ...
    , Gc_LIF_C_DFLT ...
    , tau_rise_S_DFLT, tau_fall_S_DFLT, Es_DFLT, Gs_DFLT ...
    , tau_rise_D_DFLT, tau_fall_D_DFLT, tau_rise_F_DFLT, tau_fall_F_DFLT ...
    , Dn_DFLT ...
    , [], [] ...
    , SD_Out_DFLT, ARP_DFLT, 0, TMAX_DFLT, SD_In_DFLT ... 
    , 1, 1, 1 } ...
    , varargin{ : } );

if ~ismember( model, { 'LIF', 'Inap', 'I-cell', 'E-cell', 'LIF_C' } )
    error( 'unrecognized cellular model' )
end
switch model
    case 'LIF'
        modeln                  = 0;
    case 'Inap'
        modeln                  = 1;
    case 'I-cell'
        modeln                  = 2;
    case 'E-cell'
        modeln                  = 3;
    case 'LIF_C'
        modeln                  = 4;
end
        
if ~ismember( synmodel, { 'null', 'DF', 'D', 'F' } )
    error( 'unrecognized synaptic plasticity model' )
end

%-------------------------------------------------------------------------
% parameters
%-------------------------------------------------------------------------

% input spike train parameters
if N > 1 && length( SD_In ) == 1
    SD_In                       = repmat( SD_In, [ 1 N ] );
end
if ~isequal( size( SD_In ), [ 1 N ] )
    error( 'SD_In: input size mismatch' )
end

% synaptic parameters
if N > 1 && length( tau_rise_S ) == 1
    tau_rise_S                  = repmat( tau_rise_S, [ 1 N ] );
end
if N > 1 && length( tau_fall_S ) == 1
    tau_fall_S                  = repmat( tau_fall_S, [ 1 N ] );
end
if N > 1 && length( Es ) == 1
    Es                          = repmat( Es, [ 1 N ] );
end
if ~isequal( size( tau_rise_S ), [ 1 N ] )
    error( 'tau_rise_S: input size mismatch' )
end
if ~isequal( size( tau_fall_S ), [ 1 N ] )
    error( 'tau_fall_S: input size mismatch' )
end
if ~isequal( size( Es ), [ 1 N ] )
    error( 'Es: input size mismatch' )
end

% synaptic plasticity parameters
switch synmodel
    case 'DF'
        synp                    = 1;
        F0                      = 0;
    case 'D'
        synp                    = 1;
        F0                      = 1;
        tau_rise_F              = inf;
        tau_fall_F              = inf;
    case 'F'
        synp                    = 1;
        F0                      = 0;
        tau_rise_D              = inf;
        tau_fall_D              = inf;
    otherwise
        synp                    = 0;
        tau_rise_F              = inf;
        tau_fall_F              = inf;
        tau_rise_D              = inf;
        tau_fall_D              = inf;
end
if N > 1 && length( tau_rise_D ) == 1
    tau_rise_D                  = repmat( tau_rise_D, [ 1 N ] );
end
if N > 1 && length( tau_fall_D ) == 1
    tau_fall_D                  = repmat( tau_fall_D, [ 1 N ] );
end
if N > 1 && length( tau_rise_F ) == 1
    tau_rise_F                  = repmat( tau_rise_F, [ 1 N ] );
end
if N > 1 && length( tau_fall_F ) == 1
    tau_fall_F                  = repmat( tau_fall_F, [ 1 N ] );
end
if ~isequal( size( tau_rise_D ), [ 1 N ] )
    error( 'tau_rise_D: input size mismatch' )
end

if ~isequal( size( tau_fall_D ), [ 1 N ] )
    error( 'tau_fall_D: input size mismatch' )
end

if ~isequal( size( tau_rise_F ), [ 1 N ] )
    error( 'tau_rise_F: input size mismatch' )
end

if ~isequal( size( tau_fall_F ), [ 1 N ] )
    error( 'tau_fall_F: input size mismatch' )
end

% connectivity parameters
if N > 1 && length( Gs ) == 1
    Gs                          = repmat( Gs, [ 1 N ] );
end
if ~isequal( size( Gs ), [ 1 N ] )
    error( 'Gs: input size mismatch' )
end
if any( Gs( : ) )
    doS                         = 1;
else
    doS                         = 0;
end

% membrane potential parameters - post synaptic cell
Ibias_model                     = 0;
switch modeln
    
    case 1                                                                  % 'Inap'
        
        nvars                   = 2;

        Gl                      = Gl_Inap;                                  % leak conductance
        El                      = El_Inap;                                  % leak reversal potential
        
        Gp                      = Gp_Inap;
        Ena                     = 55;                                       % [mV]

        Gh                      = Gh_Inap;
        Eh                      = -20;                                      % [mV]
        rtaunap                 = 100;                                      % [ms]; fixed (not voltage dependent)
        
        Vth                     = Vth_Inap;                                 % spike threshold
        Vreset                  = Vreset_Inap;                           	% post-spike reset
        rreset                  = rreset_Inap;                              % post-spike reset of h-current 
        Vpeak                   = Vpeak_Inap;                             	% spike peak (for plotting the spike manually)

        Ibias_model             = -1.8;                                     % [muA/cm^2]

        if dt > MAX_DT_INAP
            fprintf( 1, 'Inap will be unstable with dt=%0.3g, reducing to %0.3g\n', dt, MAX_DT_INAP )
            dt                  = MAX_DT_INAP;
        end
    
    case { 0, 4 }                                                           % 'LIF', 'LIF_C'

        nvars                   = 1;

        Gl                      = Gl_LIF;                               	% leak conductance
        El                      = El_LIF;                                   % leak reversal potential
        
        Vth                     = Vth_LIF;                                  % spike threshold
        Vreset                  = Vreset_LIF;                           	% post-spike reset
        Vpeak                   = Vpeak_LIF;                             	% spike peak (for plotting the spike manually)

        if modeln == 4                                                      % 'LIF_C'
            Gc                  = Gc_LIF_C;                                 % Ca conductance
            Ec                  = 100;                                      % [mV]; Ca reversal potential 
            tau_act_c           = 50;                                       % [ms]; Ca rise time constant
            tau_inact_c         = 5;                                        % [ms]; Ca fall time constant
            tau_dec_Nc          = 70;                                       % [ms]; Nc fall time constant
            Nc_reset            = 0.1;                                      % reset value of Nc (gating variable)
        end
        
        if dt > MAX_DT_LIF
            fprintf( 1, 'LIF will be unstable with dt=%0.3g, reducing to %0.3g\n', dt, MAX_DT_LIF )
            dt                  = MAX_DT_LIF;
        end
        
    case 2                                                                  % 'I-cell'

        nvars                   = 3;

        Gl                      = 0.1;                                      % [mS/cm^2]; leak conductance
        El                      = -65;                                      % [mV]; leak reversal potential

        Gna                     = 35;                                       % [mS/cm^2]; Sodium conductance
        Ena                     = 55;                                       % [mV]; Sodium reversal potential

        Gk                      = 9;                                        % [mS/cm^2]; Potassium conductance
        Ek                      = -90;                                      % [mV]; Potassium reversal potential

        Ibias_model             = -2.0;                                     % [muA/cm^2]; bias current (~exact F/noise curve of the LIF model: spikes start @ Dn=2.2, 12 spks/s @ Dn=4)
        
        phi                     = 5;                                        % time scale for activation variables
        
        if Dn > 5 && dt > MAX_DT_I_CELL_WITH_NOISE
            fprintf( 1, 'I-cell will be unstable with dt=%0.3g, reducing to %0.3g\n', dt, MAX_DT_I_CELL_WITH_NOISE )
            dt                  = MAX_DT_I_CELL_WITH_NOISE;
        end
        
    case 3                                                                  % 'E-cell'
        
        nvars                   = 4;
        
        Gl                      = 0.1;                                      % [mS/cm^2]; leak conductance
        El                      = -67;                                      % [mV]; leak reversal potential

        Gna                     = 100;                                      % [mS/cm^2]; Sodium conductance
        Ena                     = 50;                                       % [mV]; Sodium reversal potential

        Gk                      = 80;                                       % [mS/cm^2]; Potassium conductance
        Ek                      = -100;                                     % [mV]; Potassium reversal potential

        Gh                      = 0.485;                                    % [mS/cm^2]; h-current conductance 
        Eh                      = -33;                                      % [mV]; h-current reversal potential
        
        Ibias_model             = -2.7;
        
        if dt > MAX_DT_E_CELL
            fprintf( 1, 'E-cell will be unstable with dt=%0.3g, reducing to %0.3g\n', dt, MAX_DT_E_CELL )
            dt                  = MAX_DT_E_CELL;
        end
        
end % modeln

if length( Ibias ) ~= 1
    Ibias                       = Ibias_model;
end

%-------------------------------------------------------------------------
% input handling
%-------------------------------------------------------------------------
% spike trains
nspikes_pre                     = zeros( N, 1 );
tmax_pre                        = zeros( N, 1 );
for i                           = 1 : N
    st{ i }                     = st{ i }( : );
    nspikes_pre( i )            = length( st{ i } );
    if nspikes_pre( i ) 
        tmax_pre( i )           = max( st{ i } );
    end
end
Tmax                            = max( max( tmax_pre( : ).' + SD_In ), Tmax );

% time and external input vectors
nt                              = Tmax / dt;
nGB                             = ( 2 * synp * N + N + nvars ) * nt * 4 / 2^30;         	% single (4 bytes/element)
if nGB > maxGB
    fprintf( 1, 'RAM required: %0.2g GB; this is more than the maximal allowed (%0.2g GB); exiting\n', nGB, maxGB )
    return
end

t                               = ( dt : dt : Tmax )';                      % time vector
nt                              = length( t );
dtSD_In                         = ceil( SD_In / dt );                       % used for analogization of all input spikes
dtSD_Out                        = ceil( SD_Out / dt );                      % LIF model only
dtARP                           = ceil( ARP / dt );                         % LIF model only
Iext                            = Iext( : );
if numel( Iext ) ~= nt
    Iext                        = [];
    if verbose
        fprintf( 1, 'Iext ignored; using only bias current (Ibias: %0.2g)\n', Ibias )
    end
end
if isempty( Iext )
    Iext                        = ones( nt, 1 ) * Ibias;
else
    Iext                        = Iext + Ibias;
end

% membrane potentials of presynaptic cells
Vpre                            = ones( nt, N ) * El_pre;
for i                           = 1 : N
    if isempty( st{ i } )
        continue
    end
    st{ i }                     = ceil( st{ i } / dt );
    idx                         = st{ i };
    idx                         = bsxfun( @plus, idx, ( 1 : dtSD_In( i ) ) - 1 );
    Vpre( idx, i )              = Vpeak_pre;
end

% feedforward net, so do not consider coupling if no input spikes
if nspikes_pre == 0
    doS                         = 0;
end

%-------------------------------------------------------------------------
% initialization:
%-------------------------------------------------------------------------
% variables
V                               = zeros( nt, 1 );
S                               = zeros( nt, N );
if synp
    D                           = zeros( nt, N );
    F                           = zeros( nt, N );
end
switch modeln
    case 1                                                                  % 'Inap'
        r                       = zeros( nt, 1 );
    case 2                                                                  % 'I-cell'
        h                       = zeros( nt, 1 );
        n                       = zeros( nt, 1 );
    case 3                                                                  % 'E-cell'
        h                       = zeros( nt, 1 );
        n                       = zeros( nt, 1 );
        r                       = zeros( nt, 1 );
    case 4                                                                  % 'LIF_C'
        c                       = zeros( nt, 1 );
        Nc                      = zeros( nt, 1 );
end
        
% initial conditions:
V( 1 )                          = El;                                       % steady state
if doS
    S( 1, : )                   = 0;                                        % no synaptic conductance
    if synp
        D( 1, : )              	= 1;                                        % no depression
        F( 1, : )           	= F0;                                       %
    end
end

% LIF_C model: start with a spike
if modeln == 4 && Gc > 0                                                    % 'LIF_C'
    V( 1 )                      = Vth + 1;
end

% E-cell model: compute I.C. from a noiseless run
if calcIC
    if modeln == 3                                                          % 'E-cell'
        [ ~, V0, ~, ~, avar0 ]  = NR_syntransmit( {}, 'model', model ...
            , 'tau_rise_S', tau_rise_S( 1 ), 'tau_fall_S', tau_fall_S( 1 ), 'Es', Es( 1 ), 'Gs', Gs( 1 ) ...
            , 'Dn', 0, 'calcIC', 0, 'graphics', 0, 'Tmax', 100 );
        V( 1 )                  = V0( end );
        h( 1 )                  = avar0( end, 1 );
        n( 1 )                  = avar0( end, 2 );
        r( 1 )                  = avar0( end, 3 );
    end
else
    if modeln == 3                                                          % 'E-cell'
        V( 1 )                  = -65.559;
        h( 1 )                  = 0.9931;
        n( 1 )                  = 0.0473;
        r( 1 )                  = 0.1874;
    end
end

% noise:
wn                              = randn( nt, 1 );

%-------------------------------------------------------------------------
% simulation:
%-------------------------------------------------------------------------

newspk                          = 0;                                        % LIF model

for i                           = 1 : ( nt - 1 )
    
    % slope at i
    k1v                         = Iext( i ) - Gl * ( V( i ) - El );         % external and leak currents
    if modeln > 0
        if modeln == 1                                                      % 'Inap'
            k1v                 = k1v - Gp * pinfnap( V( i ) ) * ( V( i ) - Ena );
            k1v                 = k1v - Gh * r( i ) * ( V( i ) - Eh );
            k1r                 = ( rinfnap( V( i ) ) - r( i ) ) / rtaunap;
        elseif modeln == 2                                                  %'I-cell'
            k1v                 = k1v - Gna * h( i ) * minfwb( V( i ) )^3 * ( V( i ) - Ena);
            k1v                 = k1v - Gk * n( i )^4 * ( V( i ) - Ek );
            k1h                 = phi * ( hinfwb( V( i ) ) - h( i ) ) / htauwb( V( i ) );
            k1n                 = phi * ( ninfwb( V( i ) ) - n( i ) ) / ntauwb( V( i ) );
        elseif modeln == 3                                                  % 'E-cell'
            k1v                 = k1v - Gna * h( i ) * minfek( V( i ) )^3 * ( V( i ) - Ena );
            k1v                 = k1v - Gk * n( i )^4 * ( V( i ) - Ek );
            k1v                 = k1v - Gh * r( i ) * ( V( i ) - Eh );
            k1h                 = ( hinfek( V( i ) ) - h( i ) ) / htauek( V( i ) );
            k1n                 = ( ninfek( V( i ) ) - n( i ) ) / ntauek( V( i ) );
            k1r                 = ( rinfek( V( i ) ) - r( i ) ) / rtauek( V( i ) );
        elseif modeln == 4                                                  % 'LIF_C'
            k1v                 = k1v - Gc * c( i ) * ( V ( i ) - Ec );    
            k1c                 = Nc( i ) * ( 1 - c( i ) ) / tau_act_c - c( i ) / tau_inact_c;
            k1nc                = -Nc( i ) / tau_dec_Nc;
        end
    end
    if doS
        if synp
            SDF                	= S( i, : ) .* D( i, : ) .* F( i, : );
        else
            SDF                 = S( i, : );
        end
        k1v                     = k1v - sum( Gs .* SDF .* ( V( i ) - Es ) );      % assume linear summation of multiple synaptic currents
    end
    k1v                         = k1v / C;                                  % capacitance
    if doS
        k1s                  	= H( Vpre( i, : ) ) .* ( 1 - S( i, : ) ) ./ tau_rise_S - S( i, : ) ./ tau_fall_S;
        if synp
            k1d              	= -H( Vpre( i, : ) ) .* D( i, : ) ./ tau_rise_D + ( 1 - D( i, : ) ) ./ tau_fall_D;
            k1f               	= H( Vpre( i, : ) ) .* ( 1 - F( i, : ) ) ./ tau_rise_F - F( i, : ) ./ tau_fall_F;
        end
    end
    
    % value at endpoint ( i + 1 ), based on slope at i
    av1                         = V( i ) + k1v * dt;
    av1                         = av1 + Dn * dt ^ 0.5 * wn( i );
    if doS
        as1                     = S( i, : ) + k1s * dt;
        if synp
            ad1             	= D( i, : ) + k1d * dt;
            af1              	= F( i, : ) + k1f * dt;
        end
    end
    if modeln > 0
        if modeln == 1                                                      % 'Inap'
            ar1                 = r( i ) + k1r * dt;
        elseif modeln == 2                                                  % 'I-cell'
            ah1                 = h( i ) + k1h * dt;
            an1                 = n( i ) + k1n * dt;
        elseif modeln == 3                                                  % 'E-cell'
            ah1                 = h( i ) + k1h * dt;
            an1                 = n( i ) + k1n * dt;
            ar1                 = r( i ) + k1r * dt;
        elseif modeln == 4                                                  % 'LIF_C'
            ac1                 = c( i ) + k1c * dt;
            anc1                = Nc( i ) + k1nc * dt;
        end
    end
    
    % slope at ( i + 1 )
    k2v                         = Iext( i + 1 ) - Gl * ( av1 - El );    	% external and leak currents
    if modeln > 0
        if modeln == 1                                                      % 'Inap'
            k2v                 = k2v - Gp * pinfnap( av1 ) * ( av1 - Ena );
            k2v                 = k2v - Gh * ar1 * ( av1 - Eh );
            k2r                 = ( rinfnap( av1 ) - ar1 ) / rtaunap;
        elseif modeln == 2                                                  % 'I-cell'
            k2v                 = k2v - Gna * ah1 * minfwb( av1 )^3 * ( av1 - Ena );
            k2v                 = k2v - Gk * an1^4 * ( av1 - Ek );
            k2h                 = phi * ( hinfwb( av1 ) - ah1 ) / htauwb( av1 );
            k2n                 = phi * ( ninfwb( av1 ) - an1 ) / ntauwb( av1 );
        elseif modeln == 3                                                  % 'E-cell'
            k2v                 = k2v - Gna * ah1 * minfek( av1 )^3 * ( av1 - Ena );
            k2v                 = k2v - Gk * an1^4 * ( av1 - Ek );
            k2v                 = k2v - Gh * ar1 * ( av1 - Eh );
            k2h                 = ( hinfek( av1 ) - ah1 ) / htauek( av1 );
            k2n                 = ( ninfek( av1 ) - an1 ) / ntauek( av1 );
            k2r                 = ( rinfek( av1 ) - ar1 ) / rtauek( av1 );
        elseif modeln == 4                                                  % 'LIF_C'
            k2v                 = k2v - Gc * ac1 * ( av1 - Ec );
            k2c                 = anc1 * ( 1 - ac1 ) / tau_act_c - ac1 / tau_inact_c;
            k2nc                = -anc1 / tau_dec_Nc;
        end
    end
    if doS
        if synp
            SDF                	= as1 .* ad1 .* af1;
        else
            SDF                 = as1;
        end
        k2v                     = k2v - sum( Gs .* SDF .* ( av1 - Es ) );  	% assume linear summation of multiple synaptic currents
        % should save as1, ad1, af1, and av1 over a buffer of Tsyn samples back
    end
    k2v                         = k2v / C;                               	% capacitance
    if doS
        k2s                  	= H( Vpre( i, : ) ) .* ( 1 - as1 ) ./ tau_rise_S - as1 ./ tau_fall_S;
        if synp
            k2d               	= -H( Vpre( i, : ) ) .* ad1 ./ tau_rise_D + ( 1 - ad1 ) ./ tau_fall_D;
            k2f                	= H( Vpre( i, : ) ) .* ( 1 - af1 ) ./ tau_rise_F - af1 ./ tau_fall_F;
        end
    end
    
    % value at endpoint ( i + 1 ), based on average slope
    V( i + 1 )                  = V( i ) + ( k1v + k2v ) * dt / 2;
    V( i + 1 )                  = V( i + 1 ) + Dn * dt ^ 0.5 * wn( i );
    if doS
        S( i + 1, : )         	= S( i, : ) + ( k1s + k2s ) * dt / 2;
        if synp
            D( i + 1, : )      	= D( i, : ) + ( k1d + k2d ) * dt / 2;
            F( i + 1, : )     	= F( i, : ) + ( k1f + k2f ) * dt / 2;
        end
    end
    if modeln == 0 || modeln == 1 || modeln == 4                            % 'Inap', 'LIF', 'LIF_C'
        if modeln == 1                                                      % 'Inap'
            r( i + 1 )          = r( i ) + ( k1r + k2r ) * dt / 2;
        end
        if modeln == 4                                                      % 'LIF_C'
            c( i + 1 )          = c( i ) + ( k1c + k2c ) * dt / 2;
            Nc( i + 1 )         = Nc( i ) + ( k1nc + k2nc ) * dt / 2;
        end
        if V( i + 1 ) > Vth && newspk == 0                                  % spike now
            V( i + 1 )          = Vpeak;
            newspk              = 1;
            if modeln == 4                                                  % 'LIF_C'
                Nc( i + 1 )     = Nc_reset;
            end
        elseif newspk >= 1 && newspk < dtSD_Out                             % same spike
            V( i + 1 )          = Vpeak;
            newspk              = newspk + 1;
        elseif newspk >= 1 && newspk < dtARP                                % refractory period
            V( i + 1 )          = Vreset;
            if modeln == 1                                                  % 'Inap'
                r( i + 1 )      = rreset;
            end
            newspk              = newspk + 1;
        elseif newspk >= 1 && newspk >= dtSD_Out && newspk >= dtARP         % back to dynamics
            V( i + 1 )          = Vreset;
            if modeln == 1                                                  % 'Inap'
                r( i + 1 )      = rreset;
            end
            newspk              = 0;
        end
    elseif modeln == 2                                                      % 'I-cell'
        h( i + 1 )              = h( i ) + ( k1h + k2h ) * dt / 2;
        n( i + 1 )              = n( i ) + ( k1n + k2n ) * dt / 2;
    elseif modeln == 3                                                      % 'E-cell'
        h( i + 1 )              = h( i ) + ( k1h + k2h ) * dt / 2;
        n( i + 1 )              = n( i ) + ( k1n + k2n ) * dt / 2;
        r( i + 1 )              = r( i ) + ( k1r + k2r ) * dt / 2;
    end
    
end


%-------------------------------------------------------------------------
% prepare output
%-------------------------------------------------------------------------

% organize activation variables
switch model
    case 'Inap'
        avar                    = r;
    case 'LIF'
        avar                    = [];
    case 'LIF_C' 
        avar                    = [ c Nc ];
    case 'I-cell'
        avar                    = [ h n ];
    case 'E-cell'
        avar                    = [ h n r ];
end

% digitize spikes
pst1                            = find( V > 0 );                            % threshold crossings 
mat                             = local_parse( pst1 );
st1                             = ceil( mean( mat, 2 ) ) * dt;              % spike times [ms]

%-------------------------------------------------------------------------
% graphics
%-------------------------------------------------------------------------
if ~graphics
    return
end

%-------------------------------------------------------------------------
% plot

xlims                           = [ 0 max( t ) ];

figure
subplot( 3, 1, 1 )
plot_raster( st, t );
ylabel( 'Unit number' )
title( sprintf( '%d spikes (%0.2g spks/s)', sum( nspikes_pre ), sum( nspikes_pre ) / Tmax * 1000 ) )

subplot( 3, 1, 2 )
plotTraces( t, S, -1, 1, [ NaN NaN ] );
ylabel( 'S' )
ylim( [ 0 N ] )

subplot( 3, 1, 3 )
plot( t, V, '-b' )
ylabel( 'V [mV]' )
xlabel( 't [ms]' )
title( sprintf( '%d spikes (%0.2g spks/s)', length( st1 ), length( st1 ) / Tmax * 1000 ) )

for spi                         = 1 : 3
    subplot( 3, 1, spi )
    set( gca, 'tickdir', 'out', 'box', 'off' );
    xlim( xlims )
end

return % NR_syntransmit

%------------------------------------------------------------------------
% mat = local_parse( vec )
% parse a vector to sets of points, each of consecutive values
%------------------------------------------------------------------------
function mat = local_parse( vec )

if isempty( vec )
    mat                         = [];
    return
end

vec                             = vec( : );
n                               = length( vec );
if n == 1
    mat                         = [ vec vec ];
    return
end

st                              = zeros( n, 1 );
et                              = zeros( n, 1 );
j                               = 1;
st( 1 )                         = vec( 1 );
for i                           = 1 : ( length( vec ) - 1 )
    if ( vec( i ) + 1 ) ~= vec( i + 1 ) && vec( i ) ~= vec( i + 1 )
        et( j )                 = vec( i );
        j                       = j + 1;
        st( j )                 = vec( i + 1 );
    end
end
et( j )                         = vec( i + 1 );
mat                             = [ st( 1 : j ) et( 1 : j ) ];
    
return % local_parse

% EOF