% NR_einet      network of E and I cells: LIF/Inap, all-to-all chemical connectivity
%
% call          [ st, Ve, Vi, See, Sei, Sii, Sie ] = NR_einet( N, Iext )
%
% gets          N               {10}; number of units
%                                   supported formats:
%                                   (1) [ N ]
%                                   (2) [ Ne Ni ]
%               Iext            {0}; [muA/cm^2]; external input and time vectors
%                                   supported formats:
%                                   (1) 3-column matrix: [ Ie Ii t ]
%                                   (2) (N+1)-column matrix: [ [ Ne columns ] [ Ni columns ] t ]
% 
% optional arguments (name/value pairs)
% 
%               fE              {0.8}   fraction of E-cells
%                                           ignored if Ne, Ni defined explicitly 
%               fEprune         {0}     fraction of EE synapses to
%                                           prune (set Gee to 0)
%               Tmax            {2000}  [ms]    simulation time 
%                                           ignored if t defined explicitly via Iext
%               autapseE        {0}     flag: 0 := no autapses between E cells (all diagonal elements of Gee set to zero)
%               autapseI        {1}     flag: 0 := no autapses between I cells (all diagonal elements of Gii set to zero)
%               maxGB           {2}     [GB]    maximal amount of RAM to be used
%               thin            {1}             uses linear (in N_cells) RAM, do not save every sample for the synaptic (quadratic) variables 
%               fast            {1}             uses linear (in N_cells) computations, uses the same time constants 
%                                                   for all synapses made by a single presynaptic neuron
%               graphics        {0}     flag
%               verbose         {0}     flag
%
%           external inputs: 
% 
%               Dne             {2}     [mV] SD of membrane potential noise (G_n = 1 mS/cm^2)
%               Dni             {2}     [mV] SD of membrane potential noise (G_n = 1 mS/cm^2)
%               Ibias_e         {0}     [muA/cm^2]; additive bias current
%               Ibias_i         {0}     [muA/cm^2]; additive bias current
%
%           LIF (E- or I-cell) dynamics parameters:
%
%               C_e, C_i        {1}     [muF/cm^2]; E/I cell capacitance
%               Gl_e, Gl_i      {0.5}   [mS/cm^2];  E/I cell leak conductance
%               El_e, El_i      {-60}   [mV];       E/I cell leak reversal potential
%
%               Vth_e, Vth_i    {-50}   [mV];       E/I cell spike threshold
%               Vreset_e/i      {-70}   [mV];       E/I cell post-spike reset
%               Vpeak_e/i       {50}    [mV];       E/I cell spike peak
%
%           I_Na,p+I_h (E-cell) dynamics parameters:
%
%               Inap         	{0}                 include I_Na,p and I-h in all E cells
%               Gp_e         	{0.1}   [mS/cm^2];  E cell persistent sodium conductance
%               Ena_e           {55}    [mV];       E cell persistent sodium reversal potential 
%               Gh_e          	{1}     [mS/cm^2];  E cell h-current conductance
%               Eh_e            {-20}   [mV];       E cell h-current reversal potential 
%               rtaunap_e       {100}   [ms];       E cell h-current time constant (not voltage dependent)
%               rreset_e        {0}     [mV];       E cell h-current gating variable post-spike reset
%               
%           spike dynamics parameters:
%
%               ARP_e           {1}     [ms]        duration of E-cell ARP 
%               ARP_i           {1}     [ms]        duration of I-cell ARP
%               SD_e            {0.1}   [ms]        duration of every E spike 
%               SD_i            {0.1}   [ms]        duration of every I spike 
%
%           all spiking model parameters can be set at the level of individual units. 
%               if that is desired, the size of the array should correspond to Ne/Ni. 
%               For instance, C_e can be either a scalar or an Ne x 1 vector. 
%
%           connectivity parameters - chemical synapses:
%
%               tau_rise_Se     {0.1};  [ms];       [ N Ne ];   for all E-synapses (ee, ie)
%               tau_fall_Se     {3.0};  [ms];       [ N Ne ];   for all E-synapses
%               Ese             {0};    [mV];       [ N Ne ];   for all E-synapses (ee, ie)
%               Gee             {0.12}; [mS/cm^2];  [ Ne Ne ];  for all E-to-E synapses
%               Gie             {0.12}; [mS/cm^2];  [ Ni Ne ];  for all E-to-I synapses
%               tau_rise_Si     {0.3};  [ms];       [ N Ni ];   for all I-synapses (ii, ei)
%               tau_fall_Si     {9.0};  [ms];       [ N Ni ];   for all I-synapses
%               Esi             {-80};  [mV];       [ N Ni ];   for all I-synapses (ii, ei)
%               Gii             {0.05}; [mS/cm^2];  [ Ni Ni ];  for all I-to-I synapses
%               Gei             {0.30}; [mS/cm^2];  [ Ne Ni ];  for all I-to-E synapses
%
%           connectivity parameters - plasticity of chemical synapses:
%
%               tau_rise_De     {0.1}   [ms]        [ N Ne ];   for all E-synapses (ee, ie); rise time for synaptic depression 
%               tau_fall_De     {100}   [ms]        [ N Ne ];   for all E-synapses (ee, ie); fall time for synaptic depression 
%               tau_rise_Fe     {0.2}   [ms]        [ N Ne ];   for all E-synapses (ee, ie); rise time for synaptic facilitation 
%               tau_fall_Fe     {300}   [ms]        [ N Ne ];   for all E-synapses (ee, ie); fall time for synaptic facilitation
%               tau_rise_Di     {0.1}   [ms]        [ N Ni ];   for all I-synapses (ii, ei); rise time for synaptic depression 
%               tau_fall_Di     {100}   [ms]        [ N Ni ];   for all I-synapses (ii, ei); fall time for synaptic depression 
%               tau_rise_Fi     {0.2}   [ms]        [ N Ni ];   for all I-synapses (ii, ei); rise time for synaptic facilitation 
%               tau_fall_Fi     {300}   [ms]        [ N Ni ];   for all I-synapses (ii, ei); fall time for synaptic facilitation
%
%           all connectivity parameters can be set at the level of
%               individual synapses. if that is desired, the size of the
%               array should correspond to Ne/Ni. For instance, Gie can be
%               either a scalar or an Ni x Ne matrix. 
%
% returns       st              spike trains: sparse array, [ nt x N ]
%               Ve              membrane potentials, [ nt x Ne ]
%               Vi              membrane potentials, [ nt x Ni ]
%               t               time vector [ms]
%               See             synaptic activation (0-1), [ Ne Ne ]
%               Sei             synaptic activation (0-1), [ Ne Ni ]
%               Sii             synaptic activation (0-1), [ Ni Ni ]
%               Sie             synaptic activation (0-1), [ Ni Ne ]
%
% notation      Gee             for a matrix 
%                                   Gee = [ G11 G12; G21 G22 ]
%                               G21 is the conductance from E-cell 1 to E-cell 2
%
% calls         ParseArgPairs
%               plot_raster, plotTraces
%
% see also      NR_syntransmit

%-----------------------------------------------------------------------
% Cellular model: LIF with noise and refractory period
%
% C dV/dt   = I_ext - gL*(V-E_L) - gS*S*(V-Es) + gN*eta
%
% if V > Vth, then V = Vreset (for SD samples, then ARP)
%
%-----------------------------------------------------------------------
% E-cell alternative: with I_Na,p, I_h, noise, and spiking mechanism (2D + reset)
%
% C dV/dt = I_ext - gL*(V-E_L) - gS*S*(V-Es) + gN*eta
%           - g_P * pinf( V ) * ( V - E_Na ) 
%           - g_h * r * ( V - E_h )
% dr/dt = ( rinf( V )- r ) / rtau
%
% if V > Vth, then V = Vreset; r = rreset (for spike duration, then ARP)
%
%-----------------------------------------------------------------------
% Synaptic model: exponential rise and fall
% Borgers et al., 2012, PLOS Comput. Biol.
% 
% dS/dt     = H(Vpre) * (1-S)/tau_rise_S - S/tau_fall_S
% H( V )    = (1+tanh(V - Vth)/4)/2
%
%-----------------------------------------------------------------------
% Synaptic plasticity model:
% Stark et al, 2022, PLOS Comput. Biol.
%
% dS/dt     = H(Vpre) * (1-S)/tau_rise_S - S/tau_fall_S
% dD/dt     = -H(Vpre) * D/tau_rise_D - (1-D)/tau_fall_D
% dF/dt     = H(Vpre) * (1-F)/tau_rise_F - F/tau_fall_F
% H(V)      = (1+tanh(V - Vth)/4)/2;
% 
% C dV/dt = I_ext - gL*(V-E_L) - gS*S*D*F*(V-Es) + eta
% 
% possibly different parameters (tau_rise_D/F, tau_fall_D/F) for different synapses
%
%-----------------------------------------------------------------------
% Network connectivity:
%
% E-E, E-I, I-I, and I-E
% autapses dis/allowed
% optional random pruning of E-E synapses
%-----------------------------------------------------------------------

% 07-sep-20 ES

% last update
% 30-jun-22

function [ st, Ve, Vi, t, See, Sei, Sii, Sie ] = NR_einet( N, Iext, varargin )

% [C] = muF/cm^2
% [I] = muA/cm^2
% [G] = mS/cm^2
% [t] = ms
% [V] = mV

%-------------------------------------------------------------------------
% constants
%-------------------------------------------------------------------------
% time steps required for stability
MAX_DT_LIF                      = 0.1;                                      % [ms]

% membrane potential parameters - pre synaptic cell
Vth_pre                         = 0;                                        % [mV]; sets inflection point of presynaptic sigmoid

% graphics
colors                          = [ 1 0 0; 0 0 0.7 ];                       % E-color, I-color

%-------------------------------------------------------------------------
% default parameters
%-------------------------------------------------------------------------

% membrane potential parameters - same for all units (E,I)
C_DFLT                          = 1.0;                                      % [muF/cm^2]; capacitance
Gl_LIF_DFLT                     = 0.5;                                      % [mS/cm^2]; leak conductance
El_LIF_DFLT                     = -60;                                      % [mV]; leak reversal potential

Vth_LIF_DFLT                    = -50;                                      % [mV]; spike threshold
Vreset_LIF_DFLT                 = -70;                                      % [mV]; post-spike reset
Vpeak_LIF_DFLT                  = 50;                                       % [mV]; spike peak (for plotting the spike manually)

% Inap (E-cell) parameters
Inap_DFLT                       = 0;                                        % no Na,p and h currents in E-cells
Gl_Inap_DFLT                    = 0.1;                                      % [mS/cm^2]; leak conductance
El_Inap_DFLT                    = -65;                                      % [mV]; leak reversal potential
Gp_DFLT                         = 0.1;                                      % [mS/cm^2]; persistent sodium conductance
Ena_DFLT                        = 55;                                       % [mV]; persistent sodium reversal potential 
Gh_DFLT                         = 1;                                        % [mS/cm^2]; h-current conductance
Eh_DFLT                         = -20;                                      % [mV]; h-current reversal potential 
rtaunap_DFLT                    = 100;                                      % [ms]; fixed (not voltage dependent)
rreset_DFLT                     = 0;                                        % h-current gating variable post-spike reset

% timing parameters
dt_DFLT                         = 0.1;                                      % [ms]
Tmax_DFLT                       = 2000;                                     % [ms]

% default noise level
Dne_DFLT                         = 2;                                       % [mV] 
Dni_DFLT                         = 2;                                       % [mV]

% default synaptic parameters
tau_rise_Se_DFLT                = 0.1;                                      % [ms]; AMPA
tau_fall_Se_DFLT                = 3.0;                                      % [ms]; AMPA
tau_rise_Si_DFLT                = 0.3;                                      % [ms]; GABA_A
tau_fall_Si_DFLT                = 9.0;                                      % [ms]; GABA_A

Ese_DFLT                        = 0.0;                                      % [mV]; AMPA
Esi_DFLT                        = -80;                                      % [mV]; GABA_A

% default chemical synaptic connectivity parameters
Gee_DFLT                        = 0.12;                                     % [mS/cm^2]; from E to E
Gie_DFLT                        = 0.12;                                     % [mS/cm^2]; from E to I
Gii_DFLT                        = 0.05;                                     % [mS/cm^2]; from I to I
Gei_DFLT                        = 0.30;                                     % [mS/cm^2]; from I to E

% default chemical synaptic plasticity parameters
tau_rise_De_DFLT                = 0.1;                                      % [ms]; depression rise time, on E
tau_fall_De_DFLT                = 100;                                      % [ms]; depression decay time, on E
tau_rise_Fe_DFLT                = 0.2;
tau_fall_Fe_DFLT                = 300;
tau_rise_Di_DFLT                = 0.1;                                      % [ms]; depression rise time, on I
tau_fall_Di_DFLT                = 100;                                      % [ms]; depression decay time, on 
tau_rise_Fi_DFLT                = 0.2;
tau_fall_Fi_DFLT                = 300;

% default external current
Ibias_e_DFLT                    = 0;                                        % [muA/cm^2]
Ibias_i_DFLT                    = 0;                                        % [muA/cm^2]

% ARP, spike duration
ARP_e_DFLT                      = 1;                                        % [ms]; E-cell absolute refractory period
ARP_i_DFLT                      = 1;                                        % [ms]; I-cell absolute refractory period
SD_e_DFLT                       = 0.1;                                      % [ms]; duration of each E spike
SD_i_DFLT                       = 0.1;                                      % [ms]; duration of each I spike

%-------------------------------------------------------------------------
% functions
%-------------------------------------------------------------------------
% presynaptic sigmoid
H                               = @(v) ( 1 + tanh( ( v - Vth_pre ) / 4 ) ) / 2;

% Inap 
pinfnap                         = @(v) 1.0 ./ ( 1 + exp( -( v + 38 ) / 6.5 ) );
rinfnap                         = @(v) 1.0 ./ ( 1 + exp( ( v + 79.2 ) / 9.78 ) ); 

%-------------------------------------------------------------------------
% arguments
%-------------------------------------------------------------------------
nargs                           = nargin;

if nargs < 1 || isempty( N )
    N                           = 10;
end
if length( N ) > 2 || ~isequal( N, round( N ) )
    error( 'N must be either one or two integers' )
end
if nargs < 2 || isempty( Iext )
    Iext                        = 0;
end

% other arguments
[ fE, fEprune, dt, synmodel ...
    , C_e, Gl_e, El_e, Vth_e, Vreset_e, Vpeak_e ...
    , C_i, Gl_i, El_i, Vth_i, Vreset_i, Vpeak_i ...
    , Inap, Gp_e, Ena_e, Gh_e, Eh_e, rtaunap_e, rreset_e ...
    , tau_rise_Se, tau_fall_Se, Ese, Gee, Gie ...
    , tau_rise_Si, tau_fall_Si, Esi, Gii, Gei ...
    , tau_rise_De, tau_fall_De, tau_rise_Fe, tau_fall_Fe ...
    , tau_rise_Di, tau_fall_Di, tau_rise_Fi, tau_fall_Fi ...
    , Dne, Dni ...
    , Tmax, Ibias_e, Ibias_i ... 
    , ARP_e, ARP_i, SD_e, SD_i ...
    , autapseE, autapseI ...
    , fast, thin, maxGB, verbose, graphics ]    = ParseArgPairs(...
    { 'fE', 'fEprune', 'dt', 'synmodel' ...
    , 'C_e', 'Gl_e', 'El_e', 'Vth_e', 'Vreset_e', 'Vpeak_e' ...
    , 'C_i', 'Gl_i', 'El_i', 'Vth_i', 'Vreset_i', 'Vpeak_i' ...
    , 'Inap', 'Gp_e', 'Ena_e', 'Gh_e', 'Eh_e', 'rtaunap_e', 'rreset_e' ...
    , 'tau_rise_Se', 'tau_fall_Se', 'Ese', 'Gee', 'Gie' ...
    , 'tau_rise_Si', 'tau_fall_Si', 'Esi', 'Gii', 'Gei' ...
    , 'tau_rise_De', 'tau_fall_De', 'tau_rise_Fe', 'tau_fall_Fe' ...
    , 'tau_rise_Di', 'tau_fall_Di', 'tau_rise_Fi', 'tau_fall_Fi' ...
    , 'Dne', 'Dni'...
    , 'Tmax', 'Ibias_e', 'Ibias_i' ... 
    , 'ARP_e', 'ARP_i', 'SD_e', 'SD_i' ...
    , 'autapseE', 'autapseI' ...
    , 'fast', 'thin', 'maxGB', 'verbose', 'graphics' }...
    , { 0.8, 0, dt_DFLT, 'null' ...
    , C_DFLT, Gl_LIF_DFLT, El_LIF_DFLT, Vth_LIF_DFLT, Vreset_LIF_DFLT, Vpeak_LIF_DFLT ...
    , C_DFLT, Gl_LIF_DFLT, El_LIF_DFLT, Vth_LIF_DFLT, Vreset_LIF_DFLT, Vpeak_LIF_DFLT ...
    , Inap_DFLT, Gp_DFLT, Ena_DFLT, Gh_DFLT, Eh_DFLT, rtaunap_DFLT, rreset_DFLT ...
    , tau_rise_Se_DFLT, tau_fall_Se_DFLT, Ese_DFLT, Gee_DFLT, Gie_DFLT ...
    , tau_rise_Si_DFLT, tau_fall_Si_DFLT, Esi_DFLT, Gii_DFLT, Gei_DFLT ...
    , tau_rise_De_DFLT, tau_fall_De_DFLT, tau_rise_Fe_DFLT, tau_fall_Fe_DFLT ...
    , tau_rise_Di_DFLT, tau_fall_Di_DFLT, tau_rise_Fi_DFLT, tau_fall_Fi_DFLT ...
    , Dne_DFLT, Dni_DFLT ...
    , Tmax_DFLT, Ibias_e_DFLT, Ibias_i_DFLT ... 
    , ARP_e_DFLT, ARP_i_DFLT, SD_e_DFLT, SD_i_DFLT ...
    , 0, 1 ...
    , 1, 1, 2, 0, 0 } ...
    , varargin{ : } );

if ~ismember( synmodel, { 'null', 'DF', 'D', 'F' } )
    error( 'unrecognized synaptic plasticity model' )
end
if numel( Inap ) ~= 1 || ~ismember( Inap, [ 0 1 ] )
    error( 'Select either Inap or not (LIF) for E-cells' )
end
    
% modify default Gl and El for E cells if Inap model
if Inap
    Gl_e_by_DFLT                = 1;
    for i                       = 1 : length( varargin )
        if isequal( varargin{ i }, 'Gl_e' )
            Gl_e_by_DFLT        = 0;
            break
        end
    end
    if Gl_e_by_DFLT
        Gl_e                    = Gl_Inap_DFLT;
    end
    
    El_e_by_DFLT                = 1;
    for i                       = 1 : length( varargin )
        if isequal( varargin{ i }, 'El_e' )
            El_e_by_DFLT        = 0;
            break
        end
    end
    if El_e_by_DFLT
        El_e                    = El_Inap_DFLT;
    end
end

if fEprune < 0 || fEprune > 1
    error( 'fEprune must be between 0 and 1' )
end

%-------------------------------------------------------------------------
% parameters
%-------------------------------------------------------------------------

% number of E,I units
if length( N ) == 2
    Ne                      	= N( 1 );
    Ni                          = N( 2 );
    N                           = Ne + Ni;
else
    Ne                          = round( N * fE );
    Ni                          = N - Ne;
end

% E-cell parameters - may be different for every unit
ue                              = ones( 1, Ne );
if numel( Dne ) == 1
    Dne                         = Dne * ue;
end
if numel( Ibias_e ) == 1
    Ibias_e                   	= Ibias_e * ue;
end
if numel( C_e ) == 1
    C_e                         = C_e * ue;
end
if numel( Gl_e ) == 1
    Gl_e                        = Gl_e * ue;
end
if numel( El_e ) == 1
    El_e                        = El_e * ue;
end
if numel( Vth_e ) == 1
    Vth_e                    	= Vth_e * ue;
end
if numel( Vreset_e ) == 1
    Vreset_e                    = Vreset_e * ue;
end
if numel( Vpeak_e ) == 1
    Vpeak_e                     = Vpeak_e * ue;
end
if numel( ARP_e ) == 1
    ARP_e                      	= ARP_e * ue;
end
if numel( SD_e ) == 1
    SD_e                     	= SD_e * ue;
end
if ~isequal( size( C_e ), [ 1 Ne ] )
    error( 'C_e: input size mismatch' )
end
if ~isequal( size( Gl_e ), [ 1 Ne ] )
    error( 'Gl_e: input size mismatch' )
end
if ~isequal( size( El_e ), [ 1 Ne ] )
    error( 'El_e: input size mismatch' )
end
if ~isequal( size( Vth_e ), [ 1 Ne ] )
    error( 'Vth_e: input size mismatch' )
end
if ~isequal( size( Vreset_e ), [ 1 Ne ] )
    error( 'Vreset_e: input size mismatch' )
end
if ~isequal( size( Vpeak_e ), [ 1 Ne ] )
    error( 'Vpeak_e: input size mismatch' )
end
if ~isequal( size( Vreset_e ), [ 1 Ne ] )
    error( 'Vreset_e: input size mismatch' )
end
if ~isequal( size( ARP_e ), [ 1 Ne ] )
    error( 'ARP_e: input size mismatch' )
end
if ~isequal( size( SD_e ), [ 1 Ne ] )
    error( 'SD_e: input size mismatch' )
end

if Inap
    if numel( Gp_e ) == 1
        Gp_e                 	= Gp_e * ue;
    end
    if numel( Ena_e ) == 1
        Ena_e                 	= Ena_e * ue;
    end
    if numel( Gh_e ) == 1
        Gh_e                  	= Gh_e * ue;
    end
    if numel( Eh_e ) == 1
        Eh_e                	= Eh_e * ue;
    end
    if numel( rtaunap_e ) == 1
        rtaunap_e               = rtaunap_e * ue;
    end
    if numel( rreset_e ) == 1
        rreset_e                = rreset_e * ue;
    end
    if ~isequal( size( Gp_e ), [ 1 Ne ] )
        error( 'Gp_e: input size mismatch' )
    end
    if ~isequal( size( Ena_e ), [ 1 Ne ] )
        error( 'Ena_e: input size mismatch' )
    end
    if ~isequal( size( Gh_e ), [ 1 Ne ] )
        error( 'Gh_e: input size mismatch' )
    end
    if ~isequal( size( Eh_e ), [ 1 Ne ] )
        error( 'Eh_e: input size mismatch' )
    end
    if ~isequal( size( rtaunap_e ), [ 1 Ne ] )
        error( 'rtaunap_e: input size mismatch' )
    end
    if ~isequal( size( rreset_e ), [ 1 Ne ] )
        error( 'rreset_e: input size mismatch' )
    end
end

% I-cell parameters - may be different for every unit
ui                              = ones( 1, Ni );
if numel( Dni ) == 1
    Dni                         = Dni * ui;
end
if numel( Ibias_i ) == 1
    Ibias_i                   	= Ibias_i * ui;
end
if numel( C_i ) == 1
    C_i                         = C_i * ui;
end
if numel( Gl_i ) == 1
    Gl_i                        = Gl_i * ui;
end
if numel( El_i ) == 1
    El_i                        = El_i * ui;
end
if numel( Vth_i ) == 1
    Vth_i                    	= Vth_i * ui;
end
if numel( Vreset_i ) == 1
    Vreset_i                    = Vreset_i * ui;
end
if numel( Vpeak_i ) == 1
    Vpeak_i                     = Vpeak_i * ui;
end
if numel( ARP_i ) == 1
    ARP_i                    	= ARP_i * ui;
end
if numel( SD_i ) == 1
    SD_i                     	= SD_i * ui;
end
if ~isequal( size( C_i ), [ 1 Ni ] )
    error( 'C_i: input size mismatch' )
end
if ~isequal( size( Gl_i ), [ 1 Ni ] )
    error( 'Gl_i: input size mismatch' )
end
if ~isequal( size( El_i ), [ 1 Ni ] )
    error( 'El_i: input size mismatch' )
end
if ~isequal( size( Vth_i ), [ 1 Ni ] )
    error( 'Vth_i: input size mismatch' )
end
if ~isequal( size( Vreset_i ), [ 1 Ni ] )
    error( 'Vreset_i: input size mismatch' )
end
if ~isequal( size( Vpeak_i ), [ 1 Ni ] )
    error( 'Vpeak_i: input size mismatch' )
end
if ~isequal( size( ARP_i ), [ 1 Ni ] )
    error( 'ARP_i: input size mismatch' )
end
if ~isequal( size( SD_i ), [ 1 Ni ] )
    error( 'SD_i: input size mismatch' )
end

% (chemical) synaptic parameters - may be different for every synapse
if numel( tau_rise_Se ) == 1
    tau_rise_Se                 = repmat( tau_rise_Se, [ N Ne ] );          % dim1: unit (Ne,Ni); dim2: presynaptic source(Ne)
end
if numel( tau_fall_Se ) == 1
    tau_fall_Se                 = repmat( tau_fall_Se, [ N Ne ] );
end
if numel( Ese ) == 1
    Ese                         = repmat( Ese, [ N Ne ] );                  % every row is the Es of E-synapses on one unit
end
if numel( tau_rise_Si ) == 1
    tau_rise_Si                 = repmat( tau_rise_Si, [ N Ni ] );          % dim1: unit (Ne,Ni); dim2: presynaptic source(Ni)
end
if numel( tau_fall_Si ) == 1
    tau_fall_Si                 = repmat( tau_fall_Si, [ N Ni ] );
end
if numel( Esi ) == 1
    Esi                         = repmat( Esi, [ N Ni ] );
end
if ~isequal( size( tau_rise_Se ), [ N Ne ] )
    error( 'tau_rise_Se: input size mismatch' )
end
if ~isequal( size( tau_fall_Se ), [ N Ne ] )
    error( 'tau_fall_Se: input size mismatch' )
end
if ~isequal( size( Ese ), [ N Ne ] )
    error( 'Ese: input size mismatch' )
end
if ~isequal( size( tau_rise_Si ), [ N Ni ] )
    error( 'tau_rise_Si: input size mismatch' )
end
if ~isequal( size( tau_fall_Si ), [ N Ni ] )
    error( 'tau_fall_Si: input size mismatch' )
end
if ~isequal( size( Esi ), [ N Ni ] )
    error( 'Esi: input size mismatch' )
end

% connectivity parameters - chemical synapses - may be different for every synapse
if numel( Gee ) == 1
    Gee                         = repmat( Gee, [ Ne Ne ] );
end
if numel( Gie ) == 1
    Gie                         = repmat( Gie, [ Ni Ne ] );
end
if numel( Gii ) == 1
    Gii                         = repmat( Gii, [ Ni Ni ] );
end
if numel( Gei ) == 1
    Gei                         = repmat( Gei, [ Ne Ni ] );
end
if ~isequal( size( Gee ), [ Ne Ne ] )
    error( 'Gee: input size mismatch' )
end
if ~isequal( size( Gie ), [ Ni Ne ] )
    error( 'Gie: input size mismatch' )
end
if ~isequal( size( Gii ), [ Ni Ni ] )
    error( 'Gii: input size mismatch' )
end
if ~isequal( size( Gei ), [ Ne Ni ] )
    error( 'Gei: input size mismatch' )
end
if any( Gee( : ) < 0 )
    error( 'Gee: all values must be positive' )
end
if any( Gie( : ) < 0 )
    error( 'Gie: all values must be positive' )
end
if any( Gii( : ) < 0 )
    error( 'Gii: all values must be positive' )
end
if any( Gei( : ) < 0 )
    error( 'Gei: all values must be positive' )
end
doSee                           = 0;
doSei                           = 0;
doSii                           = 0;
doSie                           = 0;
if any( [ Gee( : ); Gei( : ); Gii( : ); Gie( : ) ] )
    if any( Gee( : ) )
        doSee                   = 1;
    end
    if any( Gei( : ) )
        doSei                   = 1;
    end
    if any( Gii( : ) )
        doSii                   = 1;
    end
    if any( Gie( : ) )
        doSie                   = 1;
    end
end

% connectivity parameters - chemical synaptic plasticity parameters
switch synmodel
    case 'DF'
        synp                    = 1;
        F0                      = 0;
    case 'D'
        synp                    = 1;
        F0                      = 1;
        tau_rise_Fe             = inf;
        tau_fall_Fe             = inf;
        tau_rise_Fi             = inf;
        tau_fall_Fi             = inf;
    case 'F'
        synp                    = 1;
        F0                      = 0;
        tau_rise_De             = inf;
        tau_fall_De             = inf;
        tau_rise_Di             = inf;
        tau_fall_Di             = inf;
    otherwise
        synp                    = 0;
        tau_rise_Fe             = inf;
        tau_fall_Fe             = inf;
        tau_rise_De             = inf;
        tau_fall_De             = inf;
        tau_rise_Fi             = inf;
        tau_fall_Fi             = inf;
        tau_rise_Di             = inf;
        tau_fall_Di             = inf;
end
if numel( tau_rise_De ) == 1
    tau_rise_De                 = repmat( tau_rise_De, [ N Ne ] );          % dim1: unit (Ne,Ni); dim2: presynaptic source(Ne)
end
if numel( tau_fall_De ) == 1
    tau_fall_De                 = repmat( tau_fall_De, [ N Ne ] );
end
if numel( tau_rise_Fe ) == 1
    tau_rise_Fe                 = repmat( tau_rise_Fe, [ N Ne ] );
end
if numel( tau_fall_Fe ) == 1
    tau_fall_Fe                 = repmat( tau_fall_Fe, [ N Ne ] );
end
if numel( tau_rise_Di ) == 1
    tau_rise_Di                 = repmat( tau_rise_Di, [ N Ni ] );          % dim1: unit (Ne,Ni); dim2: presynaptic source(Ni)
end
if numel( tau_fall_Di ) == 1
    tau_fall_Di                 = repmat( tau_fall_Di, [ N Ni ] );
end
if numel( tau_rise_Fi ) == 1
    tau_rise_Fi                 = repmat( tau_rise_Fi, [ N Ni ] );
end
if numel( tau_fall_Fi ) == 1
    tau_fall_Fi                 = repmat( tau_fall_Fi, [ N Ni ] );
end
if ~isequal( size( tau_rise_De ), [ N Ne ] )
    error( 'tau_rise_De: input size mismatch' )
end
if ~isequal( size( tau_fall_De ), [ N Ne ] )
    error( 'tau_fall_De: input size mismatch' )
end
if ~isequal( size( tau_rise_Fe ), [ N Ne ] )
    error( 'tau_rise_Fe: input size mismatch' )
end
if ~isequal( size( tau_fall_Fe ), [ N Ne ] )
    error( 'tau_fall_Fe: input size mismatch' )
end
if ~isequal( size( tau_rise_Di ), [ N Ni ] )
    error( 'tau_rise_Di: input size mismatch' )
end
if ~isequal( size( tau_fall_Di ), [ N Ni ] )
    error( 'tau_fall_Di: input size mismatch' )
end
if ~isequal( size( tau_rise_Fi ), [ N Ni ] )
    error( 'tau_rise_Fi: input size mismatch' )
end
if ~isequal( size( tau_fall_Fi ), [ N Ni ] )
    error( 'tau_fall_Fi: input size mismatch' )
end

% handle autapses/auto-coupling
if ~autapseE
    didx                        = 1 : ( Ne + 1 ) : Ne^2;
    Gee( didx )                 = 0;
end
if ~autapseI
    didx                        = 1 : ( Ni + 1 ) : Ni^2;
    Gii( didx )                 = 0;
end

% prune EE synapses randomly
if fEprune > 0
    nee                         = numel( Gee );
    nprune                      = round( nee * fEprune );
    [ ~, pidx ]                 = sort( rand( nee, 1 ) );                   % randomly choose nprune connections
    Gee( pidx( 1 : nprune ) )   = 0;                                        % actually prune
end

%-------------------------------------------------------------------------
% input handling
%-------------------------------------------------------------------------

% time and external input vectors
if size( Iext, 2 ) == 3 && size( Iext, 1 ) > 1
    nt                          = size( Iext, 1 );
    Ie                          = Iext( :, 1 ) * ones( 1, Ne ) + ones( nt, 1 ) * Ibias_e;
    Ii                          = Iext( :, 2 ) * ones( 1, Ni ) + ones( nt, 1 ) * Ibias_i;
    t                           = Iext( :, 3 );
    dt                          = t( 2 ) - t( 1 );
elseif size( Iext, 2 ) == ( N + 1 ) && size( Iext, 1 ) > 1
    nt                          = size( Iext, 1 );
    Ie                          = Iext( :, 1 : Ne ) + ones( nt, 1 ) * Ibias_e;
    Ii                          = Iext( :, ( Ne + 1 ) : N ) + ones( nt, 1 ) * Ibias_i;
    t                           = Iext( :, N + 1 );
    dt                          = t( 2 ) - t( 1 );
else
    nt                          = Tmax / dt;
    if thin
        nGB                 	= N * nt * 4 / 2^30;
    else
        
        nGB                 	= ( N^2+N ) * nt * 4 / 2^30;                % single (4 bytes/element)
    end
    if nGB > maxGB
        fprintf( 1, 'RAM required: %0.2g GB; this is more than the maximal allowed (%0.2g GB); exiting\n', nGB, maxGB )
        return
    end
    t                           = ( dt : dt : Tmax )';                      % time vector, [ms]
    nt                          = length( t );
    Ie                          = ones( nt, 1 ) * Ibias_e;                  % each column - input to one E-unit
    Ii                          = ones( nt, 1 ) * Ibias_i;                  % each column - input to one I-unit
    if verbose
        fprintf( 1, 'Iext ignored; using only bias currents (Ibias_e: %s; Ibias_i: %s)\n', num2str( Ibias_e ), num2str( Ibias_i ) )
    end
end

if dt > MAX_DT_LIF
    fprintf( 1, 'LIF will be unstable with dt=%0.3g; exiting\n', dt, MAX_DT_LIF )
    return
end

% dt-dependent parameters
dtARPE                          = ceil( ARP_e / dt );
dtARPI                          = ceil( ARP_i / dt );
dtSDE                           = ceil( SD_e / dt );
dtSDI                           = ceil( SD_i / dt );

%-------------------------------------------------------------------------
% checks:
%-------------------------------------------------------------------------

if dt > dt_DFLT
    fprintf( 1, 'Specified step size: %0.2g ms; this is larger than the minimum required to guarantee stability (%0.2g ms); exiting\n', dt, dt_DFLT )
    return
end

if thin
    nGB                         = N * nt * 4 / 2^30;
else
    nGB                         = ( N^2+N ) * nt * 4 / 2^30;                % single (4 bytes/element)
end

if nGB > maxGB
    fprintf( 1, 'RAM required: %0.2g GB; this is more than the maximal allowed (%0.2g GB); exiting\n', nGB, maxGB )
    return
end

%-------------------------------------------------------------------------
% initialization:
%-------------------------------------------------------------------------
    
% variables:
Ve                              = zeros( nt, Ne, 'single' );                % dim1: time; dim2: unit
Vi                              = zeros( nt, Ni, 'single' );
if Inap
    re                          = zeros( nt, Ne, 'single' );
end

if thin
    if Inap
    end
    See                         = zeros( Ne, Ne, 1, 'single' );
    Sei                         = zeros( Ne, Ni, 1, 'single' );
    Sii                         = zeros( Ni, Ni, 1, 'single' );
    Sie                         = zeros( Ni, Ne, 1, 'single' );
    if synp
        Dee                     = zeros( Ne, Ne, 1, 'single' );
        Dei                     = zeros( Ne, Ni, 1, 'single' );
        Dii                     = zeros( Ni, Ni, 1, 'single' );
        Die                     = zeros( Ni, Ne, 1, 'single' );
        Fee                     = zeros( Ne, Ne, 1, 'single' );
        Fei                     = zeros( Ne, Ni, 1, 'single' );
        Fii                     = zeros( Ni, Ni, 1, 'single' );
        Fie                     = zeros( Ni, Ne, 1, 'single' );
    end
else
    See                         = zeros( Ne, Ne, nt, 'single' );            % dim1: unit; dim2: presynaptic source; dim3: time
    Sei                         = zeros( Ne, Ni, nt, 'single' );
    Sii                         = zeros( Ni, Ni, nt, 'single' );
    Sie                         = zeros( Ni, Ne, nt, 'single' );
    if synp
        Dee                     = zeros( Ne, Ne, nt, 'single' );
        Dei                     = zeros( Ne, Ni, nt, 'single' );
        Dii                     = zeros( Ni, Ni, nt, 'single' );
        Die                     = zeros( Ni, Ne, nt, 'single' );
        Fee                     = zeros( Ne, Ne, nt, 'single' );
        Fei                     = zeros( Ne, Ni, nt, 'single' );
        Fii                     = zeros( Ni, Ni, nt, 'single' );
        Fie                     = zeros( Ni, Ne, nt, 'single' );
    end
end

% initial conditions:
Ve( 1, : )                      = El_e;                                   	% steady state
Vi( 1, : )                      = El_i;                                     % steady state
if Inap
    re( 1, : )                  = 0;
end
See( :, :, 1 )                  = 0;                                        % no synaptic conductance
Sei( :, :, 1 )                  = 0;                                        % no synaptic conductance
Sii( :, :, 1 )                  = 0;                                        % no synaptic conductance
Sie( :, :, 1 )                  = 0;                                        % no synaptic conductance
if synp
    Dee( :, :, 1 )              = 1;                                        % no depression
    Dei( :, :, 1 )            	= 1;                                        % no depression
    Dii( :, :, 1 )              = 1;                                        % no depression
    Die( :, :, 1 )              = 1;                                        % no depression
    Fee( :, :, 1 )              = F0;                                       %
    Fei( :, :, 1 )              = F0;                                       %
    Fii( :, :, 1 )              = F0;                                       %
    Fie( :, :, 1 )              = F0;                                       %
end

% additive noise:
wne                             = randn( nt, Ne );
wni                             = randn( nt, Ni );

% spike counter
NSe                             = zeros( 1, Ne );                           % counter: new spike in E-cell
NSi                             = zeros( 1, Ni );                           % counter: new spike in I-cell

varnames                        = whos; 
nGB_used                        = sum( [ varnames.bytes ] ) / 2^30;
if verbose
    fprintf( 1, '\n%s: RAM allocated: %0.2g GB; this is less than the maximal RAM allowed (%0.2g GB); continuing\n'...
        , upper( mfilename  ), nGB_used, maxGB )
end

%-------------------------------------------------------------------------
% simulation:
%-------------------------------------------------------------------------
for i                           = 1 : ( nt - 1 )
    
    if thin
        k                       = 1;
        kp                      = 1;
    else
        k                       = i;
        kp                      = i + 1;
    end
    
    % slope at i
    if Ne > 0
        k1ve                    = Ie( i, : );                               % row vector, Ne elements
        k1ve                    = k1ve - Gl_e .* ( Ve( i, : ) - El_e );     % row vector, Ne elements
        if Inap
            k1ve                = k1ve - Gp_e .* pinfnap( Ve( i, : ) ) .* ( Ve( i, : ) - Ena_e );
            k1ve                = k1ve - Gh_e .* re( i, : ) .* ( Ve( i, : ) - Eh_e );
            k1re                = ( rinfnap( Ve( i, : ) ) - re( i, : ) ) ./ rtaunap_e;
        end
        if doSee
            if fast
                if synp
                    SDF         = ones( Ne, 1 ) * ( See( 1, :, k ) .* Dee( 1, :, k ) .* Fee( 1, :, k ) );
                else
                    SDF       	= ones( Ne, 1 ) * See( 1, :, k );
                end
            else
                if synp
                    SDF         = See( :, :, k ) .* Dee( :, :, k ) .* Fee( :, :, k );
                else
                    SDF       	= See( :, :, k );
                end
            end
            i_syn               = ( Gee .* SDF ) .* ( Ve( i, : )' * ones( 1, Ne ) - Ese( 1 : Ne, : ) );
            k1ve                = k1ve - sum( i_syn, 2 )';
        end
        if Ni > 0
            if doSei
                if fast
                    if synp
                        SDF  	= ones( Ne, 1 ) * ( Sei( 1, :, k ) .* Dei( 1, :, k ) .* Fei( 1, :, k ) );
                    else
                        SDF    	= ones( Ne, 1 ) * Sei( 1, :, k );
                    end
                else
                    if synp
                        SDF  	= Sei( :, :, k ) .* Dei( :, :, k ) .* Fei( :, :, k );
                    else
                        SDF    	= Sei( :, :, k );
                    end
                end
                i_syn           = ( Gei .* SDF ) .* ( Ve( i, : )' * ones( 1, Ni ) - Esi( 1 : Ne, : ) );
                k1ve            = k1ve - sum( i_syn, 2 )';
            end
        end
        k1ve                	= k1ve ./ C_e;
        if doSee
            if fast
                k1see         	= H( Ve( i, : ) ) .* ( 1 - See( 1, :, k ) ) ./ tau_rise_Se( 1 ) - See( 1, :, k ) ./ tau_fall_Se( 1 );                               % 1 by Ne vector
                if synp
                    k1dee    	= -H( Ve( i, : ) ) .* Dee( 1, :, k ) ./ tau_rise_De( 1 ) + ( 1 - Dee( 1, :, k ) ) ./ tau_fall_De( 1 );
                    k1fee     	= H( Ve( i, : ) ) .* ( 1 - Fee( 1, :, k ) ) ./ tau_rise_Fe( 1 ) - Fee( 1, :, k ) ./ tau_fall_Fe( 1 );
                end
            else
                k1see         	= ( ones( Ne, 1 ) * H( Ve( i, : ) ) ) .* ( 1 - See( :, :, k ) ) ./ tau_rise_Se( 1 : Ne, : ) - See( :, :, k ) ./ tau_fall_Se( 1 : Ne, : ); % Ne by Ne matrix
                if synp
                    k1dee    	= ( ones( Ne, 1 ) * -H( Ve( i, : ) ) ) .* Dee( :, :, k ) ./ tau_rise_De( 1 : Ne, : ) + ( 1 - Dee( :, :, k ) ) ./ tau_fall_De( 1 : Ne, : );
                    k1fee     	= ( ones( Ne, 1 ) * H( Ve( i, : ) ) ) .* ( 1 - Fee( :, :, k ) ) ./ tau_rise_Fe( 1 : Ne, : ) - Fee( :, :, k ) ./ tau_fall_Fe( 1 : Ne, : );
                end
            end
        end
        if doSei
            if Ni > 0
                if fast
                    k1sei       = H( Vi( i, : ) ) .* ( 1 - Sei( 1, :, k ) ) ./ tau_rise_Si( 1 ) - Sei( 1, :, k ) ./ tau_fall_Si( 1 );                               % 1 by Ni vector
                    if synp
                        k1dei	= -H( Vi( i, : ) ) .* Dei( 1, :, k ) ./ tau_rise_Di( 1 ) + ( 1 - Dei( 1, :, k ) ) ./ tau_fall_Di( 1 );
                        k1fei   = H( Vi( i, : ) ) .* ( 1 - Fei( 1, :, k ) ) ./ tau_rise_Fi( 1 ) - Fei( 1, :, k ) ./ tau_fall_Fi( 1 );
                    end
                else
                    k1sei       = ( ones( Ne, 1 ) * H( Vi( i, : ) ) ) .* ( 1 - Sei( :, :, k ) ) ./ tau_rise_Si( 1 : Ne, : ) - Sei( :, :, k ) ./ tau_fall_Si( 1 : Ne, : ); % Ne by Ni matrix
                    if synp
                        k1dei	= ( ones( Ne, 1 ) * -H( Vi( i, : ) ) ) .* Dei( :, :, k ) ./ tau_rise_Di( 1 : Ne, : ) + ( 1 - Dei( :, :, k ) ) ./ tau_fall_Di( 1 : Ne, : );
                        k1fei   = ( ones( Ne, 1 ) * H( Vi( i, : ) ) ) .* ( 1 - Fei( :, :, k ) ) ./ tau_rise_Fi( 1 : Ne, : ) - Fei( :, :, k ) ./ tau_fall_Fi( 1 : Ne, : );
                    end
                end
            end
        end
        
    end
    
    if Ni > 0
        k1vi                    = Ii( i, : ) - Gl_i .* ( Vi( i, : ) - El_i );    % row vector, Ni elements
        if Ne > 0
            if doSie
                if fast
                    if synp
                        SDF   	= ones( Ni, 1 ) * ( Sie( 1, :, k ) .* Die( 1, :, k ) .* Fie( 1, :, k ) );
                    else
                        SDF   	= ones( Ni, 1 ) * Sie( 1, :, k );
                    end
                else
                    if synp
                        SDF   	= Sie( :, :, k ) .* Die( :, :, k ) .* Fie( :, :, k );
                    else
                        SDF   	= Sie( :, :, k );
                    end
                end
                i_syn           = ( Gie .* SDF ) .* ( Vi( i, : )' * ones( 1, Ne ) - Ese( ( Ne + 1 ) : N, : ) );
                k1vi            = k1vi - sum( i_syn, 2 )';
            end
        end
        if doSii
            if fast
                if synp
                    SDF     	= ones( Ni, 1 ) * ( Sii( 1, :, k ) .* Dii( 1, :, k ) .* Fii( 1, :, k ) );
                else
                    SDF        	= ones( Ni, 1 ) * Sii( 1, :, k );
                end
            else
                if synp
                    SDF     	= Sii( :, :, k ) .* Dii( :, :, k ) .* Fii( :, :, k );
                else
                    SDF        	= Sii( :, :, k );
                end
            end
            i_syn               = ( Gii .* SDF ) .* ( Vi( i, : )' * ones( 1, Ni ) - Esi( ( Ne + 1 ) : N, : ) );
            k1vi                = k1vi - sum( i_syn, 2 )';
        end
        k1vi                   	= k1vi ./ C_i;
        if doSie
            if Ne > 0
                if fast
                    k1sie       = H( Ve( i, : ) ) .* ( 1 - Sie( 1, :, k ) ) ./ tau_rise_Se( 1 ) - Sie( 1, :, k ) ./ tau_fall_Se( 1 );                                                      % 1 by Ne vector
                    if synp
                        k1die   =  -H( Ve( i, : ) ) .* Die( 1, :, k ) ./ tau_rise_De( 1 ) + ( 1 - Die( 1, :, k ) ) ./ tau_fall_De( 1 );
                        k1fie   = H( Ve( i, : ) ) .* ( 1 - Fie( 1, :, k ) ) ./ tau_rise_Fe( 1 ) - Fie( 1, :, k ) ./ tau_fall_Fe( 1 );
                    end
                else
                    k1sie       = ( ones( Ni, 1 ) * H( Ve( i, : ) ) ) .* ( 1 - Sie( :, :, k ) ) ./ tau_rise_Se( ( Ne + 1 ) : N, : ) - Sie( :, :, k ) ./ tau_fall_Se( ( Ne + 1 ) : N, : ); % Ni by Ne matrix
                    if synp
                        k1die   = ( ones( Ni, 1 ) * -H( Ve( i, : ) ) ) .* Die( :, :, k ) ./ tau_rise_De( ( Ne + 1 ) : N, : ) + ( 1 - Die( :, :, k ) ) ./ tau_fall_De( ( Ne + 1 ) : N, : );
                        k1fie   = ( ones( Ni, 1 ) * H( Ve( i, : ) ) ) .* ( 1 - Fie( :, :, k ) ) ./ tau_rise_Fe( ( Ne + 1 ) : N, : ) - Fie( :, :, k ) ./ tau_fall_Fe( ( Ne + 1 ) : N, : );
                    end
                end
            end
        end
        if doSii
            if fast
                k1sii           = H( Vi( i, : ) ) .* ( 1 - Sii( 1, :, k ) ) ./ tau_rise_Si( 1 ) - Sii( 1, :, k ) ./ tau_fall_Si( 1 );                                                       % 1 by Ni vector
                if synp
                    k1dii       = -H( Vi( i, : ) ) .* Dii( 1, :, k ) ./ tau_rise_Di( 1 ) + ( 1 - Dii( 1, :, k ) ) ./ tau_fall_Di( 1 );
                    k1fii       = H( Vi( i, : ) ) .* ( 1 - Fii( 1, :, k ) ) ./ tau_rise_Fi( 1 ) - Fii( 1, :, k ) ./ tau_fall_Fi( 1 );
                end
            else
                k1sii           = ( ones( Ni, 1 ) * H( Vi( i, : ) ) ) .* ( 1 - Sii( :, :, k ) ) ./ tau_rise_Si( ( Ne + 1 ) : N, : ) - Sii( :, :, k ) ./ tau_fall_Si( ( Ne + 1 ) : N, : ); % Ni by Ni matrix
                if synp
                    k1dii       = ( ones( Ni, 1 ) * -H( Vi( i, : ) ) ) .* Dii( :, :, k ) ./ tau_rise_Di( ( Ne + 1 ) : N, : ) + ( 1 - Dii( :, :, k ) ) ./ tau_fall_Di( ( Ne + 1 ) : N, : );
                    k1fii       = ( ones( Ni, 1 ) * H( Vi( i, : ) ) ) .* ( 1 - Fii( :, :, k ) ) ./ tau_rise_Fi( ( Ne + 1 ) : N, : ) - Fii( :, :, k ) ./ tau_fall_Fi( ( Ne + 1 ) : N, : );
                end
            end
        end
    end
    
    % value at endpoint ( i + 1 ), based on slope at i
    if Ne > 0
        a1ve                    = Ve( i, : ) + k1ve * dt;
        a1ve                    = a1ve + dt ^ 0.5 * Dne .* wne( i, : );
        if doSee
            if fast
                a1see         	= See( 1, :, k ) + k1see * dt;
                if synp
                    ad1ee      	= Dee( 1, :, k ) + k1dee * dt;
                    af1ee     	= Fee( 1, :, k ) + k1fee * dt;
                end
            else
                a1see         	= See( :, :, k ) + k1see * dt;
                if synp
                    ad1ee      	= Dee( :, :, k ) + k1dee * dt;
                    af1ee     	= Fee( :, :, k ) + k1fee * dt;
                end
            end
        end
        if doSei
            if Ni > 0
                if fast
                    a1sei    	= Sei( 1, :, k ) + k1sei * dt;
                    if synp
                        ad1ei  	= Dei( 1, :, k ) + k1dei * dt;
                        af1ei 	= Fei( 1, :, k ) + k1fei * dt;
                    end
                else
                    a1sei    	= Sei( :, :, k ) + k1sei * dt;
                    if synp
                        ad1ei  	= Dei( :, :, k ) + k1dei * dt;
                        af1ei 	= Fei( :, :, k ) + k1fei * dt;
                    end
                end
            end
        end
        if Inap
            ar1e                = re( i, : ) + k1re * dt;
        end
    end
    
    if Ni > 0
        a1vi                 	= Vi( i, : ) + k1vi * dt;
        a1vi                 	= a1vi + dt ^ 0.5 * Dni .* wni( i, : );
        if doSie
            if Ne > 0
                if fast
                    a1sie       = Sie( 1, :, k ) + k1sie * dt;
                    if synp
                        ad1ie 	= Die( 1, :, k ) + k1die * dt;
                        af1ie 	= Fie( 1, :, k ) + k1fie * dt;
                    end
                else
                    a1sie       = Sie( :, :, k ) + k1sie * dt;
                    if synp
                        ad1ie 	= Die( :, :, k ) + k1die * dt;
                        af1ie 	= Fie( :, :, k ) + k1fie * dt;
                    end
                end
            end
        end
        if doSii
            if fast
                a1sii         	= Sii( 1, :, k ) + k1sii * dt;
                if synp
                    ad1ii      	= Dii( 1, :, k ) + k1dii * dt;
                    af1ii      	= Fii( 1, :, k ) + k1fii * dt;
                end
            else
                a1sii         	= Sii( :, :, k ) + k1sii * dt;
                if synp
                    ad1ii      	= Dii( :, :, k ) + k1dii * dt;
                    af1ii      	= Fii( :, :, k ) + k1fii * dt;
                end
            end
        end
    end
   
    % slope at ( i + 1 )
    if Ne > 0
        k2ve                    = Ie( i + 1, : );
        k2ve                    = k2ve - Gl_e .* ( a1ve - El_e );
        if Inap
            k2ve                = k2ve - Gp_e .* pinfnap( a1ve ) .* ( a1ve - Ena_e );
            k2ve                = k2ve - Gh_e .* ar1e .* ( a1ve - Eh_e );
            k2re                = ( rinfnap( a1ve ) - ar1e ) ./ rtaunap_e;
        end
        if doSee
            if synp
                SDF             = a1see .* ad1ee .* af1ee;
            else
                SDF             = a1see;
            end
            i_syn               = ( Gee .* SDF ) .* ( a1ve' * ones( 1, Ne ) - Ese( 1 : Ne, : ) );
            k2ve                = k2ve - sum( i_syn, 2 )';
        end
        if Ni > 0
            if doSei
                if synp
                    SDF         = a1sei .* ad1ei .* af1ei;
                else
                    SDF        	= a1sei;
                end
                i_syn           = ( Gei .* SDF ) .* ( a1ve' * ones( 1, Ni ) - Esi( 1 : Ne, : ) );
                k2ve            = k2ve - sum( i_syn, 2 )';
            end
        end
        k2ve                   	= k2ve ./ C_e;
        if doSee
            if fast
                k2see           = H( a1ve ) .* ( 1 - a1see ) ./ tau_rise_Se( 1 ) - a1see ./ tau_fall_Se( 1 );
                if synp
                    k2dee       = -H( a1ve ) .* ad1ee ./ tau_rise_De( 1 ) + ( 1 - ad1ee ) ./ tau_fall_De( 1 );
                    k2fee       = H( a1ve ) .* ( 1 - af1ee ) ./ tau_rise_Fe( 1 ) - af1ee ./ tau_fall_Fe( 1 );
                end
            else
                k2see           = ( ones( Ne, 1 ) * H( a1ve ) ) .* ( 1 - a1see ) ./ tau_rise_Se( 1 : Ne, : ) - a1see ./ tau_fall_Se( 1 : Ne, : );
                if synp
                    k2dee       = ( ones( Ne, 1 ) * -H( a1ve ) ) .* ad1ee ./ tau_rise_De( 1 : Ne, : ) + ( 1 - ad1ee ) ./ tau_fall_De( 1 : Ne, : );
                    k2fee       = ( ones( Ne, 1 ) * H( a1ve ) ) .* ( 1 - af1ee ) ./ tau_rise_Fe( 1 : Ne, : ) - af1ee ./ tau_fall_Fe( 1 : Ne, : );
                end
            end
        end
        if doSei
            if Ni > 0
                if fast
                    k2sei       = H( a1vi ) .* ( 1 - a1sei ) ./ tau_rise_Si( 1 ) - a1sei ./ tau_fall_Si( 1 );
                    if synp
                        k2dei   = -H( a1vi ) .* ad1ei ./ tau_rise_Di( 1 ) + ( 1 - ad1ei ) ./ tau_fall_Di( 1 );
                        k2fei   = H( a1vi ) .* ( 1 - af1ei ) ./ tau_rise_Fi( 1 ) - af1ei ./ tau_fall_Fi( 1 );
                    end
                else
                    k2sei       = ( ones( Ne, 1 ) * H( a1vi ) ) .* ( 1 - a1sei ) ./ tau_rise_Si( 1 : Ne, : ) - a1sei ./ tau_fall_Si( 1 : Ne, : );
                    if synp
                        k2dei   = ( ones( Ne, 1 ) * -H( a1vi ) ) .* ad1ei ./ tau_rise_Di( 1 : Ne, : ) + ( 1 - ad1ei ) ./ tau_fall_Di( 1 : Ne, : );
                        k2fei   = ( ones( Ne, 1 ) * H( a1vi ) ) .* ( 1 - af1ei ) ./ tau_rise_Fi( 1 : Ne, : ) - af1ei ./ tau_fall_Fi( 1 : Ne, : );
                    end
                end
            end
        end
    end
    
    if Ni > 0
        k2vi                    = Ii( i + 1, : ) - Gl_i .* ( a1vi - El_i );
        if Ne > 0
            if doSie
                if synp
                    SDF     	= a1sie .* ad1ie .* af1ie;
                else
                    SDF        	= a1sie;
                end
                i_syn           = ( Gie .* SDF ) .* ( a1vi' * ones( 1, Ne ) - Ese( ( Ne + 1 ) : N, : ) );
                k2vi            = k2vi - sum( i_syn, 2 )';
            end
        end
        if doSii
            if synp
                SDF             = a1sii .* ad1ii .* af1ii;
            else
                SDF             = a1sii;
            end
            i_syn               = ( Gii .* SDF ) .* ( a1vi' * ones( 1, Ni ) - Esi( ( Ne + 1 ) : N, : ) );
            k2vi                = k2vi - sum( i_syn, 2 )';
        end
        k2vi                  	= k2vi ./ C_i;
        if doSie
            if Ne > 0
                if fast
                    k2sie     	= H( a1ve ) .* ( 1 - a1sie ) ./ tau_rise_Se( 1 ) - a1sie ./ tau_fall_Se( 1 );
                    if synp
                        k2die 	= -H( a1ve ) .* ad1ie ./ tau_rise_De( 1 ) + ( 1 - ad1ie ) ./ tau_fall_De( 1 );
                        k2fie  	= H( a1ve ) .* ( 1 - af1ie ) ./ tau_rise_Fe( 1 ) - af1ie ./ tau_fall_Fe( 1 );
                    end
                else
                    k2sie     	= ( ones( Ni, 1 ) * H( a1ve ) ) .* ( 1 - a1sie ) ./ tau_rise_Se( ( Ne + 1 ) : N, : ) - a1sie ./ tau_fall_Se( ( Ne + 1 ) : N, : );
                    if synp
                        k2die 	= ( ones( Ni, 1 ) * -H( a1ve ) ) .* ad1ie ./ tau_rise_De( ( Ne + 1 ) : N, : ) + ( 1 - ad1ie ) ./ tau_fall_De( ( Ne + 1 ) : N, : );
                        k2fie  	= ( ones( Ni, 1 ) * H( a1ve ) ) .* ( 1 - af1ie ) ./ tau_rise_Fe( ( Ne + 1 ) : N, : ) - af1ie ./ tau_fall_Fe( ( Ne + 1 ) : N, : );
                    end
                end
            end
        end
        if doSii
            if fast
                k2sii        	= H( a1vi ) .* ( 1 - a1sii ) ./ tau_rise_Si( 1 ) - a1sii ./ tau_fall_Si( 1 );
                if synp
                    k2dii      	= -H( a1vi ) .* ad1ii ./ tau_rise_Di( 1 ) + ( 1 - ad1ii ) ./ tau_fall_Di( 1 );
                    k2fii     	= H( a1vi ) .* ( 1 - af1ii ) ./ tau_rise_Fi( 1 ) - af1ii ./ tau_fall_Fi( 1 );
                end
            else
                k2sii        	= ( ones( Ni, 1 ) * H( a1vi ) ) .* ( 1 - a1sii ) ./ tau_rise_Si( ( Ne + 1 ) : N, : ) - a1sii ./ tau_fall_Si( ( Ne + 1 ) : N, : );
                if synp
                    k2dii      	= ( ones( Ni, 1 ) * -H( a1vi ) ) .* ad1ii ./ tau_rise_Di( ( Ne + 1 ) : N, : ) + ( 1 - ad1ii ) ./ tau_fall_Di( ( Ne + 1 ) : N, : );
                    k2fii     	= ( ones( Ni, 1 ) * H( a1vi ) ) .* ( 1 - af1ii ) ./ tau_rise_Fi( ( Ne + 1 ) : N, : ) - af1ii ./ tau_fall_Fi( ( Ne + 1 ) : N, : );
                end
            end
        end
    end
    
    % value at endpoint ( i + 1 ), based on average slope
    if Ne > 0
        Ve( i + 1, : )      	= Ve( i, : ) + ( k1ve + k2ve ) * dt / 2;
        Ve( i + 1, : )       	= Ve( i + 1, : ) + dt ^ 0.5 * Dne .* wne( i, : );      
        if Inap
            re( i + 1, : )      = re( i, : ) + ( k1re + k2re ) * dt / 2;
        end
        if doSee
            % buffer the synaptic currents
            if fast
                See( 1, :, kp ) = See( 1, :, k ) + ( k1see + k2see ) * dt / 2;
                if synp
                    Dee( 1, :, kp )	= Dee( 1, :, k ) + ( k1dee + k2dee ) * dt / 2;
                    Fee( 1, :, kp ) = Fee( 1, :, k ) + ( k1fee + k2fee ) * dt / 2;
                end
            else
                See( :, :, kp ) = See( :, :, k ) + ( k1see + k2see ) * dt / 2;
                if synp
                    Dee( :, :, kp )	= Dee( :, :, k ) + ( k1dee + k2dee ) * dt / 2;
                    Fee( :, :, kp ) = Fee( :, :, k ) + ( k1fee + k2fee ) * dt / 2;
                end
            end
        end
        if doSei
            if Ni > 0
                % buffer the synaptic currents
                if fast
                    Sei( 1, :, kp ) = Sei( 1, :, k ) + ( k1sei + k2sei ) * dt / 2;
                    if synp
                        Dei( 1, :, kp )	= Dei( 1, :, k ) + ( k1dei + k2dei ) * dt / 2;
                        Fei( 1, :, kp )	= Fei( 1, :, k ) + ( k1fei + k2fei ) * dt / 2;
                    end
                else
                    Sei( :, :, kp ) = Sei( :, :, k ) + ( k1sei + k2sei ) * dt / 2;
                    if synp
                        Dei( :, :, kp )	= Dei( :, :, k ) + ( k1dei + k2dei ) * dt / 2;
                        Fei( :, :, kp )	= Fei( :, :, k ) + ( k1fei + k2fei ) * dt / 2;
                    end
                end
            end
        end
    end
    
    if Ni > 0
        Vi( i + 1, : )          = Vi( i, : ) + ( k1vi + k2vi ) * dt / 2;
        Vi( i + 1, : )       	= Vi( i + 1, : ) + dt ^ 0.5 * Dni .* wni( i, : );
        if doSie
            if Ne > 0
                % buffer the synaptic currents
                if fast
                    Sie( 1, :, kp ) = Sie( 1, :, k ) + ( k1sie + k2sie ) * dt / 2;
                    if synp
                        Die( 1, :, kp ) = Die( 1, :, k ) + ( k1die + k2die ) * dt / 2;
                        Fie( 1, :, kp ) = Fie( 1, :, k ) + ( k1fie + k2fie ) * dt / 2;
                    end
                else
                    Sie( :, :, kp ) = Sie( :, :, k ) + ( k1sie + k2sie ) * dt / 2;
                    if synp
                        Die( :, :, kp ) = Die( :, :, k ) + ( k1die + k2die ) * dt / 2;
                        Fie( :, :, kp ) = Fie( :, :, k ) + ( k1fie + k2fie ) * dt / 2;
                    end
                end
            end
        end
        if doSii
            % buffer the synaptic currents
            if fast
                Sii( 1, :, kp ) = Sii( 1, :, k ) + ( k1sii + k2sii ) * dt / 2;
                if synp
                    Dii( 1, :, kp ) = Dii( 1, :, k ) + ( k1dii + k2dii ) * dt / 2;
                    Fii( 1, :, kp ) = Fii( 1, :, k ) + ( k1fii + k2fii ) * dt / 2;
                end
            else
                Sii( :, :, kp ) = Sii( :, :, k ) + ( k1sii + k2sii ) * dt / 2;
                if synp
                    Dii( :, :, kp ) = Dii( :, :, k ) + ( k1dii + k2dii ) * dt / 2;
                    Fii( :, :, kp ) = Fii( :, :, k ) + ( k1fii + k2fii ) * dt / 2;
                end
            end
        end
    end

    % decide upon spiking
    for j                    	= 1 : Ne
        if Ve( i + 1, j ) > Vth_e( j ) && NSe( j ) == 0             	% spike now
            Ve( i + 1, j )      = Vpeak_e( j );
            NSe( j )            = 1;
        elseif NSe( j ) >= 1 && NSe( j ) < dtSDE( j )               	% same spike
            Ve( i + 1, j )      = Vpeak_e( j );
            NSe( j )            = NSe( j ) + 1;
        elseif NSe( j ) >= 1 && NSe( j ) < dtARPE( j )              	% refractory period
            Ve( i + 1, j )      = Vreset_e( j );
            if Inap
                re( i + 1, j )  = rreset_e( j );
            end
            NSe( j )            = NSe( j ) + 1;
        elseif NSe( j ) >= 1 && NSe( j ) >= dtSDE( j ) && NSe( j ) >= dtARPE( j )     % back to dynamics
            Ve( i + 1, j )      = Vreset_e( j );
            if Inap
                re( i + 1, j )  = rreset_e( j );
            end
            NSe( j )            = 0;
        end
    end
    
    for j                       = 1 : Ni
        if Vi( i + 1, j ) > Vth_i( j ) && NSi( j ) == 0                 % spike now
            Vi( i + 1, j )      = Vpeak_i( j );
            NSi( j )            = 1;
        elseif NSi( j ) >= 1 && NSi( j ) < dtSDI( j )                   % same spike
            Vi( i + 1, j )      = Vpeak_i( j );
            NSi( j )            = NSi( j ) + 1;
        elseif NSi( j ) >= 1 && NSi( j ) < dtARPI( j )                  % refractory period
            Vi( i + 1, j )      = Vreset_i( j );
            NSi( j )            = NSi( j ) + 1;
        elseif NSi( j ) >= 1 && NSi( j ) >= dtSDI( j ) && NSi( j ) >= dtARPI( j )     % back to dynamics
            Vi( i + 1, j )      = Vreset_i( j );
            NSi( j )            = 0;
        end
    end
      
end

%-------------------------------------------------------------------------
% prepare output
%-------------------------------------------------------------------------

% digitize spikes (algorithm suitable also for wider spikes)
st                              = sparse( nt, N );
for j                           = 1 : Ne
    pst1                        = find( Ve( :, j ) > 0 );                	% threshold crossings
    mat                         = local_parse( pst1 );
    sidx                        = ceil( mean( mat, 2 ) );
    st( sidx, j )               = 1;
end
for j                           = 1 : Ni
    pst1                        = find( Vi( :, j ) > 0 );
    mat                         = local_parse( pst1 );
    sidx                        = ceil( mean( mat, 2 ) );
    st( sidx, Ne + j )          = 1;
end

%-------------------------------------------------------------------------
% graphics
%-------------------------------------------------------------------------
if ~graphics
    return
end

figure
subplot( 3, 1, 1 )
spc                             = mean( [ max( Ie ) - min( Ie ) max( Ii ) - min( Ii ) ] );
[ ~, ph ]                       = plotTraces( t, [ Ie Ii ], -spc, 1, [ NaN NaN ] );
set( ph( 1 : Ne ), 'color', colors( 1, : ) )
set( ph( Ne + 1 : N ), 'color', colors( 2, : ) )
ylims = ylim;
set( gca, 'ylim', ylims + [ -1 1 ] )
ylabel( 'I [\muA/cm^2]' )
set( gca, 'tickdir', 'out', 'box', 'off' );

subplot( 3, 1, 2 )
spc                             = [ Vpeak_e - Vreset_e Vpeak_i - Vreset_i ];
[ ~, ph ]                       = plotTraces( t, [ Ve Vi ], -spc, 1, [ NaN NaN ] );
set( ph( 1 : Ne ), 'color', colors( 1, : ) )
set( ph( Ne + 1 : N ), 'color', colors( 2, : ) )
ylabel( 'V [mV]' )
set( gca, 'tickdir', 'out', 'box', 'off' );

subplot( 3, 1, 3 )
if Ne > 0
    rh                          = plot_raster( st( :, 1 : Ne ), t );
    set( rh, 'color', colors( 1, : ) )
    hold on
end
if Ni > 0
    rh                          = plot_raster( st( :, ( Ne + 1 ) : N ), t, Ne + 1 );
    set( rh, 'color', colors( 2, : ) )
end
ylim( [ 0 N + 1 ] )
set( gca, 'tickdir', 'out', 'box', 'off' );
xlabel( 't [ms]' )

return % NR_einet

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