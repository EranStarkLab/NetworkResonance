% NR_sinusoids_to_cmodel        sinusoid input to single cell/network models
%
% call              [ cohs, fo, fbins, frate_vec, s, sei, z, zphs, zfvals ] = NR_sinusoids_to_cmodel( model )
%
% gets              model           {'LIF'}             'Inap', 'LIF_C', 'Inap_syn', and 'I-cell'
%        
% optional arguments (given as name/value pairs)
%
%                   Tlong           {3}     [s]         duration of every 'trial'
%                   Dn              {0}     [mV]        noise SD
%                   fvals           {1:40}  [Hz]        frequencies of input current
%                                                           if a two-element vector, will generate a chirp between the two frequencies
%                                                           otherwise will generate a series of sines at the discrete frequencies
%                   nreps           {1}                 number of trials at each frequency
%                   Amp             {}      [muA/cm^2]	'LIF': 0.115    'LIF_C': 8  'LIF_syn': 1    'Inap': 0.15    'Inap_syn': 0.15    'I-cell': 1.8;
%                   Ibias           {}      [muA/cm^2]	'LIF': 0.9      'LIF_C': -3 'LIF_syn': 0    'Inap': -1.85   'Inap_syn': -1.9    'I-cell': -2.2;
%                   prepad          {1}     [s]         pad with 1 s of zero current (to allow model transients to settle)
%            
%           biophysical parameters - all models:
%
%                   C               {1}     [muF/cm^2]  membrane capacitance (all models)
%
%           biophysical parameters - LIF: 
%
%                   Gl_LIF          {0.1}   [mS/cm^2]   leak conductance
%                   El_LIF          {-60}   [mV]        leak reversal potential 
%                   Vth_LIF         {-50}   [mV]        spiking threshold
%                   Vreset_LIF      {-60}   [mV]        post-spike reset
%
%           biophysical parameters - LIF with voltage-dependent Calcium dynamics: 
%
%                   Gl_LIF_C        {0.5}   [mS/cm^2]   leak conductance
%                   El_LIF_C        {-60}   [mV]        leak reversal potential
%                   Vth_LIF_C       {-50}   [mV]        spiking threshold 
%                   Vreset_LIF_C    {-70}   [mV]        post-spike reset 
%                   Gc_LIF_C        {0.08}  [mS/cm^2]   Calcium conductance
%
%           biophysical parameters - LIF with short-term synaptic plasticity: 
%
%                   Gl_LIF_syn      {0.1}   [mS/cm^2]   leak conductance
%                   El_LIF_syn      {-65}   [mV]        leak reversal potential
%                   Vth_LIF_syn     {-62}   [mV]        spiking threshold 
%                   Vreset_LIF_syn  {-70}   [mV]        post-spike reset 
%                   Gs_LIF_syn     	{0.2}   [mS/cm^2]   synaptic conductance
%                   Es_LIF_syn  	{0}     [mV]        synaptic reversal potential 
%                   synmodel_LIF_syn                    {'DF'}, 'D', 'F', or 'null'
%                   short_spike_LIF_syn                 {1}: SD_In, 0.1 ms; D/F rise times, 0.1/0.2 ms
%                                                        0 : SD_In, 1 ms; D/F rise times, 1/2 ms
%                   SD_Out_LIF_syn  {0.1}   [ms]        SD_Out
%
%           biophysical parameters - Inap: 
%
%                   Gp_Inap         {0.1}   [mS/cm^2]   persistent sodium conductance
%                   Gh_Inap         {1}     [mS/cm^2]   h-current conductance
%                   Vth_Inap        {-50}   [mV]        threshold (IF) for spiking
% 
%           biophysical parameters - E/I network (E are 'LIF' or 'Inap' only; I are always LIF)
% 
%                   einet_source    {'E'}   []          'E', 'I', or 'EI'
%                   N_cells         {0}                 number of cells; can be a two-element vector [ E_cells I_cells ]
%                                                           or the total count (argument to einet)
%                   einet_maxGB     {2}     [GB]        'LIF' net requires [ 0.08 7.5 745 ] GB of RAM for 10, 100, or 1000 cells
%                   einet_thin      {1}                 thin uses linear RAM, does not save synaptic variable 
%                   einet_fast      {1}                 fast uses linera computations, does not support different time constants for different synapses
%                   Dn_e            {}      [mV]        default values are: 'LIF': Dn; 'Inap': 0.0125
%                   Dn_i            {}      [mV]        default values are: 'LIF': Dn; 'Inap': 3
%                   Ibias_e         {}      [muA/cm^2]	default values are: 'LIF': Ibias; 'Inap': Ibias
%                   Ibias_i         {}      [muA/cm^2]	default values are: 'LIF': 0; 'Inap': -1
%                   Gie             {}      [mS/cm^2]   default values are: 'LIF': 0.01 or 0.1 (above/below 10 cells)
%                                                                           'Inap': 0.05 or 0.12 (above/below 10 cells)
%                   Gei             {0}     [mS/cm^2]   I-to-E synaptic conductance (all pairs)
%                   Gee             {0}     [mS/cm^2]   E-to-E synaptic conductance (all pairs)
%                   Gii             {0.05}  [mS/cm^2]   I-to-I synaptic conductance (all pairs)
%
%           analysis arguments:
%
%                   removeSpikes    {1}                 flag; replace spikes with Vth before computing impedance
%                   blackOutDuration {1}    [s]         remove first 1 s of output (to remove for model transients; see prepad)
%                   M               {1}     [1/Hz]      spectral resolution (argument to st_coherence and st_fingerprint)
%                   nFFT            {NaN}               spectral resolution (argument to st_coherence and st_fingerprint)
%                   dflag           {0}                 detrending before spectral estimation (argument to st_coherence)
%                   aFs             {1250}  [Hz]        downsampling before spectral estimation (argument to st_coherence and st_fingerprint)
%
%           graphics:
%
%                   graphics        { [1 1 1 ] }      	of the single-cell simulation: [ summary, rate, coherence ] 
%                   graphics_traces {0}                 plot traces (NR_syntransmit) of every repetition/trial (single-cell simulation)  
%                   graphics_ei     {[1 1]}             of the network simulation: [ simulation (einet), coherence summary ]
%
% returns           cohs                    coherence (st_coherence output)
%                   fo          [Hz]        frequencies coherence is evaluated on (st_coherence output)
%                   fbins       [Hz]        frequencies firing rate is evaluated on (st_fingerprint)
%                   frates_vec  [spks/s]    firing rate 
%
%                   s                       summary structure for single-cell (see code) 
%                   sei                     summary structure for network (see code)
%
%                   z           [MOhm/mm^2] impedance of subthreshold
%                   zphs        [rad]       phase of subthreshold
%                   zfvals      [Hz]        frequencies of subtreshold
% 
%                   fig         handles to figures corresponding to 'graphics'
%                   fig_ei      handles to figures corresponding to 'graphics_ei'
%
% calls             ParseArgPairs
%                   NR_einet, NR_syntransmit
%                   st_coherence, st_fingerprint

% 17-feb-21 ES

% last update
% 01-jul-22

function [ cohs, fo, fbins, frate_vec, s, sei, z, zphs, zfvals, fig, fig_ei ] = NR_sinusoids_to_cmodel( model, varargin )

% [C] = muF/cm^2
% [I] = muA/cm^2
% [G] = mS/cm^2
% [t] = ms
% [V] = mV

%------------------------------------------------------------------------
% constants
%------------------------------------------------------------------------
Fs                              = 1000/0.1;                                 % [Hz]

%------------------------------------------------------------------------
% defaults
%------------------------------------------------------------------------
model_DFLT                      = 'LIF';

% simulation
Tlong_DFLT                      = 3;                                        % [s]
Dn_DFLT                         = 0;                                        % [mV]
fvals_DFLT                      = 1 : 1 : 40;                               % [Hz]
nreps_DFLT                      = 1;
prepad_DFLT                     = 1;                                        % [s]

% general parameters
C_DFLT                          = 1;                                        % [muF/cm^2]

%---------------------------------
% LIF parameters
Gl_LIF_DFLT                     = 0.1;                                      % [mS/cm^2]
El_LIF_DFLT                     = -60;                                      % [mV]
Vth_LIF_DFLT                  	= -50;                                      % [mV]
Vreset_LIF_DFLT               	= -60;                                      % [mV]

% LIF input
Amp_LIF                         = 0.115;                                    % [muA/cm^2]
Ibias_LIF                       = 0.90;                                     % [muA/cm^2]

%---------------------------------
% LIF_C parameters
Gl_LIF_C_DFLT                   = 0.5;                                      % [mS/cm^2]
El_LIF_C_DFLT                   = -60;                                      % [mV]
Vth_LIF_C_DFLT                  = -50;                                      % [mV]
Vreset_LIF_C_DFLT               = -70;                                      % [mV]
Gc_LIF_C_DFLT                   = 0.08;                                     % [mS/cm^2]

% LIF_C input
Amp_LIF_C                       = 8;                                        % [muA/cm^2]
Ibias_LIF_C                     = -3;                                       % [muA/cm^2]

%---------------------------------
% LIF_syn parameters:

Gl_LIF_syn_DFLT              	= 0.1;                                      % [mS/cm^2], leak conductance
El_LIF_syn_DFLT                 = -65;                                      % [mV], leak reversal potential

Vth_LIF_syn_DFLT                = -62;                                      % [mV]
Vreset_LIF_syn_DFLT            	= -70;                                      % [mV]

Es_LIF_syn_DFLT                 = 0;                                        % [mV], synaptic reversal potential 
Gs_LIF_syn_DFLT                 = 0.2;                                      % [mS/cm^2], synaptic conductance for subthreshold
synmodel_LIF_syn_DFLT          	= 'DF';

% LIF_syn input:
short_spike_LIF_syn_DFLT        = 1;                                        % logical; short: 0.1 ms; not short: 1 ms
SD_In_LIF_syn_short             = 0.1;                                      % [ms]
tau_rise_D_LIF_syn_short        = [];                                       % [ms]; 0.1 for D/DF; inf for F/null
tau_rise_F_LIF_syn_short        = [];                                       % [ms]; 0.2 for F/DF; inf for D/null
SD_In_LIF_syn_long              = 1;                                        % [ms]
tau_rise_D_LIF_syn_long         = 1;                                        % [ms]
tau_rise_F_LIF_syn_long         = 2;                                        % [ms]
Amp_LIF_syn                     = 1;                                        % [muA/cm^2]
Ibias_LIF_syn                   = 0;                                        % [muA/cm^2]
SD_Out_LIF_syn_DFLT             = 0.1;                                      % [ms]

%---------------------------------
% Inap parameters
Gp_Inap_DFLT                    = 0.1;                                      % [mS/cm^2]
Gh_Inap_DFLT                    = 1;                                        % [mS/cm^2]
Vth_Inap_DFLT                   = -50;                                      % [mV]

% Inap input
Amp_Inap_sines                  = 0.115;                                    % [muA/cm^2]
Amp_Inap_chirp                  = 0.15;                                     % [muA/cm^2]
Ibias_Inap                      = -1.85;                                    % [muA/cm^2]

%---------------------------------
% Inap_syn parameters:
Gs_Inap_syn_DFLT                = 0.02;                                     % [mS/cm^2]

% Inap_syn input:
Amp_Inap_syn_sines              = 0.115;                                    % [muA/cm^2]
Amp_Inap_syn_chirp              = 0.15;                                     % [muA/cm^2]
Ibias_Inap_syn                  = -1.9;                                     % [muA/cm^2]
        
%---------------------------------
% I-cell input: 
Amp_I_cell                      = 1.8;                                      % [muA/cm^2]
Ibias_I_cell                    = -2.2;                                     % [muA/cm^2]

%---------------------------------
% NR_einet:
einet_source_DFLT               = 'E';
N_cells_DFLT                    = 0;
einet_maxGB_DFLT                = 2;                                        % [GB]
einet_thin_DFLT                 = 1;
einet_fast_DFLT                 = 1;
Dn_e_Inap_DFLT                  = 0.0125;                                   % [mV]
Dn_i_Inap_DFLT                  = 3;                                        % [mV]
Ibias_i_LIF_DFLT                = 0;                                        % [muA/cm^2]
Ibias_i_Inap_DFLT               = -1;                                       % [muA/cm^2]
Gei_DFLT                        = 0;                                        % [mS/cm^2]
Gee_DFLT                        = 0;                                        % [mS/cm^2]
Gii_DFLT                        = 0.05;                                     % [mS/cm^2]

% analysis parameters
blackOutDuration_DFLT           = 1;                                        % [s]
M_DFLT                          = 1;                                        % [1/Hz]
mtNW_DFLT                       = 3;
nFFT_DFLT                       = NaN;                                      % frequency steps will be Fs / nFFT
dflag_DFLT                      = [];                                       % [], no removal; 0, 'constant'; 1, 'linear'
dflag_spk_DFLT                  = 0;
aFs_DFLT                        = 1250;                                     % [Hz]; increase (up to Fs) if input contains frequencies above 40 Hz

%------------------------------------------------------------------------
% arguments
%------------------------------------------------------------------------
nargs                           = nargin;
if nargs < 1 || isempty( model )
    model                       = model_DFLT;
end
if ~ismember( model, { 'LIF', 'Inap', 'LIF_C', 'I-cell', 'Inap_syn', 'LIF_syn' } )
    error( 'unsupported model %s', model )
end
ca                              = varargin;
dflag_arg                       = NaN;
for i = 1 : length( ca )
    if isequal( lower( ca{ i } ), 'dflag' )
        dflag_arg               = ca{ i + 1 };
        break
    end
end
[ Tlong, Dn, fvals, nreps ...
    , C ...
    , Gl_LIF,   El_LIF,   Vth_LIF,   Vreset_LIF ...
    , Gl_LIF_C, El_LIF_C, Vth_LIF_C, Vreset_LIF_C, Gc_LIF_C ...
    , Gl_LIF_syn, El_LIF_syn, Vth_LIF_syn, Vreset_LIF_syn, Gs_LIF_syn, Es_LIF_syn ...
    , short_spike_LIF_syn, synmodel_LIF_syn, SD_Out_LIF_syn ...
    , Gp_Inap, Gh_Inap, Vth_Inap ...
    , Gs_Inap_syn ...
    , Amp, Ibias, prepad ...
    , einet_source, N_cells, einet_maxGB, einet_thin, einet_fast ...
    , Dn_e, Dn_i, Ibias_e, Ibias_i, Gie, Gei, Gee, Gii ...
    , removeSpikes, blackOutDuration ...
    , M, mtNW, nFFT, aFs ...
    , doFingerprint, doCoherence ...
    , verbose, graphics_traces, graphics, graphics_ei ...
    ]                           = ParseArgPairs(...
    { 'Tlong', 'Dn', 'fvals', 'nreps' ...
    , 'C' ...
    , 'Gl_LIF',   'El_LIF',   'Vth_LIF',   'Vreset_LIF' ...
    , 'Gl_LIF_C', 'El_LIF_C', 'Vth_LIF_C', 'Vreset_LIF_C', 'Gc_LIF_C' ...
    , 'Gl_LIF_syn', 'El_LIF_syn', 'Vth_LIF_syn', 'Vreset_LIF_syn', 'Gs_LIF_syn', 'Es_LIF_syn' ...
    , 'short_spike_LIF_syn', 'synmodel_LIF_syn', 'SD_Out_LIF_syn' ...
    , 'Gp_Inap', 'Gh_Inap', 'Vth_Inap' ...
    , 'Gs_Inap_syn' ...
    , 'Amp', 'Ibias', 'prepad' ...
    , 'einet_source', 'N_cells', 'einet_maxGB', 'einet_thin', 'einet_fast' ...
    , 'Dn_e', 'Dn_i', 'Ibias_e', 'Ibias_i', 'Gie', 'Gei', 'Gee', 'Gii' ...
    , 'removeSpikes', 'blackOutDuration' ...
    , 'M', 'mtNW', 'nFFT', 'aFs' ...
    , 'doFingerprint', 'doCoherence' ...
    , 'verbose', 'graphics_traces', 'graphics', 'graphics_ei' }...
    , { Tlong_DFLT, Dn_DFLT, fvals_DFLT, nreps_DFLT ...
    , C_DFLT ...
    , Gl_LIF_DFLT,   El_LIF_DFLT,   Vth_LIF_DFLT,   Vreset_LIF_DFLT ...
    , Gl_LIF_C_DFLT, El_LIF_C_DFLT, Vth_LIF_C_DFLT, Vreset_LIF_C_DFLT, Gc_LIF_C_DFLT ...
    , Gl_LIF_syn_DFLT, El_LIF_syn_DFLT, Vth_LIF_syn_DFLT, Vreset_LIF_syn_DFLT, Gs_LIF_syn_DFLT, Es_LIF_syn_DFLT ...
    , short_spike_LIF_syn_DFLT, synmodel_LIF_syn_DFLT, SD_Out_LIF_syn_DFLT ...
    , Gp_Inap_DFLT, Gh_Inap_DFLT, Vth_Inap_DFLT ...
    , Gs_Inap_syn_DFLT ...
    , [], [], prepad_DFLT ...
    , einet_source_DFLT, N_cells_DFLT, einet_maxGB_DFLT, einet_thin_DFLT, einet_fast_DFLT ....
    , [], [], [], [], [], Gei_DFLT, Gee_DFLT, Gii_DFLT ...
    , 1, blackOutDuration_DFLT ...
    , M_DFLT, mtNW_DFLT, nFFT_DFLT, aFs_DFLT ...
    , 1, 1 ...
    , 1, 0, 1, 1 } ...
    , varargin{ : } );

% supported synaptic models
if ~ismember( synmodel_LIF_syn, { 'null', 'DF', 'F', 'D' } )
    error( 'unsupported synmodel %s', synmodel )
end

% spike input
if isnan( dflag_arg )
    if ismember( model, { 'Inap_syn', 'LIF_syn' } )
        dflag                   = dflag_spk_DFLT;
    else
        dflag               	= dflag_DFLT;
    end
else
    dflag                       = dflag_arg;
end
if ismember( model, { 'Inap_syn', 'LIF_syn' } )
    spkInput                    = 1;
else
    spkInput                    = 0;
end

% graphics arguments
nfigs                           = 3;
if length( graphics ) == 1
    graphics                    = graphics * ones( nfigs, 1 );
end
ng                              = length( graphics );
if ng < nfigs 
    graphics                    = [ graphics( : ); zeros( nfigs - ng, 1 ) ];
end
fig                             = zeros( nfigs, 1 );

nfigs_ei                        = 2;
if length( graphics_ei ) == 1
    graphics_ei                 = graphics_ei * ones( nfigs_ei, 1 );
end
ng                              = length( graphics_ei );
if ng < nfigs_ei 
    graphics_ei                 = [ graphics_ei( : ); zeros( nfigs_ei - ng, 1 ) ];
end
fig_ei                          = zeros( nfigs_ei, 1 );

% einet bias, noise, and connectivity
if isempty( Ibias_e )
    switch model
        case 'LIF'
            Ibias_e             = Ibias;
        case 'Inap'
            Ibias_e             = Ibias;
    end
end
if isempty( Ibias_i )
    switch model
        case 'LIF'
            Ibias_i             = Ibias_i_LIF_DFLT;
        case 'Inap'
            Ibias_i             = Ibias_i_Inap_DFLT;
    end
end

if isempty( Dn_e )
    switch model
        case 'LIF'
            Dn_e                = Dn;
        case 'Inap'
            Dn_e                = Dn_e_Inap_DFLT;
    end
end
if isempty( Dn_i )
    switch model
        case 'LIF'
            Dn_i                = Dn;
        case 'Inap'
            Dn_i                = Dn_i_Inap_DFLT;
    end
end

if isempty( Gie )
    switch model
        case 'LIF'
            if sum( N_cells ) >= 10
                Gie           	= 0.01;
            else
                Gie            	= 0.1;
            end
        case 'Inap'
            if sum( N_cells ) >= 10
                Gie            	= 0.05;
            else
                Gie            	= 0.12;
            end
    end
end

% set input parameters and generate chirps
dt                              = 1 / Fs;                                   % OK for LIF and for I-cell with Dn < 5
nprepad                         = length( dt : dt : prepad );
T                               = Tlong + prepad;
t0                              = ( dt : dt : Tlong )';
t                               = ( dt : dt : T )';
nfvals                          = length( fvals );
if nfvals == 2
    sig                         = 'chirp';
    f0                          = fvals( 1 );
    f1                          = fvals( 2 );
    x0                          = cos( pi + 2 * pi * f0 * t0 + pi * ( f1 - f0 ) * t0.^2 / Tlong );
    x0                          = ( x0 + 1 ) / 2;
    x0                          = [ zeros( nprepad, 1 ); x0 ];
    nfvals                      = nreps;
else
    sig                         = 'sines';
end
fROI                            = [ min( fvals ) max( fvals ) ];            % [Hz] determine frequency resolution
if fROI( 1 ) == 1
    fROI( 1 )                   = 0;
end
fvals                           = fvals( : );
aFs                             = min( aFs, Fs );

%------------------------------------------------------------------------
% part 1: go over trials (frequencies, or chirp repetitions) and run single-cell simulation
%------------------------------------------------------------------------

%----------------------------------------------------------------------
% part 1.1: choose the proper parameters
%----------------------------------------------------------------------

switch model
    case { 'Inap' }
        if isempty( Amp )
            switch sig
                case 'chirp'
                    Amp         = Amp_Inap_chirp;
                case 'sines'
                    Amp         = Amp_Inap_sines;
            end
        end
        if isempty( Ibias )
            Ibias           	= Ibias_Inap;
        end
        Vth                     = Vth_Inap;
    case 'Inap_syn'
        if isempty( Amp )
            switch sig
                case 'chirp'
                    Amp         = Amp_Inap_syn_chirp;
                case 'sines'
                    Amp         = Amp_Inap_syn_sines;
            end
        end
        if isempty( Ibias )
            Ibias           	= Ibias_Inap_syn;
        end
        Vth                     = Vth_Inap;
    case 'LIF'
        if isempty( Amp )
            Amp                 = Amp_LIF;
        end
        if isempty( Ibias )
            Ibias               = Ibias_LIF;
        end
        Vth                     = Vth_LIF;
    case 'LIF_C'
        if isempty( Amp )
            Amp                 = Amp_LIF_C;
        end
        if isempty( Ibias )
            Ibias               = Ibias_LIF_C;
        end
        Vth                     = Vth_LIF_C;
    case 'LIF_syn'
        if isempty( Amp )
            Amp                 = Amp_LIF_syn;
        end
        if isempty( Ibias )
            Ibias               = Ibias_LIF_syn;
        end
        Vth                     = Vth_LIF_syn;
        if short_spike_LIF_syn
            SD_In_LIF_syn       = SD_In_LIF_syn_short;
            tau_rise_D_LIF_syn  = tau_rise_D_LIF_syn_short;
            tau_rise_F_LIF_syn  = tau_rise_F_LIF_syn_short;
        else
            SD_In_LIF_syn       = SD_In_LIF_syn_long;
            tau_rise_D_LIF_syn  = tau_rise_D_LIF_syn_long;
            tau_rise_F_LIF_syn  = tau_rise_F_LIF_syn_long;
        end
   
    case 'I-cell'
        if isempty( Amp )
            Amp              	= Amp_I_cell;
        end
        if isempty( Ibias )
            Ibias            	= Ibias_I_cell;
        end
        Vth                     = Vth_LIF;
    otherwise
        error( 'unsupported model' )
end

% initialize output
Nspk                            = NaN( nfvals, 1 );
Ncyc                            = NaN( nfvals, 1 );
z                               = NaN( nfvals, 1 );
zphs                            = NaN( nfvals, 1 );
zfvals                          = NaN( nfvals, 1 );
stmat                           = cell( nfvals, 1 );
stmat_in                        = cell( nfvals, 1 );
kidx                            = t > blackOutDuration;
nt                              = sum( kidx );
t                               = t( kidx ) - blackOutDuration;             % [s]
Iemat                           = NaN( nt, nfvals );
Vmat                            = NaN( nt, nfvals );
model_str                       = model;
tfvals                          = fvals;
mvals                           = NaN( length( tfvals ), 1 );
if spkInput
    switch sig
        case 'sines'
            tfvals              = fvals;
            mvals           	= NaN( length( tfvals ), 1 );
        case 'chirp'
            tfvals              = max( min( f0, f1 ), 1 ) : 1 : max( f0, f1 );
            fedges           	= tfvals( 1 ) - 0.5 : 1 : tfvals( end ) + 0.5;
            mvals             	= NaN( length( tfvals ), 1 );
    end
    if isequal( model, 'LIF_syn' )
        model_str               = sprintf( '%s, %s', model, synmodel_LIF_syn );
    end
end

%----------------------------------------------------------------------
% part 1.2: run the actual single-cell simulations
%----------------------------------------------------------------------
if verbose
    fprintf( 1, 'Simulating single-cell' )
    if isequal( sig, 'sines' )
        fprintf( 1, '; f [Hz] = ' )
    else
        fprintf( 1, ' with [ f0 f1 ] = [ %0.2g %0.2g ] Hz... ', f0, f1 )
    end
end

for i                           = 1 : nfvals
    
    % scale the input by the amplitude
    if isequal( sig, 'sines' )
        f                       = fvals( i );
        if verbose
            fprintf( 1, '%0.3g ', f );
        end
        Ie                      = [ zeros( nprepad, 1 ); ( sin( 2 * pi * f * t0 - pi / 2 ) + 1 ) / 2 * Amp ];
    else
        f                       = ( f0 + f1 ) / 2;
        Ie                      = x0 * Amp;
    end
    
    % run the single-cell simulation
    st_in                       = [];
    switch model
        
        case 'Inap'
            [ st, V, tt  ]      = NR_syntransmit( [], 'model', 'Inap', 'Iext', Ie, 'Tmax', T * 1000, 'Dn', Dn ...
                , 'graphics', graphics_traces, 'Ibias', Ibias ...
                , 'Gp_Inap', Gp_Inap, 'Gh_Inap', Gh_Inap, 'Vth_Inap', Vth_Inap, 'C', C );
            
        case 'LIF'
            [ st, V, tt  ]      = NR_syntransmit( [],  'model', 'LIF', 'Iext', Ie, 'Tmax', T * 1000, 'Dn', Dn ...
                , 'graphics', graphics_traces, 'ARP', 0, 'SD_Out', 1, 'Ibias', Ibias ...
                , 'Gl_LIF', Gl_LIF, 'El_LIF', El_LIF, 'Vreset_LIF', Vreset_LIF, 'Vth_LIF', Vth_LIF );
            
        case 'LIF_C'
            [ st, V, tt  ]      = NR_syntransmit( [],  'model', 'LIF_C', 'Iext', Ie, 'Tmax', T * 1000, 'Dn', Dn ...
                , 'graphics', graphics_traces, 'ARP', 0, 'SD_Out', 1, 'Ibias', Ibias ...
                , 'Gl_LIF', Gl_LIF_C, 'El_LIF', El_LIF_C, 'Vreset_LIF', Vreset_LIF_C, 'Vth_LIF', Vth_LIF_C ...
                , 'Gc_LIF_C', Gc_LIF_C );
            
        case 'I-cell'
            [ st, V, tt  ]      = NR_syntransmit( [], 'model', 'I-cell', 'Iext', Ie, 'Tmax', T * 1000, 'Dn', Dn ...
                , 'graphics', graphics_traces, 'Ibias', Ibias );
            
        case { 'LIF_syn', 'Inap_syn' }
            tmp                 = local_parse( find( ( Ie ) / max( Ie ) > 0.95 ) );
            st_in               = cell( 1 );
            st_in{ 1 }          = ceil( mean( tmp, 2 ) * dt * 1000 );       % one deterministic train, no phase offset
            if isequal( model, 'LIF_syn' )
                [ st, V, tt ]   = NR_syntransmit( st_in, 'model', 'LIF', 'Iext', [], 'Tmax', T * 1000, 'Dn', Dn ...
                    , 'synmodel', synmodel_LIF_syn, 'tau_rise_D', tau_rise_D_LIF_syn, 'tau_rise_F', tau_rise_F_LIF_syn ...
                    , 'Gs', Gs_LIF_syn, 'Es', Es_LIF_syn ...
                    , 'Gl_LIF', Gl_LIF_syn, 'El_LIF', El_LIF_syn, 'Vreset_LIF', Vreset_LIF_syn, 'Vth_LIF', Vth_LIF_syn ...
                    , 'SD_In', SD_In_LIF_syn, 'SD_Out', SD_Out_LIF_syn ...
                    , 'Ibias', Ibias, 'graphics', graphics_traces, 'verbose', 0 );
            else
                [ st, V, tt  ]  = NR_syntransmit( st_in, 'model', 'Inap', 'Iext', [], 'Tmax', T * 1000, 'Dn', Dn ...
                    , 'Gp_Inap', Gp_Inap, 'Gh_Inap', Gh_Inap, 'Vth_Inap', Vth_Inap, 'C', C ...
                    , 'Gs', Gs_Inap_syn, 'graphics', graphics_traces, 'Ibias', Ibias, 'SD_In', 1, 'verbose', 0 );
            end
            
    end

    % remove transients
    if blackOutDuration ~= 0
        Ie                      = Ie( kidx );
        V                       = V( kidx );
        tt                      = tt( kidx ) - blackOutDuration * 1000;     % [ms]
        if ~isempty( st )
            st                  = st - blackOutDuration * 1000;             % [ms]
            st( st < 0 )        = [];
        end
        if spkInput
            st_in{ 1 }          = st_in{ 1 } - blackOutDuration * 1000;     % [ms]
        end
    end
    
    % accumulate current input and Vm output
    Iemat( :, i )               = Ie;                                       % accumulation of multiple trials/frequencies
    Vmat( :, i )                = V;

    % rudimentary Vm analysis (compute impedance and phase):
    Vz                          = V;
    if removeSpikes
        Vz( Vz > Vth )          = Vth;
    end
    pers                        = local_parse( find( Ie >= ( max( Ie ) - eps ) ) ); % prevent two consecutive nearly-identical values
    pdurs                       = diff( pers, [], 2 ) + 1;
    if any( pdurs > 1 )
        for j                   = 1 : size( pers, 1 )
            idx                 = ( pers( j, 1 ) + 1 ) : pers( j, 2 );
            pval                = Ie( pers( j, 1 ) );
            Ie( idx )           = pval - 2 * eps;
        end
    end
    [ z( i ), zphs( i ) ]       = local_calc_z( tt / 1000, Ie, Vz );        % impedance mV / (muA/cm^2) = KOhm * cm^2
    zfvals( i )                 = f;
    
    % rudimentary subthreshold post-synaptic potential analysis
    if spkInput && isequal( sig, 'sines' )
        
        % analysis of sub-threshold - pooled across all cycles (ignore first cyc0_sth - 1 cycles):
        st0                     = ceil( st_in{ 1 } * Fs / 1000 );
        cyc0_sth                = 3;
        st0                     = st0( cyc0_sth : length( st0 ) );
        nspks                   = length( st0 );
        vals                    = NaN * ones( nspks, 1 );
        for k                   = 1 : nspks
            sti                 = st0( k );
            if k < length( st0 )
                eti             = st0( k + 1 ) - 1;
            else
                eti             = length( V );
            end
            if Es_LIF_syn == 0
                vals( k )       = max( V( sti : eti  ) );
            else
                vals( k )       = min( V( sti : eti  ) );
            end
        end
        mvals( i )              = mean( vals );
        
    elseif spkInput && isequal( sig, 'chirp' )
        
        % analysis of sub-threshold - cycle by cycle:
        % for each input spike - the instantaneous frequency is the mean of
        % the pre-ISI and the post-ISI.
        st0                     = st_in{ 1 };
        nst                     = length( st0 );
        ds                      = diff( [ 0; st0; T * 1000 ] );
        fs                      = 1000 ./ ( ( ds( 1 : nst ) + ds( 2 : nst + 1 ) ) / 2 );
        fs( nst )               = [];
        nst                     = nst - 1;                                	% ignore the last spike
        % compute the min/max after each spike
        vals                = NaN * ones( nst, 1 );
        for k                   = 1 : nst
            if k == nst
                kidx            = round( st0( k ) * Fs / 1000 ) : length( V );
            else
                kidx            = round( st0( k ) * Fs / 1000 ) : round( st0( k + 1 ) * Fs / 1000 );
            end
            vk                  = V( kidx );
            if Es_LIF_syn == 0
                vals( k )       = max( vk );
            else
                vals( k )       = min( vk );
            end
        end
        % now re-bin the frequencies
        [ ~, bidx ]             = histc( fs, fedges );
        ub                      = unique( bidx( bidx ~= 0 ) );
        for k                   = 1 : length( ub )
            mvals( ub( k ) )    = mean( vals( bidx == ub( k ) ) );
        end
    end % synaptic
    
    
    % rudimentary spiking analysis (spike and cycle counts):
    Nspk( i )                   = length( st );                             % number of spikes overall
    ncycles                     = T * f;
    Ncyc( i )                   = Nspk( i ) / ncycles;                      % number of spikes/cycle
    stmat{ i }                  = st;
    if spkInput
        stmat_in{ i }           = st_in{ 1 };
    end
    
end
z                               = z / 10;                                   % [KOhm * cm^2] -> [MOhm * mm^2]
if verbose
    fprintf( 1, 'done single-cell!\n' )
end

%------------------------------------------------------------------------
% part 2.1: graphical summary of the simulation
%------------------------------------------------------------------------
if graphics( 1 )

    fig( 1 )                    = figure;
    subplot( 3, 1, 1 )
    switch sig
        case 'sines'
            plot( fvals, Nspk / T, '.-b' )
            ylabel( 'Spikes/s' )
        case 'chirp'
            plot( tt / 1000, V, 'b' );
            ylabel( 'V [mV]' )
    end
    set( gca, 'tickdir', 'out', 'box', 'off' )
    tstr                        = sprintf( '%s: Amplitude, bias: %0.3g, %0.3g \\muA/cm^2; Noise=%0.2g mV; T=%d s' ...
        , model_str, Amp, Ibias, Dn, ceil( T ) );
    local_fig_title( tstr )

    subplot( 3, 1, 2 )
    switch sig
        case 'sines'
            plot( zfvals, z, '.-b' )
            ylabel( 'Impedance [M\Omega mm^2]' )
            xlabel( 'Frequency [Hz]' )
        case 'chirp'
            plot( tt / 1000, Ie + Ibias, 'r' )
            ylabel( 'I' )
            xlabel( 'Time [s]' )
    end
    set( gca, 'tickdir', 'out', 'box', 'off' )

    if isequal( sig, 'sines' )
        subplot( 3, 1, 3 )
        uwphs                   = unwrap( zphs );
        if max( uwphs ) > 2 * pi
            uwphs               = uwphs - 2 * pi;
        end
        plot( zfvals, uwphs, '.-b' )
        line( xlim, [ 0 0 ], 'color', [ 0 0 0 ], 'linestyle', '--' );
        ylim( [ -1 1 ] * pi / 2 )
        set( gca, 'ytick', ( -1 : 0.5 : 1 ) * pi / 2 )
        ylabel( 'Phase [rad]' )
        xlabel( 'Frequency [Hz]' )
        line( xlim, [ 1 1 ] * 2 * pi, 'color', [ 0 0 0 ], 'linestyle', '--' );
        set( gca, 'tickdir', 'out', 'box', 'off' )
    end
    
end

%------------------------------------------------------------------------
% part 2.3: firing rate and coherence analyses
%------------------------------------------------------------------------

% prepare format
stSamples                       = cell( size( stmat ) );
ntrials                         = nfvals;
for i                           = 1 : ntrials
    stSamples{ i }              = round( stmat{ i } * Fs / 1000 );
    cil                         = stSamples{ i } == 0;
    stSamples{ i }( cil )       = 1;
end

Ievec                           = Iemat( : );
Vvec                            = Vmat( : );
mat                             = sparse( size( Iemat, 1 ), size( Iemat, 2 ) );
for i                           = 1 : ntrials
    mat( stSamples{ i }, i )    = 1;
end
stvec                           = find( mat( : ) );

mat_in                          = sparse( size( Iemat, 1 ), size( Iemat, 2 ) );
for i                           = 1 : ntrials
    tmp                         = round( stmat_in{ i } * Fs / 1000 );
    cil                         = tmp == 0;
    tmp( cil )                  = 1;
    mat_in( tmp, i )            = 1;
end
stvec_in                        = find( mat_in( : ) );

% analyze as a single trial
if doFingerprint
    [ frate_vec, ~, ~, fbins ]  = st_fingerprint( stvec, Ievec, 'spkFs', Fs, 'xFs', Fs ...
        , 'fROI', fROI, 'Fs', aFs, 'graphics', graphics( 2 ), 'M', M, 'nFFT', nFFT );
    if graphics( 2 )
        local_fig_title( model_str );
        fig( 2 )                = gcf;
    end
else
    frate_vec                   = [];
    fbins                       = [];
end

if doCoherence
    if spkInput
        [ cohs, phs, fo ...
            , pvals, fidx ]     = st_coherence( stvec, stvec_in, 'spkFs', Fs, 'xFs', Fs, 'xSpk', 1, 'ntmax', length( Ievec ) ...
            , 'fROI', fROI, 'Fs', aFs, 'graphics', graphics( 3 ), 'M', M, 'mtNW', mtNW, 'nFFT', nFFT, 'dflag', dflag );
    else
        [ cohs, phs, fo ...
            , pvals, fidx ]     = st_coherence( stvec, Ievec, 'spkFs', Fs, 'xFs', Fs ...
            , 'fROI', fROI, 'Fs', aFs, 'graphics', graphics( 3 ), 'M', M, 'mtNW', mtNW, 'nFFT', nFFT, 'dflag', dflag );
    end
    if graphics( 3 )
        fig( 3 )                = gcf;
        local_fig_title( model_str );
    end
else
    cohs                        = [];
    phs                         = [];
    fo                          = [];
    pvals                       = [];
    fidx                        = [];
end

% organize output structure
s.frate_vec                     = frate_vec;
s.fbins                         = fbins;
s.cohs                          = cohs;
s.phs                           = phs;
s.fo                            = fo;
s.pvals                         = pvals;
s.fidx                          = fidx;

s.tfvals                        = tfvals;                                   % subthreshold (PSP)
s.mvals                         = mvals;                                    % subthreshold (PSP)

s.z                             = z;                                        % subthreshold (sines only)
s.zphs                          = zphs;                                     % subthreshold (sines only)
s.zfvals                        = zfvals;

s.Fs                            = Fs;                                       % [Hz]
s.Iemat                         = Iemat;                                    % [muA/cm^2]
s.Vmat                          = Vmat;                                     % [mV]
s.t                             = t( : );                                   % [s]
s.mat_in                        = mat_in;                                   % input spike times
s.mat                           = mat;                                      % output spike times 

%------------------------------------------------------------------------
% part 3: use the same input as to the single cell (part 1), apply to an NR_einet
%------------------------------------------------------------------------
if ~sum( N_cells ) > 0 || ~ismember( model, { 'LIF', 'Inap' } ) || ~isequal( sig, 'chirp' )
    sei                         = [];
    return
end
fprintf( 1, 'Simulating %d-cell network...', sum( N_cells ) )

% determine which cells receive the sinusoid input
ze                              = zeros( size( Iemat( : ) ) );
switch einet_source
    case 'E'
        Iext                    = [ Iemat( : ) ze t( : ) * 1000 ];
    case 'I'
        Iext                    = [ ze Iemat( : ) t( : ) * 1000 ];
    case 'EI'
        Iext                    = [ Iemat( : ) Iemat( : ) t( : ) * 1000 ];
end

% run the network simulation
switch model
    case 'LIF'
        [ sr, Ve, Vi, tei ]     = NR_einet( N_cells, Iext, 'Ibias_e', Ibias_e, 'Ibias_i', Ibias_i, 'Dne', Dn_e, 'Dni', Dn_i ...
            , 'Gie', Gie, 'Gei', Gei, 'Gee', Gee, 'Gii', Gii, 'Tmax',  T * 1000 ...
            , 'Gl_e', 0.1, 'Gl_i', 0.1, 'Vreset_e', -60, 'Vreset_i', -60, 'SD_e', 1, 'ARP_e', 0, 'SD_i', 1, 'ARP_i', 0 ...
            , 'maxGB', einet_maxGB, 'thin', einet_thin, 'fast', einet_fast ...
            , 'graphics', graphics_ei( 1 ) );
        
    case 'Inap'
        [ sr, Ve, Vi, tei ]     = NR_einet( N_cells, Iext, 'Ibias_e', Ibias_e, 'Ibias_i', Ibias_i, 'Dne', Dn_e, 'Dni', Dn_i ...
            , 'Gie', Gie, 'Gei', Gei, 'Gee', Gee, 'Gii', Gii, 'Tmax',  T * 1000 ...
            , 'Gl_e', 0.1, 'Gl_i', 0.1, 'Vreset_e', -70, 'Vreset_i', -60, 'El_e', -65, 'El_i', [] ...
            , 'SD_e', 1, 'ARP_e', 1, 'Inap', 1 ...
            , 'maxGB', einet_maxGB, 'thin', einet_thin, 'fast', einet_fast ...
            , 'graphics', graphics_ei( 1 ) );
end
if graphics_ei( 1 )
    fig_ei( 1 )                 = gcf;
end
fprintf( 1, 'done einet!\n' )

% analyze the coherence
cohs_ei                         = NaN( size( cohs, 1 ), sum( N_cells ) );
phs_ei                          = NaN( size( cohs, 1 ), sum( N_cells )  );
max_coh_frq_ei                  = NaN( 1, sum( N_cells ) );
max_coh_val_ei                  = NaN( 1, sum( N_cells ) );
for j                           = 1 : sum( N_cells )
    if sum( sr( :, j ) ) == 0
        continue
    end
    [ coh_j, phs_j ]            = st_coherence( find( sr( :, j ) ) , Iemat( : ), 'spkFs', Fs, 'xFs', Fs ...
        , 'fROI', fROI, 'Fs', aFs, 'graphics', 0, 'M', M, 'mtNW', mtNW, 'nFFT', nFFT, 'dflag', dflag );
    cohs_ei( :, j )             = coh_j;
    phs_ei( :, j )              = phs_j;
    [ ~, maxidx ]               = max( cohs_ei( :, j ) );
    max_coh_frq_ei( j )         = fo( maxidx );
    max_coh_val_ei( j )         = cohs_ei( maxidx, j );
end

% summarize the einet
Ne                              = size( Ve, 2 );
Ni                              = size( Vi, 2 );
isE                             = logical( [ ones( Ne, 1 ); zeros( Ni, 1 ) ] );
if Ne > 0
    mcoh_e                      = mean( cohs_ei( :, isE == 1 ), 2 );
    mphs_e                      = local_circ_mean( phs_ei( :, isE == 1 )' )';
    mphs_e                      = unwrap( mphs_e );
    mphs_e                      = mod( mphs_e + pi, 2 * pi ) - pi;
else
    mcoh_e                      = NaN( size( cohs, 1 ), 1 );
    mphs_e                      = NaN( size( cohs, 1 ), 1 );
end
if Ni > 0
    mcoh_i                      = mean( cohs_ei( :, isE == 0 ), 2 );
    mphs_i                      = local_circ_mean( phs_ei( :, isE == 0 )' )';
    mphs_i                      = unwrap( mphs_i );
    mphs_i                      = mod( mphs_i + pi, 2 * pi ) - pi;
else
    mcoh_i                      = NaN( size( cohs, 1 ), 1 );
    mphs_i                      = NaN( size( cohs, 1 ), 1 );
end

% organize output structure
sei.N_cells                     = [ Ne Ni ];
sei.Dn_ei                       = [ Dn_e Dn_i ];
sei.Gie                         = Gie;
sei.Gei                         = Gei;
sei.Gee                         = Gee;
sei.Gii                         = Gii;
sei.t                           = tei;                                      % [ms]
sei.Ie                          = Iext( :, 1 );
sei.Ii                          = Iext( :, 2 );
sei.Ve                          = Ve;
sei.Vi                          = Vi;
sei.sr                          = sr;
sei.cohs_ei                     = cohs_ei;
sei.phs_ei                      = phs_ei;
sei.mcoh_e                      = mcoh_e;
sei.mcoh_i                      = mcoh_i;
sei.mphs_e                      = mphs_e;
sei.mphs_i                      = mphs_i;
sei.max_coh_frq_ei              = max_coh_frq_ei;
sei.max_coh_val_ei              = max_coh_val_ei;
sei.frates                      = full( sum( sr ) ) / t( end );

% plot
if graphics_ei( 2 )
    
    fig_ei( 2 )                 = figure;
    
    subplot( 2, 1, 1 )
    hold on
    if Ne > 0
        ph                      = plot( fo, cohs_ei( :, isE == 1 ), 'r' );
        set( ph, 'color', [ 1 0.5 0.5 ] )
        ph                      = plot( fo, mcoh_e, 'r' );
        set( ph, 'linewidth', 2 )
    end
    if Ni > 0
        ph                      = plot( fo, cohs_ei( :, isE == 0 ), 'b' );
        set( ph, 'color', [ 0.5 0.5 1 ] )
        ph                      = plot( fo, mcoh_i, 'b' );
        set( ph, 'linewidth', 2 )
    end
    set( ph, 'linewidth', 2 )
    xlim( fROI )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Coherence' )
    title( sprintf( 'nE = %d (Dn\\_e=%0.3g); nI = %d (Dn\\_i=%0.3g)', Ne, Dn_e, Ni, Dn_i ) )
    
    % phase
    subplot( 2, 1, 2 )
    hold on
    uwphs                       = unwrap( phs_ei );
    uwphs                       = mod( uwphs + pi, 2 * pi ) - pi;
    if Ne > 0
        ph                      = plot( fo, uwphs( :, isE == 1 ), 'r' );
        set( ph, 'color', [ 1 0.5 0.5 ] )
        ph                      = plot( fo, mphs_e, 'r' );
        set( ph, 'linewidth', 2 )
    end
    if Ni > 0
        ph                      = plot( fo, uwphs( :, isE == 0 ), 'b' );
        set( ph, 'color', [ 0.5 0.5 1 ] )
        ph                      = plot( fo, mphs_i, 'b' );
        set( ph, 'linewidth', 2 )
    end
    xlim( fROI )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Phase [rad]' )
    line( xlim, [ 0 0 ], 'color', [ 0 0 0 ], 'linestyle', '--' );
    ylim( [ -pi pi ] )
    xlabel( 'Frequency [Hz]' )
    title( sprintf( 'Source=%s; Gie=%0.3g; Gei=%0.3g', einet_source, Gie, Gei ) )
    
end % graphics

return % NR_sinusoids_to_cmodel

%------------------------------------------------------------------------
% [ z, phs ] = local_calc_z( t, I, V )
% impedance and phase during stationary periodic input
%------------------------------------------------------------------------
function [ z, phs ] = local_calc_z( t, I, V )

% constants
dead_time                       = 1;                                     % [s]

% initialize output
z                               = NaN;
phs                             = NaN;

% arguments
t                               = t( : );
I                               = I( : );
V                               = V( : );
if ~isequal( length( t ), length( I ), length( V ) )
    error( 'input size mistmatch' )
end
    
% ignore onset transients
zidx                            = t > dead_time;
tz                              = t( zidx );
vz                              = V( zidx );
iz                              = I( zidx );
if isempty( tz )
    return
end

% find all local extrema
[ iext_idx, iext_val, iext_typ ] = local_find_local_extrema( iz );         	% this will fail if two consecutive samples have same value
pidx                            = find( iext_typ == 1 );
npks                            = length( pidx );
if npks < 2 
    maxi                        = max( iz );
    iext_idx1                   = find( iz == maxi );
    iext_val1                   = maxi * ones( size( iext_idx1 ) );
    iext_typ1                   = ones( size( iext_idx1 ) );
    mini                        = min( iz );
    iext_idx2                   = find( iz == mini );
    iext_val2                   = mini * ones( size( iext_idx2 ) );
    iext_typ2                   = -ones( size( iext_idx2 ) );
    iext_idx                    = [ iext_idx1; iext_idx2 ];
    iext_val                    = [ iext_val1; iext_val2 ];
    iext_typ                    = [ iext_typ1; iext_typ2 ];
    [ ~, sidx ]                 = sort( iext_idx );
    iext_idx                    = iext_idx( sidx );
    iext_val                    = iext_val( sidx );
    iext_typ                    = iext_typ( sidx );
    pidx                        = find( iext_typ == 1 );
    npks                        = length( pidx );
end
if npks < 2
    return
end

% go over segments between local maxima
zi                              = NaN( npks, 1 );
ai                              = NaN( npks, 1 );
slct                            = 2 : npks;
for i                           = slct
    if iext_typ( pidx( i ) - 1 ) ~= -1
        continue
    end
    si                          = iext_idx( pidx( i - 1 ) );                % time of previous peak
    ei                          = iext_idx( pidx( i ) );                    % time of present peak
    maxi                        = iext_val( pidx( i ) );                    % value at present peak
    mini                        = iext_val( pidx( i ) - 1 );                % value at the trough before the present peak
    di                          = maxi - mini;                              % input swing
    ti                          = ei - si;                                  % [samples] cycle duration
    v                           = vz( si : ei );
    [ maxv, tv ]                = max( v );
    minv                        = min( v );
    dv                          = maxv - minv;                              % output swing
    zi( i )                     = dv / di;                                  % |output|/|input|
    ai( i )                     = ( tv - 1 ) / ti * 2 * pi;                 % phase of output peak
end

% average over all segments
z                               = nanmean( zi( slct ) );
phs                             = local_circ_mean( ai( slct ) );

return % local_calc_z

%------------------------------------------------------------------------
% phi = local_circ_mean( t )
% compute mean direction
%------------------------------------------------------------------------
function phi = local_circ_mean( t )

% trigonometric functions
nans                            = isnan( t );
n                               = sum( ~nans, 1 );
x                               = cos( t );
y                               = sin( t );

% compute direction
sumx                            = nansum( x, 1 );
sumy                            = nansum( y, 1 );
C                               = sumx ./ n;
S                               = sumy ./ n;
phi                             = mod( atan2( S, C ), 2 * pi );

return % local_circ_mean

%------------------------------------------------------------------------
% th = local_fig_title( tstr )
% title for the figure
%------------------------------------------------------------------------
function th = local_fig_title( tstr )

idx                             = strfind( tstr, '_' ); 
idx                             = [ idx length( tstr ) + 1 ];
str                             = tstr( 1 : idx( 1 ) - 1 );
for i                           = 1 : ( length( idx ) - 1 )
    str                         = sprintf( '%s%s%s', str, '\_' ...
        , tstr( ( idx( i ) + 1 ) : ( idx( i + 1 ) - 1 ) ) ); 
end

ah0                             = gca;
axes( 'position', [ 0.5 0.95 0.01 0.01 ] );
axis off;
th                              = text( 0.5, 0.5, str, 'Fontsize', 12 );
set( th, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' )
subplot( ah0 );

return % local_fig_title

%------------------------------------------------------------------------
% [ idx, vals, etype ] = local_find_local_extrema( x )
% detect all local extrema in a matrix
%------------------------------------------------------------------------
function [ idx, vals, etype ] = local_find_local_extrema( x )

if isempty( x )
    return
end
sx                              = size( x );
if length( sx ) == 2 && sum( sx == 1 ) >= 1 && sx( 2 ) > 1
    x                           = x';
end
if length( sx ) > 2
    error( 'input size mismatch: x must be a vector or a matrix' )
end
m                               = sx( 1 );
n                               = sx( 2 );

% compute second 'derivative'
d2                              = diff( sign( diff( x ) ) );

% identify extrema
[ rowMin, colMin ]              = find( d2 > 1 );
[ rowMax, colMax ]              = find( d2 < -1 );
etype                           = [ -1 * ones( length( rowMin ), 1 ); ones( length( rowMax ), 1 ) ];
mat                             = [ [ [ rowMin colMin ]; [ rowMax colMax ] ] etype ];
smat                            = sortrows( mat, [ 2 1 ] );
row                             = smat( :, 1 );
col                             = smat( :, 2 );
etype                           = smat( :, 3 );
row                             = row + 1;

% organize output
if n == 1
    idx                         = row;
else
    idx                         = [ row col ];
end
vals                            = x( row + ( col - 1 ) * m );

return % local_find_local_extrema

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