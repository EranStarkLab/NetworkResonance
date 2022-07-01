% NR_sinusoids_to_cmodel_run       wrapper for sinusoids_to_cmodel
%
% call                      [ cohs, frates, fo ] = NR_sinusoids_to_cmodel_run( model, sig, param_name )
%                           [ ..., figs, pvals ] = NR_sinusoids_to_cmodel_run
%
% gets                      model           {'LIF'}, 'LIF_C', 'LIF_syn', 'Inap', 'Inap_syn'
%                           sig             {'chirp'}, 'sines'
%                           param_name      {'Gl'}, 'Vr', 'El', 'Dn', 'Amp', 'Ibias'
%                                           'Gc', 'Gs', 'Es', 'Vth'
% 
% optional arguments (given as name/value pairs)
%
%                           param_vals      {[]}, vector of parameter values
%                           output          {'spiking'}, 'subthreshoild'
%                           nFFT            {NaN} determine spectral resolution
%
% does
% (1) according to the selected model, signal, and parameter, determines
% the grid search parameters. default values can be overriden using param_vals. 
% (2) call sinusoids_to_cmodel for each, saves the coherence and firing rates only
% (3) (optionally) saves the data
% (4) (optionally) plots and saves the plot
%
% calls                     ParseArgPairs
%                           NR_sinusoids_to_cmodel
%
% see also                  st_coherence, st_fingerprint

% 28-jul-21 ES

% last update
% 01-jul-22

function [ cohs, frates, fo, fig, pvals ] = NR_sinusoids_to_cmodel_run( model, sig, param_name, varargin )

%------------------------------------------------------------------------
% default values
%------------------------------------------------------------------------
Dn_DFLT                     	= 0;

El_DFLT                      	= -60;
Gl_DFLT                     	= 0.1;
Vr_DFLT                         = -60;
Vth_DFLT                        = -50;

Gc_DFLT                         = 0;

Gs_DFLT                         = 0.175;
Es_DFLT                         = 0;

synmodel_DFLT                   = 'DF';
dflag_DFLT                      = '';

%------------------------------------------------------------------------
% arguments
%------------------------------------------------------------------------
nargs                           = nargin;
if nargs < 1 || isempty( model )
    model                       = 'LIF';
end
if ~ismember( model, { 'Inap', 'Inap_syn', 'LIF', 'LIF_C', 'LIF_syn' } )
    error( 'unsupported model' )
end
if nargs < 2 || isempty( sig )
    sig                         = 'chirp';
end
if ~ismember( sig, { 'chirp', 'sines' } )
    error( 'unsupported sig' )
end
if nargs < 3 || isempty( param_name )
    param_name                  = 'Gl';
end
[ param_vals, output, nFFT ...
    , Dn, El, Gl, Vr, Vth, Gc, Gs, Es ...
    , synmodel, dflag ...
    , param_fvals, param_Tlong, param_ibias ...
    , SD_Out ...
    , graphics ...
    ]                           = ParseArgPairs(...
    { 'param_vals', 'output', 'nFFT' ...
    , 'Dn', 'El', 'Gl', 'Vr', 'Vth', 'Gc', 'Gs', 'Es' ...
    , 'synmodel', 'dflag' ...
    , 'param_fvals', 'param_Tlong', 'param_ibias' ...
    , 'SD_Out' ...
    , 'graphics' }...
    , { [], 'spiking', NaN ...
    , Dn_DFLT, El_DFLT, Gl_DFLT, Vr_DFLT, Vth_DFLT, Gc_DFLT, Gs_DFLT, Es_DFLT ...
    , synmodel_DFLT, dflag_DFLT ...
    , [], [], [] ...
    , [] ...
    , 1 }...
    , varargin{ : } );

%------------------------------------------------------------------------
% prepare parameters
%------------------------------------------------------------------------

% Tlong, Amp
switch sig
    case 'chirp'
        fvals                   = [ 0 40 ];                                 % [Hz]
        nfz                     = [];
        Tlong                   = 20;                                       % [s]
        switch model
            case { 'Inap', 'Inap_syn' }
                Amp             = 0.15;
            case 'LIF'
                Amp             = 0.115;
            case 'LIF_C'
                Amp             = 8;
            case 'LIF_syn'
                Amp             = 1;
        end
    case 'sines'
        fvals                   = 1 : 40; 
        nfz                     = 40;
        Tlong                   = 3;                                        % [s]
        switch model
            case { 'Inap', 'Inap_syn' }
                Amp             = 0.15;
            case 'LIF'
                Amp             = 0.115;
            case 'LIF_C'
                Amp             = 8;
            case 'LIF_syn'
                Amp             = 1;
        end
end
if ~isempty( param_fvals )
    fvals                       = param_fvals( : );
end
if ~isempty( param_Tlong )
    Tlong                       = param_Tlong( : );
end

% grid
if isnan( nFFT )
    nfc                         = 513;
    nfr                         = 33;
    nfp                         = 41;
else
    nfc                         = ceil( nFFT / 2 ) + 1;
    nfr                         = max( fvals ) + 1;
    nfp                         = sum( ( 0 : nFFT / 2 )' * 1250 / nFFT <= max( fvals ) );
end

% Ibias
switch model
    case 'Inap'
        Ibias                   = -1.85;
    case 'Inap_syn'
        Ibias                   = -1.9;
    case 'LIF'
        Ibias                   = 0.9;
    case 'LIF_C'
        Ibias                   = -3;
    case 'LIF_syn'
        Ibias                   = 1.2;
end
if ~isempty( param_ibias )
    Ibias                       = param_ibias;
end

% summary of default values
if isequal( model, 'LIF_syn' )
    mstr                        = sprintf( '%s, %s', model, synmodel );
else
    mstr                        = model;
end
tstr                            = sprintf( '%s, %s; Amp=%0.3g; I_{bias}=%0.3g; Dn=%0.3g; Gl=%0.3g; El=%0.3g; Vth=%0.3g; Vreset=%0.3g' ...
    , mstr, sig, Amp, Ibias, Dn, Gl, El, Vth, Vr );

% varied values
switch param_name
    case 'Gl'
        Vals                    = 0.08 : 0.001 : 0.103;
        ValStr                  = 'Gl [mS/cm^2]';
        Vcan                    = Gl;
    case 'El'
        Vals                    = -70 : 1 : -51;                            % Vth is -50 mV
        ValStr                  = 'El [mV]';
        Vcan                    = El;
    case 'Vr'
        Vals                    = -70 : 1 : -51;                        	% Vth is -50 mV
        ValStr                  = 'Vreset [mV]';
        Vcan                    = El;
    case 'Dn'
        Vals                    = 0 : 0.005 : 0.3;
        ValStr                  = 'Dn [mV]';
        Vcan                    = Dn;
    case 'Amp'
        switch output
            case 'spiking'
                Vals            = 0.085 : 0.0025 : 0.3;
            case 'subthreshold'
                Vals            = 0.001 : 0.005 : 0.1;
        end
        ValStr                  = 'Amplitude [\muA/cm^2]';
        Vcan                    = Amp;
    case 'Ibias'
        switch model
            case 'Inap'
                Vals            = -2 : 0.05 : -1;                           % default: -1.85
            case 'Inap_syn'
                Vals            = -2 : 0.005 : -1.8;                        % default: -1.9
            case 'LIF'
                Vals            = 0.5 : 0.05 : 1.5;                         % default: 0.9
            case 'LIF_C'
                Vals            = -3.5 : 0.05 : -2.5;                       % default: -3
            case 'LIF_syn'
                Vals            = 1.05 : 0.01 : 1.25;                       % default: 1.2
        end
        ValStr                  = 'Ibias [\muA/cm^2]';
        Vcan                    = Ibias;
    case 'Gc'
        Vals                    = 0.04 : 0.001 : 0.12;
        ValStr                  = 'Gc [mS/cm^2]';
        Vcan                    = Gc;
    case 'Gs'
        Vals                    = 0.1 : 0.005 : 0.3;
        ValStr                  = 'Gs [mS/cm^2]';
        Vcan                    = Gs;
    case 'Es'
        Vals                    = -80 : 5 : 0;                              % Vth is -50 mV
        ValStr                  = 'Es [mV]';
        Vcan                    = Es;
    case 'Vth'
        Vals                    = -59 : 1 : -40;                        	% Vr is -60 mV
        ValStr                  = 'Vth [mV]';
        Vcan                    = Vth;
end
if ~isempty( param_vals )
    Vals                        = param_vals( : );
end
nVals                           = length( Vals );

%------------------------------------------------------------------------
% simulate
%------------------------------------------------------------------------
cohs                            = NaN( nfc, nVals );
frates                          = NaN( nfr, nVals );
zmat                            = NaN( nfz, nVals );
zphs                            = NaN( nfz, nVals );
pvals                           = NaN( nfp, nVals );
for i                           = 1 : nVals
    val                         = Vals( i );
    fprintf( 1, '%s=%0.3g...\n', param_name, val )
    switch param_name
        case 'Gl'
            Gl                  = val;
        case 'El'
            El                  = val;
        case 'Vr'
            Vr                  = val;
        case 'Dn'
            Dn                  = val;
        case 'Ibias'
            Ibias               = val;
        case 'Amp'
            Amp                 = val;
        case 'Gc'
            Gc                  = val;
        case 'Gs'
            Gs                  = val;
        case 'Es'
            Es                  = val;
        case 'Vth'
            Vth                 = val;
    end
    switch model
        case 'LIF_C'
            [ coh1, fo, fbins, frate1, s, ~ ...
                , z1, zphs1, zfvals ]   = NR_sinusoids_to_cmodel( model, 'Tlong', 1 + Tlong ...
                , 'fvals', fvals, 'Dn', Dn, 'Ibias', Ibias, 'Amp', Amp ...
                , 'graphics', [ 0 0 0 ], 'nFFT', nFFT ...
                , 'prepad', 0, 'blackOutDuration', 1 ...
                , 'Gl_LIF_C', Gl, 'El_LIF_C', El, 'Vreset_LIF_C', Vr, 'Gc_LIF_C', Gc ...
                , 'dflag', dflag );
            
        case 'LIF_syn'
            [ coh1, fo, fbins, frate1, s, ~ ...
                , z1, zphs1, zfvals ]   = NR_sinusoids_to_cmodel( model, 'Tlong', Tlong ...
                , 'fvals', fvals, 'Dn', Dn, 'Ibias', Ibias, 'Amp', Amp ...
                , 'graphics', [ 0 0 0 ], 'nFFT', nFFT ...
                , 'prepad', 0, 'blackOutDuration', 0 ...
                , 'Gs_LIF_syn', Gs, 'Es_LIF_syn', Es, 'synmodel_LIF_syn', synmodel ...
                , 'SD_Out_LIF_syn', SD_Out ...
                , 'Gl_LIF_syn', Gl, 'El_LIF_syn', El, 'Vth_LIF_syn', Vth, 'Vreset_LIF_syn', Vr ...
                , 'dflag', dflag );
            
        otherwise
            [ coh1, fo, fbins, frate1, s, ~ ...
                , z1, zphs1, zfvals ]   = NR_sinusoids_to_cmodel( model, 'Tlong', Tlong ...
                , 'fvals', fvals, 'Dn', Dn, 'Ibias', Ibias, 'Amp', Amp ...
                , 'graphics', [ 0 0 0 ], 'nFFT', nFFT ...
                , 'Gl_LIF', Gl, 'El_LIF', El, 'Vreset_LIF', Vr ...
                , 'dflag', dflag );
            
    end
    cohs( :, i )                = coh1;
    frates( :, i )              = frate1;
    zmat( :, i )                = z1;
    zphs( :, i )                = zphs1;
    pvals( :, i )               = s.pvals( 1 : nfp, 3 );                    % [shuffling HF time-reversal]
end

%------------------------------------------------------------------------
% graphics
%------------------------------------------------------------------------
if ~graphics
    return
end

fig                             = figure;
colormap( myjet )
fidx                            = fo <= fvals( end );

subplot( 2, 2, 1 )
switch output
    case 'spiking'
        imagesc( fo( fidx ), 1 : nVals, cohs( fidx, : )' )
        ststr                   = 'Coherence';
    case 'subthreshold'
        imagesc( zfvals, 1 : nVals, zmat' )
        ststr                   = 'Impedance [M\Omega mm^2]';
end
axis xy
clims                           = get( gca, 'clim' );
set( gca, 'clim', [ 0 clims( 2 ) ] )
ch                              = colorbar;
set( ch, 'tickdir', 'out', 'box', 'off' )
title( ch, ststr )
set( ch, 'tickdir', 'out', 'box', 'off', 'FontName', 'Arial' )

xlabel( 'Frequency [Hz]' )
ylabel( ValStr )
set( gca, 'tickdir', 'out', 'box', 'off', 'FontName', 'Arial' )
ticks                           = get( gca, 'ytick' );
ticks                           = [ 1 ticks + 1 ];
if ticks( end ) > nVals
    ticks( end )                = [];
end
set( gca, 'ytick', ticks, 'YTickLabel', Vals( ticks ) );
[ ~, minidx ]                   = min( abs( Vals - Vcan ) );
line( xlim, [ 1 1 ] * minidx, 'color', [ 1 1 1 ], 'linestyle', '--' );

subplot( 2, 2, 2 )
imagesc( fbins, 1 : nVals, frates' )
axis xy
clims                           = get( gca, 'clim' );
set( gca, 'clim', [ 0 clims( 2 ) ] )
ch                              = colorbar;
set( ch, 'tickdir', 'out', 'box', 'off' )
title( ch, 'Firing rate [spks/s]' )
set( ch, 'tickdir', 'out', 'box', 'off', 'FontName', 'Arial' )
xlabel( 'Frequency [Hz]' )
ylabel( ValStr )
set( gca, 'tickdir', 'out', 'box', 'off', 'FontName', 'Arial' )
ticks                           = get( gca, 'ytick' );
ticks                           = [ 1 ticks + 1 ];
if ticks( end ) > nVals
    ticks( end )                = [];
end
set( gca, 'ytick', ticks, 'YTickLabel', Vals( ticks ) );
[ ~, minidx ]                   = min( abs( Vals - Vcan ) );
line( xlim, [ 1 1 ] * minidx, 'color', [ 1 1 1 ], 'linestyle', '--' );

local_fig_title( tstr );

return % NR_sinusoids_to_cmodel_run

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

% EOF
