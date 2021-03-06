% NR_make_Fig5      graphics for Figure 5 in Stark, Levi, Rotstein, 2022
%
% This script creates the sub-figures that compose Figure 5 in 
% Stark, Levi, Rotstein, 2022, PLOS Computational Biology

% Figure 5: Resonance generated at the level of postsynaptic potentials can be inherited to the network level
% (A-E) Synaptic and spiking resonance in a model of synaptic depression and facilitation 
% (H-I) Network resonance in a model of synaptic plasticity (facilitation + depression)
% (J-L) Network resonance in a model of synaptic plasticity
%
% Calls     NR_sinusoids_to_cmodel, NR_sinusoids_to_cmodel_run, NR_syntransmit
%           st_coherence, 
%           ParseArgPairs, plot_raster
%            
%
% 08-aug-21 ES

% last update
% 01-jul-22 AL

function NR_make_Fig5( savef, pstr, outdir, varargin )

% constants
prefix                          = 'resonance_network';
renderer_name                   = 'painters';
resize                          = '-bestfit';
color_output                    = [ 0 0 0 ];                                % output (Vm/spikes): black 
colors_layers_mean          	= [ 1 0 0; 0 0 0.7 ];                       % layer 1: "PYR"; layer 2: "INT"
colors_layers                   = [ 1 0.5 0.5; 0.7 0.7 1 ];

colors_synp                     = [ 1 0 0; 0 0 0; 0 0 1; 0 0.7 0 ];
synmodels_sth                   = { 'null', 'DF', 'F', 'D' };

% arguments
nargs                           = nargin;
if nargs < 1 || isempty( savef )
    savef                       = 1;
end
if nargs < 2 || isempty( pstr )
    pstr                        = '-dpdf';
end
if nargs < 3 || isempty( outdir )
    outdir                      = pwd;
end
datadir                         = outdir;

% arguments
[ Overwrite, graphics_each ...
    , graphics_all ]         	= ParseArgPairs(...
    { 'Overwrite', 'graphics_each' ...
    , 'graphics_all' }...
    , { -2, 0, 1 } ...
    , varargin{ : } );

% ----------------------------------------------------------------------
% (A-E) Synaptic and spiking resonance in a model of synaptic depression and facilitation 

% general parameters:
spike_length                = 1;
spike_width                 = 2;
spike_color                 = [ 1 1 1 ] * 0.5;
nmodels                     = length( synmodels_sth );

Gs_LIF_syn                  = 0.175;
Es_LIF_syn                  = 0;
Ibias_LIF_syn               = 1.3;

Tlong_sth                   = 20;
fvals_sth                   = ( 1 : 40 );
if length( fvals_sth ) == 2
    sig                     = 'chirp';
else
    sig                     = 'sines';
end

dflag                       = 0;

%---------------------------------------------------------
% subthreshold parameters:

Dn_syn                      = 0;
Vth_LIF_syn                 = 0;

%---------------------------------------------------------
% subthreshold - examples:

% what to show in the traces examples
fvals_slct                  = [ 1 10 40 100 ];                              % Inap - 4, 8, 16
nfvals                      = length( fvals_slct );
ncycles                     = 10;
synmodel                    = 'DF';

% 2D firing rate map (rate vs. frequency and phase)
fvals_all                   = [ fvals_sth( : ); fvals_slct( : ) ];
fvals_all                   = min( fvals_all ) : max( fvals_all );
[ ~, ~, ~, ~, s ]           = NR_sinusoids_to_cmodel( 'LIF_syn', 'Tlong', Tlong_sth, 'fvals', fvals_all ...
    , 'Dn', Dn_syn, 'Gs_LIF_syn', Gs_LIF_syn, 'Es_LIF_syn', Es_LIF_syn ...
    , 'nreps', 1, 'graphics', [ 0 0 0 0 ], 'synmodel_LIF_syn', synmodel ...
    , 'Ibias', Ibias_LIF_syn, 'Vth_LIF_syn', Vth_LIF_syn ...
    , 'dflag', dflag ...
    );

% Steady state EPSP
fig9( 1 )                   = figure;
subplot( 2, 2, 1 ) 
ph                          = plot( s.tfvals, s.mvals, '-b' );
cj                          = ismember( synmodels_sth, synmodel );
set( ph, 'color', colors_synp( cj, : ), 'linewidth', 2 );
set( gca, 'tickdir', 'out', 'box', 'off' )
xlabel( 'Presynaptic spike rate [spks/s]' )
ylabel( 'Steady-state EPSP [mV]' )

% plot a few selected examples
fig9( 2 )                   = figure;
ylims                       = NaN( nfvals, 2 );
xlimsec                     = [ -0.01 0.03 ];
for j                       = 1 : nfvals
    f                       = fvals_slct( j );
    [ ~, fidx ]             = min( abs( s.zfvals - f ) );
    
    subplot( nfvals, 2, 2 * ( j - 1 ) + 1 )
    tlim                    = [ 0 ncycles / f ];
    tidx                    = s.t >= tlim(1) & s.t <= tlim(2);
    t                       = s.t( tidx );
    y                       = s.Vmat( tidx, fidx );
    ph                      = plot( t, y, 'b' );
    set( ph, 'color', color_output );
    axis tight
    xlim( [ 0 max( xlim ) ] )
    ylims( j, : )           = ylim;
    title( sprintf( '%0.3g Hz', f ) )
    ylabel( 'V [mV]' )
    if j == nfvals
        xlabel( 'Time [s]' )
    end

    subplot( nfvals, 2, 2 * ( j - 1 ) + 2 )
    st                      = find( s.mat_in( :, fidx ) );
    t0                      = s.t( st( ncycles ) );
    tlim                    = t0 + xlimsec;
    tidx                    = s.t >= tlim(1) & s.t <= tlim(2);
    t                       = s.t( tidx ) - t0;
    y                       = s.Vmat( tidx, fidx );
    ph                      = plot( t * 1000, y, 'b' );
    set( ph, 'color', color_output );
    axis tight
    title( sprintf( '%0.3g Hz', f ) )
    ylabel( 'V [mV]' )
    if j == nfvals
        xlabel( 'Time [ms]' )
    end
    
end
ylims                       = [ min( ylims( :, 1 ) ) max( ylims( :, 2 ) ) ] + [ -0.1 0.5 ];
for j                       = 1 : nfvals
    subplot( nfvals, 2, 2 * ( j - 1 ) + 1 )
    f                       = fvals_slct( j );
    [ ~, fidx ]             = min( abs( s.zfvals - f ) );
    tlim                    = [ 0 ncycles / f ];
    tidx                    = s.t >= tlim(1) & s.t <= tlim(2);
    st                      = s.mat_in( tidx, fidx );
    plot_raster( { s.t( find( st ) ) }, [], ylims( 1 ) - 0.5, spike_length, spike_width, spike_color );
    set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims + [ -1 0 ] )
    
    subplot( nfvals, 2, 2 * ( j - 1 ) + 2 )
    st                      = find( s.mat_in( :, fidx ) );
    t0                      = s.t( st( ncycles ) );
    tlim                    = t0 + xlimsec;
    tidx                    = s.t >= tlim(1) & s.t <= tlim(2);
    
    st                      = s.mat_in( tidx, fidx );
    str                     = ( s.t( find( st ) ) + xlimsec( 1 ) ) * 1000;      % [ms]
    plot_raster( { str }, [], ylims( 1 ) - 0.5, spike_length, spike_width, spike_color );
    set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims + [ -1 0 ] )
    xlim( xlimsec * 1000 )
    
end

%---------------------------------------------------------
% subthreshold - 4 models:
clear s
for j                       = 1 : nmodels
    [ ~, ~, ~, ~, sj ]      = NR_sinusoids_to_cmodel( 'LIF_syn', 'Tlong', Tlong_sth, 'fvals', fvals_sth ...
        , 'Dn', Dn_syn, 'Gs_LIF_syn', Gs_LIF_syn, 'Es_LIF_syn', Es_LIF_syn ...
        , 'nreps', 1, 'graphics', [ 0 0 0 0 ], 'synmodel_LIF_syn', synmodels_sth{ j } ...
        , 'Ibias', Ibias_LIF_syn, 'Vth_LIF_syn', Vth_LIF_syn ...
        , 'dflag', dflag ...
        );
    s( j ) = sj;
end
fvals                       = s( 1 ).tfvals;
mvals                       = NaN( length( fvals ), nmodels );
for j                       = 1 : nmodels
    mvals( :, j )           = s( j ).mvals;
end

% plot
fig9( 3 )                	= figure;

subplot( 2, 2, 1 )
ph                        	= plot( fvals, mvals );
xlabel( 'Input rate [spikes/s]' )
ylabel( 'Steady-state EPSP [mV]' )
for j                     	= 1 : nmodels
    synmodel              	= synmodels_sth{ j };
    cj                     	= ismember( synmodels_sth, synmodel );
    set( ph( j ), 'color', colors_synp( cj, : ), 'linewidth', 2 );
end
legend( ph, synmodels_sth, 'Location', 'best' );
set( gca, 'tickdir', 'out', 'box', 'off', 'FontSize', 12 );
limx(1)                     = min (fvals );
limx(2)                     = max (fvals );
xlim([limx(1) limx(2)])
title( sprintf( '%s; T=%0.2g s, D=%0.2g; Gs=%0.2g; Es=%0.2g', sig, Tlong_sth, Dn_syn, Gs_LIF_syn, Es_LIF_syn ) )

subplot( 2, 2, 2 )
ph                       	= plot( fvals, local_scale( mvals,'linear' ) );
xlabel( 'Input rate [spikes/s]' )
ylabel( 'Steady-state EPSP (scaled)' )
for j                     	= 1 : nmodels
    synmodel             	= synmodels_sth{ j };
    cj                      = ismember( synmodels_sth, synmodel );
    set( ph( j ), 'color', colors_synp( cj, : ), 'linewidth', 2 );
end    
legend( ph, synmodels_sth, 'Location', 'best' );
set( gca, 'tickdir', 'out', 'box', 'off', 'FontSize', 12 );
limx(1)                     = min (fvals );
limx(2)                     = max (fvals );
xlim([limx(1) limx(2)])
title( sprintf( '%s; T=%0.2g s, D=%0.2g; Gs=%0.2g; Es=%0.2g', sig, Tlong_sth, Dn_syn, Gs_LIF_syn, Es_LIF_syn ) )

%---------------------------------------------------------
% spiking parameters:

Dn_syn                      = 0.05;
Vth_LIF_syn                 = -50;

%---------------------------------------------------------
% spiking example:

% what to show in the traces examples
fvals_slct                  = [ 1 10 20 40 ];                        	% Inap - 4, 8, 16
nfvals                      = length( fvals_slct );
ncycles                     = 10;
synmodel                    = 'DF';

% 2D firing rate map (rate vs. frequency and phase)
fvals_all                   = [ fvals_sth( : ); fvals_slct( : ) ];
fvals_all                   = min( fvals_all ) : max( fvals_all );
[ ~, ~, ~, ~, s ]      = NR_sinusoids_to_cmodel( 'LIF_syn', 'Tlong', Tlong_sth, 'fvals', fvals_all ...
    , 'Dn', Dn_syn, 'Gs_LIF_syn', Gs_LIF_syn, 'Es_LIF_syn', Es_LIF_syn ...
    , 'nreps', 1, 'graphics', [ 0 1 0 0 ], 'synmodel_LIF_syn', synmodel ...
    , 'Ibias', Ibias_LIF_syn, 'Vth_LIF_syn', Vth_LIF_syn ...
    , 'dflag', dflag ...
    );

fig9( 4 )                   = gcf;

% firing rate and coherence
fig9( 5 )                   = figure;
subplot( 2, 2, 1 ) % firing rate
ph                          = plot( s.fbins, s.frate_vec, '.-b' );
set( ph, 'color', color_output )
ylims                       = ylim;
ylims                       = [ 0 max( ylims( : ) ) ];
set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims )
ylabel( 'Firing rate [spks/s]' )

subplot( 2, 2, 3 ) % coherence
ph                          = plot( s.fo( s.fidx ), s.cohs( s.fidx ), '.-k' );
set( ph, 'color', color_output )
ylims                       = ylim;
ylims                       = [ 0 max( ylims( : ) ) ];
set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims )
ylabel( 'Coherence' )
xlabel( 'Frequency [Hz]' )

ylims                       = NaN( nfvals, 2 );
for j                       = 1 : nfvals
    f                       = fvals_slct( j );
    [ ~, fidx ]             = min( abs( s.zfvals - f ) );
    
    subplot( nfvals, 2, 2 * ( j - 1 ) + 2 )
    tlim                    = [ 0 ncycles / f ];
    tidx                    = s.t >= tlim(1) & s.t <= tlim(2);
    t                       = s.t( tidx );
    y                       = s.Vmat( tidx, fidx );
    ph                      = plot( t, y, 'b' );
    set( ph, 'color', color_output );
    axis tight
    xlim( [ 0 max( xlim ) ] )
    ylims( j, : )           = ylim;
    title( sprintf( '%0.3g Hz', f ) )
    ylabel( 'V [mV]' )
    if j == nfvals
        xlabel( 'Time [s]' )
    end

end

ylims                       = [ min( ylims( :, 1 ) ) max( ylims( :, 2 ) ) ] + [ -0.1 0.5 ];

for j                       = 1 : nfvals
    subplot( nfvals, 2, 2 * ( j - 1 ) + 2 )
    f                       = fvals_slct( j );
    [ ~, fidx ]             = min( abs( s.zfvals - f ) );
    tlim                    = [ 0 ncycles / f ];
    tidx                    = s.t >= tlim(1) & s.t <= tlim(2);
    st                      = s.mat_in( tidx, fidx );
    plot_raster( { s.t( find( st ) ) }, [], ylims( 1 ) + diff( ylims ) / 4, diff( ylims ) / 2, spike_width, spike_color );
    set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims + [ -1 0 ] )
end

%---------------------------------------------------------
% spiking - 4 models:
clear s
for j                       = 1 : nmodels
    [ ~, ~, ~, ~, sj ]      = NR_sinusoids_to_cmodel( 'LIF_syn', 'Tlong', Tlong_sth, 'fvals', fvals_sth ...
        , 'Dn', Dn_syn, 'Gs_LIF_syn', Gs_LIF_syn, 'Es_LIF_syn', Es_LIF_syn ...
        , 'nreps', 1, 'graphics', [ 0 0 0 0 ], 'synmodel_LIF_syn', synmodels_sth{ j } ...
        , 'Ibias', Ibias_LIF_syn, 'Vth_LIF_syn', Vth_LIF_syn ...
        , 'dflag', dflag ...
        );
    s( j ) = sj;
end
fidx                        = s( 1 ).fidx;
fvals                       = s( 1 ).fo;
cohs                        = NaN( length( fvals ), nmodels );
phs                         = NaN( length( fvals ), nmodels );
for j                       = 1 : nmodels
    cohs( :, j )            = s( j ).cohs;
    phs( :, j )             = s( j ).phs;
end

% plot
fig9( 6 )                	= figure;

subplot( 2, 2, 1 )
ph                        	= plot( fvals, cohs );
xlabel( 'Input rate [spikes/s]' )
ylabel( 'Coherence' )
for j                     	= 1 : nmodels
    synmodel              	= synmodels_sth{ j };
    cj                     	= ismember( synmodels_sth, synmodel );
    set( ph( j ), 'color', colors_synp( cj, : ), 'linewidth', 2 );
end
legend( ph, synmodels_sth, 'Location', 'best' );
set( gca, 'tickdir', 'out', 'box', 'off', 'FontSize', 12 );
limx(1)                     = min (fvals(fidx) );
limx(2)                     = max (fvals(fidx) );
xlim([limx(1) limx(2)])
title( sprintf( '%s; T=%0.2g s, D=%0.2g; Gs=%0.2g; Es=%0.2g', sig, Tlong_sth, Dn_syn, Gs_LIF_syn, Es_LIF_syn ) )

subplot( 2, 2, 2 )
ph                       	= plot( fvals, local_scale( cohs,'linear' ) );
xlabel( 'Input rate [spikes/s]' )
ylabel( 'Coherence (scaled)' )
for j                     	= 1 : nmodels
    synmodel             	= synmodels_sth{ j };
    cj                      = ismember( synmodels_sth, synmodel );
    set( ph( j ), 'color', colors_synp( cj, : ), 'linewidth', 2 );
end    
legend( ph, synmodels_sth, 'Location', 'best' );
set( gca, 'tickdir', 'out', 'box', 'off', 'FontSize', 12 );
limx(1)                     = min (fvals(fidx) );
limx(2)                     = max (fvals(fidx) );
xlim([limx(1) limx(2)])
title( sprintf( '%s; T=%0.2g s, D=%0.2g; Gs=%0.2g; Es=%0.2g', sig, Tlong_sth, Dn_syn, Gs_LIF_syn, Es_LIF_syn ) )

subplot( 2, 2, 3 )
uwphs                       = unwrap( phs );
ph                        	= plot( fvals, uwphs );
limx(1)                     = min (fvals(fidx) );
limx(2)                     = max (fvals(fidx) );
xlim([limx(1) limx(2)])
xlabel( 'Input rate [spikes/s]' )
ylabel( 'Phase [rad]' )
for j                     	= 1 : nmodels
    synmodel              	= synmodels_sth{ j };
    cj                     	= ismember( synmodels_sth, synmodel );
    set( ph( j ), 'color', colors_synp( cj, : ), 'linewidth', 2 );
end
local_alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
ylim( [ -pi pi ] )
legend( ph, synmodels_sth, 'Location', 'best' );
set( gca, 'tickdir', 'out', 'box', 'off', 'FontSize', 12 );
title( sprintf( '%s; T=%0.2g s, D=%0.2g; Gs=%0.2g; Es=%0.2g', sig, Tlong_sth, Dn_syn, Gs_LIF_syn, Es_LIF_syn ) )

%-----
fig                         = fig9( fig9 ~= 0 );
if savef
    for i                   = 1 : length( fig )
        figname                 = [ outdir filesep prefix '_FIG5A-E_part' num2str( i ) ];
        figure( fig( i ) );
        figi              	= gcf;
        figi.Renderer      	= renderer_name;
        pause( 0.2 )
        if isequal( pstr, '-dpdf' )
            rsz          	= resize;
        else
            rsz          	= '';
        end
        print( figi, pstr, figname, rsz )
    end
end
    
% ----------------------------------------------------------------------
% (H-I) Network resonance in a model of synaptic plasticity
% ----------------------------------------------------------------------
% (1) set up the parameters

% input frequencies
Tlong_sth                       = 20;
fvals_sth                       = [ 0 40 ];

% first-to-second layer:
Es_LIF_syn                      = 0;
Vth_LIF_syn                     = -50;
synmodel                        = 'DF';

ntrials                         = 50;

Gs1                             = 0.2;
Ibias_spk1                      = 1.2;
Dn_spk1                         = 0.25;
Gs2                             = 0.12;

% second layer
Ibias_spk2                      = 0;
Dn_spk2                         = 0;

SD_Out_LIF_syn                  = 0.1;                                      % output of second layer, input to third layer
                                                                       
% analysis parameters
dflag                           = 0;
Fs                              = 10000;
cFs                             = 1250;
fROI                            = [ 0 40 ];
M                               = 1;
mtNW                            = 3;
nFFT                            = NaN;

fig5H                           = zeros( 3, 1 );

%---------------------------------------------------------
% (2) run the simulation and analyze the coherence + phase
% initialize output
st1                             = cell( 1, ntrials );
nfc                             = 513;
cohs1                           = NaN( nfc, ntrials );
phs1                            = NaN( nfc, ntrials );

% simulate first and second layers
fprintf( 1, 'Simulating second layer: ' )
for i                           = 1 : ntrials
    
    fprintf( 1, '%d ', i )
    % simulate a single train
    [ ~, ~, ~, ~, s ]           = NR_sinusoids_to_cmodel( 'LIF_syn', 'Tlong', Tlong_sth, 'fvals', fvals_sth ...
        , 'Dn', Dn_spk1, 'Gs_LIF_syn', Gs1, 'Es_LIF_syn', Es_LIF_syn ...
        , 'nreps', 1, 'graphics', [ 0 0 0 0 ], 'verbose', 0, 'synmodel_LIF_syn', synmodel ...
        , 'Ibias', Ibias_spk1, 'Vth_LIF_syn', Vth_LIF_syn ...
        , 'SD_Out_LIF_syn', SD_Out_LIF_syn ...
        );

    % analog stimulus
    Ievec                       = s.Iemat( : );

    % input spike train
    stvec_in                    = find( s.mat_in( : ) );
    
    % output spike train
    stvec                       = find( s.mat( : ) );
    st1{ i }                    = stvec;
    
    % coherence
    [ cohs, phs, fo ]           = st_coherence( stvec, stvec_in ...
        , 'spkFs', Fs, 'xFs', Fs, 'xSpk', 1, 'ntmax', length( Ievec ) ...
        , 'fROI', fROI, 'Fs', cFs, 'graphics', 0, 'M', M ...
        , 'mtNW', mtNW, 'nFFT', nFFT, 'dflag', dflag );
    cohs1( :, i )               = cohs;
    phs1( :, i )                = phs;
    
end
fprintf( 1, 'done!\n' )
st0                             = stvec_in;

% simulate third layer
fprintf( 1, 'Simulating third layer: ' )
st2                             = NR_syntransmit( st1, 'Gs', Gs2, 'Ibias', Ibias_spk2, 'Dn', Dn_spk2, 'SD_In', SD_Out_LIF_syn );
fig5H( 1 )                      = gcf;
[ cohs2, phs2 ]                 = st_coherence( st2, st0, 'spkFs', Fs, 'xFs', Fs, 'xSpk', 1, 'ntmax', length( Ievec ) ...
    , 'fROI', fROI, 'Fs', cFs, 'graphics', 0, 'M', M, 'mtNW', mtNW, 'nFFT', nFFT, 'dflag', dflag );

%---------------------------------------------------------
% (3) plot the results
fig5H( 2 )                      = figure;
subplot( 2, 1, 1 )
plot_raster( st1 ); 
title( sprintf( 'n=%d, Dn=%0.2g mV, Ibias=%0.2g \\muA/cm^2; Gs=%0.2g \\muS/cm^2', ntrials, Dn_spk1, Ibias_spk1, Gs2 ) )
ylims                           = ylim;
plot_raster( { st0 }, [], 10, 20, [], [ 1 0.7 0.7 ] );
ylim( ylims )
plot_raster( st1 );

subplot( 2, 1, 2 )
plot_raster( { st0 }, [], 2, 4, [], [ 1 0.7 0.7 ] );
plot_raster( { st2 }, [], 2, 1, [], [ 0 0 1 ] );
ylim( [ 0 3 ] )

fig5H( 3 )                      = figure;
fidx                            = fo >= fROI(1) & fo <= fROI(2);       %inrange( fo, fROI );
subplot( 2, 1, 1 )
ph1                             = plot( fo( fidx ), cohs1( fidx, : ), 'b' );
set( ph1, 'color', colors_layers( 1, : ) )
title( sprintf( 'n=%d, Dn=%0.2g mV, Ibias=%0.2g \\muA/cm^2; Gs=%0.2g \\muS/cm^2', ntrials, Dn_spk1, Ibias_spk1, Gs2 ) )
hold on
ph1m                            = plot( fo( fidx ), mean( cohs1( fidx, : ), 2 ), 'b' );
set( ph1m, 'color', colors_layers_mean( 1, : ), 'linewidth', 2 )
ph2                             = plot( fo( fidx ), cohs2( fidx, : ), 'r' );
set( ph2, 'color', colors_layers_mean( 2, : ), 'linewidth', 2 )
set( gca, 'tickdir', 'out', 'box', 'off' );
xlabel( 'Frequency [Hz]' )
ylabel( 'Coherence' )
lh                              = legend( [ ph1( 1 ) ph2 ], { 'Pre', 'Post' }, 'Location', 'northeast' );
set( lh, 'box', 'off' )

subplot( 2, 1, 2 )
ph1                             = plot( fo( fidx ), phs1( fidx, : ), 'b' );
set( ph1, 'color', colors_layers( 1, : ) )
hold on
uwphs1                          = unwrap( local_circ_mean( phs1' )' );
uwphs1                          = mod( uwphs1 + pi, 2 * pi ) - pi;
uwphs2                          = unwrap( phs2 );
uwphs2                          = mod( uwphs2 + pi, 2 * pi ) - pi;
ph1m                            = plot( fo( fidx ), uwphs1( fidx ), 'b' );
set( ph1m, 'color', colors_layers_mean( 1, : ), 'linewidth', 2 )
ph2                             = plot( fo( fidx ), uwphs2( fidx ), 'r' );
set( ph2, 'color', colors_layers_mean( 2, : ), 'linewidth', 2 )
set( gca, 'tickdir', 'out', 'box', 'off' );
ylim( [ -1 1 ] * pi/2 )
set( gca, 'ytick', ( -1 : 0.25 : 1 ) * pi / 2 )
xlabel( 'Frequency [Hz]' )
ylabel( 'Phase [rad]' )
local_alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
lh                              = legend( [ ph1( 1 ) ph2 ], { 'Pre', 'Post' }, 'Location', 'southeast' );
set( lh, 'box', 'off' )

%-----
fig                             = fig5H( fig5H ~= 0 );
if savef
    for i                       = 1 : length( fig )
        figname                 = [ outdir filesep prefix '_FIG5H-I_part' num2str( i ) ];
        figure( fig( i ) );
        figi                    = gcf;
        figi.Renderer           = renderer_name;
        pause( 0.2 )
        if isequal( pstr, '-dpdf' )
            rsz                 = resize;
        else
            rsz                 = '';
        end
        print( figi, pstr, figname, rsz )
    end
end

%---------------------------------------------------------
% (J-L) Network resonance in a model of synaptic plasticity
%---------------------------------------------------------

% (1) set up the general parameters

% input frequencies
Tlong_sth                       = 20;
fvals_sth                       = [ 0 40 ];

% first-to-second layer:
Es_LIF_syn                      = 0;
Vth_LIF_syn                     = -50;
synmodel                        = 'DF';

ntrials                         = 50;                                   % number of layer2 neurons
Gs1                             = 0.2;
ref_var                         = 'DnL';
Ibias1                          = 1.2;                                  % [muA/cm^2]

Dn_spk1_list                    = 0 : 0.025 : 2;
Ibias_spk1_list                 = Ibias1 * ones( 1, length( Dn_spk1_list ) );
ref_val                         = Dn_spk1_list;
ref_val_label                   = 'Dn [mV]';
str                             = sprintf( '%s_Ibias_%0.2g', ref_var, Ibias1 );
ref_var_str                     = sprintf( '%s#%s', str( 1 : ( strfind( str, '.' ) - 1 ) ) ...
    , str( ( strfind( str, '.' ) + 1 ) : length( str ) ) );

%---------------------------------------------------------
% (2) compute
nruns                           = length( Dn_spk1_list );
fig5J                           = zeros( 3, nruns );

fname                           = sprintf( '%s/NR_fig5_%s.mat', datadir, ref_var_str );
if Overwrite < 0 && exist( fname, 'file' )
    
    load( fname, '-mat', 'fo', 'fidx', 'cohs1_all', 'phs1_all', 'nst1_all' ...
        , 'cohs2_all', 'phs2_all', 'nst2_all', 'uwphs1_all', 'uwphs2_all' ...
        , 'ref_val', 'ref_val_label', 'Tlong_sth' )
    
else
    % compute   
    Gs2                         = 0.12;
    
    % second layer
    Ibias_spk2                  = 0;                                        % [muA/cm^2]
    Dn_spk2                     = 0;                                        % [mV]
    
    SD_Out_LIF_syn              = 0.1;                                      % output of second layer, input to third layer

    % analysis parameters
    dflag                   	= 0;
    Fs                          = 10000;
    cFs                       	= 1250;
    fROI                        = [ 0 40 ];
    M                           = 1;
    mtNW                        = 3;
    nFFT                        = NaN;
    nfc                         = 513;
    
    % set up the output fields
    cohs1_all                   = NaN( nfc, ntrials, nruns );
    phs1_all                    = NaN( nfc, ntrials, nruns );
    nst1_all                    = NaN( 1, ntrials, nruns );
    cohs2_all                   = NaN( nfc, 1, nruns );
    phs2_all                    = NaN( nfc, 1, nruns );
    nst2_all                    = NaN( 1, 1, nruns );
    uwphs1_all                  = NaN( nfc, 1, nruns );
    uwphs2_all                  = NaN( nfc, 1, nruns );
    
    % run a specific simulation and analyze the coherence + phase
    for rn                      = 1 : nruns
        
        % determine the specific parameters
        Ibias_spk1              = Ibias_spk1_list( rn );
        Dn_spk1                 = Dn_spk1_list( rn );
        
        % initialize output
        st1                     = cell( 1, ntrials );
        cohs1                   = NaN( nfc, ntrials );
        phs1                    = NaN( nfc, ntrials );
        nst1                   	= NaN( 1, ntrials );
        
        % simulate first and second layers
        fprintf( 1, '%d) Simulating second layer (Dn: %0.3g; Ibias: %0.3g):\t', rn, Dn_spk1, Ibias_spk1 )
        for i                   = 1 : ntrials
            
            fprintf( 1, '%d ', i )
            % simulate a single train
            [ ~, ~, ~, ~, s ]   = NR_sinusoids_to_cmodel( 'LIF_syn' ...
                , 'Tlong', Tlong_sth, 'fvals', fvals_sth ...
                , 'Dn', Dn_spk1, 'Gs_LIF_syn', Gs1, 'Es_LIF_syn', Es_LIF_syn ...
                , 'nreps', 1, 'graphics', [ 0 0 0 0 ], 'verbose', 0, 'synmodel_LIF_syn', synmodel ...
                , 'Ibias', Ibias_spk1, 'Vth_LIF_syn', Vth_LIF_syn ...
                , 'SD_Out_LIF_syn', SD_Out_LIF_syn ...
                );
            
            % analog stimulus
            Ievec               = s.Iemat( : );
            
            % input spike train
            stvec_in            = find( s.mat_in( : ) );
            
            % output spike train
            stvec               = find( s.mat( : ) );
            st1{ i }            = stvec;
            nst1( i )           = length( stvec );
            
            % coherence
            [ cohs, phs, fo ]   = st_coherence( stvec, stvec_in ...
                , 'spkFs', Fs, 'xFs', Fs, 'xSpk', 1, 'ntmax', length( Ievec ) ...
                , 'fROI', fROI, 'Fs', cFs, 'graphics', 0, 'M', M ...
                , 'mtNW', mtNW, 'nFFT', nFFT, 'dflag', dflag );
            cohs1( :, i )       = cohs;
            phs1( :, i )        = phs;
            
        end
        fprintf( 1, 'done!\n' )
        st0                     = stvec_in;
        
        % simulate third layer
        fprintf( 1, '%d) Simulating third layer (Dn: %0.3g; Ibias: %0.3g):\t' ...
            , rn, Dn_spk2, Ibias_spk2 )
        st2                     = NR_syntransmit( st1, 'Gs', Gs2 ...
            , 'Ibias', Ibias_spk2, 'Dn', Dn_spk2, 'SD_In', SD_Out_LIF_syn ...
            , 'graphics', graphics_each, 'verbose', 0 );
        nst2                    = length( st2 );
        if graphics_each
            fig5J( 1, rn )      = gcf;
        end
        [ cohs2, phs2 ]         = st_coherence( st2, st0, 'spkFs', Fs ...
            , 'xFs', Fs, 'xSpk', 1, 'ntmax', length( Ievec ) ...
            , 'fROI', fROI, 'Fs', cFs, 'graphics', 0, 'M', M ...
            , 'mtNW', mtNW, 'nFFT', nFFT, 'dflag', dflag );
        fprintf( 1, 'done!\n' )
        
        % unwrap the phases
        uwphs1                  = unwrap( local_circ_mean( phs1' )' );
        uwphs1                  = mod( uwphs1 + pi, 2 * pi ) - pi;
        uwphs2                  = unwrap( phs2 );
        uwphs2                  = mod( uwphs2 + pi, 2 * pi ) - pi;
        
        %---------------------------------------------------------
        % summarize the results from the specific run
        fidx                    = fo >= fROI( 1 ) & fo <= fROI( 2 );
        cohs1_all( :, :, rn )   = cohs1;
        phs1_all( :, :, rn )    = phs1;
        nst1_all( :, :, rn )    = nst1;
        cohs2_all( :, :, rn )   = cohs2;
        phs2_all( :, :, rn )    = phs2;
        nst2_all( :, :, rn )    = nst2;
        
        uwphs1_all( :, :, rn )  = uwphs1;
        uwphs2_all( :, :, rn )  = uwphs2;
        
        if ~graphics_each
            continue
        end
        
        %---------------------------------------------------------
        % plot the single-run results
        fig5J( 2, rn )          = figure;
        subplot( 2, 1, 1 )
        plot_raster( st1 );
        tstr                    = sprintf( 'n=%d, Dn=%0.2g mV, Ibias=%0.2g \\muA/cm^2; Gs=%0.2g \\muS/cm^2' ...
            , ntrials, Dn_spk1, Ibias_spk1, Gs2 );
        title( tstr )
        ylims                	= ylim;
        plot_raster( { st0 }, [], 10, 20, [], [ 1 0.7 0.7 ] );
        ylim( ylims )
        plot_raster( st1 );
        
        subplot( 2, 1, 2 )
        plot_raster( { st0 }, [], 2, 4, [], [ 1 0.7 0.7 ] );
        plot_raster( { st2 }, [], 2, 1, [], [ 0 0 1 ] );
        ylim( [ 0 3 ] )
        
        fig5J( 3, rn )          = figure;
        subplot( 2, 1, 1 );
        ph1                     = plot( fo( fidx ), cohs1( fidx, : ), 'b' );
        set( ph1, 'color', colors_layers( 1, : ) )
        tstr                    = sprintf( 'n=%d, Dn=%0.2g mV, Ibias=%0.2g \\muA/cm^2; Gs=%0.2g \\muS/cm^2' ...
            , ntrials, Dn_spk1, Ibias_spk1, Gs2 );
        title( tstr )
        hold on
        ph1m                    = plot( fo( fidx ), mean( cohs1( fidx, : ), 2 ), 'b' );
        set( ph1m, 'color', colors_layers_mean( 1, : ), 'linewidth', 2 )
        ph2                     = plot( fo( fidx ), cohs2( fidx, : ), 'r' );
        set( ph2, 'color', colors_layers_mean( 2, : ), 'linewidth', 2 )
        set( gca, 'tickdir', 'out', 'box', 'off' );
        xlabel( 'Frequency [Hz]' )
        ylabel( 'Coherence' )
        lh                      = legend( [ ph1( 1 ) ph2 ] ...
            , { 'Pre', 'Post' }, 'Location', 'northeast' );
        set( lh, 'box', 'off' )
        
        subplot( 2, 1, 2 );
        ph1                  	= plot( fo( fidx ), phs1( fidx, : ), 'b' );
        set( ph1, 'color', colors_layers( 1, : ) )
        hold on
        ph1m                    = plot( fo( fidx ), uwphs1( fidx ), 'b' );
        set( ph1m, 'color', colors_layers_mean( 1, : ), 'linewidth', 2 )
        ph2                     = plot( fo( fidx ), uwphs2( fidx ), 'r' );
        set( ph2, 'color', colors_layers_mean( 2, : ), 'linewidth', 2 )
        set( gca, 'tickdir', 'out', 'box', 'off' );
        ylim( [ -1 1 ] * pi/2 )
        set( gca, 'ytick', ( -1 : 0.25 : 1 ) * pi / 2 )
        xlabel( 'Frequency [Hz]' )
        ylabel( 'Phase [rad]' )
        local_alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
        lh                      = legend( [ ph1( 1 ) ph2 ] ...
            , { 'Pre', 'Post' }, 'Location', 'southeast' );
        set( lh, 'box', 'off' )
        
    end % rn
    
    fprintf( 1, 'done all!\n' )
    
end

%---------------------------------------------------------
% (3) save the results
if Overwrite >= 1 || ( Overwrite ~= -1 && ~exist( fname, 'file' ) )
    save( fname, 'fo', 'fidx', 'cohs1_all', 'phs1_all', 'nst1_all' ...
        , 'cohs2_all', 'phs2_all', 'nst2_all', 'uwphs1_all', 'uwphs2_all' ...
        , 'ref_val', 'ref_val_label', 'Tlong_sth', '-v6' )
end

%---------------------------------------------------------
% (4) plot global summary:

% mean and SEM
cohs1_list                      = permute( mean( cohs1_all( fidx, :, : ), 2 ), [ 1 3 2 ] );
cohs2_list                      = permute( cohs2_all( fidx, :, : ), [ 1 3 2 ] );
nst1_list                       = permute( mean( nst1_all, 2 ), [ 3 1 2 ] );
nst2_list                       = permute( nst2_all, [ 3 1 2 ] );
uphs1_list                      = permute( uwphs1_all( fidx, :, : ), [ 1 3 2 ] );
uphs2_list                      = permute( uwphs2_all( fidx, :, : ), [ 1 3 2 ] );
nst1_sem_list                   = permute( local_calc_sem( nst1_all, 2 ), [ 3 1 2 ] ); 

% maximal coherence and f_res
fvals                           = fo( fidx );
fres1                           = NaN( nruns, ntrials );
fcoh1                           = NaN( nruns, ntrials );
fres2                           = NaN( nruns, 1 );
fcoh2                           = NaN( nruns, 1 );
for rn                          = 1 : nruns
    [ maxval, maxidx ]          = max( cohs1_all( fidx, :, rn ) );
    vidx                        = ~isnan( maxval );
    fres1( rn, vidx )           = fvals( maxidx( ~isnan( maxval ) ) );
    fcoh1( rn, vidx )           = maxval( ~isnan( maxval ) );
    
    [ maxval, maxidx ]          = max( cohs2_all( fidx, :, rn ) );
    vidx                        = ~isnan( maxval );
    fres2( rn, vidx )           = fvals( maxidx( ~isnan( maxval ) ) );
    fcoh2( rn, vidx )           = maxval( ~isnan( maxval ) );
end

if graphics_each
    fig5J                       = fig5J( : );
else
    fig5J                       = [];
end

if graphics_all
    
    fig                         = figure;
    colormap( myjet )
    
    % coherence
    subplot( 2, 3, 1 )
    imagesc( fo( fidx ), ref_val, cohs1_list' )
    axis xy
    tstr                        = sprintf( 'Layer 2 [%0.2g %0.2g]' ...
        , min( cohs1_list( : ) ), max( cohs1_list( : ) ) );
    th                          = title( tstr );
    set( th, 'color', colors_layers_mean( 1, : ) );
    set( gca, 'tickdir', 'out', 'box', 'off' );
    xlabel( 'Presynaptic spike rate [spks/s]' )
    ylabel( ref_val_label )
    
    subplot( 2, 3, 2 )
    imagesc( fo( fidx ), ref_val, cohs2_list' ), axis xy
    tstr                        = sprintf( 'Layer 3 [%0.2g %0.2g]' ...
        , min( cohs2_list( : ) ), max( cohs2_list( : ) ) );
    th                          = title( tstr );
    set( th, 'color', colors_layers_mean( 2, : ) );
    set( gca, 'tickdir', 'out', 'box', 'off' );
    xlabel( 'Presynaptic spike rate [spks/s]' )
    ylabel( ref_val_label )
    
    % firing rate
    subplot( 2, 3, 4 )
    hold on
    
    mu                          = nst1_list / Tlong_sth;
    se                          = nst1_sem_list / Tlong_sth;
    local_patch_band( ref_val, mu, se, colors_layers_mean( 1, : ), colors_layers( 1, : ) );
    ph2                         = plot( ref_val, nst2_list / Tlong_sth, '-' );
    set( ph2, 'color', colors_layers_mean( 2, : ) )
    xlabel( ref_val_label )
    ylabel( 'FR [spk/s]' )
    ylim( [ 0 max( ylim ) ] )
    
    subplot( 2, 3, 5 )
    hold on
    ph1                         = plot( fo( fidx ), uphs1_list, 'b' );
    ph2                         = plot( fo( fidx ), uphs2_list, 'r' );
    set( ph1, 'color', colors_layers_mean( 1, : ) )
    set( ph2, 'color', colors_layers_mean( 2, : ) )
    ylim( [ -1 1 ] * pi/2 )
    set( gca, 'ytick', ( -1 : 0.25 : 1 ) * pi / 2 )
    xlabel( 'Presynaptic spike rate [spks/s]' )
    ylabel( 'Phase [rad]' )
    local_alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    
    % add the f_res and the coherence peak as a function of ref_val
    f_res                       = [ nanmean( fres1, 2 ) nanmean( fres2, 2 ) ];
    f_res_sem                   = [ local_calc_sem( fres1, 2 ) local_calc_sem( fres2, 2 ) ];
    f_coh                       = [ nanmean( fcoh1, 2 ) nanmean( fcoh2, 2 ) ];
    f_coh_sem                   = [ local_calc_sem( fcoh1, 2 ) local_calc_sem( fcoh2, 2 ) ];

    subplot( 2, 3, 3 )
    vidx                        = ~isnan( f_res );
    local_patch_band( ref_val( vidx( :, 1 ) ), f_res( vidx( :, 1 ), 1 ) ...
        , f_res_sem( vidx( :, 1 ), 1 ), colors_layers_mean( 1, : ), colors_layers( 1, : ) );
    hold on
    local_patch_band( ref_val( vidx( :, 2 ) ), f_res( vidx( :, 2 ), 2 ) ...
        , f_res_sem( vidx( :, 2 ), 2 ), colors_layers_mean( 2, : ), colors_layers( 2, : ) );
    xlabel( 'Noise level [mV]' )
    ylabel( 'Resonant frequency [Hz]' )
    ylim( [ 0 max( ylim ) ] )
    
    subplot( 2, 3, 6 )
    vidx                        = ~isnan( f_coh );
    local_patch_band( ref_val( vidx( :, 1 ) ), f_coh( vidx( :, 1 ), 1 ) ...
        , f_coh_sem( vidx( :, 1 ), 1 ), colors_layers_mean( 1, : ), colors_layers( 1, : ) );
    hold on
    local_patch_band( ref_val( vidx( :, 2 ) ), f_coh( vidx( :, 2 ), 2 ) ...
        , f_coh_sem( vidx( :, 2 ), 2 ), colors_layers_mean( 2, : ), colors_layers( 2, : ) );
    xlabel( 'Noise level [mV]' )
    ylabel( 'Coherence' )
    ylim( [ 0 max( ylim ) ] )
    
    [ ~, midx ]                 = max( f_coh( :, 1 ) ); 
    title( sprintf( 'L2: Dn=%0.2g mV', ref_val( midx ) ) )
    local_alines( ref_val( midx ), 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );    

    subplot( 2, 3, 3 )
    title( sprintf( 'L2: Dn=%0.2g mV, f_{SR}=%0.2g Hz', ref_val( midx ), f_res( midx, 1 ) ) )
    local_alines( ref_val( midx ), 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );    

    subplot( 2, 3, 4 )
    local_alines( ref_val( midx ), 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );    

    subplot( 2, 3, 1 )
    local_alines( f_res( midx, 1 ), 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );    
    local_alines( ref_val( midx ), 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );    
    
    for i                       = 1 : 6
        subplot( 2, 3, i )
        set( gca, 'tickdir', 'out', 'box', 'off', 'FontName', 'Arial' )
        axis square
    end
    
end

fig5J                           = [ fig5J; fig ];

%---------------------------------------------------------
% (5) save global summary plot:
fig                             = fig5J( fig5J ~= 0 );
if savef
    for i                       = 1 : length( fig )
        figname                 = [ outdir filesep prefix '_FIG5J-L_part' num2str( i ) ];
        figure( fig( i ) );
        figi                    = gcf;
        figi.Renderer           = renderer_name;
       
        pause( 0.2 )
        if isequal( pstr, '-dpdf' )
            rsz                 = resize;
        else
            rsz                 = '';
        end
        print( figi, pstr, figname, rsz )
    end
end

return % NR_make_Fig5

%------------------------------------------------------------------------
function y = local_calc_sem( x, dim )

if isempty( x )
    y                           = NaN;
    return
end

if all( isnan( x( : ) ) )
    if nargin < 2 || isempty( dim )
        dim                     = 1;
    end
    y                           = nanmean( x, dim );
    return 
end 
if ~exist( 'dim', 'var' ) || isempty( dim )
    if any( size( x ) == 1 )
        x                       = x( : );
    end
    dim                         = 1;
end

y                               = nanstd( x, [], dim ) ./ sqrt( sum( ~isnan( x ), dim ) );

return % local_calc_sem

%------------------------------------------------------------------------
function z = local_scale( mat, scaling )

bias                        = NaN;
if ismember( scaling, { 'log', 'db' } ) && isnan( bias )
    if sum( mat( : ) > 0 )
        mat( mat <= 0 )     = min( mat( mat > 0 ) );
    end
    bias                    = 0;
end
switch scaling
    case 'log'
        z                   = log10( abs( mat ) + bias );
    case 'db'
        z                   = 20 * log10( abs( mat ) + bias );
    case 'linear'
        z                   = mat;        
end

return % local_scale

%------------------------------------------------------------------------
function ph = local_patch_band( x, y, s, c1, c2, sy, nanflag )

ph                          = [];

nargs = nargin;
if nargs < 2 || ~isequal( numel( x ), numel( y ) ) || isempty( x )
    return
else
    x                       = x( : );
    y                       = y( : );
end
if isempty( s )
    s                       = zeros( size( x ) );
else
    s                       = s( : );
end
if nargs < 4 || isempty( c1 ) || length( c1 ) ~= 3
    c1                      = [ 0 0 1 ];
else
    c1                      = c1( : )';
end
if nargs < 5 || isempty( c2 ) || length( c2 ) ~= 3
    c2                      = c1;
    c2( c2 == 0 )           = 0.3;
else
    c2                      = c2( : )';
end
if nargs < 6 || isempty( sy )
    sy                      = zeros( size( y ) );
else
    sy                      = sy( : );
end
if nargs < 7 || isempty( nanflag )
    nanflag                 = 0;
end

if nanflag
    nans                    = isnan( x ) | isnan( y );
    x( nans )               = [];
    y( nans )               = [];
    s( nans )               = [];
    sy( nans )              = [];
end
if isempty( x )
    return
end

h0 = ishold;
if ~h0
    hold on
end
ph( 2 )                     = patch( [ x + sy; flipud( x - sy ) ], [ y + s; flipud( y - s ) ], c2 );
set( ph( 2 ), 'edgecolor', c2 )
ph( 1 )                     = line( x, y, 'color', c1, 'linewidth', 2 );
if ~h0
    hold off
end

return % local_patch_band

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
function lh = local_alines( x, ax, varargin )

if nargin < 2 || isempty( ax ), ax = 'x'; end
ax                              = lower( ax( 1 ) );

x                               = x( : );
nx                              = length( x );
ah                              = gca;
if ax == 'x'
    X                           = [ x x NaN * ones( nx, 1 ) ]';
    X                           = X( : );
    Y                           = repmat( [ get( ah, 'ylim' ) NaN ]', nx, 1 );
elseif ax == 'y'
    Y                           = [ x x NaN * ones( nx, 1) ]';
    Y                           = Y( : );
    X                           = repmat( [ get( ah, 'xlim' ) NaN ]', nx, 1 );
else
    return
end
lh                              = line( X, Y, varargin{:} );

return % local_alines

%EOF