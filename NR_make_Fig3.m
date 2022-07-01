% NR_make_Fig3      graphics for Figure 3 in Stark, Levi, Rotstein, 2022
%
% This script creates the sub-figures that compose Figure 3 in 
% Stark, Levi, Rotstein, 2022, PLOS Computational Biology
%
% Figure 3: Resonance can be generated directly at the spiking level
% (A-E) Spiking resonance in the LIF model
% (F-H) Spiking resonance without subthreshold resonance in a model with voltage-dependent Calcium dynamics
%
% Calls      NR_sinusoids_to_cmodel,NR_sinusoids_to_cmodel_run.
%
% 08-aug-21 ES

% last update
% 01-jul-22 AL

function NR_make_Fig3( savef, pstr, outdir)

% constants
prefix                          = 'resonance_network';
renderer_name                   = 'painters';
resize                          = '-bestfit';
color_input                     = [ 0 0 0.7 ];                              % input (current): blue
color_output                    = [ 0 0 0 ];                                % output (Vm/spikes): black 
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
% ----------------------------------------------------------------------
%  (A-E) Spiking resonance in the LIF model
% ----------------------------------------------------------------------
    
% how to run the simulations
model                       = 'LIF';
Dn                          = 0;
Ibias                       = 0.9;
Amp                         = 0.115;                                    % supraTH
fvals                       = 1 : 40;
Tlong                       = 3;                                        % [s]

% what to show in the traces examples
fvals_slct                  = [ 4 12 16 ];                              % Inap - 4, 8, 16
nfvals                      = length( fvals_slct );
tlim                        = [ 0 1 ];

% analysis parameters
dflag                       = '';                                       % with 0, not as clear

% 2D firing rate map (rate vs. frequency and phase)
[ ~, ~, ~, ~, s ]           = NR_sinusoids_to_cmodel( model, 'Tlong', Tlong, 'fvals', fvals ...
    , 'Dn', Dn, 'nreps', 1, 'Amp', Amp, 'Ibias', Ibias, 'graphics', [ 0 1 0 0 ], 'nFFT', 1250, 'dflag', dflag );
fig6( 1 )                   = gcf;

% firing rate and coherence 
fig6( 2 )                   = figure;
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

% add a few selected examples (~identical to code in fignums(3), different scaling of input)
tidx                        = s.t >= tlim(1) & s.t <= tlim(2);       %inrange( s.t, tlim );
t                           = s.t( tidx );
ylims                       = NaN( nfvals, 2 );
for j                       = 1 : nfvals
    f                       = fvals_slct( j );
    [ ~, fidx ]             = min( abs( s.zfvals - f ) );
    y                       = s.Vmat( tidx, fidx );
    subplot( nfvals, 2, 2 * ( j - 1 ) + 2 )
    ph                      = plot( t * 1000, y, 'b' );
    set( ph, 'color', color_output );
    axis tight
    xlim( [ 0 max( xlim ) ] )
    ylims( j, : )           = ylim;
    title( sprintf( '%0.3g \\muA/cm^2; %0.3g Hz', Amp, f ) )
    ylabel( 'V [mV]' )
    if j == nfvals
        xlabel( 'Time [ms]' )
    end
end
ylims                       = [ min( ylims( :, 1 ) ) max( ylims( :, 2 ) ) ] + [ -1 1 ] * 0.1;
for j                       = 1 : nfvals
    subplot( nfvals, 2, 2 * ( j - 1 ) + 2 )
    set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims + [ -1 0 ] )
    f                       = fvals_slct( j );
    [ ~, fidx ]             = min( abs( s.zfvals - f ) );
    x                       = s.Iemat( tidx, fidx );
    hold on
    ph                      = plot( t * 1000, local_scale( x ,'linear') * 10 + ylims( 1 ) + 1, 'b' );
    set( ph, 'color', color_input );
end

%--------------------------------------------------------------------
fig6( 8 )                   = figure;

% data:
fin                         = s.fo( s.fidx );
Z                           = s.cohs( s.fidx );
Phi                         = s.phs( s.fidx );
pvals                       = s.pvals( :, 2 );
pTH                         = 0.05;
vidx                        = pvals <= pTH;
Phi( isnan( Z ) )           = NaN;

subplot( 3, 3, 1 ) % firing rate
ph                          = plot( s.fbins, s.frate_vec, '.-b' );
set( ph, 'color', color_output )
ylims                       = ylim;
ylims                       = [ 0 max( ylims( : ) ) ];
set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims )
ylabel( 'Firing rate [spks/s]' )

subplot( 3, 3, 4 ) % coherence
ph                          = plot( fin, Z, '.-k' );
set( ph, 'color', color_output )
hold on
ph                          = plot( fin( vidx ), Z( vidx ), '.-r' ); 
set( ph, 'color', [ 1 0 0 ] )
ylims                       = ylim;
ylims                       = [ 0 max( ylims( : ) ) ];
set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims )
ylabel( 'Coherence' )
%axis square
set( gca, 'xtick', 0 : 10 : 40 )
[ ~, midx ]                 = max( Z ); 
fres                        = fin( midx );
alines( fres, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
title( sprintf( '%s, Amp=%0.2g; f_{res} = %0.2g Hz', model, Amp, fres ) )

subplot( 3, 3, 7 ) % phase
ylims                       = [ -1 1 ] * pi;
dy                          = ylims( 1 ) : diff( ylims ) / 8 : ylims( 2 );
uwp                         = unwrap( Phi );
if max( uwp( vidx ) ) > 2 * pi
    uwp                     = uwp - 2 * pi;
end
ph                          = plot( fin, uwp, '-' );
set( ph, 'color', [ 0 0 0 ], 'linewidth', 0.5 );
hold on
ph                          = plot( fin( vidx ), uwp( vidx ), '.-' ); 
set( ph, 'color', [ 1 0 0 ], 'linewidth', 2 );
xlim( [ 0 40 ] )
ylim( ylims * 1.1 )
line( xlim, [ 0 0 ], 'color', [ 0 0 0 ], 'linestyle', '--' )
line( xlim, [ 1 1 ] * ylims( 1 ), 'color', [ 0 0 0 ], 'linestyle', '--' )
line( xlim, [ 1 1 ] * ylims( 2 ), 'color', [ 0 0 0 ], 'linestyle', '--' )
set( gca, 'tickdir', 'out', 'box', 'off' )
set( gca, 'ytick', dy, 'yticklabel', round( dy * 100 ) / 100 )
xlabel( 'Frequency [Hz]' )
ylabel( 'Phase [rad]' )
set( gca, 'xtick', 0 : 10 : 40 )
mat                         = [ uwp( 1 : end - 1 ) < 0 uwp( 2 : end ) > 0 vidx( 1 : end - 1 ) ];
fphs                        = fin( all( mat, 2 ) ) + diff( fin( 1 : 2 ) ) / 2;
alines( fphs, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
title( sprintf( 'f_{phs} = %0.2g Hz', fphs ) )

% add a few selected examples (~identical to code in fignums(3), different scaling of input)
tidx                        = s.t >= tlim(1) & s.t <= tlim(2);      %inrange( s.t, tlim );
t                           = s.t( tidx );
ylims                       = NaN( nfvals, 2 );
for j                       = 1 : nfvals
    f                       = fvals_slct( j );
    [ ~, fidx ]             = min( abs( s.zfvals - f ) );
    y                       = s.Vmat( tidx, fidx );
    subplot( nfvals, 2, 2 * ( j - 1 ) + 2 )
    ph                      = plot( t * 1000, y, 'b' );
    set( ph, 'color', color_output );
    axis tight
    xlim( [ 0 max( xlim ) ] )
    ylims( j, : )           = ylim;
    title( sprintf( '%0.3g \\muA/cm^2; %0.3g Hz', Amp, f ) )
    ylabel( 'V [mV]' )
    if j == nfvals
        xlabel( 'Time [ms]' )
    end
end
ylims                       = [ min( ylims( :, 1 ) ) max( ylims( :, 2 ) ) ] + [ -1 1 ] * 0.1;
for j                       = 1 : nfvals
    subplot( nfvals, 2, 2 * ( j - 1 ) + 2 )
    set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims + [ -1 0 ] )
    f                       = fvals_slct( j );
    [ ~, fidx ]             = min( abs( s.zfvals - f ) );
    x                       = s.Iemat( tidx, fidx );
    hold on
    ph                      = plot( t * 1000, local_scale( x ,'linear') * 10 + ylims( 1 ) + 1, 'b' );
    set( ph, 'color', color_input );
end

%-----
fig                         = fig6;
if savef
    for i                   = 1 : length( fig )
        figname = [ outdir filesep prefix '_FIG3A-E_part' num2str( i ) ];
        try
            figure( fig( i ) );
        catch
            continue
        end
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

%---------------------------------------------------------
% 2D plots with Noise, Amp, and Ibias levels
sig                         = 'sines';
[ ~, ~, ~, figs ]           = NR_sinusoids_to_cmodel_run( model, sig, 'Dn', 'nFFT', 1250, 'dflag', dflag );
fig6( 3 )                   = figs;
[ ~, ~, ~, figs ]           = NR_sinusoids_to_cmodel_run( model, sig, 'Amp', 'nFFT', 1250, 'dflag', dflag );
fig6( 4 )                   = figs;
[ ~, ~, ~, figs ]           = NR_sinusoids_to_cmodel_run( model, sig, 'Ibias', 'nFFT', 1250, 'dflag', dflag );
fig6( 5 )                   = figs;

%---------------------------------------------------------
% subthreshold response
Amp                         = 0.05;
[ ~, ~, ~, ~, s ]           = NR_sinusoids_to_cmodel( model, 'Tlong', Tlong, 'fvals', fvals ...
    , 'Dn', Dn, 'nreps', 1, 'Amp', Amp, 'Ibias', Ibias, 'graphics', 0, 'nFFT', 1250 );

fin                         = s.zfvals;
Z                           = s.z;
Phi                         = s.zphs;

fig6( 6 )                   = figure;

subplot( 2, 2, 1 )
ph                          = plot( fin, Z );                           % MOhm * mm^2
set( ph, 'color', [ 0 0 0 ], 'linewidth', 2 );
ylim( [ 0 1.05 ] )
xlim( [ 0 40 ] )
set( gca, 'tickdir', 'out', 'box', 'off' )
ylabel( 'Impedance [M\Omega*mm^2]' )
axis square
set( gca, 'xtick', 0 : 10 : 40 )

subplot( 2, 2, 3 )
Phi( isnan( Z ) )           = NaN;
uwp                         = unwrap( Phi );
if max( uwp ) > 2 * pi
    uwp                     = uwp - 2 * pi;
end
ph                          = plot( fin, uwp );
set( ph, 'color', [ 0 0 0 ], 'linewidth', 2 );
xlim( [ 0 40 ] )
ylim( [ -1 1 ] * pi/2 * 1.1 )
line( xlim, [ 0 0 ], 'color', [ 0 0 0 ], 'linestyle', '--' )
line( xlim, [ 1 1 ] * pi/2, 'color', [ 0 0 0 ], 'linestyle', '--' )
line( xlim, [ 1 1 ] * -pi/2, 'color', [ 0 0 0 ], 'linestyle', '--' )
set( gca, 'tickdir', 'out', 'box', 'off' )
set( gca, 'ytick', ( -1 : 0.25 : 1 ) * pi / 2, 'yticklabel', round( ( -1 : 0.25 : 1 ) * pi / 2 * 100 ) / 100 )
xlabel( 'Frequency [Hz]' )
ylabel( 'Phase [rad]' )
axis square
set( gca, 'xtick', 0 : 10 : 40 )

% add a few selected examples (~identical to code in fignums(3), different scaling of input)
tidx                        = s.t >= tlim(1) & s.t <= tlim(2);       %inrange( s.t, tlim );
t                           = s.t( tidx );
ylims                       = NaN( nfvals, 2 );
for j                       = 1 : nfvals
    f                       = fvals_slct( j );
    [ ~, fidx ]             = min( abs( s.zfvals - f ) );
    y                       = s.Vmat( tidx, fidx );
    subplot( nfvals, 2, 2 * ( j - 1 ) + 2 )
    ph                      = plot( t * 1000, y, 'b' );
    set( ph, 'color', color_output );
    axis tight
    xlim( [ 0 max( xlim ) ] )
    ylims( j, : )           = ylim;
    title( sprintf( '%0.3g \\muA/cm^2; %0.3g Hz', Amp, f ) )
    ylabel( 'V [mV]' )
    if j == nfvals
        xlabel( 'Time [ms]' )
    end
end
ylims                       = [ min( ylims( :, 1 ) ) max( ylims( :, 2 ) ) ] + [ -1 1 ] * 0.1;
for j                       = 1 : nfvals
    subplot( nfvals, 2, 2 * ( j - 1 ) + 2 )
    set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims + [ -1 0 ] )
    f                       = fvals_slct( j );
    [ ~, fidx ]             = min( abs( s.zfvals - f ) );
    x                       = s.Iemat( tidx, fidx );
    hold on
    ph                      = plot( t * 1000, local_scale( x,'linear' ) / 2 + ylims( 1 ) - 0.9, 'b' );
    set( ph, 'color', color_input );
end



%---------------------------------------------------------
% 2D plots with extended Amp range
sig                         = 'sines';
amp_vals                    = 0 : 0.0025 : 0.3; % extended range
[ ~, ~, ~ ...
    , afig, ~ ]             = NR_sinusoids_to_cmodel_run( model, sig, 'Amp','nFFT', 1250 ...
    , 'param_vals', amp_vals, 'dflag', dflag );
fig6( 7 )                   = afig;

%-----
fig                         = fig6;
if savef
    for i                   = 1 : length( fig )
        figname             = [ outdir filesep prefix '_FIG3A-E_part' num2str( i ) ];
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
% (F-H) Spiking resonance without subthreshold resonance in a model with voltage-dependent Calcium dynamics
% ----------------------------------------------------------------------

% how to run the simulations
model                       = 'LIF_C';
Dn                          = 0.001;                                    % [mV]
Ibias                       = -3;
Amp                         = 8;
Gc_spk                      = 0.08;                                     % [mS/cm^2]
Gc_spk0                     = 0;                                        % [mS/cm^2]
fvals                       = 1 : 40;
prepad                      = 0;                                        % [s]
Tlong                       = 4;                                        % [s]
blackOutDuration            = 1;                                        % [s]
Gl_LIF_C                    = 0.5;

% analysis
dflag                       = '';

% what to show in the traces examples
fvals_slct                  = [ 1 8 16 ];
nfvals                      = length( fvals_slct );
tlim                        = [ 1 2 ];

%---------------------------------------------------------------------
% 2D firing rate map (rate vs. frequency and phase)
[ ~, ~, ~, ~, s0 ]       	= NR_sinusoids_to_cmodel( model, 'Tlong', Tlong, 'fvals', fvals ...
    , 'Dn', Dn, 'nreps', 1, 'Amp', Amp, 'Ibias', Ibias, 'graphics', [ 0 0 0 0 ], 'nFFT', 1250 ...
    , 'Gc_LIF_C', Gc_spk0, 'Gl_LIF_C', Gl_LIF_C, 'dflag', dflag ...
    , 'prepad', prepad, 'blackOutDuration', blackOutDuration );

[ ~, ~, ~, ~, s ]           = NR_sinusoids_to_cmodel( model, 'Tlong', Tlong, 'fvals', fvals ...
    , 'Dn', Dn, 'nreps', 1, 'Amp', Amp, 'Ibias', Ibias, 'graphics', [ 0 1 1 0 ], 'nFFT', 1250 ...
    , 'Gc_LIF_C', Gc_spk, 'Gl_LIF_C', Gl_LIF_C, 'dflag', dflag ...
    , 'prepad', prepad, 'blackOutDuration', blackOutDuration );
fig8( 1 )                   = gcf;

% firing rate and coherence 
color_nocal                 = colors_synp( ismember( synmodels_sth, 'null' ), : );

fig8( 2 )                   = figure;
subplot( 2, 2, 1 ) % firing rate
ph                          = plot( s.fbins, s.frate_vec, '.-b', s0.fbins, s0.frate_vec, '.-r' );
set( ph( 1 ), 'color', color_output )
set( ph( 2 ), 'color', color_nocal )
ylims                       = ylim;
ylims                       = [ 0 max( ylims( : ) ) ];
set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims )
ylabel( 'Firing rate [spks/s]' )

subplot( 2, 2, 3 ) % coherence
ph                          = plot( s.fo( s.fidx ), s.cohs( s.fidx ), '.-k', s0.fo( s.fidx ), s0.cohs( s.fidx ), '.-r' );
set( ph( 1 ), 'color', color_output )
set( ph( 2 ), 'color', color_nocal )
ylims                       = ylim;
ylims                       = [ 0 max( ylims( : ) ) ];
set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims )
ylabel( 'Coherence' )
xlabel( 'Frequency [Hz]' )

% add a few selected examples (~identical to code in fignums(3), different scaling of input)
tidx                        = s.t >= tlim(1) & s.t <= tlim(2);        %inrange( s.t, tlim );
t                           = s.t( tidx ) - tlim( 1 );
ylims                       = NaN( nfvals, 2 );
for j                       = 1 : nfvals
    f                       = fvals_slct( j );
    [ ~, fidx ]             = min( abs( s.zfvals - f ) );
    y                       = s.Vmat( tidx, fidx );
    subplot( nfvals, 2, 2 * ( j - 1 ) + 2 )
    ph                      = plot( t * 1000, y, 'b' );
    set( ph, 'color', color_output, 'linewidth', 2 );
    axis tight
    xlim( [ 0 max( xlim ) ] )
    ylims( j, : )           = ylim;
    title( sprintf( '%0.3g \\muA/cm^2; %0.3g Hz', Amp, f ) )

    hold on
    y0                      = s0.Vmat( tidx, fidx );
    ph                      = plot( t * 1000, y0, '-r' );
    set( ph, 'color', color_nocal );
    
    ylabel( 'V [mV]' )
    if j == nfvals
        xlabel( 'Time [ms]' )
    end
end
ylims                       = [ min( ylims( :, 1 ) ) max( ylims( :, 2 ) ) ] + [ -1 1 ] * 0.1;
for j                       = 1 : nfvals
    subplot( nfvals, 2, 2 * ( j - 1 ) + 2 )
    set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims + [ -1 0 ] )
    f                       = fvals_slct( j );
    [ ~, fidx ]             = min( abs( s.zfvals - f ) );
    x                       = s.Iemat( tidx, fidx );
    hold on
    ph                      = plot( t * 1000, local_scale( x ,'linear') * 10 + ylims( 1 ) + 1, 'b' );
    set( ph, 'color', color_input );
end

%---------------------------------------------------------------------
% 2D plots with Gc levels
Gl_LIF_C                    = 0.5; 
El_LIF_C                    = -60; 
Vreset_LIF_C                = -70;
sig                         = 'sines';
[ ~, ~, ~, figs ]           = NR_sinusoids_to_cmodel_run( model, sig, 'Gc' ...
    , 'Gl', Gl_LIF_C, 'El', El_LIF_C, 'Vr', Vreset_LIF_C, 'Gc', 0.08, 'Dn', Dn ...
    , 'nFFT', 1250 );
fig8( 3 )                   = figs;

%---------------------------------------------------------------------
% subthreshold (impedance) profiles
Gl_LIF_C                    = 0.5;
[ ~, ~, ~, ~, s0 ]       	= NR_sinusoids_to_cmodel( model, 'Tlong', Tlong, 'fvals', fvals ...
    , 'Dn', Dn, 'nreps', 1, 'Amp', Amp, 'Ibias', Ibias, 'graphics', [ 0 0 0 0 ], 'nFFT', 1250 ...
    , 'Gc_LIF_C', Gc_spk0, 'Vth_LIF_C', 0, 'Gl_LIF_C', Gl_LIF_C ...
    , 'prepad', prepad, 'blackOutDuration', blackOutDuration );


[ ~, ~, ~, ~, s1 ]        	= NR_sinusoids_to_cmodel( model, 'Tlong', Tlong, 'fvals', fvals ...
    , 'Dn', Dn, 'nreps', 1, 'Amp', Amp, 'Ibias', Ibias, 'graphics', [ 0 0 0 0 ], 'nFFT', 1250 ...
    , 'Gc_LIF_C', Gc_spk, 'Vth_LIF_C', 0, 'Gl_LIF_C', Gl_LIF_C ...
    , 'prepad', prepad, 'blackOutDuration', blackOutDuration );

fig8( 4 )                   = figure;

for i = 1 : 2
    
    switch i
        case 1
            s               = s1;
            acolor          = color_output;
            astyle          = '-';
        case 2
            s               = s0;
            acolor          = color_nocal;
            astyle          = '--';
    end
    fin                     = s.zfvals;
    Z                       = s.z;
    Phi                     = s.zphs;
    
    subplot( 2, 2, 1 )
    hold on
    ph                      = plot( fin, Z );                           % MOhm * mm^2
    set( ph, 'color', acolor, 'linewidth', 2, 'linestyle', astyle );
    xlim( [ 0 40 ] )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Impedance [M\Omega*mm^2]' )
    axis square
    set( gca, 'xtick', 0 : 10 : 40 )
    
    subplot( 2, 2, 3 )
    hold on
    Phi( isnan( Z ) )       = NaN;
    uwp                     = unwrap( Phi );
    if max( uwp ) > 2 * pi
        uwp                 = uwp - 2 * pi;
    end
    ph = plot( fin, uwp );
    set( ph, 'color', acolor, 'linewidth', 2, 'linestyle', astyle );
    xlim( [ 0 40 ] )
    ylim( [ -1 1 ] * pi/2 * 1.1 )
    line( xlim, [ 0 0 ], 'color', [ 0 0 0 ], 'linestyle', '--' )
    line( xlim, [ 1 1 ] * pi/2, 'color', [ 0 0 0 ], 'linestyle', '--' )
    line( xlim, [ 1 1 ] * -pi/2, 'color', [ 0 0 0 ], 'linestyle', '--' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    set( gca, 'ytick', ( -1 : 0.25 : 1 ) * pi / 2, 'yticklabel', round( ( -1 : 0.25 : 1 ) * pi / 2 * 100 ) / 100 )
    xlabel( 'Frequency [Hz]' )
    ylabel( 'Phase [rad]' )
    axis square
    set( gca, 'xtick', 0 : 10 : 40 )
    
    % selected examples
    tidx                    = s.t >= tlim(1) & s.t <= tlim(2);              %inrange( s.t, tlim );
    t                       = s.t( tidx ) - tlim( 1 );
    ylims                   = NaN( nfvals, 2 );
    for j                   = 1 : nfvals
        f                   = fvals_slct( j );
        [ ~, fidx ]         = min( abs( s.zfvals - f ) );
        y                   = s.Vmat( tidx, fidx );
        subplot( nfvals, 2, 2 * ( j - 1 ) + 2 )
        hold on
        ph                  = plot( t * 1000, y, 'b' );
        set( ph, 'color', acolor, 'linestyle', astyle );
        
        axis tight
        xlim( [ 0 max( xlim ) ] )
        ylims( j, : )       = ylim;
        title( sprintf( '%0.3g \\muA/cm^2; %0.3g Hz', Amp, f ) )
        ylabel( 'V [mV]' )
        if j == nfvals
            xlabel( 'Time [ms]' )
        end
    end
    ylims                   = [ min( ylims( :, 1 ) ) max( ylims( :, 2 ) ) ] + [ -1 1 ] * 0.1;
    for j                   = 1 : nfvals
        subplot( nfvals, 2, 2 * ( j - 1 ) + 2 )
        set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims + [ -1 0 ] )
        f                   = fvals_slct( j );
        [ ~, fidx ]         = min( abs( s.zfvals - f ) );
        x                   = s.Iemat( tidx, fidx );
        ph                  = plot( t * 1000, local_scale( x,'linear' ) / 2 + ylims( 1 ) - 0.75, 'b' );
        set( ph, 'color', color_input );
    end
    
end 

%-----
fig                         = fig8;
if savef
    for i                   = 1 : length( fig )
        figname = [ outdir filesep prefix '_FIG3F-H_part' num2str( i ) ];
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

return

%------------------------------------------------------------------------
function z = local_scale( mat, scaling )
bias                        = NaN;          % illegal hack for preventing visualizing 0's as NaNs
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
return

%EoF