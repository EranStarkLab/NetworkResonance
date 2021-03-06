% NR_make_Fig2      graphics for Figure 2 in Stark, Levi, Rotstein, 2022
%
% This script creates the sub-figures that compose Figure 2 in 
% Stark, Levi, Rotstein, 2022, PLOS Computational Biology
%
% Figure 2: Resonance generated at the level of membrane potential fluctuations can be inherited to the network level
% (A)   subTH resonance in Inap/h model                       
% (B-E) Subthreshold to spiking resonance in the Inap/h model            
% (F-H) Subthreshold to network resonance in the Inap/h model     
%
% Calls     NR_sinusoids_to_cmodel, NR_sinusoids_to_cmodel_run, plotTraces, st_fingerprint
%
% 08-aug-21 ES

% last update
% 01-jul-22 AL

function NR_make_Fig2( savef, pstr, outdir )

% constants
prefix                          = 'resonance_network';
renderer_name                   = 'painters';
resize                          = '-bestfit';
color_input                     = [ 0 0 0.7 ];                              % input (current): blue
color_output                    = [ 0 0 0 ];                                % output (Vm/spikes): black 
colors_EI                       = [ 106 27 154; 46 125 50 ] / 255;          % PYR, INT
colors_EI_input                 = [ 156 77 204; 96 173 94 ] / 255;          % PYR, INT

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
% (A) resonance in Inap/h model
% ----------------------------------------------------------------------

% how to run the simulations
model                           = 'Inap';
Dn                              = 0;
Ibias                           = -1.85;
Amp                             = 0.05;                                     % subTH
fvals                           = [ 0.1 0.5 : 0.5 : 40 ];                   % [Hz]
Tlong_DFLT                      = 3;                                        % [s]
Tlong_f                         = 3 * ceil( 1 / min( fvals ) );             % at least 3 full cycles
Tlong                           = max( Tlong_DFLT, Tlong_f );

% what to show in the traces examples
fvals_slct                      = [ 4 8 16 ];
nfvals                          = length( fvals_slct );
tlim                            = [ 1 2.002 ];

fig2A                           = zeros( 4, 1 );

for i                           = 1 : 4
    
    switch i 
        case 1 % LPF
            Gp_Inap             = 0;
            Gh_Inap             = 0;
            C                   = 1;
        case 2 % HPF
            Gp_Inap             = 0;
            Gh_Inap             = 1;
            C                   = 0.1;
        case 3 % BPF
            Gp_Inap             = 0;
            Gh_Inap             = 1;
            C                   = 1;
        case 4 % sharper BPF
            Gp_Inap             = 0.1;
            Gh_Inap             = 1;
            C                   = 1;
    end
    
    [ ~, ~, ~, ~, s ]           = NR_sinusoids_to_cmodel( model, 'Tlong', Tlong, 'fvals', fvals ...
        , 'Dn', Dn, 'nreps', 1, 'Amp', Amp, 'Ibias', Ibias, 'graphics', 0 ...
        , 'nFFT', 1250 ...
        , 'Gp_Inap', Gp_Inap,'Gh_Inap', Gh_Inap, 'C', C );      
    
    fin                         = s.zfvals; 
    Z                           = s.z; 
    Phi                         = s.zphs;

    fig2A( i )                  = figure;
    
    subplot( 2, 2, 1 )
    ph                          = plot( fin, Z );                           % MOhm * mm^2
    set( ph, 'color', [ 0 0 0 ], 'linewidth', 2 );
    if i == 4
        ylim( [ 0 3.05 ] )
    else
        ylim( [ 0 1.05 ] )
    end
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
    ph = plot( fin, uwp );
    set( ph, 'color', [ 0 0 0 ], 'linewidth', 2 );
    xlim( [ 0 40 ] )
    ylim( [ -1 1 ] * pi/2 * 1.1 )
    line( xlim, [ 0 0 ], 'color', [ 0 0 0 ], 'linestyle', '--' )
    line( xlim, [ 1 1 ] * pi/2, 'color', [ 0 0 0 ], 'linestyle', '--' )
    line( xlim, [ 1 1 ] * -pi/2, 'color', [ 0 0 0 ], 'linestyle', '--' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    set( gca, 'ytick', ( -1 : 0.25 : 1 ) * pi / 2 ...
        , 'yticklabel', round( ( -1 : 0.25 : 1 ) * pi / 2 * 100 ) / 100 )
    xlabel( 'Frequency [Hz]' )
    ylabel( 'Phase [rad]' )
    axis square
    set( gca, 'xtick', 0 : 10 : 40 )
    
    % add a few selected examples
    tidx                        = s.t >= tlim(1) & s.t <= tlim(2); 
    t                           = s.t( tidx ) - tlim( 1 );
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
        ylims( j, : )       = ylim;
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
        ph                      = plot( t * 1000, local_scale( x ,'linear') / 2 + ylims( 1 ) - 0.75, 'b' );
        set( ph, 'color', color_input );
    end
    
end

%-----
fig                             = fig2A;
if savef
    for i                       = 1 : length( fig )
        figname                 = [ outdir filesep prefix '_FIG2A_part' num2str( i ) ];
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


% ----------------------------------------------------------------------
% (B-E): Subthreshold to spiking resonance in the Inap/h model (including examples)

% show examples with several discrete sines (2, 6, 12 Hz) at one Amp level
% show firing rate, coherence, phase for the same data
% then show 2D plot with Amp, and 2D plot with noise level (coherence only for both)
% ----------------------------------------------------------------------

% how to run the simulations
model                           = 'Inap';
Dn                              = 0;
Ibias                           = -1.85;
Amp                             = 0.15;                                     % subTH
fvals                           = 1 : 40;
Tlong                           = 3;                                        % [s]

% what to show in the traces examples
fvals_slct                      = [ 4 8 16 ];
nfvals                          = length( fvals_slct );
tlim                            = [ 0 1 ];
fig2B                           = zeros( 10, 1 );

% 2D firing rate map (rate vs. frequency and phase)
[ ~, ~, ~, ~, s ]               = NR_sinusoids_to_cmodel( model, 'Tlong', Tlong, 'fvals', fvals ...
    , 'Dn', Dn, 'nreps', 1, 'Amp', Amp, 'Ibias', Ibias, 'graphics', [ 0 1 0 0 ], 'nFFT', 1250 );
fig2B( 1 )                      = gcf;

% firing rate and coherence 
fig2B( 2 )                      = figure;

subplot( 2, 2, 1 ) % firing rate
ph                              = plot( s.fbins, s.frate_vec, '.-b' );
set( ph, 'color', color_output )
ylims                           = ylim;
ylims                           = [ 0 max( ylims( : ) ) ];
set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims )
ylabel( 'Firing rate [spks/s]' )

subplot( 2, 2, 3 ) % coherence
ph                              = plot( s.fo( s.fidx ), s.cohs( s.fidx ), '.-k' );
set( ph, 'color', color_output )
ylims                           = ylim;
ylims                           = [ 0 max( ylims( : ) ) ];
set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims )
ylabel( 'Coherence' )
xlabel( 'Frequency [Hz]' )

% add a few selected examples (~identical to code in fignums(3), different scaling of input)
tidx                            = s.t >= tlim( 1 ) & s.t <= tlim( 2 );
t                               = s.t( tidx );
ylims                           = NaN( nfvals, 2 );
for j                           = 1 : nfvals
    f                           = fvals_slct( j );
    [ ~, fidx ]                 = min( abs( s.zfvals - f ) );
    y                           = s.Vmat( tidx, fidx );
    subplot( nfvals, 2, 2 * ( j - 1 ) + 2 )
    ph                          = plot( t * 1000, y, 'b' );
    set( ph, 'color', color_output );
    axis tight
    xlim( [ 0 max( xlim ) ] )
    ylims( j, : )               = ylim;
    title( sprintf( '%0.3g \\muA/cm^2; %0.3g Hz', Amp, f ) )
    ylabel( 'V [mV]' )
    if j == nfvals
        xlabel( 'Time [ms]' )
    end
end
ylims                           = [ min( ylims( :, 1 ) ) max( ylims( :, 2 ) ) ] + [ -1 1 ] * 0.1;
for j                           = 1 : nfvals
    subplot( nfvals, 2, 2 * ( j - 1 ) + 2 )
    set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims + [ -1 0 ] )
    f                           = fvals_slct( j );
    [ ~, fidx ]                 = min( abs( s.zfvals - f ) );
    x                           = s.Iemat( tidx, fidx );
    hold on
    ph                          = plot( t * 1000 ...
        , local_scale( x, 'linear' ) * 10 + ylims( 1 ) + 1, 'b' );
    set( ph, 'color', color_input );
end

%--------------------------------------------------------------------
% phase of coherence
Amp                             = 0.17;
[ ~, ~, ~, ~, s ]               = NR_sinusoids_to_cmodel( model, 'Tlong', Tlong, 'fvals', fvals ...
    , 'Dn', Dn, 'nreps', 1, 'Amp', Amp, 'Ibias', Ibias, 'graphics', [ 0 1 0 0 ], 'nFFT', 1250 );
fig2B( 3 )                      = gcf;

fig2B( 4 )                      = figure;

% data:
fin                             = s.fo( s.fidx );
Z                               = s.cohs( s.fidx );
Phi                             = s.phs( s.fidx );
pvals                           = s.pvals( :, 2 );
pTH                             = 0.05;
vidx                            = pvals <= pTH;
Phi( isnan( Z ) )               = NaN;

subplot( 3, 3, 1 ) % firing rate
ph                              = plot( s.fbins, s.frate_vec, '.-b' );
set( ph, 'color', color_output )
ylims                           = ylim;
ylims                           = [ 0 max( ylims( : ) ) ];
set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims )
ylabel( 'Firing rate [spks/s]' )

subplot( 3, 3, 4 ) % coherence
ph                              = plot( fin, Z, '.-k' );
set( ph, 'color', color_output )
hold on
ph                              = plot( fin( vidx ), Z( vidx ), '.-r' ); 
set( ph, 'color', [ 1 0 0 ])
ylims                           = ylim;
ylims                           = [ 0 max( ylims( : ) ) ];
set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims )
ylabel( 'Coherence' )
%axis square
set( gca, 'xtick', 0 : 10 : 40 )
[ ~, midx ]                     = max( Z ); 
fres                            = fin( midx );
line( [ 1 1 ] * fres, ylim, 'color', [ 0 0 0 ], 'linestyle', '--' );
title( sprintf( '%s, Amp=%0.2g; f_{res} = %0.2g Hz', model, Amp, fres ) )

subplot( 3, 3, 7 ) % phase
ylims                           = [ -1 1 ] * pi;
dy                              = ylims( 1 ) : diff( ylims ) / 8 : ylims( 2 );
uwp                             = unwrap( Phi );
if max( uwp ) > 2 * pi
    uwp                         = uwp - 2 * pi;
end
ph                              = plot( fin, uwp, '-' );
set( ph, 'color', [ 0 0 0 ], 'linewidth', 0.5 );
hold on
ph                              = plot( fin( vidx ), uwp( vidx ), '.-' ); 
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
mat                             = [ uwp( 1 : end - 1 ) < 0 uwp( 2 : end ) > 0 vidx( 1 : end - 1 ) ];
fphs                            = fin( all( mat, 2 ) ) + diff( fin( 1 : 2 ) ) / 2;
if isempty( fphs )
    fphs                        = NaN;
else
    line( [ 1 1 ] * fphs, ylim, 'color', [ 0 0 0 ], 'linestyle', '--' );
end
title( sprintf( 'f_{phs} = %0.2g Hz', fphs ) )

% add a few selected examples 
tidx                            = s.t >= tlim( 1 ) & s.t <= tlim( 2 );
t                               = s.t( tidx );
ylims                           = NaN( nfvals, 2 );
for j                           = 1 : nfvals
    f                           = fvals_slct( j );
    [ ~, fidx ]                 = min( abs( s.zfvals - f ) );
    y                           = s.Vmat( tidx, fidx );
    subplot( nfvals, 2, 2 * ( j - 1 ) + 2 )
    ph                          = plot( t * 1000, y, 'b' );
    set( ph, 'color', color_output );
    axis tight
    xlim( [ 0 max( xlim ) ] )
    ylims( j, : )               = ylim;
    title( sprintf( '%0.3g \\muA/cm^2; %0.3g Hz', Amp, f ) )
    ylabel( 'V [mV]' )
    if j == nfvals
        xlabel( 'Time [ms]' )
    end
end
ylims                           = [ min( ylims( :, 1 ) ) max( ylims( :, 2 ) ) ] + [ -1 1 ] * 0.1;
for j                           = 1 : nfvals
    subplot( nfvals, 2, 2 * ( j - 1 ) + 2 )
    set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims + [ -1 0 ] )
    f                           = fvals_slct( j );
    [ ~, fidx ]                 = min( abs( s.zfvals - f ) );
    x                           = s.Iemat( tidx, fidx );
    hold on
    ph                          = plot( t * 1000, local_scale( x ,'linear') * 10 + ylims( 1 ) + 1, 'b' );
    set( ph, 'color', color_input );
end


%--------------------------------------------------------------------
% 2D plots with Noise, Amp, and Ibias levels
sig                             = 'sines';
[ ~, ~, ~, figs ]               = NR_sinusoids_to_cmodel_run( model, sig, 'Dn' ...
    ,'nFFT', 1250, 'dflag', 0 );
fig2B( 5 )                      = figs;
[ ~, ~, ~, figs ]               = NR_sinusoids_to_cmodel_run( model, sig, 'Amp' ...
    , 'nFFT', 1250, 'dflag', 0 );
fig2B( 6 )                      = figs;
[ ~, ~, ~, figs ]               = NR_sinusoids_to_cmodel_run( model, sig, 'Ibias' ...
    , 'nFFT', 1250, 'dflag', 0 );
fig2B( 6 )                      = figs;

%--------------------------------------------------------------------
% check noise:
dn_vals                         = 0 : 0.005 : 2;
[ ~, ~, ~, afig, ~ ]            = NR_sinusoids_to_cmodel_run( model, sig ...
    , 'Dn', 'nFFT', 1250, 'param_vals', dn_vals, 'dflag', 0 );
fig2B( 7 : 8 )                  = afig;

amp_vals                        = 0 : 0.005 : 1;
[ ~, ~, ~, afig, ~ ]            = NR_sinusoids_to_cmodel_run( model, sig ...
    , 'Amp', 'nFFT', 1250, 'param_vals', amp_vals, 'dflag', 0 );
fig2B( 9 : 10 )                 = afig;

%-----
fig                             = fig2B;
if savef
    for i                       = 1 : length( fig )
        figname                 = [ outdir filesep prefix '_FIG2B-E_part' num2str( i ) ];
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
    
% ----------------------------------------------------------------------
% F-H: Subthreshold to network resonance in the Inap/h model (including examples)
% ----------------------------------------------------------------------

% how to run the simulations
model                           = 'Inap';
Dn                              = 0.0125;
Dn_i                            = 3;
Ibias                           = -1.85;
Amp                             = 0.14125;
Tlong                           = 20;
fvals                           = [ 0 40 ];

% how to plot the results
gr                              = 0;
fROI                            = [ 0 40 ];

% fingerprint analysis params
fROI_stfp                       = [ 0 40 ];
Fs_stfp                         = 10000;
aFs_stfp                        = 1250;
M_stfp                          = 1;
nFFT_stfp                       = NaN;

% summary params
nfigs                           = 6;
fig2F                           = zeros( nfigs, 1 );
ylims2                          = zeros( nfigs, 2 );
fig2G                           = zeros( nfigs, 1 );
fig2H                           = zeros( nfigs, 1 );

for i                           = 1 : nfigs
    switch i
        case 1
            Dn_e                = Dn;
            Amp_e               = Amp;
            Amp_i               = 0;
            [ ~, ~, ~, ~, s, sei ]  = NR_sinusoids_to_cmodel( model, 'Tlong', Tlong ...
                , 'fvals', fvals, 'Dn', Dn, 'nFFT', 1250 ...
                , 'nreps', 1, 'Amp', Amp_e, 'Ibias', Ibias, 'graphics', gr, 'N_cells', [ 1 1 ] ...
                , 'Dn_e', Dn_e, 'Dn_i', Dn_i, 'Gie', 1 );
        case 2
            Dn_e                = Dn;
            Amp_e               = Amp;
            Amp_i               = 0;
            [ ~, ~, ~, ~, s, sei ] 	= NR_sinusoids_to_cmodel( model, 'Tlong', Tlong ...
                , 'fvals', fvals, 'Dn', Dn, 'nFFT', 1250 ...
                , 'nreps', 1, 'Amp', Amp_e, 'Ibias', Ibias, 'graphics', gr, 'N_cells', 20 ...
                , 'Dn_e', Dn_e, 'Dn_i', Dn_i, 'einet_source', 'E' );
        case { 3, 4, 5, 6, 7 }
            Dn_e                = Dn;
            Amp_e               = 0;
            switch i
                case 3
                    Amp_i       = Amp;
                case 4
                    Amp_i       = Amp * 4;
                case 5
                    Amp_i       = Amp * 8;
                case 6
                    Amp_i       = Amp * 16;
                case 7
                    Amp_i       = Amp * 32;
            end
            [ ~, ~, ~, ~, s, sei ] 	= NR_sinusoids_to_cmodel( model, 'Tlong', Tlong ...
                , 'fvals', fvals, 'Dn', Dn_e, 'nFFT', 1250 ...
                , 'nreps', 1, 'Amp', Amp_i, 'Ibias', Ibias, 'graphics', gr, 'N_cells', 20 ...
                , 'Dn_e', Dn_e, 'Dn_i', Dn_i, 'einet_source', 'I' );
    end
    mfr_e                       = mean( sei.frates( 1 : sei.N_cells( 1 ) ) );
    mfr_i                       = mean( sei.frates( sei.N_cells( 1 ) + 1 : sum( sei.N_cells ) ) );
    close( gcf )
    close( gcf )
    
    % add a few selected examples
    Ne                          = sei.N_cells( 1 );
    Ni                          = sei.N_cells( 2 );
    N                           = Ne + Ni;
    t                           = sei.t / 1000;                             % [ms]
    Ve                          = sei.Ve;
    Vi                          = sei.Vi;
    Ie                          = sei.Ie;
    Ii                          = sei.Ii;
    isE                         = true( N, 1 );
    isE( Ne + 1 : N )           = 0;
    
    fig2F( i )                  = figure;
    
    sh1                         = axes( 'position', [ 0.13	0.35	0.775	0.55 ] );
    sh2                         = axes( 'position', [ 0.13	0.11	0.775	0.2 ] );
    
    subplot( sh1 )
    % plot the traces
    spc                         = [ ones( 1, Ne ) * 120 ones( 1, Ni ) * 110 ] + 40;
    [ ~, ph ]                   = plotTraces( t, [ Ve Vi ], -spc, 1, [ NaN NaN ] );
    set( ph( 1 : Ne ), 'color', colors_EI( 1, : ) )
    set( ph( Ne + 1 : N ), 'color', colors_EI( 2, : ) )
    set( gca, 'tickdir', 'out', 'box', 'off' );
    local_calibration( [ 1 100 ], [ 0.9 0.1 ] * complex( 0, 1 ), sh1, { 's', 'mV' } );
    % add the current inputs
    ylims0                      = ylim;
    f                           = max( spc );
    ylims                       = [ ylims0( 1 ) - 3.5 * f ylims0( 2 ) + 0.5 * f ];
    ylim( ylims )
    hold on
    levE                        = ylims0( 1 ) - 3 * f;
    levI                        = ylims0( 1 ) - 1.5 * f;
    den                         = max( [ Ie; Ii ] );
    ph                          = plot( t, Ie / den * f + levE, '-r', t, Ii / den * f + levI, '-b' );
    set( ph( 1 ), 'color', colors_EI_input( 1, : ) )
    set( ph( 2 ), 'color', colors_EI_input( 2, : ) )
    axis off
    
    % plot the coherences
    subplot( sh2 )
    hold on
    if Ne > 0
        ph                      = plot( s.fo, sei.cohs_ei( :, isE == 1 ), 'r' );
        set( ph, 'color', colors_EI_input( 1, : ) )
        ph                      = plot( s.fo, sei.mcoh_e, 'r' );
        set( ph, 'color', colors_EI( 1, : ) )
        set( ph, 'linewidth', 2 )
    end
    if Ni > 0
        ph                      = plot( s.fo, sei.cohs_ei( :, isE == 0 ), 'b' );
        set( ph, 'color', colors_EI_input( 2, : ) )
        ph                      = plot( s.fo, sei.mcoh_i, 'b' );
        set( ph, 'color', colors_EI( 2, : ) )
        set( ph, 'linewidth', 2 )
    end
    xlim( fROI )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Coherence' )
    xlabel( 'Frequency [Hz]' )
    tstr                        = sprintf( 'Amp E/I: %0.3g/%0.3g; Dn E/I: %0.3g/%0.3g mV; E/I: %0.3g/%0.3g spikes/s'...
        , Amp_e, Amp_i, Dn_e, Dn_i, mfr_e, mfr_i );
    title( tstr )
    
    % gather the coherence values
    ylims2( i, : )              = ylim;
    
    % plot fingerprint for an E-cell and for an I-cell
    eidx                        = 1;
    iidx                        = sei.N_cells( 1 ) + 1;
    switch i
        case { 1, 2 }
            Ievec               = sei.Ie( : );
        case { 3, 4, 5, 6, 7 }
            Ievec               = sei.Ii( : );
    end
    stvec                       = find( sei.sr( :, eidx ) );
    st_fingerprint( stvec, Ievec, 'spkFs', Fs_stfp, 'xFs', Fs_stfp ...
        , 'fROI', fROI_stfp, 'Fs', aFs_stfp, 'graphics', 1, 'M', M_stfp, 'nFFT', nFFT_stfp );
    fig2G( i )                  = gcf;

    stvec                       = find( sei.sr( :, iidx ) );
    st_fingerprint( stvec, Ievec, 'spkFs', Fs_stfp, 'xFs', Fs_stfp ...
        , 'fROI', fROI_stfp, 'Fs', aFs_stfp, 'graphics', 1, 'M', M_stfp, 'nFFT', nFFT_stfp );
    fig2H( i )                  = gcf;
    
end

ylims2                          = [ 0 max( ylims2( :, 2 ) ) ];
for i                           = 1 : nfigs
    figure( fig2F( i ) )
    ylim( ylims2 )
end

%-----
fig                             = [ fig2F( : ); fig2G( : ); fig2H( : ) ];
if savef
    for i                       = 1 : length( fig )
        figname                 = [ outdir filesep prefix '_FIG2F-H_part' num2str( i ) ];
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

return % NR_make_Fig2

%------------------------------------------------------------------------
function [ lh, th ] = local_calibration( len, xy, h, s, varargin )

% constants
loc0                        = [ 0.1 0.1 ]; % default xy in relative coordinates

% arguments
lh                          = [];
th                          = [];
nargs                       = nargin;
if nargs < 1 || isempty( len )
    len                     = -[ 0.1 0.1 ];
end
if length( len( : ) ) ~= 2
    fprintf( '%s: input size mismatch (len)', upper( mfilename ) )
    return
end
len                         = double( len );
if nargs < 2 || isempty( xy )
    xy                      = [];
end
xy                          = double( xy );
if nargs < 3 || isempty( h )
    h                       = gca;
end
if nargs < 4 || isempty( s ) || ~isa( s, 'cell' ) || length( s ) ~= 2
    if ~exist( 's', 'var' ) || ~isequal( s, 0 )
        s                   = { '' '' };
    end
end
if ~isempty( varargin )
    lineparams              = varargin;
else
    lineparams              = { 'color', [ 0 0 0 ], 'linewidth', 2 };
end

% get limits, xy
try
    subplot( h )
catch
    fprintf( '%s: input type mismatch (h)', upper( mfilename ) )
    return
end
xlims                       = xlim;
ylims                       = ylim;
dx                          = diff( xlims );
dy                          = diff( ylims );
if isempty( xy )
    xy                      = [ min( xlims ) + loc0( 1 ) * dx min( ylims ) + loc0( 2 ) * dy ];
end
if length( xy( : ) ) ~= 2
    fprintf( '%s: input size mismatch (xy)', upper( mfilename ) )
    return
end
if all( imag( xy ) > 0 )
    f                       = abs( imag( xy ) );
    xy                      = double( [ min( xlims ) min( ylims ) ] + [ dx dy ] .* f( : ).' );
end

% get len
dxy0                        = [ dx dy ];
dxy( len >= 0 )             = len( len >= 0 );
dxy( len < 0 & len >= -1 )  = dxy0( len < 0 & len >= -1 ) .* abs( len( len < 0 & len >= -1 ) );
dxy( len < -1 )             = NaN;
if sum( isnan( dxy ) )
    fprintf( '%s: input type mismatch (len)', upper( mfilename ) )
    return
end

% draw the line
x                           = xy( 1 ) + [ 0 0 dxy( 1 ) ];
y                           = xy( 2 ) + [ dxy( 2 ) 0 0 ];
lh                          = line( x, y, lineparams{ : } );

% add the text
if isequal( s, 0 )
    str1                    = '';
    str2                    = '';
else
    str1                    = sprintf( '%0.3g %s', dxy( 1 ), s{ 1 } );
    str2                    = sprintf( '%0.3g %s', dxy( 2 ), s{ 2 } );
end
if dxy( 1 ) == 0
    dxy( 1 )                = dx * 0.1;
    str1                    = '';
end
if dxy( 2 ) == 0
    dxy( 2 )                = dy * 0.1;
    str2                    = '';
end
th( 1 )                     = text( xy( 1 ) + dxy( 1 ) / 2, xy( 2 ) - dxy( 2 ) / 2, str1 );
set( th( 1 ), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'rotation', 0 )
th( 2 )                     = text( xy( 1 ) - dxy( 1 ) / 2, xy( 2 ) + dxy( 2 ) / 2, str2 );
set( th( 2 ), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'rotation', 90 )

return % local_calibration

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

% EOF
