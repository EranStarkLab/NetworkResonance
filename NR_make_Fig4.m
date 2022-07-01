% NR_make_Fig4      graphics for Figure 4 in Stark, Levi, Rotstein, 2022
%
% This script creates the sub-figures that compose Figure 4 in 
% Stark, Levi, Rotstein, 2022, PLOS Computational Biology

% Figure 4: Resonance generated at the spiking level can be inherited to the network level
% (A-E). Spiking to network resonance in the LIF model
%
% Calls      NR_sinusoids_to_cmodel, plotTraces, st_fingerprint, ParseArgPairs
%
% 08-aug-21 ES

% last update
% 01-jul-22 AL

function NR_make_Fig4( savef, pstr, outdir, varargin)

% constants
prefix                          = 'resonance_network';
renderer_name                   = 'painters';
resize                          = '-bestfit';

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
[ rerun_analysis_fig4 ...
    ]                           = ParseArgPairs(...
    { 'rerun_analysis_fig4' }...
    , { 0 } ...
    , varargin{ : } );

% set plotting parameters:
datadir                         = outdir;

% ----------------------------------------------------------------------
% (A-E) Spiking to network resonance in the LIF model (including examples)
% ----------------------------------------------------------------------

% how to run the simulations
model                           = 'LIF';
Dn_es                           = [ 0.02 0.04 0.08 0.16 0.24 ];
Ibias                           = 0.9;
Amp                             = 0.115;                                    % supraTH
Dn_i                            = 2;

Tlong                           = 20;
fvals                           = [ 0 40 ];                                 % Hz

gr                              = 0;
fROI                            = [ 0 40 ];

% fingerprint analysis params
fROI_stfp                       = [ 0 40 ];
Fs_stfp                         = 10000;
aFs_stfp                        = 1250;
M_stfp                          = 1;
nFFT_stfp                       = NaN;

% summary graphics
nfigs                           = 11;
fig4A                           = zeros( nfigs, 1 ); % fig7 -> fig4A
fig4B                           = zeros( nfigs, 1 ); % fig7E -> fig4B
fig4C                           = zeros( nfigs, 1 ); % fig7I -> fig4C
ylims2                          = zeros( nfigs, 2 );

for i                           = 1 : nfigs
    switch i
        case 1
            Dn_e                = Dn_es( 1 );
            Amp_e               = Amp;
            Amp_i               = 0;
            [ ~, ~, ~, ~, s, sei ]           = NR_sinusoids_to_cmodel( model, 'Tlong', Tlong ...
                , 'fvals', fvals, 'Dn', Dn_e, 'nFFT', 1250 ...
                , 'nreps', 1, 'Amp', Amp_e, 'Ibias', Ibias ...
                , 'graphics', gr, 'N_cells', [ 1 1 ], 'Dn_e', Dn_e ...
                , 'Dn_i', Dn_i, 'Gie', 1 );
            
        case { 2, 3, 4, 5, 6 }
            Dn_e                = Dn_es( i - 1 );
            Amp_e               = Amp;
            Amp_i               = 0;
            [ ~, ~, ~, ~, s, sei ]           = NR_sinusoids_to_cmodel( model, 'Tlong', Tlong ...
                , 'fvals', fvals, 'Dn', Dn_e, 'nFFT', 1250 ...
                , 'nreps', 1, 'Amp', Amp_e, 'Ibias', Ibias ...
                , 'graphics', gr, 'N_cells', 20, 'Dn_e', Dn_e ...
                , 'Dn_i', Dn_i, 'einet_source', 'E' );
            
        case { 7, 8, 9, 10, 11 }
            Dn_e                = Dn_es( 1 );
            Amp_e               = 0;
            switch i
                case 7
                    Amp_i       = Amp;
                case 8
                    Amp_i       = Amp * 4;
                case 9
                    Amp_i       = Amp * 8;
                case 10
                    Amp_i       = Amp * 16;
                case 11
                    Amp_i       = Amp * 32;
            end
            [ ~, ~, ~, ~, s, sei ]           = NR_sinusoids_to_cmodel( model, 'Tlong', Tlong ...
                , 'fvals', fvals, 'Dn', Dn_e, 'nFFT', 1250 ...
                , 'nreps', 1, 'Amp', Amp_i, 'Ibias', Ibias ...
                , 'graphics', gr, 'N_cells', 20, 'Dn_e', Dn_e ...
                , 'Dn_i', Dn_i, 'einet_source', 'I' );
            
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
    
    fig4A( i )                  = figure;
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
    tstr = sprintf( 'Amp E/I: %0.3g/%0.3g; Dn E/I: %0.3g/%0.3g mV; E/I: %0.3g/%0.3g spikes/s' ...
        , Amp_e, Amp_i, Dn_e, Dn_i, mfr_e, mfr_i );
    title( tstr )
    
    % gather the coherence values
    ylims2( i, : )              = ylim;
    
    % plot fingerprint for an E-cell and for an I-cell
    eidx                        = 1;
    iidx                        = sei.N_cells( 1 ) + 1;
    switch i
        case { 1, 2, 3, 4, 5, 6 }
            Ievec               = sei.Ie( : );
        case { 7, 8, 9, 10, 11 }
            Ievec               = sei.Ii( : );
    end
    stvec                       = find( sei.sr( :, eidx ) );
    st_fingerprint( stvec, Ievec, 'spkFs', Fs_stfp, 'xFs', Fs_stfp ...
        , 'fROI', fROI_stfp, 'Fs', aFs_stfp, 'graphics', 1, 'M', M_stfp, 'nFFT', nFFT_stfp );
    tstr                        = sprintf( 'E; Amp E/I: %0.3g/%0.3g; Dn E/I: %0.3g/%0.3g mV' ...
        , Amp_e, Amp_i, Dn_e, Dn_i );
    local_fig_title( tstr );
    fig4B( i )                  = gcf;

    stvec                       = find( sei.sr( :, iidx ) );
    st_fingerprint( stvec, Ievec, 'spkFs', Fs_stfp, 'xFs', Fs_stfp ...
        , 'fROI', fROI_stfp, 'Fs', aFs_stfp, 'graphics', 1, 'M', M_stfp, 'nFFT', nFFT_stfp );
    tstr                        = sprintf( 'I; Amp E/I: %0.3g/%0.3g; Dn E/I: %0.3g/%0.3g mV' ...
        , Amp_e, Amp_i, Dn_e, Dn_i );
    local_fig_title( tstr );
    fig4C( i )                  = gcf;
    
    
end

ylims2                          = [ 0 max( ylims2( :, 2 ) ) ];
for i                           = 1 : nfigs
    figure( fig4A( i ) )
    ylim( ylims2 )
end

%-------------------------------------------------------------------
% network by noise level
if ~rerun_analysis_fig4
    try
        % load precomputed results results
        load( [ datadir '/NR_sinusoids_to_cmodel_LIF_chirp_einet_16E_4I_by_noise' ] ...
            , 'Dns', 's_keep', 'sei_keep', 'mcoh_e', 'mcoh_i', 'f_res', 'f_res_sem', 'f_coh', 'f_coh_sem' )
    catch
        rerun_analysis_fig4     = 1;
    end
end

if rerun_analysis_fig4
    
    % parameters
    model                       = 'LIF';
    Ibias                       = 0.9;
    Amp                         = 0.115;                                    % supraTH
    Dn_i                        = 2;
    Tlong                       = 20;
    fvals                       = [ 0 40 ];                                 % Hz
    gr                          = 0;
    gr_ei                       = [ 0 0 ];
    Dns                         = 0 : 0.01 : 0.3;
    aFs                         = 1250;
    
    % initialize
    nFFT                        = aFs;
    Fs                          = aFs;
    fROI                        = [ min( fvals ) max( fvals ) ];
    fo                          = ( 0 : nFFT / 2 )' * Fs / nFFT;            % all frequencies up to Nyquist
    fidx                        = fo >= fROI( 1 ) & fo <= fROI( 2 );        % bandpass (used for normalization and significance testing)
    fbins                       = fo( fidx );
    nfbins                      = length( fbins );
    
    nDns                        = length( Dns );
    f_res                       = NaN( nDns, 2 );
    f_res_sem                   = NaN( nDns, 2 );
    f_coh                       = NaN( nDns, 2 );
    f_coh_sem                   = NaN( nDns, 2 );
    mcoh_e                      = NaN( nfbins, nDns );
    mcoh_i                      = NaN( nfbins, nDns );
    
    % compute
    for i                       = 1 : nDns
        Dn                      = Dns( i );
        fprintf( 1, '%d/%d: Dn = %0.3g\n', i, length( Dns ), Dn )
        
        [ ~, ~, ~, ~, s, sei ]           = NR_sinusoids_to_cmodel( model, 'Tlong', Tlong ...
            , 'fvals', fvals, 'Dn', Dn, 'nFFT', nFFT, 'aFs', aFs ...
            , 'nreps', 1, 'Amp', Amp, 'Ibias', Ibias, 'graphics', gr, 'graphics_ei', gr_ei ...
            , 'N_cells', 20, 'Dn_e', Dn, 'Dn_i', Dn_i, 'einet_source', 'E' );
        
        mcoh_e( :, i )          = sei.mcoh_e( s.fidx );
        mcoh_i( :, i )          = sei.mcoh_i( s.fidx );
        
        f_res( i, 1 )           = mean( sei.max_coh_frq_ei( 1 : sei.N_cells( 1 ) ) );
        f_res( i, 2 )           = mean( sei.max_coh_frq_ei( sei.N_cells( 1 ) + 1 : sum( sei.N_cells ) ) );
        f_res_sem( i, 1 )       = local_calc_sem( sei.max_coh_frq_ei( 1 : sei.N_cells( 1 ) ) );
        f_res_sem( i, 2 )       = local_calc_sem( sei.max_coh_frq_ei( sei.N_cells( 1 ) + 1 : sum( sei.N_cells ) ) );
        f_coh( i, 1 )           = mean( sei.max_coh_val_ei( 1 : sei.N_cells( 1 ) ) );
        f_coh( i, 2 )           = mean( sei.max_coh_val_ei( sei.N_cells( 1 ) + 1 : sum( sei.N_cells ) ) );
        f_coh_sem( i, 1 )       = local_calc_sem( sei.max_coh_val_ei( 1 : sei.N_cells( 1 ) ) );
        f_coh_sem( i, 2 )       = local_calc_sem( sei.max_coh_val_ei( sei.N_cells( 1 ) + 1 : sum( sei.N_cells ) ) );
        
        fprintf( 1, '\n\n' )
    end
    s_keep                      = s;
    sei_keep                    = sei;
    
    % save the results
    save( [ datadir '/NR_sinusoids_to_cmodel_LIF_chirp_einet_16E_4I_by_noise' ] ...
        , 'Dns', 's_keep', 'sei_keep', 'mcoh_e', 'mcoh_i', 'f_res', 'f_res_sem', 'f_coh', 'f_coh_sem' )

end

% plot
fig4D                           = figure;
colormap( myjet )

subplot( 2, 2, 1 )
imagesc( s_keep.fo( s_keep.fidx ), 1 : length( Dns ), mcoh_e' )
axis xy
clim                            = get( gca, 'clim' );
set( gca, 'clim', [ 0 clim( 2 ) ] );
ah                              = gca;
ch                              = colorbar;
set( ch, 'tickdir', 'out', 'box', 'off' )
title( ch, 'Coherence' )
subplot( ah( 1 ) )
yticks                          = get( gca, 'ytick' );
yticks                          = [ 1 yticks + 1 ];
yticklabels                     = Dns;
set( gca, 'ytick', yticks, 'yticklabel', yticklabels( yticks ) )
xlabel( 'Frequency [Hz]' )
ylabel( 'Dn [mV]' )
set( gca, 'tickdir', 'out', 'box', 'off', 'FontName', 'Arial' )
title( 'E-cells' )

subplot( 2, 2, 2 )
imagesc( s_keep.fo( s_keep.fidx ), 1 : length( Dns ), mcoh_i' )
axis xy
clim                            = get( gca, 'clim' );
set( gca, 'clim', [ 0 clim( 2 ) ] );
ah                              = gca;
ch                              = colorbar;
set( ch, 'tickdir', 'out', 'box', 'off' )
title( ch, 'Coherence' )
subplot( ah( 1 ) )
xlabel( 'Frequency [Hz]' )
ylabel( 'Dn [mV]' )
yticks                          = get( gca, 'ytick' );
yticks                          = [ 1 yticks + 1 ];
yticklabels                     = Dns;
set( gca, 'ytick', yticks, 'yticklabel', yticklabels( yticks ) )
set( gca, 'tickdir', 'out', 'box', 'off', 'FontName', 'Arial' )
title( 'I-cells' )

subplot( 2, 2, 3 )
local_patch_band( Dns, f_res( :, 1 ), f_res_sem( :, 1 ), colors_EI( 1, : ), colors_EI_input( 1, : ) );
hold on
local_patch_band( Dns, f_res( :, 2 ), f_res_sem( :, 2 ), colors_EI( 2, : ), colors_EI_input( 2, : ) );
ph = plot( Dns, f_res( :, 1 ), '.-r', Dns, f_res( :, 2 ), '.-b' );
set( ph( 1 ), 'color', colors_EI( 1, : ) );
set( ph( 2 ), 'color', colors_EI( 2, : ) );
lh = legend( ph, { 'E-cells', 'I-cells' }, 'Location', 'best' );
set( lh, 'box', 'off' )
set( gca, 'tickdir', 'out', 'box', 'off', 'FontName', 'Arial' )
xlabel( 'Noise level [mV]' )
ylabel( 'Resonant frequency [Hz]' )
title( sprintf( 'EI net (N=%d E- and %d I-cells); Gie=%0.2g, Gei=%0.2g' ...
    , sei_keep( 1 ).N_cells( 1 ), sei_keep( 1 ).N_cells( 2 ), sei_keep( 1 ).Gie, sei_keep( 1 ).Gei ) )

subplot( 2, 2, 4 )
local_patch_band( Dns, f_coh( :, 1 ), f_coh_sem( :, 1 ), colors_EI( 1, : ), colors_EI_input( 1, : ) );
hold on
local_patch_band( Dns, f_coh( :, 2 ), f_coh_sem( :, 2 ), colors_EI( 2, : ), colors_EI_input( 2, : ) );
ph = plot( Dns, f_coh( :, 1 ), '.-r', Dns, f_coh( :, 2 ), '.-b' );
set( ph( 1 ), 'color', colors_EI( 1, : ) );
set( ph( 2 ), 'color', colors_EI( 2, : ) );
lh = legend( ph, { 'E-cells', 'I-cells' }, 'Location', 'best' );
set( lh, 'box', 'off' )
set( gca, 'tickdir', 'out', 'box', 'off', 'FontName', 'Arial' )
xlabel( 'Noise level [mV]' )
ylabel( 'Coherence' )
title( sprintf( 'EI net (N=%d E- and %d I-cells); Gie=%0.2g, Gei=%0.2g' ...
    , sei_keep( 1 ).N_cells( 1 ), sei_keep( 1 ).N_cells( 2 ), sei_keep( 1 ).Gie, sei_keep( 1 ).Gei ) )

%-------------------------------------------------------------------
fig                             = [ fig4A( : ); fig4B( : ); fig4C( : ); fig4D( : ) ];
if savef
    for i                       = 1 : length( fig )
        figname                 = [ outdir filesep prefix '_FIG4A-E_part' num2str( i ) ];
        figure( fig( i ) );
        figi                    = gcf;
        figi.Renderer           = renderer_name;
        pause( 0.2 )
        if isequal( pstr, '-dpdf' )
            rsz          	= resize;
        else
            rsz          	= '';
        end
        print( figi, pstr, figname, rsz )
    end
end

return % NR_make_Fig4

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
function ph = local_patch_band( x, y, s, c1, c2, sy, nanflag )

ph = [];

nargs = nargin;
if nargs < 2 || ~isequal( numel( x ), numel( y ) ) || isempty( x )
    return
else
    x                           = x( : );
    y                           = y( : );
end
if isempty( s )
    s                           = zeros( size( x ) );
else
    s = s( : );
end
if nargs < 4 || isempty( c1 ) || length( c1 ) ~= 3
    c1                          = [ 0 0 1 ];
else
    c1                          = c1( : )';
end
if nargs < 5 || isempty( c2 ) || length( c2 ) ~= 3
    c2                          = c1;
    c2( c2 == 0 )               = 0.3;
else
    c2                          = c2( : )';
end
if nargs < 6 || isempty( sy )
    sy                          = zeros( size( y ) );
else
    sy                          = sy( : );
end
if nargs < 7 || isempty( nanflag )
    nanflag                     = 0;
end

if nanflag
    nans                        = isnan( x ) | isnan( y );
    x( nans )                   = [];
    y( nans )                   = [];
    s( nans )                   = [];
    sy( nans )                  = [];
end
if isempty( x )
    return
end

h0                              = ishold;
if ~h0
    hold on
end
ph( 2 )                         = patch( [ x + sy; flipud( x - sy ) ], [ y + s; flipud( y - s ) ], c2 );
set( ph( 2 ), 'edgecolor', c2 )
ph( 1 )                         = line( x, y, 'color', c1, 'linewidth', 2 );
if ~h0
    hold off
end

return % local_patch_band

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