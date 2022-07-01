% NR_make_Fig8      graphics for Figure 8 in Stark, Levi, Rotstein, 2022
%
% This script creates the sub-figures that compose Figure 8 in 
% Stark, Levi, Rotstein, 2022, PLOS Computational Biology
%
% Figure 8:  Inhibition-induced network resonance is sharpened by presynaptic high-pass filtering
% (A-B) Gamma INT - subthreshold and spiking
% (C)   Network resonance (E,IE,IE-gamma) with chirps (single trial + 20 rasters for each)

% Calls      NR_eisim,NR_calc_z_spectral,st_coherence,st_fingerprint,NR_sinusoids_to_cmodel_images 
%
% 08-aug-21 ES

% last update
% 01-jul-22 AL

function NR_make_Fig8( savef, pstr, outdir )

% constants
prefix                          = 'resonance_network';
renderer_name                   = 'painters';
resize                          = '-bestfit';
color_output                    = [ 0 0 0 ];                                % output (Vm/spikes): black 

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

% ----------------------------------------------------------------------
% (A-B) Gamma INT - subthreshold and spiking
% ----------------------------------------------------------------------

fig14                       = [];

% parameters for coherence and firing rate calculations:
aFs                        	= 1250;                                     % st_coherence
aFs2                    	= 5000;                                     % st_fingerprint
fROI                     	= [ 0 80 ];
graphics0                	= 0;
graphics1                 	= 1;
M                        	= 1;
mtNW                     	= 3;
nFFT                       	= NaN;
dflag                     	= '';

%----------------------------------------------------------------
% plot some examples of gamma resonance (subthreshold and spiking)

%----------------------------------------------------------------
% subthreshold gamma resonance in the INT
MDL                         = 7;
MDL_params                  = [ 0.5 3.8 0 0 ];                          % [ Ain, Iapp, Gie (i to e...) D_i ]
[ t, ~, Vi_1, Ie ]    = NR_eisim( { MDL, MDL_params });
Fs                          = 1 / diff( t( 1 : 2 ) ) * 1000;

fig14( 1 )                  = figure;
subplot( 2, 1, 1 )
plot( t / 10000 * 80, Vi_1, 'b' ); 
hold on
plot( t / 10000 * 80, Ie * 2 - 71, 'k' )
set( gca, 'tickdir', 'out', 'box', 'off' ); 
ylim( [ -73 -62 ] );
xlim( fROI )
ylabel( 'Vm [mV]' )
xlabel( 'Frequency [Hz]' );
title( sprintf( 'MDL=%d; Ain=%0.3g; Iapp=%0.3g; Dn=%0.3g', MDL, MDL_params( 1 ), MDL_params( 2 ), MDL_params( 4 ) ) )

[ zI, phsI, foI, fidxI ]  	= NR_calc_z_spectral( Ie( : ), Vi_1( : ), 'xFs', Fs ...
        , 'fROI', fROI, 'Fs', aFs, 'graphics', graphics0, 'M', M, 'mtNW', mtNW );
fin                         = foI( fidxI );
Z                           = zI( fidxI );
Phi                         = phsI( fidxI );

xticks                      = 0 : 10 : fROI( 2 );

subplot( 2, 3, 4 )
ph                          = plot( fin, Z );                           % MOhm * mm^2
set( ph, 'color', [ 0 0 0 ], 'linewidth', 2 );
ylims                       = ylim;
ylim( [ 0 ylims( 2 ) ] )
xlim( fROI )
set( gca, 'tickdir', 'out', 'box', 'off' )
xlabel( 'Frequency [Hz]' )
ylabel( 'Impedance [M\Omega*mm^2]' )
axis square
set( gca, 'xtick', xticks )

subplot( 2, 3, 5 )
Phi( isnan( Z ) )           = NaN;
uwp                         = unwrap( Phi );
if max( uwp ) > 2 * pi
    uwp                     = uwp - 2 * pi;
end
ph                          = plot( fin, uwp );
set( ph, 'color', [ 0 0 0 ], 'linewidth', 2 );
xlim( fROI )
ylim( [ -1 1 ] * pi/2 * 1.1 )
line( xlim, [ 0 0 ], 'color', [ 0 0 0 ], 'linestyle', '--' )
line( xlim, [ 1 1 ] * pi/2, 'color', [ 0 0 0 ], 'linestyle', '--' )
line( xlim, [ 1 1 ] * -pi/2, 'color', [ 0 0 0 ], 'linestyle', '--' )
set( gca, 'tickdir', 'out', 'box', 'off' )
set( gca, 'ytick', ( -1 : 0.25 : 1 ) * pi / 2, 'yticklabel', round( ( -1 : 0.25 : 1 ) * pi / 2 * 100 ) / 100 )
xlabel( 'Frequency [Hz]' )
ylabel( 'Phase [rad]' )
axis square
set( gca, 'xtick', xticks )

%----------------------------------------------------------------
% spiking gamma resonance in the INT - at low Ain
MDL_params                  = [ 0.9 3.8 0 0.1 ];
[ t, ~, Vi_1, ~, Ie ]       = NR_eisim( { MDL, MDL_params } );

% plot raw:
fig14( 2 )                  = figure;
subplot( 2, 1, 1 )
plot( t / 10000 * 80, Vi_1, 'b' ); 
hold on
plot( t / 10000 * 80, Ie * 2 - 91, 'k' )
set( gca, 'tickdir', 'out', 'box', 'off' ); 
ylim( [ -100 40 ] );
xlim( fROI )
ylabel( 'Vm [mV]' )
xlabel( 'Frequency [Hz]' );
title( sprintf( 'MDL=%d; Ain=%0.3g; Iapp=%0.3g; Dn=%0.3g', MDL, MDL_params( 1 ), MDL_params( 2 ), MDL_params( 4 ) ) )

% compute coherence and firing rate
stimes                      = local_parse( find( Vi_1 > 0 ) );
st                          = round( mean( stimes, 2 ) );
[ cohI, phsI, foI, ~, fidxI ]   = st_coherence( st, Ie, 'spkFs', Fs, 'xFs', Fs ...
    , 'fROI', fROI, 'Fs', aFs, 'graphics', graphics0, 'M', M, 'mtNW', mtNW, 'nFFT', nFFT, 'dflag', dflag );
s.cohs                      = cohI;
s.fidx                      = fidxI;
s.fo                        = foI;
s.phs                       = phsI;

[ frateI, ~, ~, fbinsI ]  	= st_fingerprint( st, Ie, 'spkFs', Fs, 'xFs', Fs ...
    , 'fROI', fROI, 'Fs', aFs2, 'graphics', graphics1, 'M', M, 'nFFT', nFFT );
s.frate_vec                 = frateI;
s.fbins                     = fbinsI;
fig14( 3 )                  = gcf;

% firing rate and coherence
figure( fig14( 2 ) )
subplot( 2, 3, 4 ) % firing rate
ph                          = plot( s.fbins, s.frate_vec, '.-b' );
set( ph, 'color', color_output )
ylims                       = ylim;
ylims                       = [ 0 max( ylims( : ) ) ];
xlim( fROI )
set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims )
ylabel( 'Firing rate [spks/s]' )
xlabel( 'Frequency [Hz]' )
axis square

subplot( 2, 3, 5 ) % coherence
ph                          = plot( s.fo( s.fidx ), s.cohs( s.fidx ), '.-k' );
set( ph, 'color', color_output )
ylims                       = ylim;
ylims                       = [ 0 max( ylims( : ) ) ];
xlim( fROI )
set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims )
ylabel( 'Coherence' )
xlabel( 'Frequency [Hz]' )
axis square

%----------------------------------------------------------------
% spiking gamma resonance in the INT - at high Ain
MDL_params                  = [ 2.0 3.8 0 0.1 ];
[ t, ~, Vi_1, ~, Ie ]       = NR_eisim( { MDL, MDL_params } );

% plot raw:
fig14( 4 )                  = figure;
subplot( 2, 1, 1 )
plot( t / 10000 * 80, Vi_1, 'b' ); 
hold on
plot( t / 10000 * 80, Ie * 2 - 91, 'k' )
set( gca, 'tickdir', 'out', 'box', 'off' ); 
ylim( [ -100 40 ] );
xlim( fROI )
ylabel( 'Vm [mV]' )
xlabel( 'Frequency [Hz]' );
title( sprintf( 'MDL=%d; Ain=%0.3g; Iapp=%0.3g; Dn=%0.3g', MDL, MDL_params( 1 ), MDL_params( 2 ), MDL_params( 4 ) ) )

% compute coherence and firing rate
stimes                      = local_parse( find( Vi_1 > 0 ) );
st                          = round( mean( stimes, 2 ) );
Fs                          = 1 / diff( t( 1 : 2 ) ) * 1000;
[ cohI, phsI, foI, ~, fidxI ]   = st_coherence( st, Ie, 'spkFs', Fs, 'xFs', Fs ...
    , 'fROI', fROI, 'Fs', aFs, 'graphics', graphics0, 'M', M, 'mtNW', mtNW, 'nFFT', nFFT, 'dflag', dflag );
s.cohs                      = cohI;
s.fidx                      = fidxI;
s.fo                        = foI;
s.phs                       = phsI;

[ frateI, ~, ~, fbinsI ]  	= st_fingerprint( st, Ie, 'spkFs', Fs, 'xFs', Fs ...
    , 'fROI', fROI, 'Fs', aFs2, 'graphics', graphics1, 'M', M, 'nFFT', nFFT );
s.frate_vec                 = frateI;
s.fbins                     = fbinsI;
fig14( 5 )                  = gcf;

% firing rate and coherence
figure( fig14( 4 ) )
subplot( 2, 3, 4 ) % firing rate
ph                          = plot( s.fbins, s.frate_vec, '.-b' );
set( ph, 'color', color_output )
ylims                       = ylim;
ylims                       = [ 0 max( ylims( : ) ) ];
xlim( fROI )
set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims )
ylabel( 'Firing rate [spks/s]' )
xlabel( 'Frequency [Hz]' )
axis square

subplot( 2, 3, 5 ) % coherence
ph                          = plot( s.fo( s.fidx ), s.cohs( s.fidx ), '.-k' );
set( ph, 'color', color_output )
ylims                       = ylim;
ylims                       = [ 0 max( ylims( : ) ) ];
xlim( fROI )
set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims )
ylabel( 'Coherence' )
xlabel( 'Frequency [Hz]' )
axis square

%----------------------------------------------------------------
% grid search for Ain
% 
% parameter scaling:
% subthreshold -    start at 0.05, go in 0.05 steps to 0.75
% spikes -          start at Ain 0.8; get to DC around 3.0; go from 0.1 : 0.1 : 3.5
% overall -         start at 0.1, go to 2.1

reCompute                   = 0;
MDL                         = 7;                                        % INT gamma

% parameters for I-gamma simulations:
Iapp                        = 3.8;
Ain                         = 0.1 : 0.1 : 2.1;
nAin                        = length( Ain );
D_i                         = 0.1; 
nreps                       = 40;

% load or simply recompute (spikes only)
filename                    = sprintf( '%s/NR_eisim_MDL%d_%d_reps.mat', datadir, MDL, nreps );
if ~reCompute
    try
        L                   = load( filename, '-mat' );
        st                  = L.st;
        finput              = L.finput;
        t                   = L.t;
        nreps               = L.nreps;
        Iapp                = L.Iapp;
        Ain                 = L.Ain;
    catch
        reCompute           = 1;
    end
end
if reCompute
    st                      = cell( nAin, nreps );
    for i                   = 1 : nAin
        for j               = 1 : nreps
            fprintf( 1, 'MDL%d: Ain=%0.3g, rep %d/%d\n', MDL, Ain( i ), j, nreps )
            [ t, ~, Vi_1, ~, finput_i ] = NR_eisim( { MDL, [ Ain( i ), Iapp NaN D_i ] } );
            stimes          = local_parse( find( Vi_1 > 0 ) );
            if ~isempty( stimes )
                st{ i, j }	= round( mean( stimes, 2 ) );
            end
        end
    end
    finput                  = finput_i;
    save( filename, 'st', 'finput', 't', 'nreps', 'Iapp', 'Ain', 'D_i' );
    L.st                    = st;
    L.finput                = finput;
    L.nreps                 = nreps;
    L.Iapp                  = Iapp;
    L.Ain                   = Ain;
    L.D_i                   = D_i;
end

% parameters for coherence and firing rate calculations:
Fs                          = 1 / diff( L.t( 1 : 2 ) ) * 1000;

% analyze coherence together for every same-Ain trials:
nc                          = 513;
nf                         	= 66;
ttu                     	= 2;
cohs                        = NaN( nc, nAin );
phs                         = NaN( nc, nAin );
frates                   	= NaN( nf, nAin );
fprintf( 1, 'Computing coherences for %d Ains, %d trials/level:\n', nAin, ttu )
for i                     	= 1 : nAin
    fprintf( 1, '%0.2g ', Ain( i ) )
    stcells                	= st( i, 1 : ttu );
    Iemat                 	= finput * ones( 1, ttu);
    [ cohI, phsI, foI ]  	= st_coherence( stcells, Iemat, 'spkFs', Fs, 'xFs', Fs ...
        , 'fROI', fROI, 'Fs', aFs, 'graphics', graphics0, 'M', M, 'mtNW', mtNW, 'nFFT', nFFT, 'dflag', dflag );
    cohs( :, i )            = cohI;
    phs( :, i )             = phsI;
    
    [ frateI, ~, ~, fbinsI ]= st_fingerprint( stcells, Iemat, 'spkFs', Fs, 'xFs', Fs ...
        , 'fROI', fROI, 'Fs', aFs2, 'graphics', graphics0, 'M', M, 'nFFT', nFFT );
    frates( :, i )      	= frateI;        
end
fprintf( 'done!\n' );
fo                          = foI;
fidx                        = fo >= fROI( 1 ) & fo <= fROI( 2 );        % bandpass (used for normalization and significance testing)
cohs                     	= cohs( fidx, : );

% plot this summary:
fig14( 5 )              	= NR_sinusoids_to_cmodel_images( cohs, fo( fidx ), fROI, Ain, 'frates', frates, 'fbins', fbinsI, 'Vcan', 0.9 );

% simulate subthreshold (noiselesss) and compute impedances
if true
    Ain_sth                 = 0.05 : 0.05 : 0.7;
    nAin_sth            	= length( Ain_sth );
    zmat                  	= NaN( nc, nAin_sth );
    phsZ                 	= NaN( nc, nAin_sth );
    D_i_sth              	= 0;
    
    for i                 	= 1 : nAin_sth
        fprintf( 1, 'MDL%d: Ain=%0.3g\n', MDL, Ain_sth( i ) )
        [ t, ~, Vi_1, ~, finput_i ]  = NR_eisim( { MDL, [ Ain_sth( i ), Iapp NaN D_i_sth ] } );
        [ zI, phsI ]       	= NR_calc_z_spectral( finput_i( : ), Vi_1( : ), 'xFs', Fs ...
            , 'fROI', fROI, 'Fs', aFs, 'graphics', graphics0, 'M', M, 'mtNW', mtNW );
        zmat( :, i )        = zI;
        phsZ( :, i )      	= phsI;
    end
    
    fig14( 6 )          	= NR_sinusoids_to_cmodel_images( zmat( fidx, : ), fo( fidx ), fROI, Ain_sth, 'output', 'subthreshold', 'Vcan', 0.5 );
end

%-----
fig                         = fig14( fig14 ~= 0 );
if savef
    for i                   = 1 : length( fig )
        figname = [ outdir filesep prefix '_FIG8AB_part' num2str( i ) ];
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
% (C) Network resonance (E,IE,IE-gamma) with chirps (single trial + 20 rasters for each)
% ----------------------------------------------------------------------
    
fig8                       = [];

% parameters for coherence and firing rate calculations:
aFs                      	= 1250;                                     % st_coherence
aFs2                      	= 5000;                                     % st_fingerprint
fROI                        = [ 0 40 ];
graphics0                   = 0;
graphics1                   = 1;
M                           = 1;
mtNW                        = 3;
nFFT                        = NaN;
dflag                       = '';

% general parameters
reCompute                   = 0;
MDLs                        = [ 5 4 6 ];                                % E, IE, IE-gamma

% parameters for simulations:
k                           = 0;
nreps                       = 20;
IhFlag                      = 1;
filename                    = sprintf( '%s/NR_eisim_%d_reps_Ih%d.mat', datadir, nreps, IhFlag );
MDL4_params                 = { 4, [ 0.5 -0.5 0.4 0.1 0.1 ] }; % [ Ain, Iapp, Gie (i to e...) D_i D_e ]
MDL5_params                 = { 5, [ 0.2 -2.7 NaN 0.1 0.1 ] };
MDL6_params                 = { 6, [ 2.1  3.7 0.4 0.1 0.1 ] };

% load or simply recompute (spikes only)
if ~reCompute
    try
        L                   = load( filename, '-mat' );
    catch
        reCompute           = 1;
    end
end

if reCompute
    all_st_e                = cell( length( MDLs ), 1 );
    all_st_i                = cell( length( MDLs ), 1 );
    all_finput              = cell( length( MDLs ), 1 );
    all_t                   = cell( length( MDLs ), 1 );
end

% network resonance with chirps
for i                       = 1 : length( MDLs )
    
    st_e                    = cell( 1, nreps );
    st_i                    = cell( 1, nreps );
    for j                   = 1 : nreps
        if j                == 1
            graphics        = 1;
        else
            graphics        = 0;
        end
        if reCompute || j == 1
            
            if MDLs( i ) == 4
                MDL_call    = MDL4_params;
            elseif MDLs( i ) == 5
                MDL_call    = MDL5_params;
            elseif MDLs( i ) == 6
                MDL_call    = MDL6_params;
            end
            
            fprintf( 1, 'MDL%d: rep %d/%d\n', MDLs( i ), j, nreps )
            [ t, Ve_1, Vi_1, finput_e, finput_i ] = NR_eisim( MDL_call );

            stimes          = local_parse( find( Ve_1 > 0 ) );
            if ~isempty( stimes )
                st_e{ j }   = round( mean( stimes, 2 ) );
            end
            stimes          = local_parse( find( Vi_1 > 0 ) );
            if ~isempty( stimes )
                st_i{ j }   = round( mean( stimes, 2 ) );
            end
            
            
            if MDLs( i )    == 5
                finput      = finput_e;
            else
                finput      = finput_i;
            end
            
        end
        if j == 1
            k               = k + 1;
            fig8( 1, k )   = gcf;
        end
    end
    if ~reCompute
        st_e                = L.all_st_e{ i };
        st_i                = L.all_st_i{ i };
        finput              = L.all_finput{ i };
        t                   = L.all_t{ i };
    end
    
    % plot raster (w/ firing rate and spike timing resonance) and
    % fingerprint for I and for E cell
    
    for rc                  = 0 : 1
        
        if rc == 0
            st              = st_i;
            tpstr           = 'I-cell';
            rn              = 4;
        elseif rc == 1
            st              = st_e;
            tpstr           = 'E-cell';
            rn              = 2;
        end
        
        fig8( rn, k )  	= figure;
        
        % plot the raster
        subplot( 2, 1, 1 )
        fin                 = t / 10000 * 40;                           % 10 seconds, 40 Hz
        sh                  = plot_raster( st, fin );
        hold on
        plot( fin, ( finput / max( finput ) + 1 )/2 - 1, 'k' )
        ylim( [ -1.5 0.5 ] + [ 0 nreps ] )
        xlim( fROI )
        set( gca, 'tickdir', 'out', 'box', 'off' );
        title( sprintf( '%s; Model %d', tpstr, MDLs( i ) ) )
        
        % plot coherence and firing rate
        Ie                  = finput( : ) * ones( 1, length( st ) );
        Fs                  = 1 / diff( t( 1 : 2 ) ) * 1000;
        [ cohI, phsI, foI, ~, fidxI ]   = st_coherence( st, Ie, 'spkFs', Fs, 'xFs', Fs ...
            , 'fROI', fROI, 'Fs', aFs, 'graphics', graphics0, 'M', M, 'mtNW', mtNW, 'nFFT', nFFT, 'dflag', dflag );
        s.cohs              = cohI;
        s.fidx             	= fidxI;
        s.fo                = foI;
        s.phs               = phsI;
        
        [ frateI, ~, ~, fbinsI ]= st_fingerprint( st, Ie, 'spkFs', Fs, 'xFs', Fs ...
            , 'fROI', fROI, 'Fs', aFs2, 'graphics', graphics1, 'M', M, 'nFFT', nFFT );
        fig8( rn + 1, k ) 	= gcf;
        frateI( isnan( frateI ) ) = 0;
        s.frate_vec         = frateI;
        s.fbins             = fbinsI;
        subplot( 3, 2, 2 )
        title( sprintf( '%s; Model %d', tpstr, MDLs( i ) ) )
        
        % firing rate and coherence
        figure( fig8( rn, k ) )
        subplot( 2, 3, 4 )      % firing rate
        ph                   = plot( s.fbins, s.frate_vec, '.-b' );
        set( ph, 'color', color_output )
        ylims                = ylim;
        ylims                = [ 0 max( ylims( : ) ) ];
        set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims )
        ylabel( 'Firing rate [spks/s]' )
        xlabel( 'Frequency [Hz]' )
        axis square
        
        subplot( 2, 3, 5 )      % coherence
        ph                    = plot( s.fo( s.fidx ), s.cohs( s.fidx ), '.-k' );
        set( ph, 'color', color_output )
        ylims                 = ylim;
        ylims                 = [ 0 max( ylims( : ) ) ];
        set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims )
        ylabel( 'Coherence' )
        xlabel( 'Frequency [Hz]' )
        axis square
    end
    
    if reCompute
        all_st_e{ i }           = st_e;
        all_st_i{ i }           = st_i;
        all_finput{ i }         = finput;
        all_t{ i }              = t;
    end
    
end
if reCompute
    save( filename, 'all_st_e', 'all_st_i', 'all_finput', 'all_t' );
end    

%-----
fig                             = fig8( fig8 ~= 0 );
if savef
    for i                       = 1 : length( fig )
        figname = [ outdir filesep prefix '_FIG8C_part' num2str( i ) ];
        figure( fig( i ) );
        figi              	    = gcf;
        figi.Renderer      	    = renderer_name;
        pause( 0.2 )
        if isequal( pstr, '-dpdf' )
            rsz          	    = resize;
        else
            rsz          	    = '';
        end
        print( figi, pstr, figname, rsz )
    end
end


return

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

% EoF