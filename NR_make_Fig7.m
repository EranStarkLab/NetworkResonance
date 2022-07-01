% NR_make_Fig7      graphics for Figure 7 in Stark, Levi, Rotstein, 2022
%
% This script creates the sub-figures that compose Figure 7 in 
% Stark, Levi, Rotstein, 2022, PLOS Computational Biology

% Figure 7: Inhibition-induced network resonance can be inherited from the level of membrane potential fluctuations
% (A-B) Network resonance (E,IE,IE-gamma) with chirps (single trial + 20 rasters for each)
%
% Calls: st_coherence, st_fingerprint,NR_eisim,plot_raster   
%
% 08-aug-21 ES

% last update
% 01-jul-22 AL


function NR_make_Fig7( savef, pstr, outdir )

% constants
prefix                          = 'resonance_network';
renderer_name                   = 'painters';
resize                          = '-bestfit';
color_output                    = [ 0 0 0 ];                                % output (Vm/spikes): black 


% arguments
nargs = nargin;
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
% (A-B) Network resonance (E,IE,IE-gamma) with chirps (single trial + 20 rasters for each)
% ----------------------------------------------------------------------

fig13                       = [];

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
            fig13( 1, k )   = gcf;
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
        
        fig13( rn, k )  	= figure;
        
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
        fig13( rn + 1, k ) 	= gcf;
        frateI( isnan( frateI ) ) = 0;
        s.frate_vec         = frateI;
        s.fbins             = fbinsI;
        subplot( 3, 2, 2 )
        title( sprintf( '%s; Model %d', tpstr, MDLs( i ) ) )
        
        % firing rate and coherence
        figure( fig13( rn, k ) )
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
        ph                   = plot( s.fo( s.fidx ), s.cohs( s.fidx ), '.-k' );
        set( ph, 'color', color_output )
        ylims                = ylim;
        ylims                = [ 0 max( ylims( : ) ) ];
        set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims )
        ylabel( 'Coherence' )
        xlabel( 'Frequency [Hz]' )
        axis square
    end
    
    if reCompute
        all_st_e{ i }       = st_e;
        all_st_i{ i }       = st_i;
        all_finput{ i }     = finput;
        all_t{ i }          = t;
    end
    
end
if reCompute
    save( filename, 'all_st_e', 'all_st_i', 'all_finput', 'all_t' );
end    

%-----
fig                         = fig13( fig13 ~= 0 );
if savef
    for i                   = 1 : length( fig )
        figname = [ outdir filesep prefix '_FIG7_part' num2str( i ) ];
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