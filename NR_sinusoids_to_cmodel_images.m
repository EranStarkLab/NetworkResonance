% NR_sinusoids_to_cmodel_images    plot two matrices
%
% call                          fig = NR_sinusoids_to_cmodel_images( cohs, fo, fvals, Vals, varargin )
%
% gets                          cohs        coherences (nfreqs x nvals ) (overloaded with zmat)
%                               fo          frequencies evaluated (overloaded with zfvals)
%                               fvals       vector of frequencies (only last element is used)
%                             	Vals        vector of values (1 x nvals)
% 
% optional arguments (given as name/value pairs)
% 
%                               frates      {[]}    firing rates (nfrates x nvals) 
%                               fbins       {[]}    frequencies (nfrates)
%                               tstr        {''}    title string
%                               Vcan        {NaN}   candidate value (where to plot horizontal dashed line)
%                               ValStr      {''}    y-axis label
%                               output      {'spiking'} or 'subthreshold'
%                                                   'subthreshold' plots impedances
%
% calls                         ParseArgPairs
%                               myjet
%
% see also                      NR_sinusoids_to_cmodel

% 18-oct-21 ES

% last update
% 30-jun-22

function fig = NR_sinusoids_to_cmodel_images( cohs, fo, fvals, Vals, varargin )

% arguments
[ frates, fbins, tstr ...
    , Vcan, ValStr, output ...
    ]                           = ParseArgPairs(...
    { 'frates', 'fbins', 'tstr' ...
    , 'Vcan', 'ValStr', 'output' }...
    , { [], [], '' ...
    , NaN, '', 'spiking' }...
    , varargin{ : } );
if isnan( Vcan )
    Vcan = Vals( 1 );
end
nVals                           = length( Vals );
if ~isequal( nVals, size( cohs, 2 ) )
    error( 'input size mismatch' )
end
if ~isempty( frates )
    if ~isequal( nVals, size( frates, 2 ) )
        error( 'input size mismatch' )
    end
end
switch output
    case 'spiking'
    case 'subthreshold'
        zfvals                  = fo;
        zmat                    = cohs;
end

% cohs              coherences (nfreqs x nvals ) (overloaded with zmat)
% fo                frequencies evaluated (overloaded with zfvals)
% fvals             frequency range (only last element is used)
% nVals             number of values

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

if ~isempty( frates )
    subplot( 2, 2, 2 )
    imagesc( fbins, 1 : nVals, frates' )
    axis xy
    clims                       = get( gca, 'clim' );
    set( gca, 'clim', [ 0 clims( 2 ) ] )
    ch                          = colorbar;
    set( ch, 'tickdir', 'out', 'box', 'off' )
    title( ch, 'Firing rate [spks/s]' )
    set( ch, 'tickdir', 'out', 'box', 'off', 'FontName', 'Arial' )
    xlabel( 'Frequency [Hz]' )
    ylabel( ValStr )
    set( gca, 'tickdir', 'out', 'box', 'off', 'FontName', 'Arial' )
    ticks                     	= get( gca, 'ytick' );
    ticks                       = [ 1 ticks + 1 ];
    if ticks( end ) > nVals
        ticks( end )            = [];
    end
    set( gca, 'ytick', ticks, 'YTickLabel', Vals( ticks ) );
    [ ~, minidx ]           	= min( abs( Vals - Vcan ) );
    line( xlim, [ 1 1 ] * minidx, 'color', [ 1 1 1 ], 'linestyle', '--' );
end

if ~isempty( tstr )
    local_fig_title( tstr )
end

return % NR_sinusoids_to_cmodel_images

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
