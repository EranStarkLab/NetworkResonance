% plotTraces            in a matrix
%
% CALL                  [ py, ph, colors, ch ] = plotTraces( x, y, spacing, scaling, calib )
%
% GETS                  x             vector, m-elements
%                       y             matrix, m by n (trials)
%                       spacing       {1} scalar, multiples of max peak-to-peak of data OR absolute
%                                         number (negative)
%                       scaling       {1}; scalar (multiply each column by the same number) 
%                                         or two-element vector, e.g. [ 0 1 ] scales the matrix between 0-1
%                                         or four-element vector, also scales by first row
%                                         of the matrix so [ 0 1 0 9 ] arranges the traces to
%                                         amplitude 0-1 for all together then distributes between 0-9
%                                         (just shifting, no scaling) 
%                       calib         {[1 1]} calibration bar (original x, y units)
%
% RETURNS               py              returns the scaled values
%
% DOES                  plots
%                       to only get the scaled values, use py = plotTraces...
%                       to plot and get the scaled values, use [ py, ph, .. ] = plotTraces... 
% 
% CALLS                 myjet

% 05-oct-12 ES

% last update
% 01-jul-22

function [ py, ph, colors, ch ] = plotTraces( x, y, spacing, scaling, calib )

% constants
CALIBCOLOR                      = [ 0 0.7 0 ];
CALIBWIDTH                      = 2;

% arguments
nargs                           = nargin;
nout                            = nargout;
if nargs == 1
    y                           = x;
    x                           = 1 : size( y, 1 );
end
if isempty( x ) || isempty( y )
    error( 'missing arguments' )
end
x                               = double( x );
y                               = double( y );
if nargs < 3 || isempty( spacing )
    spacing                     = 1;
end
spacing                         = spacing( 1 );
if nargs < 4 || isempty( scaling )
    scaling                     = 1;
end
if length( scaling ) > 4
    scaling                     = scaling( 1 : 4 );
end
if nargs < 5 || isempty( calib )
    calib                       = [ 1 1 ];                                  % x, y calibration bar
end
if length( calib ) < 2
    calib                       = calib * [ 1 1 ];
end
    
x                               = x( : );
if ~ismatrix( y )
    error( 'higher than 2D array plotting not supported' )
end
[ m, n ]                        = size( y );
if length( x ) ~= m
    y                           = y';
    [ m, n ]                    = size( y );
end
if length( x ) ~= m 
    error( 'input size mismatch' )
end

% preparations
p2p                             = max( max( y ) - min( y ) );
if spacing == 0
    oset                        = zeros( 1, n );
elseif spacing > 0
    oset                        = ( 0 : p2p / n : p2p - p2p / n ) * spacing;
else
    oset                        = 0 : -spacing : ( n - 1 ) * ( -spacing );
end
if isempty( oset )
    oset                        = 0;
end
if length( scaling ) == 1
    py                          = y * scaling + ones( m, 1 ) * oset;
    cy                          = calib( 2 ) * scaling;
elseif length( scaling ) == 2
    py                          = y + ones( m, 1 ) * oset;
    py                          = ( py - min( py( : ) ) ) / (  max( py( : ) ) -  min( py( : ) ) );
    py                          = py * diff( scaling ) + scaling( 1 );
    cy                          = calib( 2 ) * diff( scaling ) + scaling( 1 );
elseif length( scaling ) == 4
    py                          = y + ones( m, 1 ) * oset;
    py                          = ( py - min( py( : ) ) ) / (  max( py( : ) ) -  min( py( : ) ) );
    py                          = py * diff( scaling( 1 : 2 ) ) + scaling( 1 );
    cy                          = calib( 2 ) * diff( scaling( 1 : 2 ) ) + scaling( 1 );
    bias                        = scaling( 3 );
    scale                       = diff( scaling( 3 : 4 ) );
    py                          = bsxfun( @minus, py, py( 1, : ) - scale * ( 1/(n+1) : 1/(n+1) : 1-1/(n+1) ) ) + bias;
else
    py                          = y + ones( m, 1 ) * oset;
    cy                          = calib( 2 );
end

% plot
if nout == 0 || nout > 1
    % traces
    map                         = myjet; 
    nmap                        = size( map, 1 );
    colors                      = map( ceil( nmap / n : nmap / n : nmap ), : );
    if ~ishold
        newplot
    end
    ph                          = plot( x, py );
    set( gca, 'tickdir', 'out', 'box', 'off', 'xlim', x( [ 1 end ] ) )
    for i                       = 1 : n
        set( ph( i ), 'color', colors( i, : ) );
    end
    ylims                       = [ min( py( : ) ) max( py( : ) ) ];
    if ylims( 1 ) ~= ylims( 2 )
        ylim( ylims )
    end
    % calibration
    if calib( 1 ) > 0 || calib( 2 ) > 0
        calib                   = [ calib( 1 ) cy ];
        xlims                   = xlim;
        x0                      = xlims( 1 ) + diff( xlims ) * 0.05;
        y0                      = min( py( 1, : ) ) + range( py( 1, : ) ) * 0.15;
        ch                      = line( x0 + [ 0 0 calib( 1 ) ], y0 + [ calib( 2 ) 0 0 ] ...
            , 'color', CALIBCOLOR, 'linewidth', CALIBWIDTH );
    end
end

return

% EOF
