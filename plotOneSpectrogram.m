% plotOneSpectrogram    plot a spectrogram/phasogram
% 
% call                  [ fig, ah, zout ] = plotOneSpectrogram( x, y, z, USF )
%
% gets                  x, y, z         time/phase; frequency; power (or any other measure)
%                       USF             {10}; up-sampling factor
%
% note                  can be called with empty x and/or y, but not with empty z
%
% returns               fig             handle to figure
%                       ah              handle to axis
%                       zout            the upsampled matrix
%
% calls                 myjet

% 13-mar-13 ES

% last update
% 30-jun-22

function [ fig, ah, zout ] = plotOneSpectrogram( x, y, z, USF )

% constants
sepColor                        = [ 1 1 1 ] * 1;                            % white
nCycles                         = 1.5;                                      % phase cycles (for phasograms)
ylab                            = 'Frequency [Hz]';

% arguments
if nargout >= 3
    out                         = 1;
else
    out                         = 0;
end
nargs                           = nargin;
if nargs < 3 || isempty( z )
    error( 'missing/empty argument' )
end
if isempty( x )
    x                           = 1 : size( z, 1 );
end
if isempty( y )
    y                           = 1 : size( z, 2 );
end
if length( x ) == size( z, 2 ) && length( y ) == size( z, 1 ) && size( z, 1 ) ~= size( z, 2 )
    z                           = permute( z, [ 2 1 3 ] );
end
if nargs < 4 || isempty( USF )
    USF                         = 10;
end

% preparations
x                               = x( : ).';
if all( x >= -pi & x <= pi ) || all( x >= 0 & x <= 2 * pi )
    istime                      = 0;
    xlab                        = 'Phase [rad]';
else
    istime                      = 1;
    xlab                        = 'Time';
end
if ~istime
    x                           = [ x - 2 * pi  x x + 2 * pi ];
end
n                               = size( z, 3 );
if n > 1
    error( 'single spectrogram only' )
end

% where to plot
newplot
ah                              = gca;
fig                             = gcf;
klab                            = 1;

% plot
if istime
    zout                        = NaN( size( z, 2 ), size( z, 1 ), n );
else
    zout                        = NaN( size( z, 2 ), 3 * size( z, 1 ), n );
end
for k                           = 1 : n
    subplot( ah( k ) )
    if istime
        zk                      = z( :, :, k )';
    else
        zk                      = repmat( z( :, :, k )', [ 1 3 ] );
    end
    [ ~, ch ]                   = contourf( x, y, zk, USF ^ 2 );
    set( ch, 'linestyle', 'none' );
    if out
        zout( :, :, k )         = zk;
    end
    title( sprintf( '%d: %0.3g', k, max( zk( : ) ) ) )
    if ~istime
        xlim( [ -1 1 ] * nCycles * pi )
        for ln                  = [ -1 0 1 ]
            line( ln * [ 1 1 ] * pi, ylim, 'color', sepColor, 'linestyle', '--', 'linewidth', 0.5 );
        end
        set( ah( k ), 'xtick', ( -2 * pi : pi : 2 * pi )...
            , 'xticklabel', round( ( -2 * pi : pi : 2 * pi ) ./ 0.01 ) * 0.01 );
    end
    set( ah( k ), 'box', 'off', 'tickdir', 'out' )
    if k == klab
        ylabel( ylab )
        if n <= 6
            xlabel( xlab )
        end
    else
        set( ah( k ), 'xticklabel', '', 'yticklabel', '' )
    end
end
for k                           = ( n + 1 ) : size( ah, 1 )
    subplot( ah( k ) )
    axis off
end

colormap( myjet )

return 

% EOF
