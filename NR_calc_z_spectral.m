% NR_calc_z_spectral        impedance and phase in the frequency domain (any input)
% 
% call                      [ z, phs, fo ] = NR_calc_z_spectral( I, V )
%                           [ ..., fidx, cohs, ccohs ] = NR_calc_z_spectral( I, V )
%
% gets                      I               [mA]    current
%                           V               [mV]    voltage
%
% optional arguments (given as name/value pairs):
%
%       sampling parameters:
%
%                           xFs         {1250}      [Hz], sampling rate of x
%                           Fs          {1250}      [Hz], intervally-used rate
%                           fROI        {[0 20]}    [Hz], frequency ROI
%                           dflagtrial  {0}         logical; remove mean from every trial (column of st)
%
%       multi-taper parameters:  
%
%                           M           {1}         determinant of spectral resolution (1/M Hz)
%                           nFFT        {NaN}       alternative manner to determine spectral resolution
%                           mtNW        {3}         MT: coherence time-freq parameter
%                           dflag       {0}         MT: no detrending (may remove mean for each segment separately)
%                                                     argument to detrend (''/[], 'constant'/0, 'linear'/1, 2, ...)
%                           nSpecSmooth {5}         number of FFT frequency bins (used only for normalized coherence, ccohs)
%       flow control:
%                           graphics    {0}         
%
% returns                   z           [Ohm]       impedance
%                           phs         [rad]       phase
%                           fo          [Hz]        frequency bins
%                           fidx        logical vector with nfidx true bits
%                                           nfidx corresponds to the number of frequency bins within fROI
%                           cohs        coherence: vector of nFFT/2+1
%                           ccohs       scaled coherence (scaled by the input spectrum within fROI)
%
% does
% (1) (optionally) removes mean for each trial separately, then concatenates
% (2) computes coherence and phase
% (3) computes impedance := sqrt( power( V ) / power( I ) )
% (4) (optionally) plots
%
% notes:
% (1) in contrast to st_coherence for field-field, here the default dflag is 0 
% (2) coherence is NOT a good metric for impedance, since it merely measures
% that consistency between the signals (linear transformations are
% irrelevant). however, the phase of the cross spectrum is relevant 
%
% calls                     ParseArgPairs
%
% see also                  st_coherence (spike-field, spike-spike)
%                           st_fingerprint (spike-field)

% 05-aug-21 ES

% last update:
% 01-jul-22

function [ z, phs, fo, fidx, cohs, ccohs ] = NR_calc_z_spectral( x, y, varargin )

%------------------------------------------------------------------------
% constants

% sampling rates and ROI
xFs_DFLT                        = 1250;                                     % [Hz]
Fs_DFLT                         = 1250;                                     % [Hz]
fROI_DFLT                       = [ 0 20 ];                                 % [Hz]

% coherence analyses
dflagtrial_DFLT                 = 0;                                        % 1: remove mean of every trial; 0: do not
M_DFLT                      	= 1;                                        % determinant of spectral resolution (1/M Hz)
nFFT_DFLT                       = NaN;                                      % frequency steps will be Fs / nFFT
mtNW_DFLT                     	= 3;                                        % MT: coherence time-freq parameter
dflag_DFLT                   	= 0;                                        % MT: 'constant' detrending (may also remove mean for each trial separately, dflagtrial)

% for normalizing the coherence
nSpecSmooth_DFLT             	= 5;                                        % number of FFT frequency bins (only for normalized coherence)

%------------------------------------------------------------------------
% arguments
nargs                           = nargin;
if nargs < 2 || isempty( x ) || isempty( y )
    to_return                   = 1;
else
    to_return                   = 0;
end
ca                              = varargin;
dflag_arg                       = NaN;
for i = 1 : length( ca )
    if isequal( lower( ca{ i } ), 'dflag' )
        dflag_arg               = ca{ i + 1 };
        break
    end
end
[ xFs, Fs ...
    , fROI ...
    , M, nFFT, mtNW, dflag, dflagtrial ...
    , nSpecSmooth ...
    , graphics ]                = ParseArgPairs(...
    { 'xFs', 'Fs' ...
    , 'fROI' ...
    , 'M', 'nFFT', 'mtNW', 'dflag', 'dflagtrial' ...
    , 'nSpecSmooth' ...
    , 'graphics' }...
    , { xFs_DFLT, Fs_DFLT ...
    , fROI_DFLT ...
    , M_DFLT, nFFT_DFLT, mtNW_DFLT, dflag_DFLT, dflagtrial_DFLT ...
    , nSpecSmooth_DFLT ...
    , 0 } ...
    , varargin{ : } );
if isempty( dflag_arg ) || ~isnan( dflag_arg )
    dflag                       = dflag_arg;
elseif isnan( dflag_arg )
    dflag                       = dflag_DFLT;
end
if round( xFs / Fs ) ~= ( xFs / Fs )
    error( 'downsampling from xFs to Fs must be by an interger ratio' )
end

% determine frequencies for spectral analysis
if isnan( nFFT )
    nFFT                        = 2 ^ floor( log2( M * Fs ) );              % ~1/M Hz resolution 
end
fo                              = ( 0 : nFFT / 2 )' * Fs / nFFT;            % all frequencies up to Nyquist

% return NaN vectors if no spikes
if to_return
    nans                        = NaN( size( fo ) );
    z                           = nans;
    phs                         = nans;
    fo                          = nans;
    fidx                        = nans;
    cohs                        = nans;
    ccohs                       = nans;
    return
end

%------------------------------------------------------------------------
% (1) preparations

% (1.1) x: downsample
if Fs ~= xFs
    x                           = resample( x, Fs, xFs );                   % downsample to Fs
    y                           = resample( y, Fs, xFs );                   % downsample to Fs
end
[ nt, ntrials ]                 = size( x );
T                               = nt / Fs;                                  % [s]
if ~isequal( [ nt, ntrials ], size( y ) )
    error( 'input size mismatch' )
end

% (1.2) x: remove mean from each column and reorganize as a vector
if dflagtrial
    x0                       	= bsxfun( @minus, x, mean( x, 1 ) );        % matrix; remove mean for each input segment separately
    y0                       	= bsxfun( @minus, y, mean( y, 1 ) );        % matrix; remove mean for each input segment separately
else
    x0                          = x;
    y0                          = y;
end
x0                              = reshape( x0, [ nt * ntrials 1 ] );        % vector
y0                              = reshape( y0, [ nt * ntrials 1 ] );        % vector
x0                              = x0 - mean( x0 );                          % remove global mean
y0                              = y0 - mean( y0 );                          % remove global mean

% (1.3) determine parameters for spectral analysis
nWindow                         = nFFT;                                     % window of same duration as nFFT
nOverlap                        = nWindow / 2;                              % overlap of 50%

fidx                            = fo >= fROI( 1 ) & fo <= fROI( 2 );        % bandpass (used for normalization and significance testing)

fwin                            = triang( nSpecSmooth );                    % window for smoothing normalization spectrum
fwin                            = fwin / sum( fwin );

%------------------------------------------------------------------------
% (2) compute coherence and impedance

% compute trial averaged spectra
yo                              = local_mtcsd1( [ x0 y0 ], nFFT, Fs, nWindow, nOverlap, mtNW, dflag );

% compute coherence, phase, and impedance
p1                              = yo( :, 1, 1 );                            % input spectrum; yo( :, 2, 1 )
p2                              = yo( :, 1, 2 );                            % output spectrum
c12                             = yo( :, 2, 2 );                            % cross specrtum
cohs                            = single( abs( c12 .^ 2 ) ./ ( p1 .* p2 ) );
phs                             = single( atan2( imag( c12 ), real( c12 ) ) );
z                               = sqrt( abs( p2 ) ) ./ sqrt( abs( p1 ) );   % impedance := |FFT( Y )| / |FFT( X )|

% correction by input spectrum withing the bandpass (fROI)
s1f                             = local_firfilt( sqrt( abs( p1 ) ), fwin );
normf                           = s1f( fidx ) ./ max( s1f( fidx ) );
ccohs                           = cohs;
ccohs( fidx, : )                = cohs( fidx, : ) ./ normf;

%------------------------------------------------------------------------
% (3) plot

if ~graphics
    return
end

fig                             = figure;

% coherence
subplot( 3, 1, 1 )
hold on
plot( fo, cohs, '.-b' )                                                     % actual coherence
xlim( fROI )
set( gca, 'tickdir', 'out', 'box', 'off' )
ylabel( 'Coherence' )
tstr                            = sprintf( '%0.2g-%0.2g Hz, %d trials, T=%d s', fROI( 1 ), fROI( 2 ), ntrials, ceil( T ) );
title( tstr )

% impedance
subplot( 3, 1, 2 )
plot( fo, z, '.-b' )
xlim( fROI )
set( gca, 'tickdir', 'out', 'box', 'off' )
ylabel( 'Impedance [M\Omega mm^2]' )
xlabel( 'Frequency [Hz]' )

% phase
subplot( 3, 1, 3 )
uwphs                           = unwrap( phs );
if max( uwphs ) > 2 * pi
    uwphs                       = uwphs - 2 * pi;
end
uwphs                           = mod( uwphs + pi, 2 * pi ) - pi;
plot( fo, uwphs, '.-b' )
xlim( fROI )
set( gca, 'tickdir', 'out', 'box', 'off' )
ylabel( 'Phase [rad]' )
line( xlim, [ 0 0 ] , 'color', [ 0 0 0 ], 'linestyle', '--' );
ylim( [ -1 1 ] * pi / 2 )

return % NR_calc_z_spectral

%------------------------------------------------------------------------
% Y = local_firfilt( x, W )
% zero-phase lag low-pass filtering of x's columns with the FIR W
%------------------------------------------------------------------------
function Y = local_firfilt( x, W )

C                               = length( W );
D                               = ceil( C / 2 ) - 1;
Y                               = filter( W, 1, [ flipud( x( 1 : C, : ) ); x; flipud( x( end - C + 1 : end, : ) ) ] );
Y                               = Y( 1 + C + D : end - C + D, : );

return % local_firfilt

%------------------------------------------------------------------------
% [ y, f ] = local_mtcsd1( x, nFFT, Fs, WinLength, nOverlap, NW, Detrend )
% cross-spectra between one signal (first column) and all others
%------------------------------------------------------------------------
function [ y, f ] = local_mtcsd1( x, nFFT, Fs, WinLength, nOverlap, NW, Detrend )

nTapers                         = 2 * NW - 1;
winstep                         = WinLength - nOverlap;
nChannels                       = size( x, 2 );
nSamples                        = size( x, 1 );

% check for column vector input
if nSamples == 1
	x                           = x';
	nSamples                    = size( x, 1 );
	nChannels                   = 1;
end

% calculate number of FFTChunks per channel
nFFTChunks                      = 1 + floor( ( ( nSamples - WinLength ) / winstep ) );

% allocate memory
y                               = complex( zeros( nFFT, 2, nChannels ) );   % output array
Periodogram                     = complex( zeros( nFFT, nTapers, nChannels ) ); % intermediate FFTs
Temp1                           = complex( zeros( nFFT, nTapers ) );
Temp2                           = complex( zeros( nFFT, nTapers ) );
Temp3                           = complex( zeros( nFFT, nTapers ) );
eJ                              = complex( zeros( nFFT, 1 ) );

% calculate Slepian sequences
Tapers                          = dpss( WinLength, NW, nTapers, 'calc' );

% compute tapered periodogram with FFT 
TaperingArray                   = repmat( Tapers, [ 1 1 nChannels ] );
for j                           = 1 : nFFTChunks
	Segment                     = x( ( j - 1 ) * winstep + ( 1 : WinLength ), : );
	if ~isempty( Detrend )
		Segment                 = detrend( Segment, Detrend );
    end
	SegmentsArray               = permute(repmat(Segment, [ 1 1 nTapers ] ), [ 1 3 2 ] );
	TaperedSegments             = TaperingArray .* SegmentsArray;

	Periodogram( :, :, : )      = fft( TaperedSegments, nFFT );

	% products 
    Temp1( :, : )               = squeeze( Periodogram( :, :, 1 ) );        % Ch1
    for Ch2                     = 1 : nChannels 
        
        Temp2( :, : )           = squeeze( Periodogram( :, :, Ch2 ) );
        Temp2                   = conj( Temp2 );
        Temp3( :, : )           = Temp1 .* Temp2;
        eJ( : )                 = sum( Temp3, 2 );
        y( :, 2, Ch2 )          = y( :, 2, Ch2 ) + eJ / ( nTapers * nFFTChunks ); % cross

        Temp3                   = squeeze( Periodogram( :, :, Ch2 ) ) .* Temp2;
        eJ                      = sum( Temp3, 2 );
        y( :, 1, Ch2 )          = y( :, 1, Ch2 ) + eJ / ( nTapers * nFFTChunks ); % auto
    end
    
end

% select subset of y, set up f array
if ~any( any( imag( x ) ) )    % x purely real
	if rem( nFFT, 2 )
		select                  = 1 : ( nFFT + 1 ) / 2;
	else
		select                  = 1 : nFFT / 2 + 1;
	end
	y                           = y( select, :, : );
else
	select                      = 1 : nFFT;
end
f                               = ( select - 1 )' * Fs / nFFT;

return % local_mtcsd1

% EOF