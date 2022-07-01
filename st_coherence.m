% st_coherence          coherence between spike trains and analog signal
%
% call                  cohs = st_coherence( st, x )
%                       [ ..., phs, fo, pvals, fidx, sigCoh, ccohs ] = ...
%
% gets                  st      [samples] spike times, sampled at spkFs. can be: 
%                                   -a single vector (one 'trial')
%                                   -a cell array of ntrials vectors
%                                   -a sparse matrix (trials in columns)
%                       x       continuous signal, sampled at Fs. can be: 
%                                   -a single vector 
%                                   -a matrix with ntrials columns
%
% thus, multiple segments are supported. different segments may use
%   different inputs (columns of x), but all must have same duration. 
%
% optional arguments (given as name/value pairs)
%
%       sampling parameters:
%
%                       spkFs       {20000}     [Hz], sampling rate of st
%                       xFs         {1250}      [Hz], sampling rate of x
%                       Fs          {1250}      [Hz], intervally-used rate
%                       fROI        {[0 20]}    [Hz], frequency ROI
%                       dflagtrial  {0}         logical; remove mean from every trial (column of st)
%
%       multi-taper parameters:  
%
%                       M           {1}         determinant of spectral resolution (1/M Hz)
%                       nFFT        {NaN}       alternative manner to determine spectral resolution
%                       mtNW        {3}         MT: coherence time-freq parameter
%                       dflag       {''}        MT: no detrending (may remove mean for each segment separately)
%                                                   argument to detrend (''/[], 'constant'/0, 'linear'/1, 2, ...)
%                       nSpecSmooth {5}         number of FFT frequency bins (used only for normalized coherence, ccohs)
%
%       significance testing parameters:
% 
%                       pTH         {0.001}     all methods:      p-value threhsold for coherence significance
%                       minBand     {5}         all methods:      [Hz], minimal width of band for TH-crossing of coherence
%                       fMaxMults   {2}         HF coherence:     [fMax multiples] lower frequency range for fHIGH
%                       nreps0      {0}         number of shuffling repetitions
%                       jmode       {0}         argument to svmix (0-shuffle; 1-spike time; 2-interval)
%                       jhwinS      {0.1}       [s], jittering half-width, in seconds (relevant only to jmode 1 and 2)
%
%       flow control:
%                       graphics    {0}         
%
% returns               cohs        coherence: vector of nFFT/2+1
%                       phs         phases
%                       fo          frequency bins
%                       pvals       3-column matrix, nfidx by 3
%                       fidx        logical vector with nfidx true bits
%                                       nfidx corresponds to the number of frequency bins within fROI
%                       sigCoh      logical vector, 1 x 3
%                                       indicates whether there is at least one continuous segment
%                                       (minBand Hz wide) of significant coherence
%                       ccohs       scaled coherence (scaled by the input spectrum within fROI)
%                     
% does
% (1) removes mean for each trial separately, then concatenates
% (2) computes coherence, phase, and scaled coherence
% (3) computes significance (three methods):
%       -shuffling (whole trial, spike-time jitter, or interval jitter)
%       -high-frequency coherence
%       -time reversal
%   significance based on random segments is not supported in this routine (see spikeWB for an implementation)
% (4) (optionally) plots
%
% notes
% (1) increasing M makes narrow effects narrower (so does reducing mtNW)
%       in st_on the other hand, if M is too large, some st_fingerprint bins 
%       (which does not frequency-smooth the input) may be empty
% (2) if all trials use nearly-identical conditions (i.e., x0 is the same
%       for all trials), then trial-specific mean removal (the default) is
%       appropriate. However, if every trial consists of e.g., different
%       frequencies, mean removal (or linear detrending) is inappropriate
%
% example
% Assume that st is a vector of spike times, measured in samples, and Iext 
% is the input current. Both simulated at dt = 0.0001 s.
% Then, spkFs = 1/0.0001 = 10000 Hz, and xFs is the same. To run the
% analysis, type: 
% >> fROI = [ 0 40 ]; Fs = 10000; 
% >> [ cohs, phs, fo ]               = st_coherence( st, Iext, 'spkFs', Fs, 'xFs', Fs, 'fROI', fROI, 'Fs', 1250 );
%
% calls                         ParseArgPairs
%
% see also                      st_fingerprint

% 25-mar-21 ES

% last update
% 01-jul-22

function [ cohs, phs, fo, pvals, fidx, sigCoh, ccohs ] = st_coherence( st, x, varargin )

% st: [ms]
% x: @ xFs

% sampling rates and ROI
spkFs_DFLT                      = 20000;                                    % [Hz]
xFs_DFLT                        = 1250;                                     % [Hz]
Fs_DFLT                         = 1250;                                     % [Hz]
fROI_DFLT                       = [ 0 20 ];                                 % [Hz]

% coherence analyses
dflagtrial_DFLT                 = 0;                                        % 1: remove mean of every trial; 0: do not
M_DFLT                      	= 1;                                        % determinant of spectral resolution (1/M Hz)
nFFT_DFLT                       = NaN;                                      % frequency steps will be Fs / nFFT
mtNW_DFLT                     	= 3;                                        % MT: coherence time-freq parameter
dflag_DFLT                   	= '';                                       % MT: no detrending (may remove mean for each segment separately, dflagtrial)

% for normalizing the coherence
nSpecSmooth_DFLT             	= 5;                                        % number of FFT frequency bins (only for normalized coherence)

% for significance testing
pTH_DFLT                     	= 0.001;                                    % all methods:      p-value threhsold for coherence significance
minBand_DFLT                   	= 5;                                        % all methods:      [Hz], minimal width of band for TH-crossing of coherence

fMaxMults_DFLT                 	= 2;                                        % HF coherence:     [fMax multiples] lower frequency range for fHIGH

nrepsShuf_DFLT                  = 0;                                        % shuffling:        number of shuffles
jmode_DFLT                      = 0;                                        % shuffling:        permutation by shuffling (1: spike time jitter; 2: interval jitter)
jhwinS_DFLT                     = 0.1;                                      % shuffling:        [s] jitter half-window (relevant only for jmode 1 or 2)


%------------------------------------------------------------------------
% arguments
nargs                           = nargin;
if nargs < 2 || isempty( st ) || isempty( x )
    to_return                   = 1;
else
    to_return                   = 0;
end
[ spkFs, xFs, Fs ...
    , xSpk, ntmax ...
    , fROI ...
    , nreps0, jmode, jhwinS ...
    , M, nFFT, mtNW, dflag, dflagtrial ...
    , nSpecSmooth, pTH, minBand, fMaxMults ...
    , graphics ]                = ParseArgPairs(...
    { 'spkFs', 'xFs', 'Fs' ...
    , 'xSpk', 'ntmax' ...
    , 'fROI' ...
    , 'nreps0', 'jmode', 'jhwinS' ...
    , 'M', 'nFFT', 'mtNW', 'dflag', 'dflagtrial' ...
    , 'nSpecSmooth', 'pTH', 'minBand', 'fMaxMults' ...
    , 'graphics' }...
    , { spkFs_DFLT, xFs_DFLT, Fs_DFLT ...
    , 0, [] ...
    , fROI_DFLT ...
    , nrepsShuf_DFLT, jmode_DFLT, jhwinS_DFLT ...
    , M_DFLT, nFFT_DFLT, mtNW_DFLT, dflag_DFLT, dflagtrial_DFLT ...
    , nSpecSmooth_DFLT, pTH_DFLT, minBand_DFLT, fMaxMults_DFLT ...
    , 0 } ...
    , varargin{ : } );
if round( xFs / Fs ) ~= ( xFs / Fs )
    error( 'downsampling from xFs to Fs must be by an interger ratio' )
end
if spkFs < Fs
    error( 'spkFs must be larger than (or equal to) Fs' )
end
jhwin                           = jhwinS * Fs;
if isnan( dflag )
    dflag                       = dflag_DFLT;
end

% determine frequencies for spectral analysis
if isnan( nFFT )
    nFFT                        = 2 ^ floor( log2( M * Fs ) );              % ~1/M Hz resolution 
end
fo                              = ( 0 : nFFT / 2 )' * Fs / nFFT;            % all frequencies up to Nyquist

% return NaN vectors if no spikes
if to_return
    nans                        = NaN( size( fo ) );
    cohs                        = nans;
    phs                         = nans;
    pvals                       = nans * ones( 1, 3 );
    fidx                        = nans;
    sigCoh                      = [];
    ccohs                       = [];
    return
end

%------------------------------------------------------------------------
% (1) preparations

if xSpk
    
    % (1.1) x: type handling (cell array of ntrials vectors (each double), or a single vector)
    ntrials                     = size( x, 2 );
    nt                          = ceil( Fs / xFs * ntmax );
    T                         	= nt / Fs;                                  % [s]
    if isa( x, 'cell' )
        if ~isvector( x )
            error( 'should be a vector of cell arrays' )
        end
        ntrials2                = length( x );
    else
        if ~isvector( x )
            error( 'every element should be a vector' )
        end
        ntrials2              	= 1;
        x                       = { x };
    end
    if ntrials2 ~= ntrials
        error( 'input size mismatch: st and x must contain same number of trials' )
    end
    
    % (1.2) x: convert x to samples @ Fs
    for i                       = 1 : ntrials
        if ~isvector( x{ i } )
            error( 'every element should be a vector' )
        end
        if ~issorted( x{ i } )
            error( 'every element should be a sorted vector of spike times' )
        end
        ns                      = length( x{ i } );
        if ns > 0
            xI                  = x{ i } / spkFs;                           % [samples @ spkFs] -> [s]
            x{ i }              = round( xI * Fs );                         % [samples @ spkFs] -> [samples @ Fs]
            cil                 = x{ i } == 0;
            x{ i }( cil )       = 1;
            if x{ i }( 1 ) < 1
                error( 'all spike times must be positive integers' )
            end
            if x{ i }( ns ) > nt
                error( 'all spike times must be smaller than %d', nt )
            end
        end
    end
    
    % (1.5) x: accumulate zero-mean spike trains for all trials
    x0mat                       = zeros( nt, ntrials );                     % matrix
    nspks                       = NaN( 1, ntrials );
    for i                       = 1 : ntrials
        fr                      = zeros( nt, 1 );
        xI                    	= x{ i };
        nspks( i )              = length( xI );
        fr( xI )                = 1;
        if dflagtrial
            x0mat( :, i )       = fr - nspks( i ) / nt;                     % remove mean for each output segment separately
        else
            x0mat( :, i )       = fr;
        end
    end
    x0                          = x0mat( : );                               % vector
    x0                          = x0 - mean( x0 );                          % remove global mean
    
    
else
    
    % (1.1) x: downsample
    if Fs ~= xFs
        x                       = resample( x, Fs, xFs );                   % downsample to Fs
    end
    [ nt, ntrials ]          	= size( x );
    T                         	= nt / Fs;                                  % [s]
    
    % (1.2) x: remove mean from each column and reorganize as a vector
    if dflagtrial
        x0                     	= bsxfun( @minus, x, mean( x, 1 ) );    % matrix; remove mean for each input segment separately
    else
        x0                     	= x;
    end
    x0                       	= reshape( x0, [ nt * ntrials 1 ] );        % vector
    x0                        	= x0 - mean( x0 );                          % remove global mean
end

% (1.3) st: type handling (cell array of ntrials vectors (each double), or a single vector)
if isa( st, 'cell' )
    if ~isvector( st )
        error( 'should be a vector of cell arrays' )
    end
    ntrials2                    = length( st );
else
    if ~isvector( st )
        error( 'every element should be a vector' )
    end
    ntrials2                    = 1;
    st                          = { st };
end
if ntrials2 ~= ntrials
    error( 'input size mismatch: st and x must contain same number of trials' )
end

% (1.4) st: convert st to samples @ Fs
for i                           = 1 : ntrials
    if ~isempty( st{ i } )
        if ~isvector( st{ i } )
            error( 'every element should be a vector' )
        end
        if ~issorted( st{ i } )
            error( 'every element should be a sorted vector of spike times' )
        end
    end
    ns                          = length( st{ i } );
    if ns > 0
        stI                     = st{ i } / spkFs;                          % [samples @ spkFs] -> [s]
        st{ i }                 = round( stI * Fs );                        % [samples @ spkFs] -> [samples @ Fs]
        cil                     = st{ i } == 0;
        st{ i }( cil )          = 1;
        if st{ i }( 1 ) < 1
            error( 'all spike times must be positive integers' )
        end
        if st{ i }( ns ) > nt 
            error( 'all spike times must be smaller than %d', nt )
        end
    end
end

% (1.5) st: accumulate zero-mean spike trains for all trials
f0mat                           = zeros( nt, ntrials );                     % matrix
nspks                           = NaN( 1, ntrials );
for i                           = 1 : ntrials
    fr                          = zeros( nt, 1 );
    stI                         = st{ i };
    nspks( i )                  = length( stI );
    fr( stI )                   = 1;
    if dflagtrial
        f0mat( :, i )           = fr - nspks( i ) / nt;                     % remove mean for each output segment separately
    else
        f0mat( :, i )           = fr;
    end
end
f0                              = f0mat( : );                               % vector
f0                              = f0 - mean( f0 );                          % remove global mean 

% (1.6) determine parameters for spectral analysis
nWindow                         = nFFT;                                     % window of same duration as nFFT
nOverlap                        = nWindow / 2;                              % overlap of 50%

fidx                            = fo >= fROI( 1 ) & fo <= fROI( 2 );        % bandpass (used for normalization and significance testing)
nfidx                           = sum( fidx );                              % number of frequency samples in fROI

fwin                            = triang( nSpecSmooth );                    % window for smoothing normalization spectrum
fwin                            = fwin / sum( fwin );

%------------------------------------------------------------------------
% (2) compute coherence 

% compute trial averaged spectra
yo                              = local_mtcsd1( [ x0 f0 ], nFFT, Fs, nWindow, nOverlap, mtNW, dflag );

% compute coherence and phase
p1                              = yo( :, 1, 1 );                            % input spectrum; yo( :, 2, 1 )
p2                              = yo( :, 1, 2 );                            % output spectrum
c12                             = yo( :, 2, 2 );                            % cross specrtum
cohs                            = single( abs( c12 .^ 2 ) ./ ( p1 .* p2 ) );
phs                             = single( atan2( imag( c12 ), real( c12 ) ) );

% correction by input spectrum withing the bandpass (fROI)
s1f                             = local_firfilt( sqrt( abs( p1 ) ), fwin );
normf                           = s1f( fidx ) ./ max( s1f( fidx ) );
ccohs                           = cohs;
ccohs( fidx, : )                = cohs( fidx, : ) ./ normf;


%------------------------------------------------------------------------
% (3) estimate significance (three methods)

% (3.0) determine parameters for significance estimation
fHIGH                           = [ fROI( 2 ) * fMaxMults Fs / 2 ];         % definition of 'high-frequency coherence'
hidx                            = fo >= fHIGH( 1 ) & fo <= fHIGH( 2 );
alfa                            = pTH / ( nFFT / 2 + 1 );
TH                              = abs( norminv( alfa / 2, 0, 1 ) );         % number of SDs corresponding to alfa (assuming normal distribution)
contHighBins                    = minBand / ( Fs / nFFT );

pvals                           = NaN( nfidx, 3 );
sigCoh                          = NaN( 1, 3 );

% (3.1) method #1 - spike time shuffling
%       -uniformly within the entire segment (default)
%       -each spike jhwin backwards or forwards (alternative)
%       -within an interval of jhwin samples (alternative)
if ismember( jmode, 0 : 2 )
    f0s                     	= zeros( nt, nreps0, ntrials );
    for i                   	= 1 : ntrials
        fr                      = zeros( nt, nreps0 );
        stI                     = st{ i };
        nspksI                  = length( stI );
        if nspksI == 0
            continue
        end
        switch jmode
            case 0
                h               = max( stI( 1 ) - 1, nt - stI( nspksI ) );
                stIs            = local_svmix( stI, h, nreps0, 0 );         % shuffle
            case { 1, 2 }
                stIs            = local_svmix( stI, jhwin, nreps0, jmode );	% jitter (spike/interval based)
                stIs            = round( stIs );
        end
        stIs( stIs > nt )       = NaN;
        ofst                    = ( ones( nspksI, 1 ) * ( 0 : ( nreps0 - 1 ) ) ) * nt;
        nnans                   = ~isnan( stIs );
        idx                     = stIs( nnans ) + ofst( nnans );
        fr( idx )               = 1;
        mfr                     = nanmean( fr, 1 );
        f0s( :, :, i )          = fr - mfr;                                 % remove mean for each output segment separately
    end
    f0s                        	= permute( f0s, [ 1 3 2 ] );                % nt x ntrials x nreps
    f0s                       	= reshape( f0s, [ nt * ntrials nreps0 ] );
    
    % compute coherence for the shuffled trains
    yo                      	= local_mtcsd2( x0 * ones( 1, nreps0 ), f0s, nFFT, Fs, nWindow, nOverlap, mtNW, dflag );
    p1s                      	= permute( yo( :, 1, : ), [ 1 3 2 ] );      % input spectrum; yo( :, 2, 1 )
    p2s                        	= permute( yo( :, 2, : ), [ 1 3 2 ] );      % output spectrum
    c12s                     	= permute( yo( :, 3, : ), [ 1 3 2 ] );      % cross specrtum
    
    cohs0                     	= single( abs( c12s .^ 2 ) ./ ( p1s .* p2s ) );
    
    % compute point-wise p-values (based on a Gaussian assumption/empirical)
    c                       	= cohs( fidx, : );
    chat                     	= cohs0( fidx, : );
    myu                     	= nanmean( chat, 2 );
    sig                      	= nanstd( chat, [], 2 );
    if nreps0 < 1000
        pvals0               	= 1 - normcdf( c, myu, sig );
    else
        pvals0                	= ( sum( repmat( c, [ 1 nreps0 ] ) <= chat, 2 ) + 1 ) / ( nreps0 + 1 );
    end
    pvals( :, 1 )               = pvals0;
    
    % determine whether there is at least one instance of a minBand of low p-values
    nchi                      	= diff( local_parse( find( pvals0 <= pTH ) ), [], 2 ) + 1;
    if max( nchi ) >= contHighBins
        sigCoh( 1 )            	= 1;
    else
        sigCoh( 1 )           	= 0;
    end
end

% (3.2) method #2 - high-frequency coherence (requires fROI with some margin from Nyquist)
if fHIGH( 1 ) <= fHIGH( 2 )
    
    % point-wise p-values (Gaussian assumption)
    c                           = cohs( fidx, : );
    ch                          = nanmean( cohs( hidx, : ), 2 );
    mu                          = mean( ch );                               % mean high-frequency coherence
    sg                          = std( ch );
    pvals( :, 2 )               = 1 - normcdf( c, mu, sg );
    
    % determine minBand of high coherence
    sidx                        = c > ( TH * mu );
    nchi                      	= diff( local_parse( find( sidx ) ), [], 2 ) + 1;
    if max( nchi ) >= contHighBins
        sigCoh( 2 )            	= 1;
    else
        sigCoh( 2 )           	= 0;
    end
    
end

% (3.3) method #3 - time reversal (always possible)
% time reverse each trial separately
f1                              = flipud( f0mat );

% compute trial averaged spectra
yo                              = local_mtcsd1( [ x0 f1( : ) ], nFFT, Fs, nWindow, nOverlap, mtNW, dflag );

% compute coherence and phase
p1                              = yo( :, 1, 1 );                            % input spectrum; yo( :, 2, 1 )
p2                              = yo( :, 1, 2 );                            % output spectrum
c12                             = yo( :, 2, 2 );                            % cross specrtum
cohs1                           = single( abs( c12 .^ 2 ) ./ ( p1 .* p2 ) );

% compute point-wise p-values (based on a Gaussian assumption)
ch                              = nanmean( cohs1, 2 );
mu                              = mean( ch );                               % mean coherence (all frequencies) for time-reversed spike trains
sg                              = std( ch );
pvals( :, 3 )                   = 1 - normcdf( c, mu, sg );

% determine minBand of high coherence
nchi                            = diff( local_parse( find( pvals( :, 3 ) <= pTH ) ), [], 2 ) + 1;
if max( nchi ) >= contHighBins
    sigCoh( 3 )                 = 1;
else
    sigCoh( 3 )                 = 0;
end

%------------------------------------------------------------------------
% (4) plot

if ~graphics
    return
end

figure

% coherence
subplot( 3, 1, 1 )
hold on
if ~exist( 'sidx', 'var' )
    sidx                        = false( nfidx, 1 );
end
plot( fo, cohs, '-b', fo( ~sidx ), cohs( ~sidx ), '.k', fo( sidx ), cohs( sidx ), '.b' ) % actual coherence
if ismember( jmode, 0 : 2 )
    plot( fo( fidx ), myu, 'k', fo( fidx ), myu + sig, '--k', fo( fidx ), myu - sig, '--k' ) % shuffled coherence
end
xlim( fROI )
set( gca, 'tickdir', 'out', 'box', 'off' )
ylabel( 'Coherence' )
title( sprintf( '%0.2g-%0.2g Hz, %d trials, T=%d s, %d spikes', fROI( 1 ), fROI( 2 ), ntrials, ceil( T ), sum( nspks ) ) )

% phase
subplot( 3, 1, 2 )
uwphs                           = unwrap( phs );
uwphs                           = mod( uwphs + pi, 2 * pi ) - pi;
plot( fo, uwphs, '-b', fo( sidx ), uwphs( sidx ), '.b' )
xlim( fROI )
set( gca, 'tickdir', 'out', 'box', 'off' )
ylabel( 'Phase [rad]' )
line( xlim, [ 0 0 ] , 'color', [ 0 0 0 ], 'linestyle', '--' );
ylim( [ -pi pi ] )

% p-values
if ismember( jmode, 0 : 2 ) && nreps0 > 1
    pvals0                      = pvals( :, 1 );                            % shuffled p-values
    tmethod                     = 'shuffling';
elseif fHIGH( 1 ) <= fHIGH( 2 )
    pvals0                      = pvals( :, 2 );                            % high-frequency p-values
    tmethod                     = 'high frequency coherence';
else
    pvals0                      = pvals( :, 3 );                            % time-reversed p-values
    tmethod                     = 'time reversal';
end
subplot( 3, 1, 3 )
pvals0( pvals0 <= 0 )           = sqrt( eps );
line( xlim, -log10( pTH ) * [ 1 1 ], 'linestyle', '--', 'color', [ 0 0 0 ] );
plot( fo( fidx ), -log10( pvals0 ), '.-b' )
set( gca, 'tickdir', 'out', 'box', 'off' )
title( sprintf( 'Significance by %s', tmethod ) )
ylabel( '-log_{10}( pval )' )
xlabel( 'Frequency [Hz]' )

return % st_coherence

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

%------------------------------------------------------------------------
% xm = local_svmix( x, h, p, j )
% mix/jitter spikes in a vector
%------------------------------------------------------------------------
function xm = local_svmix( x, h, p, j )

sz                              = size( x );
x                               = x( : );
n                               = length( x );
rng( round( rand( 1 ) * sum( 100 * clock ) ), 'v4' )

switch j
    
    case 0                                                                  % shuffle (re-distribute over duration)
        sv                      = max( x( 1 ) - h, 1 );
        ev                      = x( n ) + h;
        dur                     = ev - sv + 1;
        xm                      = floor( dur * rand( n, p ) ) + sv;
    case 1                                                                  % spike-time jitter: rectangular window, 2 * h + 1
        r                       = round( 2 * h * rand( n, p ) - h );
        xm                      = x + r;
    case 2                                                                  % interval jitter: rectangular window
        w                       = 2 * h + 1;
        wn                      = ceil( x / w );
        r                       = ceil( rand( n, p ) * w );
        xm                      = ( wn * ones( 1, p ) - 1 ) * w + r;
end

xm                              = sort( xm );
if ~isequal( sz, [ sz( 1 ) 1 ] )
    xm                          = reshape( xm, [ sz p ] );
end

return

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
    Temp1                       = squeeze( Periodogram( :, :, 1 ) );        % Ch1
    for Ch2                     = 1 : nChannels 
        
        Temp2                   = squeeze( Periodogram( :, :, Ch2 ) );
        Temp2                   = conj( Temp2 );
        Temp3                   = Temp1 .* Temp2;
        eJ                      = sum( Temp3, 2 );
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

%------------------------------------------------------------------------
% [ y, f ] = local_mtcsd2( x1, x2, nFFT, Fs, WinLength, nOverlap, NW, Detrend )
% cross-spectra between one signal (first matrix) and another (second matrix)
%------------------------------------------------------------------------
function [ y, f ] = local_mtcsd2( x1, x2, nFFT, Fs, WinLength, nOverlap, NW, Detrend )

nTapers                         = 2 * NW - 1;
winstep                         = WinLength - nOverlap;
nColumns                        = size( x1, 2 );                            % actually, ntrials (there are 2 channels, x1 and x2)
nSamples                        = size( x1, 1 );

% check for column vector input
if nSamples == 1
	x1                          = x1';
    x2                          = x2';
	nSamples                    = size( x1, 1 );
	nColumns                    = 1;
end

% calculate number of FFTChunks per channel
nFFTChunks                      = 1 + floor( ( ( nSamples - WinLength ) / winstep ) );

% allocate memory
y                               = complex( zeros( nFFT, 3, nColumns ) ); % output array
Periodogram1                    = complex( zeros( nFFT, nTapers, nColumns ) ); % intermediate FFTs
Periodogram2                    = complex( zeros( nFFT, nTapers, nColumns ) );
Temp1                           = complex( zeros( nFFT, nTapers ) );
Temp2                           = complex( zeros( nFFT, nTapers ) );
Temp3                           = complex( zeros( nFFT, nTapers ) );
eJ                              = complex( zeros( nFFT, 1 ) );

% calculate Slepian sequences
Tapers                          = dpss( WinLength, NW, nTapers, 'calc' );

% compute tapered periodogram with FFT 
TaperingArray                   = repmat( Tapers, [ 1 1 nColumns ] );
for j                           = 1 : nFFTChunks
	Segment1                    = x1( ( j - 1 ) * winstep + ( 1 : WinLength ), : );
	Segment2                    = x2( ( j - 1 ) * winstep + ( 1 : WinLength ), : );
	if ~isempty( Detrend ) 
		Segment1                = detrend( Segment1, Detrend );
		Segment2                = detrend( Segment2, Detrend );
    end
	SegmentsArray1              = permute( repmat( Segment1, [ 1 1 nTapers ] ), [ 1 3 2 ] );
	SegmentsArray2              = permute( repmat( Segment2, [ 1 1 nTapers ] ), [ 1 3 2 ] );
	TaperedSegments1            = TaperingArray .* SegmentsArray1;
	TaperedSegments2            = TaperingArray .* SegmentsArray2;

	Periodogram1                = fft( TaperedSegments1, nFFT );
	Periodogram2                = fft( TaperedSegments2, nFFT );

	% products
    for Col                     = 1 : nColumns
        Temp1                   = squeeze( Periodogram1( :, :, Col ) );
        Temp2                   = squeeze( Periodogram2( :, :, Col ) );
        Temp2                   = conj( Temp2 );
        Temp3                   = Temp1 .* Temp2;
        eJ                      = sum( Temp3, 2 );
        y( :, 3, Col )          = y( :, 3, Col ) + eJ / ( nTapers * nFFTChunks ); % cross
        
        Temp2                   = conj( Temp1 );
        Temp3                   = Temp1 .* Temp2;
        eJ                      = sum(Temp3, 2);
        y( :, 1, Col )          = y( :, 1, Col ) + eJ / ( nTapers * nFFTChunks ); % auto1
        
        Temp1                   = squeeze( Periodogram2( :, :, Col ) );
        Temp2                   = conj( Temp1 );
        Temp3                   = Temp1 .* Temp2;
        eJ                      = sum( Temp3, 2 );
        y( :, 2, Col )          = y( :, 2, Col ) + eJ / ( nTapers * nFFTChunks ); % auto2
    end
end

% select subset of y, set up f array
if ~any( any( imag( x1 ) ) ) && ~any( any( imag( x2 ) ) )    % x purely real
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

return % local_mtcsd2

% EOF
