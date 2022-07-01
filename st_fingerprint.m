% st_fingerprint        firing rate map (frequency-phase) for spike trains and analog signal
%
% call                  [ frate_vec, phase_vec, rate_map, fbins, pbins ] = st_fingerprint( st, x )
%                       [ ..., count_vec, s ] = st_fingerprint( st, x )
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
%
%       frequency-phase parameters:  
%
%                       M           {1}         determinant of spectral resolution (1/M Hz)
%                       nFFT        {NaN}       alternative manner to determine spectral resolution
%                       frqs        {NaN}       NaN of edges [Hz] of binning grid
%                       phases      {20}        number or edges [rad] of binning grid
%                       nSpecSmooth {5}         number of FFT frequency bins (used only for normalized coherence, ccohs)
%
%       convolution rate parameters:
% 
%                       sdGauss     {10}        [ms], SD of Gaussian FIR (used for frate_mean_vec, frate_std_vec, 
%                                                   frate_range_vec, and frate_bounds_vec metrics only)
%                       alfa        {0.05}      1-confidence interval (used for frate_bounds_vec metric)
%
%       flow control:
%                       graphics    {0}         
%
% returns               frate_vec   [spks/s], 1D: firing rate per frequency (spike based, no smoothing)
%                       phase_vec	[spks/s], 1D: phase per frequency (spike based, no smoothing)
%                       rate_map    [spks/s], 2D: firing rate per frequency and phase (spike based, no smoothing)
%                       fbins       [Hz] frequency grid (used for all metrics)
%                       pbins       [rad] phase grid (used only for the 2D metric)
%                                       nfidx corresponds to the number of frequency bins within fROI
%                       count_vec   [spks/cycle], 1D: number of spikes per cycle (of binned frequency)
%                       s           strcture with other metrics, including: 
%   
%                       pfhist_mat          [spks], 3D; nominator of the rate_mat
%                       time_mat            [s], 3D; denominator of the rate_mat
%                       rate_mat            [spks/s], 3D; a rate_map for every trial
%                       tfreqs              [s], 1D, time spent in every binned frequency 
%                       frate_mean_vec      [spks/s], 1D, mean firing rate per binned frequency (same as frate_vec, yet evaluated using convolution)
%                       frate_std_vec       [spks/s], 1D, SD firing rate per binned frequency (evaluated using convolution)
%                       frate_range_vec     [spks/s], 1D, firing rate range per binned frequency (evaluated using convolution)
%                       frate_bounds_vec    [spks/s], 1D, CL of firing rate per binned frequency (evaluated using convolution)
%
% does
% (1) divide data into individual cycles according to local minima
% (2) assign an instantaneous frequency to each cycle
% (3) assign a frequency and phase to every spike
% (4) bin in 2D (frequency-phase) and compute rate
% (5) bin in 1D (frequency) and compute rate and mean phase
% (6) compute convolution-based (time dependent) firing rate, then compute statistics (mean, SD, range, CL) for every frequency bin
%
% calls                 ParseArgPairs
%                       plotOneSpectrogram
%                       inranges
%
% see also              st_coherence

% 05-apr-21 ES

% last update
% 01-jul-22

function [ frate_vec, phase_vec, rate_map, fbins, pbins, count_vec, s ] = st_fingerprint( st, x, varargin )

% sampling rates and ROI
spkFs_DFLT                      = 20000;                                    % [Hz]
xFs_DFLT                        = 1250;                                     % [Hz]
Fs_DFLT                         = 1250;                                     % [Hz]
fROI_DFLT                       = [ 0 20 ];                                 % [Hz]

% spectral resolution
M_DFLT                      	= 1;                                        % determinant of spectral resolution (1/M Hz)
nFFT_DFLT                       = NaN;                                      % frequency steps will be Fs / nFFT
frqs_DFLT                       = NaN;                                      % frequency grid
phases_DFLT                     = 20;                                       % determinant of the phase resolution (2*pi/phases rad)

% temporal smoothing
sdGauss_DFLT                    = 10;                                       % [ms]; if NaN, will be the reciprocal of max( fROI )
alfa_DFLT                       = 0.05;                                     % 1-confidence bounds for convolution-based firing rate

%------------------------------------------------------------------------
% arguments
nargs                           = nargin;
if nargs < 2 || isempty( st ) || isempty( x )
    to_return                   = 1;
else
    to_return                   = 0;
end

[ spkFs, xFs, Fs, fROI ...
    , M, nFFT, frqs, phases ...
    , sdGauss, alfa ...
    , graphics ]                = ParseArgPairs(...
    { 'spkFs', 'xFs', 'Fs', 'fROI' ...
    , 'M', 'nFFT', 'frqs', 'phases' ...
    , 'sdGauss', 'alfa' ...
    , 'graphics' }...
    , { spkFs_DFLT, xFs_DFLT, Fs_DFLT, fROI_DFLT ...
    , M_DFLT, nFFT_DFLT, frqs_DFLT, phases_DFLT ...
    , sdGauss_DFLT, alfa_DFLT ...
    , 0 } ...
    , varargin{ : } );

% create frequency binnning grid
if ~isnan( frqs ) && isvector( frqs ) && all( diff( frqs ) > 0 )
    fedges                      = frqs( : );
else
    if isnan( nFFT )
        nFFT                    = 2 ^ floor( log2( M * Fs ) );              % ~1/M Hz resolution
    end
    df                          = Fs / nFFT;
    fo                          = ( 0 : nFFT / 2 )' * Fs / nFFT;            % all frequencies up to Nyquist
    fidx                        = fo >= fROI( 1 ) & fo <= fROI( 2 );        % bandpass (used for normalization and significance testing)
    fbins                       = fo( fidx );
    fedges                      = [ fbins - df / 2; fbins( end ) + df / 2 ];
end
nfBins                          = length( fedges ) - 1;
fbins                           = ( fedges( 1 : nfBins ) + fedges( 2 : nfBins + 1 ) ) / 2;

% create a phase binning grid
if length( phases ) == 1
    npBins                      = abs( round( phases ) );
    phsBinSize                  = 2 * pi / npBins;
    pedges                      = ( ( -pi + phsBinSize / 2 ) : phsBinSize : ( pi + phsBinSize / 2 ) )';
else
    pedges                      = phases;
end
if length( unique( mod( pedges( [ 1 end ] ), 2 * pi ) ) ) ~= 1
    error( 'Phase edges should cover the entire unit circle' )
end
pedges                          = pedges( : );
npBins                          = length( pedges ) - 1;
pbins                           = ( pedges( 1 : npBins ) + pedges( 2 : npBins + 1 ) ) / 2;

% return NaN vectors if no spikes
if to_return
    nans                        = NaN( size( fbins ) );
    frate_vec                   = nans;
    phase_vec                 	= nans;
    rate_map                 	= NaN( npBins, nfBins );
    count_vec                	= nans;
    s                           = [];
    return
end

%------------------------------------------------------------------------
% (1) preparations

% (1.1) x: downsample
if Fs ~= xFs
    x                           = resample( x, Fs, xFs );                   % downsample to Fs
end
[ nt, ntrials ]                 = size( x );
T                               = nt / Fs;                                  % [s]

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
nspks                           = NaN( ntrials, 1 );
for i                           = 1 : ntrials
    if ~isempty( st{ i } )
        if ~isvector( st{ i } )
            error( 'every element should be a vector' )
        end
        if ~issorted( st{ i } )
            error( 'every element should be a sorted vector of spike times' )
        end
    end
    nspks( i )                  = length( st{ i } );
    if nspks( i ) > 0
        stI                     = st{ i } / spkFs;                          % [samples @ spkFs] -> [s]
        st{ i }                 = round( stI * Fs );                        % [samples @ spkFs] -> [samples @ Fs]
        cil                     = st{ i } == 0;
        st{ i }( cil )          = 1;
        if st{ i }( 1 ) < 1
            error( 'all spike times must be positive integers' )
        end
        if st{ i }( nspks( i ) ) > nt 
            error( 'all spike times must be smaller than %d', nt )
        end
    end
end
ntot                            = sum( nspks );

% (1.5) prepare smoother
if isnan( sdGauss )
    sdGauss                     = 1000 / max( fROI );                       % [ms]
end
sigmaX                          = 2 * floor( sdGauss / 1000 * Fs / 2 ) + 1;
K                               = local_gausskernel( sigmaX, 6 * sigmaX + 1 ); 
gwin                            = K / sum( K );

%------------------------------------------------------------------------
% (2) computations

% 2D bin edges
bins2d                          = { pedges fedges };

% allocate memory for cycle-based statistics
cyc_cells                       = cell( 1, ntrials );
ncycs                           = NaN( 1, ntrials );

% allocate memory for spike-based statistics
spk_periods                     = NaN( ntot, 2 );                           % for every spike, period of the corresponding cycle
spk_freqs                       = NaN( ntot, 1 );                           % for every spike, frequency of the corresponding cycle
spk_phs                         = NaN( ntot, 1 );                           % for every spike, phase in the corresponding cycle
cs                              = cumsum( nspks );
spk_blks                        = [ [ 1; ( cs( 1 : ntrials - 1 ) + 1 ) ] cs ];

% allocate memory for summary statistics
frate_mat                       = NaN( nfBins, ntrials );                   % firing rate in the bandpass
count_mat                       = NaN( nfBins, ntrials );                   % number of spikes per cycle in the bandpass
n_mat                           = NaN( nfBins, ntrials );                   % number of spikes
phase_mat                       = NaN( nfBins, ntrials );                   % [rad] mean phase in the bandpass
phase_lock                      = NaN( nfBins, ntrials );                   % phase lock in the bandpass
phase_std                       = NaN( nfBins, ntrials );                   % [rad] SD of phase in the bandpass

frate_mean                      = NaN( nfBins, ntrials );                   % firing rate in the bandpass, estimated by convolution
frate_std                       = NaN( nfBins, ntrials );                   % SD of (convolution-based) firing rate in the bandpass
frate_min                       = NaN( nfBins, ntrials );                   % min (convolution-based) firing rate in the bandpass
frate_max                       = NaN( nfBins, ntrials );                   % max (convolution-based) firing rate in the bandpass
frate_lo                        = NaN( nfBins, ntrials );                   % low confidence limit of (convolution-based) firing rate in the bandpass
frate_hi                        = NaN( nfBins, ntrials );                   % high confidence limit of (convolution-based) firing rate in the bandpass

pfhist_mat                      = NaN( npBins, nfBins, ntrials );           % [spikes/bin]
time_mat                        = NaN( npBins, nfBins, ntrials );           % [s]

for i                           = 1 : ntrials

    xi                          = x( :, i );
    
    % (1) assume sinusoid, divide into cycles according to local minima
    lmin                        = [ 1; local_find_local_minima( xi ); nt ];
    ncycles                     = length( lmin ) - 1;
    periods                     = [ lmin( 1 : ncycles ) lmin( 2 : ncycles + 1 ) - 1 ]; % period: from minimum to just before the next minimum
    periods( ncycles, 2 )       = nt;                                       % [samples]; alternatively, may decide to remove the last (and possibly the first) cycle
    
    % (2) assign an instantaneous frequency to each cycle
    cyc_T                       = ( diff( periods, [], 2 ) + 1 ) / Fs;      % [s]; duration of each and every cycle
    cyc_freq                    = 1 ./ cyc_T;                               % [Hz]; frequency of each and every cycle
    [ ncyc_freq, cyc_bin ]      = histc( cyc_freq, fedges );                % number of cycles at each frequency bin (some cycles may be excluded since not in fedges: sum( ncyc_freq ), length( cyc_freq ) )
    ncyc_freq( nfBins )         = ncyc_freq( nfBins ) + ncyc_freq( nfBins + 1 );
    ncyc_freq( nfBins + 1 )     = [];
    
    % accumulate number of cycles per frequency over trials
    cyc_cells{ i }              = cyc_freq;                                 % [Hz]; this list includes the cycles that are out of fROI
    ncycs( i )                  = length( cyc_freq );
    
    % (3) determine the period, frequency, and phase of each spike
    % each cycle is from (trough to trough), so phase is linearly interpolated from (-pi to pi)
    % THIS IS NOT good for the first cycle of a chirp that starts at DC..
    sti                         = st{ i };                                  % [samples]
    [ ~, spk_seg ]              = inranges( sti, periods );                 % assign a segment to every spike
    spk_f                       = cyc_freq( spk_seg );                      % [Hz], assign a frequency to every spike
    spk_per                    	= periods( spk_seg, : );                    % [samples], onset and offset of the cycle of every spike (-pi to pi)
    dsti                        = sti - spk_per( :, 1 ) + 1;                % [samples], time from onset to spike
    spd                         = diff( spk_per, [], 2 ) + 1;               % [samples], period duration
    spk_phi                     = ( dsti ./ spd - 0.5 ) * 2 * pi;           % [rad], 0: peak, -pi: onset trough; pi: offset trough (this assumes symmetric cycles with uniform phase distribution)
    sidx                        = spk_phi < pedges( 1 );                    % wrap on the circle to ensure all within pedges
    spk_phi( sidx )             = spk_phi( sidx ) + 2 * pi;
    
    % (4) bin (time and spike count) in 2D
    % -this will yield the spikes/s (per frquency and phase bin) metric
    % bin time
    timf_i                      = zeros( nfBins, 1 );                       % [s]; time per binned frequency
    for j = 1 : nfBins
        timf_i( j )             = sum( cyc_T( cyc_bin == j ) );
    end
    denom_i                     = ones( npBins, 1 ) * timf_i';
    time_mat( :, :, i )         = denom_i;
    % bin spikes    
    data                        = [ spk_phi spk_f ];
    if isempty( data )
        continue
    end
    pfhist_i                    = local_bindata( data, bins2d );            % counts/bin
    pfhist_mat( :, :, i )       = pfhist_i;                                 % should all be the same: sum( histc( spk_f, fedges ) ), sum( histc( spk_phi, pedges ) ), sum( pfhist_i( : ) )
    
    % accumulate spike statistics over trials
    sidx                        = spk_blks( i, 1 ) : spk_blks( i, 2 );
    spk_periods( sidx, : )      = spk_per;
    spk_freqs( sidx, : )        = spk_f;
    spk_phs( sidx, : )          = spk_phi;

    % (4) count the number of spikes in each cycle
    [ nst, cycnum ]             = local_uhist( spk_seg );
    scount                      = zeros( ncycles, 1 );
    scount( cycnum )            = nst;                                      % should be the same: sum( nst ), sum( histc( spk_f, fedges ) )
    
    % (5) compute the firing rate in each cycle
    cyc_frate                   =  scount ./ cyc_T; 
    
    % (6) for all cycles in a given fbin, compute a weighted firing rate
    % this yields the spikes/cycle and the spikes/s (per frequency bin) metrics
    for j                       = 1 : nfBins
        idx                     = cyc_bin == j;
        v1                      = sum( cyc_frate( idx ) .* cyc_T( idx ) ) / sum( cyc_T( idx ) );
        v2                      = sum( scount( idx ) .* cyc_T( idx ) ) / sum( cyc_T( idx ) );
        frate_mat( j, i )       = v1;
        count_mat( j, i )       = v2;
        n_mat( j, i )           = sum( scount( idx ) );
    end
    
    % (7) for all spikes in a given fbin, compute a mean phase
    % -this yields the phase (per frequency bin) metric
    for j                       = 1 : nfBins
        sidx                    = ismember( spk_seg, find( cyc_bin == j ) );
        if sum( sidx ) > 0
            [ pp, rr, ss ]      = local_circ_mean( spk_phi( sidx ) );
            phase_mat( j, i )   = pp;
            phase_lock( j, i )  = rr;
            phase_std( j, i )   = ss;
        end
    end
    
    % (8) compute firing rate by convolution 
    fr                          = zeros( nt, 1 );
    for j                       = 1 : length( sti )
        fr( sti( j ) )       	= fr( sti( j ) ) + 1;
    end
    frs                         = local_firfilt( fr, gwin ) * Fs;
    for j                       = 1 : nfBins
        idx                     = cyc_bin == j;
        if sum( idx ) == 0
            continue
        end
        frs_j                   = local_getdatainranges( frs, periods( idx, : ) );
        frate_mean( j, i )      = mean( frs_j );
        frate_std( j, i )       = std( frs_j );
        frate_min( j, i )       = min( frs_j );
        frate_max( j, i )       = max( frs_j );
        lohi                    = local_bounds( frs_j, alfa );
        frate_lo( j, i )        = lohi( 1 );
        frate_hi( j, i )        = lohi( 2 );
    end
   
end

% organize accumulated cycles
totcycs                         = sum( ncycs );
cyc_freqs                       = NaN( totcycs, 1 );
cs                              = cumsum( ncycs( : ) );
cyc_blks                        = [ [ 1; ( cs( 1 : ntrials - 1 ) + 1 ) ] cs ];
for i                           = 1 : ntrials
    cidx                        = cyc_blks( i, 1 ) : cyc_blks( i, 2 );
    cyc_freqs( cidx, : )        = cyc_cells{ i };
end
if sum( Fs ./ cyc_freqs ) ~= nt * ntrials
    error( 'problem' )
end

% time spent in each frequency bin
[ h, cyc_bins ]                 = histc( cyc_freqs, fedges );
h( nfBins )                     = h( nfBins ) + h( nfBins + 1 );
h( nfBins + 1 )                 = [];
tfreqs                          = NaN( nfBins, 1 );
for i                           = 1 : nfBins
    idx                         = cyc_bins == i; 
    tfreqs( i )                 = sum( 1 ./ cyc_freqs( idx ) );             % [s], actual time spent in within-bin cycles
end

% average 1D metrics (spikes/cycle, spikes/s, and phase) over trials (or frequencies)
frate_vec                       = nanmean( frate_mat, 2 );
count_vec                       = nanmean( count_mat, 2 );
phase_vec                       = local_circ_mean( phase_mat', n_mat' )';
phase_vec( phase_vec > pi )     = phase_vec( phase_vec > pi ) - 2 * pi;     % warp back to -pi to pi
nans                            = sum( n_mat, 2 ) == 0;
phase_vec( nans )               = NaN;                                      % assign NaNs when no spikes
if ntrials == 1
    phase_vec_lock              = phase_lock;
    phase_vec_std               = phase_std;
else
    phase_vec_lock              = nansum( bsxfun( @rdivide, n_mat, nansum( n_mat, 2 ) ) .* phase_lock, 2 );
    phase_vec_std               = nansum( bsxfun( @rdivide, n_mat, nansum( n_mat, 2 ) ) .* phase_std, 2 );
end

frate_mean_vec                  = nanmean( frate_mean, 2 );
frate_std_vec                   = nanmean( frate_std, 2 );
frate_range_vec                 = nanmean( frate_max - frate_min, 2 );
frate_bounds_vec               	= nanmean( frate_hi - frate_lo, 2 );

% compute 2D metrics (phase-frequency rates): per trial and overall
rate_mat                        = pfhist_mat ./ time_mat;                   % may /0 if time_mat is null
rate_map                        = nansum( pfhist_mat, 3 ) ./ nansum( time_mat, 3 );
rate_mat                        = rate_mat * npBins;                        % [spikes/s]
rate_map                        = rate_map * npBins;                        % [spikes/s]

%------------------------------------------------------------------------
% (3) accumulate all outputs
s.frate_vec                     = frate_vec;
s.phase_vec                     = phase_vec;
s.rate_map                      = rate_map;
s.fbins                         = fbins;
s.pbins                         = pbins;
s.count_vec                     = count_vec;

s.fedges                        = fedges;
s.pedges                        = pedges;
s.pfhist_mat                    = pfhist_mat;
s.time_mat                      = time_mat;
s.rate_mat                      = rate_mat;
s.tfreqs                        = tfreqs;

s.phase_vec_lock                = phase_vec_lock;
s.phase_vec_std                 = phase_vec_std;

s.frate_mean_vec                = frate_mean_vec;
s.frate_std_vec                 = frate_std_vec;
s.frate_range_vec               = frate_range_vec;
s.frate_bounds_vec              = frate_bounds_vec;

s.frate_min                     = frate_min;
s.frate_max                     = frate_max;
s.frate_lo                      = frate_lo;
s.frate_hi                      = frate_hi;

s.spk_periods                   = spk_periods;
s.spk_freqs                     = spk_freqs;
s.spk_phs                       = spk_phs;

%------------------------------------------------------------------------
% (4) plot

if ~graphics
   return
end

figure;

subplot( 3, 2, 1 )
ph                              = plot( fbins, frate_mean_vec, 'o--k', fbins, frate_vec, '.-b', fbins, nanmean( rate_map ), ':r' );
set( ph( 1 ), 'MarkerSize', 8 )
lh                              = legend( ph, { sprintf( 'Conv. rate (SD=%0.2g ms)', sdGauss ), 'Spike count rate', '2D marginals' }, 'Location', 'best' ); 
set( lh, 'box', 'off' )
xlim( fROI )
ylim( [ 0 max( ylim ) ] )
set( gca, 'tickdir', 'out', 'box', 'off' )
ylabel( 'Firing rate [spks/s]' )
title( sprintf( '%0.2g-%0.2g Hz, %d reps, T=%d s', fROI( 1 ), fROI( 2 ), ntrials, ceil( T ) ) )

uwphs                           = unwrap( phase_vec );
uwphs                           = mod( uwphs + pi, 2 * pi ) - pi;
subplot( 3, 2, 3 )
ph                              = plot( fbins, uwphs, '.-b', fbins, uwphs + phase_vec_std, 'o--b', fbins, uwphs - phase_vec_std, 'o--b' );
set( ph( 2 : 3 ), 'MarkerSize', 6, 'color', [ 0.5 0.5 0.5 ] )
xlim( fROI )
set( gca, 'tickdir', 'out', 'box', 'off' )
ylabel( 'Phase [rad]' )
line( xlim, [ 0 0 ], 'color', [ 0 0 0 ], 'linestyle', '--' )
ylim( [ -pi pi ] )

subplot( 3, 2, 2 )
plot( fbins, frate_std_vec, '.-k' )
xlim( fROI )
ylim( [ 0 max( ylim ) ] )
set( gca, 'tickdir', 'out', 'box', 'off' )
ylabel( 'Firing rate SD [spks/s]' )

subplot( 3, 2, 4 )
ph                              = plot( fbins, frate_range_vec, '--k', fbins, frate_bounds_vec, '.-k' );
lh                              = legend( ph, { 'Range', sprintf( '%d%% CL', round( 100 - 100 * alfa ) ) }, 'Location', 'best' ); 
set( lh, 'box', 'off' )
xlim( fROI )
ylim( [ 0 max( ylim ) ] )
set( gca, 'tickdir', 'out', 'box', 'off' )
ylabel( 'Range [spks/s]' )
xlabel( 'Frequency [Hz]' )

subplot( 3, 2, 5 )
plot( fbins, count_vec, '.-r' )
xlim( fROI )
ylim( [ 0 max( ylim ) ] )
set( gca, 'tickdir', 'out', 'box', 'off' )
ylabel( 'Spikes/cycle' )
xlabel( 'Frequency [Hz]' )

subplot( 3, 2, 6 )
prate_map                       = rate_map;
nans                            = isnan( prate_map );
prate_map( nans )               = 0;
if all( prate_map( : ) == 0 )
    axis off
else
    plotOneSpectrogram( pbins, fbins, prate_map' );
    axis square
    title( sprintf( '%2.1f spks/s', max( rate_map( : ) ) ) )
end

return % st_fingerprint

%------------------------------------------------------------------------
% bdata = local_bindata( data, bsize )
% bin 2D data
%------------------------------------------------------------------------
function bdata = local_bindata( data, bsize )

% generate the edges arrays
edgesx                          = bsize{ 1 };
edgesy                          = bsize{ 2 };
edgesx                          = unique( edgesx( : ) );
edgesy                          = unique( edgesy( : ) );
edgesx( end )                   = edgesx( end ) + 100 * eps;
edgesy( end )                   = edgesy( end ) + 100 * eps;
nbinsx                          = length( edgesx ) - 1;
nbinsy                          = length( edgesy ) - 1;

% remove irrelevant data
kidx                            = data( :, 1 ) >= edgesx( 1 ) ...
    & data( :, 1 ) < edgesx( nbinsx + 1 ) ...
    & data( :, 2 ) >= edgesy( 1 ) ...
    & data( :, 2 ) < edgesy( nbinsy + 1 );
data                            = data( kidx, : );

% fast binning (can be made faster if use sort to hist too)
bdata                           = zeros( nbinsx, nbinsy );
idx                             = zeros( size( data ) );
[ sdata, sidx ]                 = sortrows( data, 1 );
if size( data, 1 ) == 1
    cn                          = [ 0; cumsum( histc( data( :, 1 ), edgesx ) )' ];
else
    cn                          = [ 0; cumsum( histc( data( :, 1 ), edgesx ) ) ];
end
for i                           = 1 : nbinsx
    yidx                        = cn( i ) + 1 : cn( i + 1 );
    [ y2, y2idx ]               = histc( sdata( yidx, 2 ), edgesy );
    bdata( i, : )               = y2( 1 : nbinsy );
    idx( sidx( yidx ), 1 )      = i;
    idx( sidx( yidx ), 2 )      = y2idx;
end

return % local_bindata

%------------------------------------------------------------------------
% b = local_bounds( x, alpha )
% lower and upper bounds of a vector
%------------------------------------------------------------------------
function b = local_bounds( x, alpha )

xs                              = sort( x );
lx                              = sum( ~isnan( x ) );
si                              = max( floor( alpha / 2 * lx ), 1 );
ei                              = min( ceil( ( 1 - alpha / 2 ) * lx ), lx );
b                               = [ xs( si )'; xs( ei )' ];

return % local_bounds

%------------------------------------------------------------------------
% [ phi, R, s ] = local_circ_mean( t, f )
% compute mean direction - matrix columns
%------------------------------------------------------------------------
function [ phi, R, s ] = local_circ_mean( t, f )

% argument handling and trigonometric functions
nargs                           = nargin;
if nargs < 2 || isempty( f )
    nans                        = isnan( t );
    n                           = sum( ~nans, 1 );
    x                           = cos( t );
    y                           = sin( t );
else
    if size( t, 2 ) == 1
        t                       = t * ones( 1, size( f, 2 ) );
    elseif size( f, 2 ) == 1
        f                       = f * ones( 1, size( t, 2 ) );
    end
    if ~isequal( size( f ), size( t ) )
        error( 'input size mismatch' )
    end
    n                           = nansum( f );
    x                           = f .* cos( t );
    y                           = f .* sin( t );
end

% compute direction
sumx                            = nansum( x, 1 );
sumy                            = nansum( y, 1 );
C                               = sumx ./ n;
S                               = sumy ./ n;
phi                             = mod( atan2( S, C ), 2 * pi );

% compute dispersion
if nargout > 1
    R                           = ( C.^2 + S.^2) .^ 0.5;
    R( R > 1 )                  = 1;
    R( R < 0 )                  = 0;
end
if nargout > 2
    s                           = ( -2 * log( R ) ) .^ 0.5;
end

return % local_circ_mean

%------------------------------------------------------------------------
% idx = local_find_local_minima( x )
% detect all local minima in a vector
%------------------------------------------------------------------------
function idx = local_find_local_minima( x )

x                               = x( : );
d2                              = diff( sign( diff( x ) ) );
idx                             = find( d2 > 1 ) + 1;

return % local_find_local_minima

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
% K = local_gausskernel( sigmaX, N )
% create a 2D gaussian
%------------------------------------------------------------------------
function K = local_gausskernel( sigmaX, N )                                 % 1D Gaussian kernel K with N samples and SD sigmaX

x                               = -( N - 1 ) / 2 : ( N - 1 ) / 2;
K                               = 1 / ( 2 * pi * sigmaX ) * exp( -( x.^2 / 2 / sigmaX^2 ) );

return

%------------------------------------------------------------------------
% out = getdatainranges( in, mat )
% according to a 2-column matrix of ranges
%------------------------------------------------------------------------
function out = local_getdatainranges( in, mat )

% argument handling
nargs                       	= nargin;
if nargs < 2 || isempty( in ) || isempty( mat )
    out                         = []; 
    return
end

% preparations
sx                              = size( in );
if length( sx ) == 2 && sum( sx == 1 ) >= 1 && sx( 2 ) > 1
    in                          = in';
end
n                               = size( in, 1 );
mat( sum( mat > n, 2 ) > 0, : ) = [];
if isempty( mat )
    out                         = [];
    return
end
m                               = size( mat, 1 );
if sum( mat( :, 1 ) <= mat( :, 2 ) ) ~= m
    sidx                        = mat( :, 1 ) > mat( :, 2 );
    t                           = mat( sidx, : );
    mat( sidx, : )              = t( :, [ 2 1 ] );
end
durs                            = diff( mat, [], 2 ) + 1;
eidx                            = cumsum( durs );
sidx                            = [ 1; eidx( 1 : m - 1 ) + 1 ];
out                             = zeros( min( eidx( m ), n ), size( in, 2 ) );

% extract
for i                           = 1 : m
    if mat( i, 2 ) > n
        continue
    end
    tidx                        = mat( i, 1 ) : mat( i, 2 );
    idxI                        = sidx( i ) : eidx( i );
    out( idxI, : )              = in( tidx, : );
end

return % local_getdatainranges

%------------------------------------------------------------------------
% [ h, u ] = local_uhist( x )
% count of the unique values in a vector
%------------------------------------------------------------------------
function [ h, u ] = local_uhist( x )

x                               = x( : )';
u                               = unique( x );
nu                              = length( u );
h                               = zeros( 1, nu );
for i                           = 1 : nu 
    idx                         = x == u( i );
    h( i )                      = sum( idx );
end

return

% EOF