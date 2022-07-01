% NR_st_coherence_demo          demonstrate the difference between coherence and firing rate metrics
%
% call                          NR_st_coherence_demo( cs )
%
% gets                          cs          1: higher rate at some frequency
%                                           2: phase locking at one frequency, uniform at others
%                                           3: strong phase locking at one frequency, weaker at others
%                                           4: one phase at one frequency, zero phase at others
%                                           5: uniform phase and rate in all frequencies
%
% optional arguments (given as name/value pairs)
%
%                               frqs    {1:40}      [Hz]    frequencies
%                               T       {10}        [s]     duration of every individual sinusoid
%                               Fs      {5000}      [Hz]    sampling frequency
%                               nspikes_generator {'fixed'} or 'poisson'
%                               rates   []          [spk/s] for sines, a vector of rate for every sinusoid
%                               myus    []                  for sines, a vector of phase for every sinusoid
%                               kappas  []                  for sines, a vector of concentration for every sinusoid
%                               BL      {20}        [spk/s] baseline rate
%                               fROI    {[0 40]}    [Hz]    argument to st_coherence and to st_fingerprint
%                               M       {1}                 argument to st_coherence and to st_fingerprint
%                               mtNW    {3}                 argument to st_coherence
%                               sdGauss {10}        [ms]    argument to st_fingerprint
%                               graphics{1}                 flag
%
% does
% (1) generate input (concatenated sinusoids)
% (2) generate spikes for every cycle (generative phenomenological model)
% (3) analyze the data using coherence, firing rates, and spike count metrics
%
% algorithmic details:
% for every cycle (of every frequency), determine:
% (1) how many spikes should there be (according to the rate and cycle
% duration), given that the rate is some rate - this could be 'fixed' or
% 'poisson'. Note that 'fixed' generates sawteeth in the rate-frequency
% map, this is due to integer rounding leading to highpass (same as LIF).
% solved using probabilistic rounding. 
% (2) what phase these spikes should appear at. here, this is done using a
% uni-modal von-Mises distribution
% 
% calls                         ParseArgPairs
%                               plot_raster, plotOneSpectrogram
%                               st_coherence, st_fingerprint

% 16-jul-21 ES

% last update 
% 30-jun-22

function [ fig, st, x ]         = NR_st_coherence_demo( cs , varargin )

%----------------------------------------------------------------------
% defaults
% input signal defaults
frqs_DFLT                       = ( 1 : 40 )';                              % [Hz], frequencies of individual sinusoids
T_DFLT                          = 10;                                       % [s], duration of every individual sinusoid
Fs_DFLT                         = 5000;

% simulataion defaults
BL_DFLT                         = 20;                                       % [spks/s]
nspikes_generator_DFLT        	= 'fixed';

% analysis defaults
fROI_DFLT                       = [ 0 40 ];                                 % [Hz]
M_DFLT                          = 1; 
mtNW_DFLT                       = 3;                                        % argument to st_coherence 
sdGauss_DFLT                    = 10;                                       % [ms]

% graphics defaults
graphics_DFLT                   = 1;

%----------------------------------------------------------------------
% arguments
nargs                           = nargin;
if nargs < 1 || isempty( cs ) || isempty( cs )
    cs                          = 1;
end
[ frqs, T, Fs ...
    , nspikes_generator ...
    , rates, myus, kappas, BL ...
    , fROI ...
    , M, mtNW, sdGauss ...
    , graphics ...
    ]                = ParseArgPairs(...
    { 'frqs', 'T', 'Fs' ...
    , 'nspikes_generator' ...
    , 'rates', 'myus', 'kappas', 'BL' ...
    , 'fROI' ...
    , 'M', 'mtNW', 'sdGauss' ...
    , 'graphics' }...
    , { frqs_DFLT, T_DFLT, Fs_DFLT ...
    , nspikes_generator_DFLT ...
    , [], [], [], BL_DFLT ...
    , fROI_DFLT ...
    , M_DFLT, mtNW_DFLT, sdGauss_DFLT ...
    , graphics_DFLT }...
    , varargin{ : } );

nspikes_generator               = lower( nspikes_generator );
if ~ismember( nspikes_generator, { 'poisson', 'fixed' } )
    error( 'nspikes_generator must be either poisson or fixed' )
end

%----------------------------------------------------------------------
% determine setup

nfrqs                           = length( frqs );

switch cs
    
    case 0
        rates                   = rates( : );
        if numel( rates ) ~= nfrqs
            if numel( rates ) ~= 1
                error( 'number of elements in rates must be as in frqs, or 1' )
            end
            rates               = rates * ones( nfrqs, 1 );
        end
        myus                    = myus( : );
        if numel( myus ) ~= nfrqs
            if numel( myus ) ~= 1
                error( 'number of elements in myus must be as in frqs, or 1' )
            end
            myus                = myus * ones( myus, 1 );
        end
        kappas                  = kappas( : );
        if numel( kappas ) ~= nfrqs
            if numel( kappas ) ~= 1
                error( 'number of elements in kappas must be as in frqs, or 1' )
            end
            kappas              = kappas * ones( nfrqs, 1 );
        end
        tit                     = 'user defined';
        
    case 1
        % first case - same rate for all, uniform phase for all, one frequency higher rate
        rates                 	= BL * ones( nfrqs, 1 );
        myus                   	= 0 * ones( nfrqs, 1 );
        kappas                 	= 0 * ones( nfrqs, 1 );
        rates( 8 : 12 )        	= BL * 2;
        tit                  	= sprintf( 'rate increase only; BL=%0.2g', BL );
        
    case 2
        % second case - same rate for all, uniform phase for all, one frequency strong locking
        rates                  	= BL * ones( nfrqs, 1 );
        myus                  	= 0 * ones( nfrqs, 1 );
        kappas                  = 0 * ones( nfrqs, 1 );
        kappas( 8 : 12 )        = 2;
        tit                     = sprintf( 'phase locking only; BL=%0.2g', BL );
        
    case 3
        % third case - same rate for all, phase locking for all, one frequency stronger locking
        rates                   = BL * ones( nfrqs, 1 );
        myus                    = 0 * ones( nfrqs, 1 );
        kappas                  = 0.5 * ones( nfrqs, 1 );
        kappas( 8 : 12 )        = 2;
        tit                     = sprintf( 'phase locking change only; BL=%0.2g', BL );
        
    case 4
        % fourth case - same rate for all, same phase locking for all, one frequency different phase
        rates                   = BL * ones( nfrqs, 1 );
        myus                  	= 0 * ones( nfrqs, 1 );
        kappas                  = 0.3 * ones( nfrqs, 1 );
        myus( 8 : 12 )          = pi * 7 / 8;
        myus( [ 7 13 ] )        = pi / 2;
        tit                     = sprintf( 'phase change only; BL=%0.2g', BL );

    case 5
        % fifth case - same rate for all, same phase locking for all, same phase for all
        rates                   = BL * ones( nfrqs, 1 );
        myus                  	= 0 * ones( nfrqs, 1 );
        kappas                  = 0 * ones( nfrqs, 1 );
        tit                     = sprintf( 'no changes; BL=%0.2g', BL );
        
    case 6
        % sixth case - same rate for all, uniform phase for all, one frequency higher rate
        rates                 	= BL * ones( nfrqs, 1 );
        myus                   	= 0 * ones( nfrqs, 1 );
        kappas                 	= 0 * ones( nfrqs, 1 );
        rates( 8 : 12 )        	= BL * 2;
        kappas( 8 : 12 )        = 2;
        tit                     = sprintf( 'rate and phase locking; BL=%0.2g', BL );
        
    otherwise
        error( 'cs %d not supported', cs )
     
end

%----------------------------------------------------------------------
% simulate trains

% initialize input and output
dt                              = 1 / Fs;
t                               = ( dt : dt : T )';
nunique                         = nfrqs;
x                               = [];
st                              = [];

for i                           = 1 : nunique                               % chirp cycle/frequency train
    
    % generate single-frequency sinusoid
    f                           = frqs( i );
    x0                          = sin( 2 * pi * f * t - pi / 2 ) / 2 + 0.5; % discrete sinusoid, amplitude 0-1
    ncycs                       = T * f;                                    % number of cycles
    k                           = i;
    
    % determine number of spikes/cycle
    exp_spks_cyc                = rates( k ) / f;
    switch nspikes_generator
        case 'poisson'
            spks_cyc            = poissrnd( exp_spks_cyc, [ ncycs 1 ] ); 	% spikes per cycle
        case 'fixed'
            spks_cyc            = ones( ncycs, 1 ) * exp_spks_cyc;
            if exp_spks_cyc ~= round( exp_spks_cyc )
                spks_cyc        = diff( [ 0; local_pround( cumsum( spks_cyc ) ) ] );
            end
    end
    nspks                       = sum( spks_cyc );
    
    % draw random phases
    phs                         = local_vmrnd( myus( k ), kappas( k ), nspks ); 	% phase of each spike
    
    % convert the phases to times
    ei                          = cumsum( spks_cyc );
    si                          = [ 1; ei( 1 : ncycs - 1 ) + 1 ];
    st0                         = NaN( nspks, 1 );
    for j                       = 1 : ncycs
        if si( j ) > ei( j )
            continue
        end
        idx                     = si( j ) : ei( j );
        phsj                    = sort( phs( idx ) + pi );                	% pi: offset from cosine (sinusoid x0 starts at 0)
        st0( idx )              = phsj / ( 2 * pi * f ) + ( j - 1 ) / f;  	% [rad] in cycle -> [s] in train (may be single-cycle train)
    end
    
    % wrap last spikes back due to the pi offsetting
    st0                         = sort( mod( st0, T ) );
    
    % accumulate over frequencies
    os                          = ( i - 1 ) * T;
    x                           = [ x; x0 ];
    st                          = [ st; st0 + os ];
    
end
st                              = ceil( st / dt );                        	% [samples] @ Fs

%----------------------------------------------------------------------
% analyze
[ cohs, ~, fo, ~, fidx ]        = st_coherence( st, x( : ), 'spkFs', Fs, 'xFs', Fs ...
    , 'fROI', fROI, 'Fs', 1250, 'graphics', 0, 'M', M, 'mtNW', mtNW );

[ frate_vec, ~, rate_map ...
    , fbins, pbins ]            = st_fingerprint( st, x( : ), 'spkFs', Fs, 'xFs', Fs ...
    , 'fROI', fROI, 'Fs', 1250, 'graphics', 0, 'M', M, 'sdGauss', sdGauss );

%----------------------------------------------------------------------
% graphics
if ~graphics
    return
end

fig                             = figure;

subplot( 3, 3, 1 )
ph                              = plot( frqs, rates, '.-b' );
set( ph, 'color', [ 0 0 0.7 ] )
ylabel( 'Rates [spks/s]' )
ylims                           = ylim;
ylim( [ 0 max( [ ylims( : ); 2.1 * BL ] ) ] )
ylims1                          = ylim;
xlabel( 'Frequency [Hz]' )
title( tit )

subplot( 3, 3, 4 )
ph                              = plot( frqs, kappas, '.-b' );
set( ph, 'color', [ 0 0 0.7 ] )
ylabel( '\kappa' )
ylims                           = ylim;
ylim( [ 0 max( [ ylims( : ); 2 ] ) ] )
xlabel( 'Frequency [Hz]' )

subplot( 3, 3, 2 )
ph                              = plot( fbins, frate_vec, '.-k' );
set( ph, 'color', [ 0 0 0 ] )
ylims                           = ylim;
ylim( [ 0 max( [ ylims( : ); 2.1 * BL ] ) ] )
ylims2                          = ylim;
ylabel( 'Firing rate [spks/s]' )
xlabel( 'Frequency [Hz]' )

subplot( 3, 3, 5 )
ph                              = plot( fo( fidx ), cohs( fidx ), '.-k' );
set( ph, 'color', [ 0 0 0 ] )
ylims                           = ylim;
ylim( [ 0 max( [ ylims( : ); 0.15 ] ) ] )
ylabel( 'Coherence' )
xlabel( 'Frequency [Hz]' )

subplot( 3, 3, 3 )
prate_map                       = rate_map;
nans                            = isnan( prate_map );
prate_map( nans )               = 0;
plotOneSpectrogram( pbins, fbins, prate_map' );
axis square
title( sprintf( '%0.2g spks/s', max( rate_map( : ) ) ) )

ylims                           = [ min( [ ylims1( 1 ) ylims2( 1 ) ] ) max( [ ylims1( 2 ) ylims2( 2 ) ] ) ];
for i                           = 1 : 5
    subplot( 3, 3, i )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    if i <= 2
        ylim( ylims )
    end
    if i == 3
        ylim( fROI )
    else
        xlim( fROI )
    end
end

% plot input and output for two selected frequencies
tplot                           = [ 1 2 ];                                  % [s]
iplot                           = [ 5 10 ];                                 % indices
for i                           = 1 : length( iplot )
    k                           = iplot( i );
    os                          = ( k - 1 ) * T;
    per                         = [ os ( os + T ) ] * Fs + [ 1 0 ];
    sidx                        = st >= per( 1 ) & st <= per( 2 );
    st0                         = st( sidx ) - per( 1 ) + 1;
    x0                          = x( per( 1 ) : per( 2 ) );
    f                           = frqs( k );
    subplot( 3, 2, 4 + i )
    ph                          = plot( ( 1 : length( x0 ) ) / Fs, x0, 'b' );
    set( ph, 'color', [ 0 0 0.7 ] );
    sh                          = plot_raster( { st0 }, t );
    set( sh, 'color', [ 0 0 0 ] );
    ylim( [ 0 1.3 ] )
    xlim( tplot )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    title( sprintf( 'F=%0.3g Hz', f ) )
end

return

%------------------------------------------------------------------------
% prx = local_pround( x )
% probabilistic rounding to nearest integer
%------------------------------------------------------------------------
function prx = local_pround( x )

sz                              = size( x );
x                               = x( : );
nx                              = length( x );
prx                             = NaN( nx, 1 );

% positive scalars:
pidx                            = x > 0;
fx                              = floor( x( pidx ) );
dx                              = x( pidx ) - fx;
rx                              = rand( sum( pidx ), 1 );
h                               = double( rx <= dx );
prx( pidx )                     = fx + h;

% negative scalars:
fx                              = ceil( x( ~pidx ) );
dx                              = -x( ~pidx ) + fx;
rx                              = rand( sum( ~pidx ), 1 );
h                               = double( rx <= dx );
prx( ~pidx )                    = fx - h;

prx                             = reshape( prx, sz );

return % local_pround

%------------------------------------------------------------------------
% vm_theta = vmrnd( myu0, kappa, n )
% von Mises distribution random number generator
%------------------------------------------------------------------------
function vm_theta = local_vmrnd( myu0, kappa, n )

tol                             = 0.001;
theta                           = 0 : tol : 2 * pi;
vm                              = 1 / ( 2 * pi * besseli( 0, kappa ) ) * exp( kappa * cos( theta - myu0 ) );

cum_vm                          = cumsum( vm ) / sum( vm );
index                           = ones( n, 1 );
rand_nums                       = rand( n, 1 );

for i                           = 1 : n
    rand_num                    = rand_nums( i );
    if min( cum_vm ) < rand_num
        index( i )              = find( cum_vm < rand_num, 1, 'last' );
    end
end
vm_theta                        = theta( index );

return % local_vmrnd

% EOF