% NR_make_Fig1      graphics for Figure 1 in Stark, Levi, Rotstein, 2022
%
% This script creates the sub-figures that compose Figure 1 in 
% Stark, Levi, Rotstein, 2022, PLOS Computational Biology
%
% Figure 1. Cycle-averaged firing rate resonance and spike timing resonance
% (A) -
% (B) Intrinsic oscillations vs. resonance                         
% (C) Cycle-averaged firing rate resonance
% (D) Spike timing resonance
%
% calls             NR_st_coherence_demo

% 08-aug-21 ES

% last update
% 01-jul-22 AL

function NR_make_Fig1( savef, pstr, outdir)

% constants
prefix                          = 'resonance_network';
renderer_name                   = 'painters';
resize                          = '-bestfit';
color_input                     = [ 0 0 0.7 ];                              % input (current): blue
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

% ----------------------------------------------------------------------
% (B) Intrinsic oscillations vs. resonance
% ----------------------------------------------------------------------
    
% Caricature - Intrinsic oscillations
dt                              = 0.001;
T                               = 2;
t                               = ( dt : dt : T )';
x                               = [ zeros( floor( T/dt/4 ), 1 ); ones( floor( T/dt/2 ), 1 ); zeros( floor( T/dt/4 ), 1 ) ];
f                               = 8; 
y                               = ( ( 1 + sin( 2 * pi * t * f - pi / 2 ) ) / 2 );
y                               = [ zeros( floor( T/dt/4 ), 1 ); y( t <= T/2 ); zeros( floor( T/dt/4 ), 1 ) ];
V                               = y;
I                               = x;

fig1( 1 )                       = figure;
subplot( 3, 2, 1 )
ph                              = plot( t, V, 'k', t, I - 1.5, 'b' ); 
set( ph( 1 ), 'color', color_output );
set( ph( 2 ), 'color', color_input );
set( gca, 'tickdir', 'out', 'box', 'off' )
ylim( [ -1.5 1 ] )

% Caricature - Resonance
fvals                           = [ 4 8 12 ];
amps                            = [ 0.5 1 0.25 ];

fig1( 2 )                       = figure;
for i                           = 1 : length( fvals )
    f                           = fvals( i );
    y                           = ( ( 1 + sin( 2 * pi * t * f - pi / 2 ) ) / 2 );
    y                           = [ zeros( floor( T/dt/4 ), 1 ); y( t <= T/2 ); zeros( floor( T/dt/4 ), 1 ) ];
    I                           = y;
    V                           = y * amps( i );
    subplot( 3, 2, 2 * i - 1 )
    ph                          = plot( t, V, 'b', t, I - 1.5, 'k' );
    ylim( [ -1.5 1 ] )
    set( ph( 1 ), 'color', color_output );
    set( ph( 2 ), 'color', color_input );
    set( gca, 'tickdir', 'out', 'box', 'off' )
end

%-----
fig                             = fig1;
if savef
    for i                       = 1 : length( fig )
        figname = [ outdir filesep prefix '_FIG1B_part' num2str( i ) ];
        figure( fig( i ) );
        figi                    = gcf;
        figi.Renderer           = renderer_name;
        pause( 0.2 )
        if isequal( pstr, '-dpdf' )
            rsz                 = resize;
        else
            rsz                 = '';
        end
        print( figi, pstr, figname, rsz )
    end
end


% ----------------------------------------------------------------------
% (C,D) Cycle-averaged firing rate resonance and spike timing resonance
% ----------------------------------------------------------------------

% Two cases: 
% (1) Firing rate increase in a specific frequency band, no phase locking
% (2) Firing rate the same in all frequencies, phase locking at a certain
% frequency band
    
fig2( 1 )                       = NR_st_coherence_demo(1);
fig2( 2 )                       = NR_st_coherence_demo(2);

%-----
fig                             = fig2;
if savef
    for i                       = 1 : length( fig )
        figname = [ outdir filesep prefix '_FIG1CD_part' num2str( i ) ];
        figure( fig( i ) );
        figi                    = gcf;
        figi.Renderer           = renderer_name;
        pause( 0.2 )
        if isequal( pstr, '-dpdf' )
            rsz                 = resize;
        else
            rsz                 = '';
        end
        print( figi, pstr, figname, rsz )
    end
end

return
    
% EOF
