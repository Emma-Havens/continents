 % THERMAL CONVECTION IN A BOX
% Free-slip mechanical boundary conditions on all boundaries
% Isothermal boundary conditions on top and bottom boundaries
% Reflecting thermal boundary conditions on lateral boundaries

clc;
clear variables;

addpath('../common_funcs')
set(0, 'DefaultFigureRenderer', 'zbuffer');

gamma               =       9;

%--------------------------------------------------------------------------
% PHYSICS (nondimensionalised values)
%--------------------------------------------------------------------------
Ra                  =      1e4;  % Rayleigh number
G                   =   [0 -1];  % Gravity vector [Gx Gy]
KAPPA               =        1;  % Thermal diffusivity
time_final          =       .5;  % Total simulation time
T_bottom            =        1;  % Temperature at the bottom boundary
T_top               =        0;  % Temperature at the top boundary

%--------------------------------------------------------------------------
% PART OF INPUT TO THERMAL AND MECHANICAL MILAMIN-SOLVERS
%--------------------------------------------------------------------------
elemType            =   'quad';  % Type of elements - quadratic
reorder             =    'amd';
method              =    'opt';
lumping             =   'true';

%--------------------------------------------------------------------------
% SPATIAL AND TEMPORAL DISCRETISATION
%--------------------------------------------------------------------------
ndim                =    2;  % Number of spatial dimensions
asp                 =    8;  % Aspect ratio of the box
L                   =    [asp*.5;.5]; % Horizontal and vertical extent of the box: x -> [-L(1) L(1)], y -> [-L(2) L(2)]
ney                 =    60; % Number of elements in vertical direction

nip                 =    4;  % Number of integration points -> nip^2
nnodel              =    9;  % Number of GCOORD per element

% Generating grid
[GCOORD, ELEM2NODE, ind] =    generate_quad_grid(asp, ney, nnodel, L);

ELEM2NODE	        =    int32(ELEM2NODE);
nel                 =    size(ELEM2NODE,2);
nnod                =    size(GCOORD,2);

MU                  =    ones(nel,nip^2);  % Viscosity
PHASES              =    ones(nel,1);

dx                  =    4*L(1)*L(2)/nel;
dt                  =    1e4*dx/(Ra/exp(-gamma)); % Initial time-stepping size
dt_fig              =    1e-4;
time_fig            =    dt_fig;
i_fig               =    1;

%--------------------------------------------------------------------------
% BOUNDARY CONDITIONS
%--------------------------------------------------------------------------
% Indeces of nodes that lie on the boundaries
t                   =    ind(:,end);
b                   =    ind(:,1);
l                   =    ind(1,:);
r                   =    ind(end,:);

% Free slip mechanical boundary conditions
bc_ind              =    [2*(l(:)-1)+1; 2*(r(:)-1)+1; 2*(b(:)-1)+2; 2*(t(:)-1)+2];
bc_val              =    zeros(size(bc_ind));

% Isothermal boundary conditions
bc_ind_therm        =    [t                  ; b                     ];
bc_val_therm        =    [T_top*ones(size(t)); T_bottom*ones(size(b))];

%--------------------------------------------------------------------------
% ALLOCATING FOR VELOCITY
%--------------------------------------------------------------------------
Vel                 =    zeros(2*nnod,1);
Vel(bc_ind)         =    bc_val;

%--------------------------------------------------------------------------
% INITIALIZING TEMPERATURE
%--------------------------------------------------------------------------
TEMP                =    .9*(T_bottom-T_top)*ones(nnod,1) - 5e-2*rand(nnod,1);
TEMP(bc_ind_therm)  =    bc_val_therm;

%--------------------------------------------------------------------------
% INTEGRATION POINTS
%--------------------------------------------------------------------------
[IP_X, IP_w]        =    ip_line(nip);
[xx, yy]            =    meshgrid(IP_X,IP_X);
IP_X                =    [xx(:) yy(:)];
NIP                 =    shp_deriv_quad_notCell(IP_X,nnodel);
X_ip                =    reshape(GCOORD(1,ELEM2NODE),size(ELEM2NODE))'*NIP;
Y_ip                =    reshape(GCOORD(2,ELEM2NODE),size(ELEM2NODE))'*NIP;

%--------------------------------------------------------------------------
% TIME LOOP
%--------------------------------------------------------------------------
name_file           =    '/Users/elvira/Documents/MATLAB_OUTPUT/THERMOMECHANICAL/gamma9_mat/';
save([name_file 'mats/info'],'GCOORD','ELEM2NODE','gamma');

time                =    0;

while time<time_final
    %----------------------------------------------------------------------
    % DIFFUSION
    %----------------------------------------------------------------------
    TEMP            =    thermal2d_transient(ELEM2NODE, PHASES, GCOORD, KAPPA, bc_ind_therm, nip, TEMP, dt, elemType, reorder, method, lumping);
    TEMP_ip         =    TEMP(ELEM2NODE)'*NIP;
    MU              =    exp(-gamma*TEMP_ip);
    
    %----------------------------------------------------------------------
    % MECHANICAL SOLVER
    %----------------------------------------------------------------------
    [Vel, Pressure] =    mechanical2d_contMU(ELEM2NODE, PHASES, GCOORD, MU, -TEMP, Ra, G, bc_ind, bc_val, nip, elemType, reorder, method);
    Vx              =    Vel(1:2:end-1)';
    Vy              =    Vel(2:2:end)';
    
    %----------------------------------------------------------------------
    % ADVECTION - SHOOTING BACK CHARACTERISTICS
    %----------------------------------------------------------------------
    % Feet of characteristics - global and local positions
    GCOORD_char     =    GCOORD - dt*[Vx; Vy];
    [U, V, IND]     =    uv_quads(L, [asp*ney; ney], GCOORD_char);
    SHP             =    shp_quad([U; V],nnodel);
    
    % Interpolating temperature at the feet of characteristics
    TEMP            =    sum(SHP.*TEMP(ELEM2NODE(:,IND)))';
    
    %----------------------------------------------------------------------
    % UPDATE TIME
    %----------------------------------------------------------------------
    time            =    time + dt;
    dt              =    5*min(dx./sqrt(Vx.^2+Vy.^2));
    
    %----------------------------------------------------------------------
    % PLOTTING
    %----------------------------------------------------------------------
    if time > time_fig
        %         figure(1);clf; hold on;
        %         ncorners        =    4;
        %         xplot           =    reshape(GCOORD(1,ELEM2NODE(1:ncorners,:)),ncorners,nel);
        %         yplot           =    reshape(GCOORD(2,ELEM2NODE(1:ncorners,:)),ncorners,nel);
        %         qplot           =    reshape(TEMP(ELEM2NODE(1:ncorners,:)),ncorners,nel);
        %         hp              =    patch(xplot,yplot,qplot);
        %         set(hp,'EdgeColor','n','FaceColor','interp');
        %         axis equal; axis off;
        %         colormap('jet'); caxis([0 1]); % Pretty colors
        %         drawnow
        %         time_fig        =    time_fig + dt_fig;
        %
        %         % save figure
        %         name_fig        =    [name_file 'ifig' num2str(i_fig)];
        %         print('-dpng','-r150',name_fig)
        %         crop([name_fig '.png']);
        
        % save data
        save([name_file 'mats/imat' num2str(i_fig)],'TEMP','Vx','Vy','time');
        
        i_fig           =    i_fig + 1;
        
    end
    disp(time)
    
    
end





