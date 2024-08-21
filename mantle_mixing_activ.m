
% Free-slip mechanical boundary conditions on all boundaries
% Reflecting thermal boundary conditions on lateral boundaries

disp('Start Sim')
% maxNumCompThreads(str2num(getenv('SLURM_NPROCS')));

clc;
clear variables;
set(0, 'DefaultFigureRenderer', 'zbuffer');

fnum                =        1;

%--------------------------------------------------------------------------
% PHYSICS (nondimensionalised values)
%--------------------------------------------------------------------------
Ra                  =      1e5;  % Rayleigh number
Bu                  =     -2;  % Compositional buoyancy ratio: rho_crust / (~200 kg/m^3); use negative values for chemically light, positive values for chemically anom
G                   =   [0 -1];  % Gravity vector [Gx Gy]
KAPPA               =   [1; 1];  % Thermal diffusivity
FUG                  =   [0; .5];  % Viscosity max/min
FUG_inc              =        .25;     % Amount viscosity increases on melt
H_activ             = [370; 330];  % Changes temperature dependence of viscosity max/min
H_activ_inc         =        0;      % Amount increases on melt
RHO                 =  [Bu; 0];  % Compositional density anomaly max/min
RHO_inc             =        0;  % Amount rho decreases on melt -2e-1
depth               =     2900;
time_final          =     2e-1;  % Total simulation time
thickness_exact     = 50/depth;  % Thickness of the compositionally anomalous layer
Z_660               =660/depth;  % Depth of 660km discontinuity
T_bottom            =        1;  % Temperature at the bottom boundary
T_top               =        0;  % Temperature at the top boundary

%--------------------------------------------------------------------------
% PART OF INPUT TO MILAMIN-SOLVER
%--------------------------------------------------------------------------
elemType            =   'quad';  % Type of elements - quadratic
reorder             =    'amd';
method              =    'opt'; %'std'
lumping             =   'true';

%--------------------------------------------------------------------------
% SPATIAL AND TEMPORAL DISCRETISATION
%--------------------------------------------------------------------------
ndim                =    2;  % Number of spatial dimensions
asp                 =    1;  % Aspect ratio of the box
L                   =    [asp*.5;.5]; % Horizontal and vertical extent of the box: x -> [-L(1) L(1)], y -> [-L(2) L(2)]

% MECHANICAL MESH
ney                 =    50; % Number of elements in vertical direction
nip                 =    4;  % Number of integration points -> nip^2
nnodel              =    9;  % Number of GCOORD per element

% Generating grid
[GCOORD, ELEM2NODE, ind] = generate_quad_grid(asp, ney, nnodel, L);

ELEM2NODE	        =    int32(ELEM2NODE);
nel                 =    size(ELEM2NODE,2);
nnod                =    size(GCOORD,2);
IND                 =    1:(asp*ney*ney);
IND                 =    reshape(IND,ney,asp*ney);

dx                  =    4*L(1)*L(2)/nel;
dt                  =    1e4*dx/Ra; % Initial time-stepping size

% THERMAL MESH - linear elements with four nodes per element; effectively splitting mechanical elements into 4;
nnodel_t            =    4;
ney_t               =    2*ney;
nip_t               =    2;
[~, ELEM2NODE_t, ind_t] = generate_quad_grid(asp, ney_t, nnodel_t, L);
nel_t               =    size(ELEM2NODE_t,2);

% adjust thickness to line up with mechanical mesh
res = L(2)*2/ney;
thickness = floor(thickness_exact/res)*res;

%--------------------------------------------------------------------------
% MARKERS
%--------------------------------------------------------------------------
nmy                 =    4*ney;
nm                  =    asp*nmy*nmy;
dm                  =    2*L(2)/nmy;
Mx                  =    linspace(-(L(1)-dm), L(1)-dm, asp*nmy);
My                  =    linspace(-(L(2)-dm), L(2)-dm, nmy);
[MX, MY]            =    ndgrid(Mx,My);
MARK.PHASES         =    2*ones(1,nm); % phase of each marker, has values 1 or 2
MARK.NODES          =    [MX(:)'; MY(:)'];
MARK.RHO            =    RHO(2)*ones(1,nm); % default mantle value
MARK.FUG             =    FUG(2)*ones(1,nm); % default mantle value
MARK.H_activ        =    H_activ(2)*ones(1,nm); % default mantle value
MARK.NUMMELT        =    zeros(1,nm);
MARK.TIMEMELT       =    zeros(1,nm);

prev_melt_bool              =    MARK.NODES(2,:) > (L(2)-thickness); % 1 if within melt region, 0 else
MARK.PHASES(prev_melt_bool) =    1;
MARK.RHO(prev_melt_bool)    =    MARK.RHO(prev_melt_bool) + RHO_inc; 
MARK.FUG(prev_melt_bool)     =    MARK.FUG(prev_melt_bool) + FUG_inc;
MARK.H_activ(prev_melt_bool)=    MARK.H_activ(prev_melt_bool) + H_activ_inc;

% figure(1); hold on;
% plot(MARK.NODES(1,:),MARK.NODES(2,:),'.k');
% plot(MARK.NODES(1,melt_bool),MARK.NODES(2,melt_bool),'om');
% axis tight

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
TEMP_nodes          =    .5*(T_bottom-T_top)*ones(nnod,1) - 5e-2*rand(nnod,1);
TEMP_nodes(bc_ind_therm)  =    bc_val_therm;

%--------------------------------------------------------------------------
% INITIALISING DENSITY FIELD
%--------------------------------------------------------------------------
% For each marker, what mechanical element it belongs to
[U, V, MARK2ELEM]     =  uv_quads(L,[asp*ney; ney], MARK.NODES);
SHP                   =  shp_quad([U; V], nnodel);

% figure(1); hold on;
% patch(reshape(GCOORD(1,ELEM2NODE(1:4,:)),4,nel),reshape(GCOORD(2,ELEM2NODE(1:4,:)),4,nel),IND(:));
% scatter(MARK.NODES(1,:),MARK.NODES(2,:),50,MARK2ELEM,'filled','MarkerEdgeColor','k')
% axis equal

% Markers to elements - eventually should transfer viscosity from markers to integration points in the same way as density
Mweights            =    accumarray(MARK2ELEM',ones(size(MARK2ELEM)),[nel 1]);
PHASEsum            =    accumarray(MARK2ELEM',MARK.PHASES,[nel 1]);
PHASES              =    round(PHASEsum./Mweights);
NUMMELTsum          =    accumarray(MARK2ELEM',MARK.NUMMELT,[nel 1]);
NUMMELT             =    round(NUMMELTsum./Mweights);
TIMEMELTsum         =    accumarray(MARK2ELEM',MARK.TIMEMELT,[nel 1]);
TIMEMELT            =    round(TIMEMELTsum./Mweights);

IND_el              =    reshape(1:nel, ney, asp*ney);

% which mechanical element each thermal element belongs to
THERM2MECH          =    zeros(ney_t,asp*ney_t);
THERM2MECH(1:2:end-1, 1:2:end-1) = IND_el; % upper lefts
THERM2MECH(2:2:end, 2:2:end)     = IND_el; % lower rights
THERM2MECH(1:2:end-1, 2:2:end)   = IND_el; % upper rights
THERM2MECH(2:2:end, 1:2:end-1)   = IND_el; % lower lefts

THERM2MECH          =    THERM2MECH(:);
PHASES_t            =    PHASES(THERM2MECH);
NUMMELT_t           =    NUMMELT(THERM2MECH);
TIMEMELT_t          =    TIMEMELT(THERM2MECH);

% figure(1); hold on;
% patch(reshape(GCOORD(1,ELEM2NODE(1:4,:)),4,nel),reshape(GCOORD(2,ELEM2NODE(1:4,:)),4,nel),PHASES);
% scatter(MARK.NODES(1,:),MARK.NODES(2,:),50,MARK.PHASES,'filled','MarkerEdgeColor','k')
% axis equal
% 
% figure(2); hold on;
% patch(reshape(GCOORD(1,ELEM2NODE_t(1:4,:)),4,nel_t),reshape(GCOORD(2,ELEM2NODE_t(1:4,:)),4,nel_t),PHASES_t);
% scatter(MARK.NODES(1,:),MARK.NODES(2,:),50,MARK.PHASES,'filled','MarkerEdgeColor','k')
% axis equal

% Markers to nodes

[U_t, V_t, MARK2ELEM_t]  =    uv_quads(L, [asp*ney_t; ney_t], MARK.NODES);
SHP_t               =    shp_quad([U_t; V_t], nnodel_t);
nodeBINS            =    ELEM2NODE_t(:,MARK2ELEM_t); % Too much memory used here and in the next line - make more efficient!
Weights_t           =    accumarray(nodeBINS(:), SHP_t(:), [nnod 1]);
markBINS_RHO        =    SHP_t.*repmat(MARK.RHO,nnodel_t,1);
RHO_nodes           =    accumarray(nodeBINS(:), markBINS_RHO(:), [nnod 1])./Weights_t;
markBINS_FUG        =    SHP_t.*repmat(MARK.FUG,nnodel_t,1);
FUG_nodes           =    accumarray(nodeBINS(:), markBINS_FUG(:), [nnod 1])./Weights_t;
markBINS_H_activ    =    SHP_t.*repmat(MARK.H_activ,nnodel_t,1);
H_activ_nodes       =    accumarray(nodeBINS(:), markBINS_H_activ(:), [nnod 1])./Weights_t;


%     figure(2); hold on;
% scatter(GCOORD(1,:),GCOORD(2,:),50,FUGPHASE_nodes,'s','filled','MarkerEdgeColor','k')
% scatter(MARK.NODES(1,:),MARK.NODES(2,:),25,MARK.FUG,'filled','MarkerEdgeColor','k')
% scatter(GCOORD(1,:),GCOORD(2,:),50,RHO_nodes,'s','filled','MarkerEdgeColor','k')
% scatter(MARK.NODES(1,:),MARK.NODES(2,:),25,MARK.RHO,'filled','MarkerEdgeColor','k')

%colorbar;
%axis equal

%--------------------------------------------------------------------------
% TEMPERATURE DEPENDANT VISCOSITY
%--------------------------------------------------------------------------

MU_z                        =    ones(nnod,1);
MU_z(GCOORD(2,:) < Z_660)   =    100 * MU_z(GCOORD(2,:) < Z_660);
MU_perturb                  =    arrayfun(@thermal_profile_perturb_activ_fug, ...
                                              TEMP_nodes, H_activ_nodes, FUG_nodes);
%MU_perturb                  =    arrayfun(@thermal_profile_perturb_activ_no_fug, TEMP_nodes, H_activ_nodes);
MU_nodes                    =    MU_z .* MU_perturb;

MU2EL                       =    reshape(MU_nodes(ELEM2NODE),nnodel,nel);
MU_el                       =    10.^(mean(log10(MU2EL),1))';

disp(min(MU_el))

% figure(fnum);clf; hold on;
%       ncorners        =    4;
%       xplot           =    reshape(GCOORD(1,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
%       yplot           =    reshape(GCOORD(2,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
%       muplot          =    log10(reshape(MU_nodes(ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t));
%       mp              =    patch(xplot,yplot,muplot);
%       axis equal; axis off;
%       colormap('jet'); colorbar; % caxis([-.2 1]);
%       drawnow

%--------------------------------------------------------------------------
% TIME LOOP
%--------------------------------------------------------------------------
time                =    0;
time_save           =    0;
dt_save             =    1e-3;          % time interval between file saves ~10Myr Depends on Ra, might want to change
num_save            =    0;

times               =    zeros(1e1,1);
MASS                =    zeros(size(times));
it                  =    1;
times(it)           =    time;
MASS(it)            =    sum(RHO_nodes);

%--------------------------------------------------------------------------
% CREATE SAVE FOLDERS
%--------------------------------------------------------------------------

if not(isfolder("output_files"))
    mkdir("output_files")
end

folder_exists = true;
num_sim = 0;
while folder_exists
    run_title = strjoin(["Bu_", string(Bu), "+", string(RHO_inc), "_Fug_", ...
        string(FUG(2)), "+", string(FUG_inc), "_H_", string(H_activ(1)), ...
        "+", string(H_activ_inc), "_thick_", string(round(thickness_exact,3)), ...
        "_asp_", string(asp), "_ney_", string(ney), "_", string(num_sim)], '');
    if isfolder(fullfile("output_files", run_title))
        num_sim = num_sim + 1;
    else
        disp(run_title)
        output_path = fullfile("output_files", run_title);
        mkdir(output_path);
        folder_exists = false;
    end
end

static_file_name = fullfile("output_files", run_title, "static_var.mat");

% save static variables
save(static_file_name, "Ra", "Bu", "FUG", "FUG_inc", "RHO", "RHO_inc", "H_activ", "H_activ_inc", ...
    "thickness_exact", "time_final", "dt_save", "asp", "ney", "L")

mkdir("output_files", run_title);
path = fullfile("output_files", run_title);
mkdir(path, "dynamic_vars")
% mkdir output_files dynamic_vars;

while time<time_final

    %----------------------------------------------------------------------
    % THERMAL DIFFUSION
    %----------------------------------------------------------------------
    TEMP_nodes      =    thermal2d_transient(ELEM2NODE_t, PHASES_t, GCOORD, KAPPA, bc_ind_therm, nip_t, TEMP_nodes, dt, elemType, reorder, method, lumping);

    %----------------------------------------------------------------------
    % MECHANICAL SOLVER
    %----------------------------------------------------------------------
    % Simple, transparent, but slow solver
    %     [Vel, Pressure] =    mechanical2d_simple(ELEM2NODE, PHASES, GCOORD, MU, RHO, Ra, G, bc_ind, bc_val, nip);

    % Optimised and efficient solver: MILAMIN
    % [Vel, Pressure] =    mechanical2d(ELEM2NODE, PHASES, GCOORD, MU, (RHO_nodes-TEMP_nodes), Ra, G, bc_ind, bc_val, nip, 'quad', 'amd', 'opt');
    [Vel, Pressure] =    mechanical2d_mu(ELEM2NODE, PHASES, GCOORD, MU_el, (RHO_nodes-TEMP_nodes), Ra, G, bc_ind, bc_val, nip, 'quad', 'amd', 'opt');

    Vx              =    Vel(1:2:end-1)';
    Vy              =    Vel(2:2:end)';

   %----------------------------------------------------------------------
    % ADVECTION: COMPOSITIONAL FIELD - MARKERS, UPDATES MARKER POSITIONS
    %----------------------------------------------------------------------
    Vx_mark         =    sum(SHP.*Vx(ELEM2NODE(:,MARK2ELEM)));
    Vy_mark         =    sum(SHP.*Vy(ELEM2NODE(:,MARK2ELEM)));
    MARK.NODES      =    MARK.NODES + dt*[Vx_mark; Vy_mark];

    %----------------------------------------------------------------------
    % RESEEDING MARKERS 
    %----------------------------------------------------------------------
    [~, ~, MARK2ELEM_t] =    uv_quads(L, [asp*ney_t; ney_t], MARK.NODES);
    Mweights_t          =    accumarray(MARK2ELEM_t',ones(size(MARK2ELEM_t)), [nel_t 1]); % counts how many markers in each thermal elemet
    if min(Mweights_t)<1

        MARK.NODES          =    MARK.NODES - dt*[Vx_mark; Vy_mark];

        MARK                =    reseed2_activ(MARK, asp, L, ney_t, nmy, ELEM2NODE_t, RHO_nodes, FUG_nodes, H_activ_nodes, PHASES_t, NUMMELT_t, TIMEMELT_t);

        % disp('RESEEDING')
        prev_melt_bool      =    MARK.NODES(2,:) > (L(2)-thickness);

        % For mechanical grid - next will be used for interpolating velocities from nodes to markers
        [U, V, MARK2ELEM]   =    uv_quads(L, [asp*ney; ney], MARK.NODES);
        SHP                 =    shp_quad([U; V], nnodel);

        Vx_mark         =    sum(SHP.*Vx(ELEM2NODE(:,MARK2ELEM)));
        Vy_mark         =    sum(SHP.*Vy(ELEM2NODE(:,MARK2ELEM)));
        MARK.NODES      =    MARK.NODES + dt*[Vx_mark; Vy_mark];

    end

    %----------------------------------------------------------------------
    % ADVECTION: THERMAL FIELD - SHOOTING BACK CHARACTERISTICS, UPDATES
    % TEMPERATURE IN THE NODES
    %----------------------------------------------------------------------
    % Feet of characteristics - global and local positions
    GCOORD_char     =    GCOORD - dt*[Vx; Vy];
    [U_char, V_char, IND_char]  =    uv_quads(L, [asp*ney; ney], GCOORD_char);
    SHP_char        =    shp_quad([U_char; V_char],nnodel);

    % Interpolating temperature at the feet of characteristics
    TEMP_nodes      =    sum(SHP_char.*TEMP_nodes(ELEM2NODE(:,IND_char)))';

    %----------------------------------------------------------------------
    % MELTING: UPDATE MARKER COMPOSITION
    %----------------------------------------------------------------------
    melt_bool               =    MARK.NODES(2,:) > (L(2)-thickness); % all the markers with phase 2
    new_melt                =    ~prev_melt_bool & melt_bool;   % 1 everywhere new melt, 0 else
    
    new_RHO                 =    new_melt & (MARK.RHO(:) < RHO(1))';
    new_FUG                  =    new_melt & (MARK.FUG(:) < FUG(1))';
    new_H_activ             =    new_melt & (MARK.H_activ(:) < H_activ(1))';
    MARK.RHO(new_RHO)           =    MARK.RHO(new_RHO) + RHO_inc;
    MARK.FUG(new_FUG)             =    MARK.FUG(new_FUG) + FUG_inc;
    MARK.H_activ(new_H_activ)   =    MARK.H_activ(new_H_activ) + H_activ_inc;
    MARK.NUMMELT            =    MARK.NUMMELT + new_melt;
    MARK.TIMEMELT(new_melt) =    MARK.TIMEMELT(new_melt) * time;
    prev_melt_bool          =    melt_bool;

    %----------------------------------------------------------------------    
    % MARKERS TO THE GRID: UPDATES COMPOSITION IN THE NODES
    % AND PHASES IN THE ELEMENTS
    %----------------------------------------------------------------------
    [U, V, MARK2ELEM]   =    uv_quads(L, [asp*ney; ney], MARK.NODES);
    SHP                 =    shp_quad([U; V], nnodel); % Won't be used until interpolating velocities from nodes to markers

    % Update PHASE in elements
    Mweights            =    accumarray(MARK2ELEM',ones(size(MARK2ELEM)),[nel 1]); % counts how many markers in each mechanical elemet
    PHASEsum            =    accumarray(MARK2ELEM',MARK.PHASES,[nel 1]);
    PHASES              =    round(PHASEsum./Mweights);
    PHASES_t            =    PHASES(THERM2MECH);
    NUMMELTsum          =    accumarray(MARK2ELEM',MARK.NUMMELT,[nel 1]);
    NUMMELT             =    round(NUMMELTsum./Mweights);
    NUMMELT_t           =    NUMMELT(THERM2MECH);
    TIMEMELTsum         =    accumarray(MARK2ELEM',MARK.TIMEMELT,[nel 1]);
    TIMEMELT            =    round(TIMEMELTsum./Mweights);
    TIMEMELT_t          =    TIMEMELT(THERM2MECH);

    % Update RHO in nodes - use linear shapefunctions to avoid overshoots
    [U_t, V_t, MARK2ELEM_t]  =    uv_quads(L, [asp*ney_t; ney_t], MARK.NODES);
    SHP_t               =    shp_quad([U_t; V_t], nnodel_t);
    nodeBINS            =    ELEM2NODE_t(:,MARK2ELEM_t); % Too much memory used here and in the next line - make more efficient!
    Weights_t           =    accumarray(nodeBINS(:), SHP_t(:), [nnod 1]);
    markBINS_RHO        =    SHP_t.*repmat(MARK.RHO,nnodel_t,1);
    RHO_nodes           =    accumarray(nodeBINS(:), markBINS_RHO(:), [nnod 1])./Weights_t;
    markBINS_FUG         =    SHP_t.*repmat(MARK.FUG,nnodel_t,1);
    FUG_nodes       =    accumarray(nodeBINS(:), markBINS_FUG(:), [nnod 1])./Weights_t;
    markBINS_H_activ    =    SHP_t.*repmat(MARK.H_activ,nnodel_t,1);
    H_activ_nodes       =    accumarray(nodeBINS(:), markBINS_H_activ(:), [nnod 1])./Weights_t;

    % Update temperature dependant viscosity in nodes and elements
    MU_perturb                  =    arrayfun(@thermal_profile_perturb_activ_fug, ...
                                              TEMP_nodes, H_activ_nodes, FUG_nodes);
    %MU_perturb                  =    arrayfun(@thermal_profile_perturb_activ_no_fug, TEMP_nodes, H_activ_nodes);
    MU_nodes                    =    MU_z .* MU_perturb;

    % clear MU_el;
    MU2EL                       =    reshape(MU_nodes(ELEM2NODE),nnodel,nel);
    MU_el                       =    10.^(mean(log10(MU2EL),1))';
    % disp('post update')
    % disp(min(MU_nodes))
    % disp(max(MU_nodes));

    
    % total_MU_z = 

    %----------------------------------------------------------------------
    % UPDATE TIME
    %----------------------------------------------------------------------
    time            =    time + dt
    time_save       =    time_save + dt;
    dt              =    4*min(dx./sqrt(Vx.^2+Vy.^2));

    it              =    it+1;
    MASS(it)        =    sum(RHO_nodes);
    times(it)       =    time;

    %----------------------------------------------------------------------
    % SAVE VARIABLES
    %----------------------------------------------------------------------
    if time_save > dt_save
        num_save    =    num_save + 1;
        file_name   =   fullfile("output_files", run_title, "dynamic_vars", ...
            strjoin(["dynamic_vars_", string(num_save), ".mat"], ''));
        save(file_name, "MARK", "TEMP_nodes", "RHO_nodes", "MU_nodes", "H_activ_nodes", "new_melt", "Vx", "Vy", "time")
        time_save   =    0;
    end

    
   %  %----------------------------------------------------------------------
   %  % PLOTTING
   %  %----------------------------------------------------------------------
      figure(fnum);clf; hold on;
      ncorners        =    4;
      xplot           =    reshape(GCOORD(1,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
      yplot           =    reshape(GCOORD(2,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
      % plot           =    reshape(TEMP_nodes(ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
      % plot          =    log10(reshape(MU_nodes(ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t));
      % plot           =    reshape(H_activ_nodes(ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t); 
      % mp              =    patch(xplot,yplot,plot);
      % set(mp,'EdgeColor','n','FaceColor','interp');

     %  plot(MARK.NODES(1,abs(MARK.RHO)>0), MARK.NODES(2,abs(MARK.RHO)>0),'.k') % plotting melt_bool markers
     %  plot(MARK.NODES(1,MARK.PHASES==1), MARK.NODES(2,MARK.PHASES==1),'.k') % plotting melt_bool markers
     scatter(MARK.NODES(1,:), MARK.NODES(2,:), [], MARK.NUMMELT(:),'Filled')  % plotting new melt

      axis equal; axis off;
      colormap('jet'); colorbar; % caxis([-.2 1]);
      % colormap('jet')
      % colorbar
      % clim([0 20])
      quiver(GCOORD(1,1:4:end),GCOORD(2,1:4:end),Vx(1:4:end),Vy(1:4:end),2,'k')
      drawnow

    % str = '/Users/emmahavens/Library/CloudStorage/OneDrive-NorthwesternUniversity/Emma_Summer_Research_Project/Emmas_Continents';
    % figure(fnum); print('-dpng','-r150',[str 'glideV_' num2str(it)]);

clear Vx Vy GCOORD_char IND_char U_char V_char SHP_char Vx_mark Vy_mark MU_perturb MU_nodes

end

% figure(fnum+100);
% plot(times,MASS./MASS(1),'-k','LineWidth',3)
% xlabel('Time'); ylabel('Mass / Initial Mass')


% MELTING ON MARKERS
% VISCOSITY - RECONSIDER PHASES?






