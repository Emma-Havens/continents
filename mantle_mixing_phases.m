
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
Bu                  =     -5e-1;  % Compositional buoyancy ratio: rho_crust / (~200 kg/m^3); use negative values for chemically light, positive values for chemically anom
G                   =   [0 -1];  % Gravity vector [Gx Gy]
KAPPA               =   [1; 1];  % Thermal diffusivity
MU                  =   [10; 1];  % Viscosity
MU_H2O              =   [1; 1];  % Amount of stiffening from dehydration
RHO                 =  [Bu; 0];  % Compositional density anomaly
depth               =     2900;
time_final          =     6e-2;  % Total simulation time
thickness           = 50/depth;  % Thickness of the compositionally anomalous layer
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
ney                 =    40; % Number of elements in vertical direction
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

%--------------------------------------------------------------------------
% MARKERS
%--------------------------------------------------------------------------
nmy                 =    5*ney;
nm                  =    asp*nmy*nmy;
dm                  =    2*L(2)/nmy;
Mx                  =    linspace(-(L(1)-dm), L(1)-dm, asp*nmy);
My                  =    linspace(-(L(2)-dm), L(2)-dm, nmy);
[MX, MY]            =    ndgrid(Mx,My);
MARK.NODES          =    [MX(:)'; MY(:)'];
MARK.PHASES         =    2*ones(1,nm); % phase of each marker, has values 1 or 2
MARK.NUMMELT        =    zeros(1,nm);
prev_melt_bool      =    MARK.NODES(2,:) > (L(2)-thickness); % 1 if within melt region, 0 else CHECK LINE
MARK.TIMEMELT       =    zeros(1,nm);

anom                =    MARK.NODES(2,:) > (L(2)-thickness); % all the markers with phase 2
MARK.PHASES(anom)   =    1;
MARK.RHO            =    RHO(MARK.PHASES)';
MARK.MU             =    MU(MARK.PHASES)' .* MU_H2O(MARK.PHASES)';

% figure(1); hold on;
% plot(MARK.NODES(1,:),MARK.NODES(2,:),'.k');
% plot(MARK.NODES(1,anom),MARK.NODES(2,anom),'om');
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
    markBINS_MU         =    SHP_t.*repmat(MARK.MU,nnodel_t,1);
    MUPHASE_nodes       =    accumarray(nodeBINS(:), markBINS_MU(:), [nnod 1])./Weights_t;


%     figure(2); hold on;
% scatter(GCOORD(1,:),GCOORD(2,:),50,MUPHASE_nodes,'s','filled','MarkerEdgeColor','k')
% scatter(MARK.NODES(1,:),MARK.NODES(2,:),25,MARK.MU,'filled','MarkerEdgeColor','k')
% scatter(GCOORD(1,:),GCOORD(2,:),50,RHO_nodes,'s','filled','MarkerEdgeColor','k')
% scatter(MARK.NODES(1,:),MARK.NODES(2,:),25,MARK.RHO,'filled','MarkerEdgeColor','k')

colorbar;
axis equal

%--------------------------------------------------------------------------
% TEMPERATURE DEPENDANT VISCOSITY
%--------------------------------------------------------------------------

MU_z                        =    ones(nnod,1);
MU_z(GCOORD(2,:) < Z_660)   =    100 * MU_z(GCOORD(2,:) < Z_660);
MU_perturb                  =    arrayfun(@thermal_profile_perturb, TEMP_nodes);
MU_nodes                    =    MU_z .* MU_perturb .* MUPHASE_nodes;

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

% MU_el                                   =    ones(nel,1);
% MU_el(GCOORD(2,ELEM2NODE(9,:)) > Z_660) =    100 * MU_el(GCOORD(2,ELEM2NODE(9,:)) > Z_660);
% MU_el_old                       =    zeros(nel,1);
% for row = 1:nel
%     MU_el_old(row)                  =    sum(MU_nodes(ELEM2NODE(:,row),1)) / nnodel;
% end


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
    run_title = strjoin(["Bu_", string(Bu), "_Mu_", ...
        string(MU(1)), "_", string(MU(2)), "_MuH2O_", string(MU_H2O(1)), ...
        "_", string(MU_H2O(2)), "_ney_", string(ney), ...
        "_", string(num_sim)], '');
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
save(static_file_name, "Ra", "Bu", "MU", "RHO", "thickness", ...
    "time_final", "dt_save", "asp", "ney", "L")

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

        MARK                =    reseed2(MARK, asp, L, ney_t, nmy, ELEM2NODE_t, RHO, MU, MU_H2O, RHO_nodes, MUPHASE_nodes, PHASES_t, NUMMELT_t, TIMEMELT_t);

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
    anom                =    MARK.NODES(2,:) > (L(2)-thickness); % all the markers with phase 2
    MARK.PHASES(anom)   =    1;
    MARK.RHO            =    RHO(MARK.PHASES)';
    MARK.MU             =    MU(MARK.PHASES)' .* MU_H2O(MARK.PHASES)';

    % melt_bool           =    MARK.NODES(2,:) > (L(2)-thickness);
    new_melt                =    ~prev_melt_bool & anom;   % 1 everywhere new melt, 0 else
    MARK.TIMEMELT(new_melt) =    MARK.TIMEMELT(new_melt) * time;
    MARK.NUMMELT            =    MARK.NUMMELT + new_melt;
    prev_melt_bool          =    anom;

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
    markBINS_MU         =    SHP_t.*repmat(MARK.MU,nnodel_t,1);
    MUPHASE_nodes       =    accumarray(nodeBINS(:), markBINS_MU(:), [nnod 1])./Weights_t;

    % Update temperature dependant viscosity in nodes and elements
    MU_perturb                  =    arrayfun(@thermal_profile_perturb, TEMP_nodes);
    MU_nodes                    =    MU_z .* MU_perturb .* MUPHASE_nodes;

    % clear MU_el;
    MU2EL                       =    reshape(MU_nodes(ELEM2NODE),nnodel,nel);
    MU_el                       =    10.^(mean(log10(MU2EL),1))';
    % disp('post update')
    % disp(min(MU_nodes))
    % disp(max(MU_nodes));

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
        save(file_name, "MARK", "TEMP_nodes", "RHO_nodes", "MU_nodes", "new_melt", "Vx", "Vy", "time")
        time_save   =    0;
    end

    
   %  %----------------------------------------------------------------------
   %  % PLOTTING
   %  %----------------------------------------------------------------------
     %  figure(fnum);clf; hold on;
     %  ncorners        =    4;
     %  xplot           =    reshape(GCOORD(1,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
     %  yplot           =    reshape(GCOORD(2,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
     %  % qplot           =    reshape(TEMP_nodes(ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
     %  muplot          =    reshape(MU_nodes(ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
     %  % hp              =    patch(xplot,yplot,qplot);    % plotting temp
     %  mp              =    patch(xplot,yplot,muplot);
     %  % set(hp,'EdgeColor','n','FaceColor','interp');
     % 
     % %  plot(MARK.NODES(1,abs(MARK.RHO)>0), MARK.NODES(2,abs(MARK.RHO)>0),'.k') % plotting anom markers
     % %  plot(MARK.NODES(1,MARK.PHASES==1), MARK.NODES(2,MARK.PHASES==1),'.k') % plotting anom markers
     % % scatter(MARK.NODES(1,:), MARK.NODES(2,:), [], MARK.NUMMELT(:),'Filled')  % plotting new melt
     % % scatter(MARK.NODES(1,new_melt==1),MARK.NODES(2,new_melt==1),'.r')
     % 
     %  axis equal; axis off;
     %  colormap('jet'); colorbar; % caxis([-.2 1]);
     %  % colormap('jet')
     %  % colorbar
     %  % clim([0 20])
     %  quiver(GCOORD(1,1:4:end),GCOORD(2,1:4:end),Vx(1:4:end),Vy(1:4:end),2,'k')
     %  drawnow

    % str = '/Users/emmahavens/Library/CloudStorage/OneDrive-NorthwesternUniversity/Emma_Summer_Research_Project/Emmas_Continents';
    % figure(fnum); print('-dpng','-r150',[str 'glideV_' num2str(it)]);

clear Vx Vy GCOORD_char IND_char U_char V_char SHP_char Vx_mark Vy_mark MU_perturb MU_nodes

end

% figure(fnum+100);
% plot(times,MASS./MASS(1),'-k','LineWidth',3)
% xlabel('Time'); ylabel('Mass / Initial Mass')


% MELTING ON MARKERS
% VISCOSITY - RECONSIDER PHASES?






