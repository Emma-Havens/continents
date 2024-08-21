% plots and analyzes output from mantle_mixing.m

% create foldwer
% LOOK AT num2string
% Ex. var = load("file.mat", "var")

% save("static_var.mat", "Ra", "Bu", "MU", "RHO", "thickness",
% "time_final", "dt_save", "asp", "ney", "L")

% save("dynamic_file.mat", "MARK", "TEMP_nodes", "RHO_nodes", "Vx", "Vy", "time")

clear variables
fnum = 1;

% sim_runs = dir("output_files_quest");
% sim_runs = sim_runs(~ismember({sim_runs(:).name},{'.','..', '.DS_Store'}));
% sim_runs_len = length(sim_runs);
% 
% for i = 1:sim_runs_len
% 
%     clear Ra Bu MU RHO thickness time_final dt_save num_save asp ney L ney_t nel_t ELEM2NODE_t ind_t

    % run_title = sim_runs(i).name;

run_title = "Bu_-2+0_Fug_1+0_H_450+20_thick_0.017_asp_1_ney_50_0";
outer_folder = "output_files_quest/2ab_Mar_2_fug_v_H_lo_res";

%
% LOAD STATIC VARIABLES
%

static_file         =   fullfile(outer_folder, run_title, "static_var.mat");
static_loaded       =   true;

try
    file_contents       =   load(static_file);
catch
    static_loaded       =   false;
    disp("unable to load static_var.mat")

    % physical param
    Ra                  =   1e5;
    Bu                  =   -2;
    FUG                  =   [1; 1];
    RHO                 =   [Bu; 0];
    thickness_exact     =   50/2900;
    time_final          =   6e-2;
    dt_save             =   1e-3;         
    
    % mesh param
    asp                 =   1;
    ney                 =   140;
    L                   =   [asp*.5;.5];
end

if static_loaded
    % physical param
    Ra                  =   file_contents.Ra;
    Bu                  =   file_contents.Bu;
    FUG                 =   file_contents.FUG;
    FUG_inc             =   file_contents.FUG_inc;
    H_activ             =   file_contents.H_activ;
    H_activ_inc         =   file_contents.H_activ_inc;
    RHO                 =   file_contents.RHO;
    RHO_inc             =   file_contents.RHO_inc;
    thickness_exact     =   file_contents.thickness_exact;
    time_final          =   file_contents.time_final;
    dt_save             =   file_contents.dt_save;         
    
    % mesh param
    asp                 =   file_contents.asp;
    ney                 =   file_contents.ney;
    L                   =   file_contents.L;
end

% calculate num_save
dynamic_folder =   fullfile(outer_folder, run_title, "dynamic_vars");
sim_runs = dir(dynamic_folder);
sim_runs = sim_runs(~ismember({sim_runs(:).name},{'.','..', '.DS_Store'}));
num_save = length(sim_runs);

% mechanical mesh vars
nnodel              =    9;  % Number of GCOORD per element

% Generating mechanical mesh
[GCOORD, ELEM2NODE, ind] = generate_quad_grid(asp, ney, nnodel, L);

% thermal mesh vars
nnodel_t            =    4;
ney_t               =    2*ney;

% Generating thermal mesh
[~, ELEM2NODE_t, ind_t] = generate_quad_grid(asp, ney_t, nnodel_t, L);

nel_t               =    size(ELEM2NODE_t,2);


for j = 1:num_save

    dynamic_file       =   fullfile(outer_folder, run_title, ...
        "dynamic_vars", strjoin(["dynamic_vars_", string(j), ".mat"], ''));
    file_contents      =   load(dynamic_file);

    %
    % LOAD DYNAMIC VARIABLES
    %
    MARK            =    file_contents.MARK;
    TEMP_nodes      =    file_contents.TEMP_nodes;
    RHO_nodes       =    file_contents.RHO_nodes;
    MU_nodes        =    file_contents.MU_nodes;
    % H_activ_nodes   =    file_contents.H_activ_nodes;
    Vx              =    file_contents.Vx;
    Vy              =    file_contents.Vy;
    time            =    file_contents.time;
    new_melt        =    file_contents.new_melt;


    nnod                =    size(GCOORD,2);
    [U_t, V_t, MARK2ELEM_t]  =    uv_quads(L, [asp*ney_t; ney_t], MARK.NODES);
    SHP_t               =    shp_quad([U_t; V_t], nnodel_t);
    nodeBINS            =    ELEM2NODE_t(:,MARK2ELEM_t); % Too much memory used here and in the next line - make more efficient!
    Weights_t           =    accumarray(nodeBINS(:), SHP_t(:), [nnod 1]);
    markBINS_NUMM       =    SHP_t.*repmat(MARK.NUMMELT,nnodel_t,1);
    NUMMELT_nodes       =    accumarray(nodeBINS(:), markBINS_NUMM(:), [nnod 1])./Weights_t;
    markBINS_H_activ    =    SHP_t.*repmat(MARK.H_activ,nnodel_t,1);
    H_activ_nodes       =    accumarray(nodeBINS(:), markBINS_H_activ(:), [nnod 1])./Weights_t;

    nel                 =    size(ELEM2NODE,2);
    dx                  =    4*L(1)*L(2)/nel;
    dt              =    4*min(dx./sqrt(Vx.^2+Vy.^2));
    % disp("dt:" + string(dt))
    % disp("Vx:" + max(Vx))
    % disp("Vy:" + max(Vy))
    % disp("MaxMU:" + max(MU_nodes))
    disp("MinMU:" + min(MU_nodes))

    % GOOD
    % dt:0.0010933
    % Vx:0.84603
    % Vy:1.0082
    % MaxMU:74.6464
    % MinMU:0.30693
    % dt:0.00074275
    % Vx:1.258
    % Vy:1.477
    % MaxMU:74.6464
    % MinMU:0.3087
    % dt:0.00052853
    % Vx:1.7776
    % Vy:2.066
    % MaxMU:74.6464
    % MinMU:0.31032
    % dt:0.00038779
    % Vx:2.4322
    % Vy:2.8084
    % MaxMU:74.6464
    % MinMU:0.31117

    % BAD
    % dt:0.00013387
    % Vx:9.4198
    % Vy:11.9523
    % MaxMU:23.6053
    % MinMU:0.068999
    % dt:9.0897e-05
    % Vx:13.6087
    % Vy:17.6024
    % MaxMU:23.6053
    % MinMU:0.069541
    % dt:6.1713e-05
    % Vx:19.5548
    % Vy:25.9263
    % MaxMU:17.7402
    % MinMU:0.022206

    % WORSE
    %   dt:0.02542
    % Vx:0.041162
    % Vy:0.051288
    %   dt:0.1375
    % Vx:0.0042333
    % Vy:0.0070281
    %   dt:0.46449
    % Vx:0.0016839
    % Vy:0.0022406
    %   dt:0.70425
    % Vx:0.0013482
    % Vy:0.0015194

    %
    % PLOT MARKERS, TEMPERATURE, NODE VELOCITY
    %

    figure(fnum);clf;
    t = tiledlayout(2,2);
    % t = tiledlayout(1,2);
    t.Padding = 'compact';
    t.TileSpacing = 'compact';

    % PLOT TEMPERATURE

    nexttile;
    hold on;
    ncorners        =    4;
    xplot           =    reshape(GCOORD(1,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
    yplot           =    reshape(GCOORD(2,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
    qplot           =    reshape(TEMP_nodes(ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
    hp              =    patch(xplot,yplot,qplot);
    set(hp,'EdgeColor','n','FaceColor','interp');

    axis equal; axis off;
    colormap('jet');  colorbar; % caxis([-.2 1]);
    %quiver(GCOORD(1,1:4:end),GCOORD(2,1:4:end),Vx(1:4:end),Vy(1:4:end),2,'k')
    hold off;

    % PLOT VISCOSITY

    nexttile;
    hold on;
    ncorners        =    4;
    xplot           =    reshape(GCOORD(1,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
    yplot           =    reshape(GCOORD(2,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
    muplot          =    log10(reshape(MU_nodes(ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t));
    mp              =    patch(xplot,yplot,muplot);
    set(mp,'EdgeColor','n','FaceColor','interp');

    axis equal; axis off;
    colormap('jet');  colorbar; % caxis([-.2 1]);
    clim([-1 3])
    %quiver(GCOORD(1,1:4:end),GCOORD(2,1:4:end),Vx(1:4:end),Vy(1:4:end),2,'k')
    hold off;

    % PLOT MARK.NUMMELT

    nexttile;
    hold on;
    % scatter(MARK.NODES(1,MARK.NUMMELT>0), MARK.NODES(2,MARK.NUMMELT>0), [], MARK.NUMMELT(MARK.NUMMELT>0), 'Filled') % plotting nummelt

    ncorners        =    4;
    xplot           =    reshape(GCOORD(1,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
    yplot           =    reshape(GCOORD(2,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
    muplot          =    reshape(NUMMELT_nodes(ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
    mp              =    patch(xplot,yplot,muplot);
    set(mp,'EdgeColor','n','FaceColor','interp');

    axis equal; axis off;
    colormap('jet')
    colorbar
    clim([0 3])
    %quiver(GCOORD(1,1:4:end),GCOORD(2,1:4:end),Vx(1:4:end),Vy(1:4:end),2,'k')
    hold off;

    % PLOT MARK.H_activ

    nexttile;
    hold on;
    scatter(MARK.NODES(1,MARK.NUMMELT>1), MARK.NODES(2,MARK.NUMMELT>1), [], MARK.NUMMELT(MARK.NUMMELT>1), 'Filled') % plotting nummelt

    ncorners        =    4;
    xplot           =    reshape(GCOORD(1,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
    yplot           =    reshape(GCOORD(2,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
    muplot          =    reshape(H_activ_nodes(ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
    %mp              =    patch(xplot,yplot,muplot);
    %set(mp,'EdgeColor','n','FaceColor','interp');

    axis equal; axis off;
    colormap('jet')
    colorbar
    clim([320 375])
    quiver(GCOORD(1,1:4:end),GCOORD(2,1:4:end),Vx(1:4:end),Vy(1:4:end),2,'k')
    hold off;
    drawnow

    % if time > 0.07
    %     disp('breakpoint');
    % end


     %str = '/Users/emmahavens/Library/CloudStorage/OneDrive-NorthwesternUniversity/Emma_Summer_Research_Project/Emmas_Continents/plotted_runs/';
     %figure(fnum); print('-dpng','-r150',[str 'Bu_-2+0_Fug_0.5+0.25_H_370+0_thick_0.017_asp_1_ney_50_0_' num2str(j) '.png']);


end

disp('done')

% end

%
% PLOT TOTAL MASS
%

% figure(fnum+100);
% plot(times,MASS./MASS(1),'-k','LineWidth',3)
% xlabel('Time'); ylabel('Mass / Initial Mass')


