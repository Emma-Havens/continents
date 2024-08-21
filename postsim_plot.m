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


run_title = "Bu_-0.5_Mu_10_1_MuH2O_1_1_ney_100_0";
outer_folder = "output_files_quest/Feb_8";

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
    Bu                  =   -5e-1;
    MU                  =   [1; 1];
    RHO                 =   [Bu; 0];
    thickness           =   10/600;
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
    MU                  =   file_contents.MU;
    % MU_H2O              =   file_contents.MU_H2O;
    RHO                 =   file_contents.RHO;
    thickness           =   file_contents.thickness;
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
    % MUPHASE_nodes   =    file_contents.MUPHASE_nodes;
    Vx              =    file_contents.Vx;
    Vy              =    file_contents.Vy;
    time            =    file_contents.time;
    new_melt        =    file_contents.new_melt;

    %
    % PLOT MARKERS, TEMPERATURE, NODE VELOCITY
    %

    % PLOT TEMPERATURE

    figure(fnum);clf; hold on;
    ncorners        =    4;
    xplot           =    reshape(GCOORD(1,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
    yplot           =    reshape(GCOORD(2,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
    qplot           =    reshape(TEMP_nodes(ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
    % muplot          =    log10(reshape(MU_nodes(ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t));
    hp              =    patch(xplot,yplot,qplot);
    % mp              =    patch(xplot,yplot,muplot);
    set(hp,'EdgeColor','n','FaceColor','interp');

   % plot(MARK.NODES(1,abs(MARK.RHO)>0), MARK.NODES(2,abs(MARK.RHO)>0),'.k') % plotting anom markers
   % plot(MARK.NODES(1,MARK.PHASES==1), MARK.NODES(2,MARK.PHASES==1),'.k') % plotting anom markers

    axis equal; axis off;
    colormap('jet');  colorbar; % caxis([-.2 1]);
    quiver(GCOORD(1,1:4:end),GCOORD(2,1:4:end),Vx(1:4:end),Vy(1:4:end),2,'k')
    drawnow

    % PLOT VISCOSITY

   %  figure(fnum);clf; hold on;
   %  ncorners        =    4;
   %  xplot           =    reshape(GCOORD(1,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
   %  yplot           =    reshape(GCOORD(2,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
   %  muplot          =    log10(reshape(MU_nodes(ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t));
   %  mp              =    patch(xplot,yplot,muplot);
   %  set(mp,'EdgeColor','n','FaceColor','interp');
   % 
   % % plot(MARK.NODES(1,abs(MARK.RHO)>0), MARK.NODES(2,abs(MARK.RHO)>0),'.k') % plotting anom markers
   % % plot(MARK.NODES(1,MARK.PHASES==1), MARK.NODES(2,MARK.PHASES==1),'.k') % plotting anom markers
   % 
   %  axis equal; axis off;
   %  colormap('jet');  colorbar; % caxis([-.2 1]);
   %  quiver(GCOORD(1,1:4:end),GCOORD(2,1:4:end),Vx(1:4:end),Vy(1:4:end),2,'k')
   %  drawnow

    % PLOT MARK.NUMMELT

   %  figure(fnum);clf; hold on;
   %  scatter(MARK.NODES(1,MARK.NUMMELT>0), MARK.NODES(2,MARK.NUMMELT>0), [], MARK.NUMMELT(MARK.NUMMELT>0), 'Filled') % plotting nummelt
   % % scatter(MARK.NODES(1,new_melt==1),MARK.NODES(2,new_melt==1),'.r')
   % 
   %  axis equal; axis off;
   %  colormap('jet')
   %  colorbar
   %  clim([0 10])
   %  quiver(GCOORD(1,1:4:end),GCOORD(2,1:4:end),Vx(1:4:end),Vy(1:4:end),2,'k')
   %  drawnow


     % str = '/Users/emmahavens/Library/CloudStorage/OneDrive-NorthwesternUniversity/Emma_Summer_Research_Project/Emmas_Continents/plotted_runs/';
     % figure(fnum); print('-dpng','-r150',[str 'Bu_-0.5_MU_10_visc' num2str(j) '.png']);


end

disp('done')

% end

%
% PLOT TOTAL MASS
%

% figure(fnum+100);
% plot(times,MASS./MASS(1),'-k','LineWidth',3)
% xlabel('Time'); ylabel('Mass / Initial Mass')


