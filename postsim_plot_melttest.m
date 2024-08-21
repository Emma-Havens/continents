% plots and analyzes output from mantle_mixing.m

% create foldwer
% LOOK AT num2string
% Ex. var = load("file.mat", "var")

% save("static_var.mat", "Ra", "Bu", "MU", "RHO", "thickness",
% "time_final", "dt_save", "asp", "ney", "L")

% save("dynamic_file.mat", "MARK", "TEMP_nodes", "RHO_nodes", "Vx", "Vy", "time")

clear variables
fnum = 1;

outer_folder = "output_files/NUMMELT_test";
sim_folders = dir(outer_folder);
sim_folders = sim_folders(~ismember({sim_folders(:).name},{'.','..', '.DS_Store'}));
sim_folders_len = length(sim_folders);

for i = 1:sim_folders_len

    disp("next folder");

    clear Ra Bu MU RHO thickness time_final dt_save num_save asp ney L ney_t nel_t ELEM2NODE_t ind_t

    run_title = sim_folders(i).name;



    % run_title = "Bu_-2_Mu_1_1_H_activ_375_320_asp_1_ney_100_3";
    
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
        % MU_inc              =   file_contents.MU_inc;
        % H_activ             =   file_contents.H_activ;
        % H_activ_inc         =   file_contents.H_activ_inc;
        RHO                 =   file_contents.RHO;
        % RHO_inc             =   file_contents.RHO_inc;
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
    nnod                =    size(GCOORD,2);
    
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
        % TEMP_nodes      =    file_contents.TEMP_nodes;
        % RHO_nodes       =    file_contents.RHO_nodes;
        % MU_nodes        =    file_contents.MU_nodes;
        % H_activ_nodes   =    file_contents.H_activ_nodes;
        Vx              =    file_contents.Vx;
        Vy              =    file_contents.Vy;
        % time            =    file_contents.time;
        % new_melt        =    file_contents.new_melt;
    
        %
        % PLOT MARKERS, TEMPERATURE, NODE VELOCITY
        %
    
        figure(fnum);clf; hold on;
          %ncorners        =    4;
          %xpatch           =    reshape(GCOORD(1,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
          %ypatch           =    reshape(GCOORD(2,ELEM2NODE_t(1:ncorners,:)),ncorners,nel_t);
          %tpatch           =    repmat(NUMMELT_t', 4, 1);
          %mp              =    patch(xpatch,ypatch,tpatch);
          %set(mp,'EdgeColor','n','FaceColor','interp');

          xplot           =    repmat([-L(1),-L(1),L(1),L(1)], 1, nnod);
          yplot           =    repelem(GCOORD(1,:),2);
          %plot(xplot(2:nnod*2+1), yplot,'k','LineWidth',.5);
          %plot(yplot, xplot(2:nnod*2+1),'k','LineWidth',.5);
          plot([-L(1), L(1), L(1), -L(1), -L(1)],L(2)*[-1, -1, 1, 1, -1],'k','LineWidth',1);
          plot([-L(1), L(1)], [L(2)-thickness, L(2)-thickness], 'k', 'LineWidth',2);
          scatter(MARK.NODES(1, MARK.NUMMELT()>0), MARK.NODES(2,MARK.NUMMELT()>0), [], MARK.NUMMELT(MARK.NUMMELT()>0), 'Filled')  % plotting new melt

          axis equal; axis off;
          colormap('jet'); colorbar; % caxis([-.2 1]);

          clim([0 max(MARK.NUMMELT)])
          quiver(GCOORD(1,1:4:4),GCOORD(2,1:4:4),Vx(1:4:4),Vy(1:4:4),0.2,'k')
          drawnow
    
    
         %str = '/Users/emmahavens/Library/CloudStorage/OneDrive-NorthwesternUniversity/Emma_Summer_Research_Project/Emmas_Continents/plotted_runs/';
         %figure(fnum); print('-dpng','-r150',[str 'Bu_-2_Mu_1_1_H_activ_375_320_asp_1_ney_100_3_' num2str(j) '.png']);
    
    
    end

end

disp('done')

% end

%
% PLOT TOTAL MASS
%

% figure(fnum+100);
% plot(times,MASS./MASS(1),'-k','LineWidth',3)
% xlabel('Time'); ylabel('Mass / Initial Mass')


