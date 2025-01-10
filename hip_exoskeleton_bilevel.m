clc; close all; clear

currentFolder = pwd;
addpath(genpath(currentFolder));

ExampleFolder= 'DataExample\Data2392';
DataFolder   = 'DataDigitalized';
%% TRC file
% You need a TRC file. To have a compatible version with Matlab, you need
% to 1) open TRC file and 2) save it as xlsx. You will open the xlsx version

% read markers
full_dir=fullfile(currentFolder,ExampleFolder,'P1S100W01_trc_converted.xlsx');

table_read=importdata(full_dir);

table_header=table_read.textdata(4,3:end);  % values set manually after checking the file
trc_header=table_header(~cellfun('isempty',table_header));
trc_data  =table_read.data(7:end,3:end);    % values set manually after checking the file
trc_time  =table_read.data(7:end,2);

numMakers  =length(trc_header); 
data_length=size(trc_data,1);

if numMakers == (size(trc_data,2))/3; else; disp('numColumns disagree') ;end % check #headers = #dataColumns 

% Convert to X horizontal, Y vertical, and Z, depth, and units into meters [m]
trc_data_conv=zeros(size(trc_data));
new_order=[3 2 1]; % first is X, second is Y, and thrid is Z
for i=1:numMakers
    ind        = ((i-1)*3 +1) + (new_order-1);
    trc_data_conv(:,((i-1)*3 +1)+[0 1 2])=trc_data(:,ind);
end

trc_data_conv=trc_data_conv/1000; % [mm] -> [m]

% Select markers you are interested for visualization purposes
markers_list = {'RASI' 'LASI' 'RPSI' 'LPSI' ... % physical markers
                'RKNE' 'LKNE' ...
                'RANK' 'LANK' 'RFMH' 'RVMH' 'RTOE' 'LFMH' 'LVMH' 'LTOE' 'RHEE' 'LHEE' ...
                'RHJC_cgm24' 'LHJC_cgm24' ... % model markers
                'RKJC_cgm24' 'LKJC_cgm24' ...
                'RAJC_cgm24' 'LAJC_cgm24' ...
                };

numMakers_list=length(markers_list);
marker_ind=zeros(1,numMakers_list);
for i=1:numMakers_list
    marker_ind(i)=find(strcmp(trc_header,markers_list(i)));
end

% Plot markers
toPlot_markerTrajectories=1;

if toPlot_markerTrajectories==1
    axis_choosen= [1 2]; % these numbers correspond to indexes aligned to X (horizontal) and Y (vertical) axes

    clf;
    subplot(1,1,1)
    for i=1:10:data_length
        for j=1:numMakers_list
            marker_sel = marker_ind(j);
            ind        = ((marker_sel-1)*3+1) + (axis_choosen-1);

            if strcmp(markers_list{j}(1),'L');    color_marker='r';                  else; color_marker='g'; end
            if contains(markers_list{j},'cgm24'); type_marker='s'; color_marker='k'; else; type_marker='o'; end 

            plot(trc_data_conv(i,ind(1)),trc_data_conv(i,ind(2)),type_marker,'MarkerSize',6,'MarkerFaceColor',color_marker,'MarkerEdgeColor','none');
            hold on;

        end
        axis([-2.3 2.7 -0.025 2.00])
        hold off;
        title(['marker trajectories t=' num2str(trc_time(i),'%1.2f')  's'])
        xlabel('horizontal (m)'); ylabel('vertical (m)'); 
        pause(0.05)
    end
end
%% Input information

% paths and folders
OutPath        = fullfile(currentFolder,'Results_bilevel');  % NEW FOR BILEVEL!

% experimental data
Misc.IKfile = {fullfile(currentFolder,ExampleFolder,'IK','IK.mot')};
Misc.IDfile = {fullfile(currentFolder,ExampleFolder,'ID','ID_solution.sto')};
model_path  =  fullfile(currentFolder,ExampleFolder,'1_scaledH_hammer2010_MIF2.osim');

% select the leg's side
Misc.side_sel      = 'r';

% input dofs - select the DOFs you want to include in the optimization
Misc.DofNames_Input={['hip_flexion_' Misc.side_sel] ['hip_adduction_' Misc.side_sel] ['hip_rotation_' Misc.side_sel] ['knee_angle_' Misc.side_sel]  ['ankle_angle_' Misc.side_sel]}; 

% ranges
time = [0.57 1.75];

% Time of toe_off helps set the limits of the timing for exoskeleton
% asssistance
Misc.toe_off_timing= 1.31;

% Original FMo values from Rajagopal might be to high to see a muscle
% activation reduction. I decreased it here
% FMo Maximum isometric force
Misc.updParams.FMo.names      ={'soleus' 'med_gas' 'lat_gas'};
Misc.updParams.FMo.values     =[4000 3000 1000]; % i.e., these are the new values

% label of muscle-tendon parameters
param_label={'FMo' 'lMo' 'lTs' 'alphao' 'vMmax' 'kpe' 'so' 'sM'};
Misc.param_label   =param_label;

% run muscle analysis
Misc.GetAnalysis = 1; % RUN ONCE IF THERE IS NEW OutPath@

% We will use previously optimized lMo, lTs, kT, and passive force-length
% curves. Run "example_fiberLength_PassCal" to have access to the optimized values
Misc.Param_read    = {'model'};
% Misc.myParams=labelParams(Results_all.kT.calibratedMRS,Results_all.params.calibratedMRS,param_label); % This enables the tuning (I will work on this later)
Misc.muscleFiberCal= 0; % This enables the tuning (I will work on this later)

% Weight terms in the objective function (otherwise it uses default values)
Misc.wEMG   =    0.00; 
Misc.wAct   =    1.00; 
Misc.wTres  = 1000.00; 
Misc.wVm    =    0.01; 

% WORKFLOW SETUP -----------------------------------------------------------
% Normal : no assistive moment
workFlow_Normal      = 1;

% setting act-force generation strength, normal ---------------------------
Misc.forGen.actFM.names       ={};
Misc.forGen.actFM.percentages =[]; % all muscles are 100% strength by default

Misc.Exo_Enable    = 0;
if workFlow_Normal== 1
    Misc.OutName= 'normal';
    [Results_normal,DatStore,Misc] = solveMuscleRedundancy_ExoCal(model_path,time,OutPath,Misc);
end

[J_normalized,Edot_normal_TS] = computeOuterLoopFunction(Misc,Results_normal); % this is my starting point, the Edot at unassisted conditions.
%%
% -------------------------------------------------------------------------
% Everything before is as in "hip_exoskeleton.m"
% The interesting part starts now....
% -------------------------------------------------------------------------
%% Bilevel formulation
% Here we will use Bayesian optimization. According to my previous study,
% Bayesian converges faster than Particle Swarm optimization and Genetic
% algorithm. I chose Bayesian for this implementation
% Here are the main components
% 
% - OUTER LEVEL
% OPT=Bayesian Optimization, AIM=min(E)
% - INNER LEVEL
% OPT=Direct collocation, AIM=min(a^2)
% - PARAMETRIZATION
% PGM: stiffness, start time, duration 
%% Loading
% run muscle analysis
Misc.GetAnalysis = 0;

Misc.model_path= model_path;
Misc.time      = time;
Misc.OutPath   = OutPath;
Misc.trc_time  = trc_time;
%% Actuator characteristics, not affected by control
% Geometry does not change based on the actuator assistive magnitude, it
% can be computed before.
% SETUP FOR EXO SIMULATION
Misc.Exo_Enable    = 1; % exo enabled
Misc.Exo_Mode="manual"; Misc.Exo_Type="active";
Misc.Exo_Dof={['hip_flexion_' Misc.side_sel]}; Misc.Exo_group=+1; % hip flexion + / hip extension -
Misc.OutName= 'hipExo_PGM';

% SETUP FOR ACTUATOR
% user characteristics
PGMparam.user.proximal_marker_reference = 'RASI';
PGMparam.user.proximal_marker_offset    = [0.0 0.0]; % [m] - respect to marker reference (I cannot be any other value)
PGMparam.user.distal_marker_reference   = 'RKNE';
PGMparam.user.distal_marker_offset      = [0.0 0.0]; % [m] - respect to marker reference (I cannot be any other value)
PGMparam.user.HJC_marker_reference      = 'RHJC_cgm24';

% instrinsic actuator properties
PGMparam.actuator.slack_length = 0.20; % [m]
PGMparam.actuator.time_lag     = 0.01; % [s] (tau)

opt_visual.flag        = 0;
opt_visual.marker_list = markers_list;
PGMparam               = DEVactuation_geometry(PGMparam,trc_header,trc_time,trc_data_conv,markers_list,opt_visual); % compute actuator: length and moment arm

%% Framework bilevel

% DEFINITION OF VARIABLES
% name      K    t_s   d  
% units     N/m  %GC  %GC
lb_list=[   1    50    5]; % lower bound
ub_list=[1500    70   25]; % upper bound, to be aware: MAX(t_s)+MAX(d)<=100%, because exo torque cannot exceed a full gait cycle (100%) 
nvars  =3;

var_containers = [];
for i = 1:nvars
    var_name = ['var', num2str(i)]; % Generate variable name dynamically
    var_i = optimizableVariable(var_name, [lb_list(i), ub_list(i)]); % v1 = optimizableVariable('var1',[lb(1) ub(1)]);
    var_containers = [var_containers, var_i]; % Append the new variable to the array
end

% FUNCTION TO BE EVALUATED, INNER OPTIMIZER IS WITHIN THIS FUNCTION
funMRS = @(x) MuscleRedundancyAndFuncAnalysis(Misc,PGMparam,{x},J_normalized); 

% SETUP AND RUN THE OUTER OPTIMIZER
MaxObjectiveEvaluations=10;
resultsBayesopt = bayesopt(funMRS,var_containers,'IsObjectiveDeterministic',true,"MaxObjectiveEvaluations",MaxObjectiveEvaluations,"UseParallel",false); %,"UseParallel",true

%% Summary of results
clc;
MinObjective=resultsBayesopt.MinObjective;
paramsAtIters=resultsBayesopt.XTrace;
valuesAtIters=resultsBayesopt.ObjectiveTrace;

bestIter  =find(valuesAtIters==MinObjective);
bestParams=paramsAtIters(bestIter,:);
bestParams_val=[bestParams.var1 bestParams.var2 bestParams.var3]; % CHECK YOUR RESULTS (1!)
range_list=ub_list-lb_list;

bestParams_per=(bestParams_val-lb_list)./range_list*100;

disp(['outer objective, metabolic cost change ' num2str(MinObjective,'%1.2f') '%'])
disp('optimal control parameters:');
disp(['stiffness= ' num2str(bestParams_val(1),'%1.0f') 'KN/m, [limits,' num2str(bestParams_per(1),'%1.1f') '%]']);
disp(['start time= ' num2str(bestParams_val(2),'%1.1f')  '%GC, [limits,' num2str(bestParams_per(2),'%1.1f') '%]']);
disp(['duration= ' num2str(bestParams_val(3),'%1.1f')  '%GC, [limits,' num2str(bestParams_per(3),'%1.1f') '%]']);
%% Get results from best iteration

% generate a torque profile based on PGM parametrization
% control profile parameters
PGMparam.torque.stiffness = bestParams_val(1); % [N/m]
PGMparam.torque.startTime = bestParams_val(2); % [gait cycle %]
PGMparam.torque.duration  = bestParams_val(3); % [gait cycle %]

opt_visual.flag       = 0;
[~,PGMmoment]  = PGMactuation_force(PGMparam,trc_time,time,opt_visual);
Misc.torqueProfile(1) = {PGMmoment}; % plot this if you want to visualize moment profile
[Result,~,Misc] = solveMuscleRedundancy_ExoCal(model_path,time,OutPath,Misc);

[J_best,Edot_bilevel_TS] = computeOuterLoopFunction(Misc,Result);
%% Plot best iteration

torque_length=length(Misc.torqueProfile{1});
gait_cycle=linspace(0,100,torque_length);
subplot(2,2,1);
plot(gait_cycle,Misc.torqueProfile{1});
set(gca,'FontSize',15);
ylabel('moment [Nm]')
title('time series')

subplot(2,2,3);
plot(gait_cycle(1:end-1),Edot_normal_TS,'k'); hold on
plot(gait_cycle(1:end-1),Edot_bilevel_TS,'r'); hold on
set(gca,'FontSize',15);
ylabel('metabolic rates (one leg) [W/kg]')
title('time series')

J_real=(mean(Edot_bilevel_TS)-mean(Edot_normal_TS))/mean(Edot_normal_TS)*100; % CHECK YOUR RESULTS (2!)

subplot(1,2,2);
x_label={'Normal' 'Bilevel'};
X = categorical(x_label);
X = reordercats(X,x_label);

data_val=[J_normalized J_best];

color_bar={'k' 'r'};

for j=1:length(x_label)
    data= data_val(j);

    hold on;
    bar(X(j),data,'FaceColor',color_bar{j})

    if j>1 % compute relative change compared to the first one
        E_normal=data_val(1);
        change_per=(data-E_normal)/E_normal*100;
        text(j,data+0.25,[num2str(change_per,'%+1.2f') ' %'],'HorizontalAlignment','center','FontSize',15)
    end
end
set(gca,'FontSize',15);
ylabel('metabolic rates (both legs) [W/kg]')
title('mean values')
%%
function [J] = MuscleRedundancyAndFuncAnalysis(Misc,PGMparam,exoParam,J_normalized)
% UNLOADING
model_path= Misc.model_path;
time      = Misc.time;
OutPath   = Misc.OutPath;
trc_time  = Misc.trc_time;

% generate a torque profile based on PGM parametrization
% control profile parameters
PGMparam.torque.stiffness = exoParam{1}.var1; % [N/m]
PGMparam.torque.startTime = exoParam{1}.var2; % [gait cycle %]
PGMparam.torque.duration  = exoParam{1}.var3; % [gait cycle %]

opt_visual.flag       = 0;
[~,PGMmoment]  = PGMactuation_force(PGMparam,trc_time,time,opt_visual);
Misc.torqueProfile(1) = {PGMmoment}; % plot this if you want to visualize moment profile

% RUN THE INNER LOOP
try
    [Result,~,Misc] = solveMuscleRedundancy_ExoCal(model_path,time,OutPath,Misc);
    [J,~]=computeOuterLoopFunction(Misc,Result);
catch ME
     J=100; % if simulation fails, the optimizer does not benefit from it
end

J=(J-J_normalized)/J_normalized*100; % as a percentage of the unassisted condition
end

function [J,Edot_normal_TS]=computeOuterLoopFunction(Misc,Results)

% window of analysis
frame_extra=5;
fSel=1+frame_extra:size(Results.MActivation.genericMRS,2)-1-frame_extra;

% number of muscles
NMuscle=length(Results.MuscleNames);

% get the mass of the subject using the function GetModelMass
modelmass = getModelMass(Misc.model_path);

% use the post processing function to compute the metabolic energy consumption
Results.E= GetMetabFromMRS(Results,Misc,modelmass);
Results.Edot.genericMRS=Results.E.genericMRS.Bargh2004.Edot;

Edot_normal      =Results.Edot.genericMRS;
Edot_normal_TS   =sum(Edot_normal)/NMuscle;
Edot_normal_mean =mean(Edot_normal_TS)*2; % 2 legs

% objectives and helper
aim_main=Edot_normal_mean;

NDof      =size(Results.RActivation.genericMRS,1);
helper_obs=sum(abs(Results.RActivation.genericMRS(:,fSel)))/NDof; % RA can be possitive (agonist) or negative (antagonist)
aim_helper=sum(helper_obs,'all');

J=aim_main+0.1*aim_helper;
end

% compute: actuator length and moment arm
function [PGMparam] = DEVactuation_geometry(PGMparam,trc_header,trc_time,trc_data_conv,markers_list,opt_visual) 

    data_length=size(trc_data_conv,1);

    % proximal marker
    ind_proximal           = find(strcmp(trc_header,PGMparam.user.proximal_marker_reference)); 
    proximal_marker_user   = trc_data_conv(:,((ind_proximal-1)*3+1)+[0 1]); % x and y
    proximal_marker_device = proximal_marker_user + PGMparam.user.proximal_marker_offset; % TO FIX - RIGID BODY ASSUMPTION NOT HERE
    
    % distal marker
    ind_distal           = find(strcmp(trc_header,PGMparam.user.distal_marker_reference)); 
    distal_marker_user   = trc_data_conv(:,((ind_distal-1)*3+1)+[0 1]); % x and y
    distal_marker_device = distal_marker_user + PGMparam.user.distal_marker_offset; % TO FIX - RIGID BODY ASSUMPTION NOT HERE
    
    ind_HJC              = find(strcmp(trc_header,PGMparam.user.HJC_marker_reference)); 
    HJC_marker_user      = trc_data_conv(:,((ind_HJC-1)*3+1)+[0 1]); % x and y

    % length
    PGM_length=zeros(data_length,1);
    for i=1:data_length
        PGM_length(i)=norm(proximal_marker_device(i,:)-distal_marker_device(i,:));
    end

    % moment arm (respect to the hip joint center)
    PGM_moment_arm = perpendicular_distances(distal_marker_device, proximal_marker_device, HJC_marker_user);

    PGMparam.geometry.length=PGM_length;
    PGMparam.geometry.moment_arm=PGM_moment_arm;

    % visualization
    if opt_visual.flag==1
        figure;

        subplot(3,1,1)
        plot(trc_time,PGM_length)
        xlim([trc_time(1) trc_time(end)])
        ylabel('PGM length [m]'); xlabel('time [s]')

        subplot(3,1,2)
        plot(trc_time,PGM_moment_arm)
        xlim([trc_time(1) trc_time(end)])
        ylabel('PGM moment arm, HJC [m]'); xlabel('time [s]')

        subplot(3,1,3)
        numMakers_list=length(markers_list);
        marker_ind=zeros(1,numMakers_list);
        for i=1:numMakers_list
            marker_ind(i)=find(strcmp(trc_header,markers_list(i)));
        end
        axis_choosen= [1 2]; % these numbers correspond to indexes aligned to X (horizontal) and Y (vertical) axes

        for i=1:10:data_length
            for j=1:numMakers_list
                marker_sel = marker_ind(j);
                ind        = ((marker_sel-1)*3+1) + (axis_choosen-1);

                if strcmp(markers_list{j}(1),'L'); color_marker='r'; else; color_marker='g'; end

                plot(trc_data_conv(i,ind(1)),trc_data_conv(i,ind(2)),'.','MarkerSize',10,'Color',color_marker);
                hold on;
            end
            plot(distal_marker_device(i,1),distal_marker_device(i,2),'.','MarkerSize',20,'Color','b');
            plot(proximal_marker_device(i,1),proximal_marker_device(i,2),'.','MarkerSize',20,'Color','b');

            axis([-2.3 2.7 -0.025 2.00])
            hold off;
            title(['marker trajectories t=' num2str(trc_time(i),'%1.2f')  's'])
            ylabel('vertical [m]'); xlabel('horizontal [m]')
            pause(0.10)
        end
    end
end

function [PGM_force,PGM_moment]        = PGMactuation_force(PGMparam,trc_time,time,opt_visual)
    ind_ti=find(trc_time==time(1));
    ind_tf=find(trc_time==time(end));

    time_series = time(1):0.01:time(end);
    data_length = length(time_series);
    gait_cycle  = linspace(0,100,length(time_series));

    % control signal generation
    control_signal=zeros(data_length,1); % defined
    ind =gait_cycle>PGMparam.torque.startTime & gait_cycle<(PGMparam.torque.startTime+PGMparam.torque.duration); % active time
    control_signal(ind)=1;

    % control signal - smoothed, 1st order adopted to mimic slow increase and
    % decrease of the pressure that fills the gel muscle actuator

    % Parameters
    tau = PGMparam.actuator.time_lag; %(s) It must be changed according to actuator properties (EXP required)
    y0  = 0;    % Initial condition

    % Time span for the solution
    tspan = time_series; % From t = 0 to t = 1
    u     = control_signal;

    u_interp = @(t) interp1(time_series, u, t, 'linear', 'extrap'); % Interpolating function for u(t)
    ode = @(t, y) (u_interp(t) - y) / tau; % first order differential eq.
    [t_int, y_int] = ode45(ode, tspan, y0); % Solve the ODE

    control_signal_ode=y_int;

    % compute forces: spring and damping

    exo_delta            = PGMparam.geometry.length(ind_ti:ind_tf)-PGMparam.actuator.slack_length;
    exo_active_delta_len = exo_delta.*control_signal;            % without ode smoothing (not used)
    exo_active_delta_len_smooth=exo_delta.*control_signal_ode; % with ode smoothing

    % exo_active_delta_vel=gradient(exo_active_delta_len_smooth)./gradient(time'); % velocity
    spring_force   =PGMparam.torque.stiffness*exo_active_delta_len_smooth; % spring force
    % damping_force  =-abs(exo_active_delta_vel*exo_damping); % damping force - compute

    actuator_force =spring_force; %+damping_force; % total force

    actuator_force_only_pos=zeros(data_length,1);
    actuator_force_only_pos(actuator_force>0)=actuator_force(actuator_force>0);

    PGM_force=actuator_force_only_pos;
    
    PGM_moment=PGMparam.geometry.moment_arm(ind_ti:ind_tf).*PGM_force;

    if opt_visual.flag ==1
        figure;

        % length
        PGMparam_length=PGMparam.geometry.length(ind_ti:ind_tf);
        PGMparam_length_min=min(PGMparam_length);
        PGMparam_length_max=max(PGMparam_length);
        PGMparam_length_tot=PGMparam_length_max - PGMparam_length_min;

        subplot(2,2,1)
        ind_actuator=find(ind);
        hold on;
        plot(time_series,PGMparam_length)
        yline(PGMparam.actuator.slack_length,'--r')
        xline(time_series(ind_actuator(1)),':k')
        xline(time_series(ind_actuator(end)),':k')
        % axis([time_series(1) time_series(end) PGMparam_length_min-0.1*PGMparam_length_tot PGMparam_length_max+0.1*PGMparam_length_tot])
        xlim([time_series(1) time_series(end)])
        xlabel('time [s]'); ylabel('length [m]');

        % moment arm
        PGMparam_length=PGMparam.geometry.moment_arm(ind_ti:ind_tf);
        PGMparam_length_min=min(PGMparam_length);
        PGMparam_length_max=max(PGMparam_length);
        PGMparam_length_tot=PGMparam_length_max - PGMparam_length_min;

        subplot(2,2,2)
        plot(time_series,PGMparam_length)
        xline(time_series(ind_actuator(1)),':k')
        xline(time_series(ind_actuator(end)),':k')
        axis([time_series(1) time_series(end) PGMparam_length_min-0.1*PGMparam_length_tot PGMparam_length_max+0.1*PGMparam_length_tot])
        xlabel('time [s]'); ylabel('moment arm [m]');

        subplot(2,2,3)
        plot(time_series,PGM_force)
        axis([time_series(1) time_series(end) 0 1.1*max(PGM_force)])
        xlabel('time [s]'); ylabel('force [N]');

        subplot(2,2,4)
        plot(time_series,PGM_moment)
        axis([time_series(1) time_series(end) 0 100])
        xlabel('time [s]'); ylabel('moment [Nm]');
    end
end

function distances = perpendicular_distances(pos_exo_UA, pos_exo_LA, marker_p)
    % Check input dimensions
    if size(pos_exo_UA, 1) ~= size(pos_exo_LA, 1) || ...
       size(pos_exo_UA, 1) ~= size(marker_p, 1)
        error('Input matrices must have the same number of rows.');
    end

    % Preallocate distance array
    distances = zeros(size(marker_p, 1), 1);

    % Compute the distances
    for i = 1:size(pos_exo_UA, 1)
        % Get points
        A = pos_exo_UA(i, :);
        B = pos_exo_LA(i, :);
        P = marker_p(i, :);
        
        % Direction vector of the line AB
        AB = B - A;
        
        % Vector from A to the marker point P
        AP = P - A;
        
        % Project point P onto line AB
        proj_length = dot(AP, AB) / norm(AB);
        proj_point = A + proj_length * (AB / norm(AB));
        
        % Calculate the distance from marker_p to the projection point
        distances(i) = norm(P - proj_point);
    end
end