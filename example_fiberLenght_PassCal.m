addpath(genpath('C:\Users\israe\Dropbox\PhDWork\Simulation\code\MuscleRedundancySolver-master\Casadi_installed'));
addpath(genpath('C:\Users\israe\Dropbox\PhDWork\Simulation\code\Exo&Calibration'));
clc; close all; clear Misc
%% Input information
currentFolder = pwd;

% paths and folders
ExampleFolder= 'DataExample';
DataFolder   = 'DataDigitalized';
OutPath      = fullfile(currentFolder,'Results_AllCalibration');  

% experimental data
model_path  =  fullfile(currentFolder,ExampleFolder,'Model.osim');
           
Misc.EMGfile= {fullfile(DataFolder,'emgFiles.mot')};
Misc.USfile = {fullfile(DataFolder,'ultrasoundFiles.mot')};

% select the leg's side
Misc.side_sel      = 'r';

% label of muscle-tendon parameters
param_label={'FMo' 'lMo' 'lTs' 'alphao' 'vMmax' 'kpe' 'so' 'sM'};
Misc.param_label   =param_label;

% run muscle analysis
Misc.GetAnalysis = 1;

% WORKFLOW SETUP
% Here we setup 3 workflows.
% Generic:     Generic parameters, No use muscle fiber lengths
% Calibration: Parameter calibration using muscle fiber lengths
% Validation:  Use of calibrated parameters, No use muscle fiber lengths
% Exoskeleton feature is disabled
Misc.Exo_Enable=0;

workFlow_generic     =1;
workFlow_calFiber    =1;
workFlow_optFiber    =1;
workFlow_calPass     =1;
workFlow_optPass_Fiber=1;

%% Muscle fiber length calibration
Misc.IKfile = {fullfile(currentFolder,ExampleFolder,'IK.mot')};
Misc.IDfile = {fullfile(currentFolder,ExampleFolder,'ID.sto')};
Misc.DofNames_Input={['hip_flexion_' Misc.side_sel] ['hip_adduction_' Misc.side_sel] ['hip_rotation_' Misc.side_sel] ['knee_angle_' Misc.side_sel]  ['ankle_angle_' Misc.side_sel]}; 
time = [0 10];

% simulation with generic params & kT =15 (equivalent to 150 N/mm with
% Rajagoal et al. original values)
if workFlow_generic==1
Misc.muscleFiberCal= 0; 
Misc.Param_read    = {'model'};
Misc.lMo_use_ward  = 0;
Misc.OutName       = 'Generic';

PF_kT=15; % plantarflexors
VA_kT=15; % vastii muscle
Misc.Set_kT_ByName = {['soleus_' Misc.side_sel],PF_kT; ['gasmed_' Misc.side_sel], PF_kT; ['gaslat_' Misc.side_sel],PF_kT;...
                      ['vasint_' Misc.side_sel],VA_kT; ['vaslat_' Misc.side_sel], VA_kT; ['vasmed_' Misc.side_sel],VA_kT};

[Results_gen,DatStore,Misc] = solveMuscleRedundancy_ExoCal(model_path,time,OutPath,Misc);
end

% simulation using digitalized muscle fiber lengths. Parameteres: lMo,
% lTs, and kT are design variables
if workFlow_calFiber==1
Misc.muscleFiberCal= 1; 
Misc.Param_read    = {'model'}; 
Misc.lMo_use_ward  = 1;
Misc.OutName       = 'ParmEst';

Misc.Set_kT_ByName ={}; % assign kT=35 for all muscles

[Results_cal,DatStore,Misc] = solveMuscleRedundancy_ExoCal(model_path,time,OutPath,Misc);
end

% simulation with previous calibrated parameters: lMo, lTs, and kT. This
% step does not use digitalized muscle fiber lengths
if workFlow_optFiber==1
Misc.muscleFiberCal= 0;    
Misc.Param_read    = {'selected'};
Misc.lMo_use_ward  = 0;
Misc.OutName       = 'Validation';

Misc.myParams=labelParams(Results_cal.kT.calibratedMRS,Results_cal.params.calibratedMRS,param_label);

[Results_val,DatStore,Misc] = solveMuscleRedundancy_ExoCal(model_path,time,OutPath,Misc); %solveMuscleRedundancy_USDigitalized_upd1
end

%% Passive joint-momement relationship calibration
ank_file =fullfile(DataPath,['TrialAnk_' Misc.side_sel '_kne0hip0_kneBen15hip0_kneBen60hip0_' model_conf{:} 'Model']);
hip_file =fullfile(DataPath,['TrialHip_' Misc.side_sel '_kneBen15ank0_kneBen60ank0_' model_conf{:} 'Model']);
kne_file =fullfile(DataPath,['TrialKne_' Misc.side_sel '_ankDor20hip0_ankPla15hip0_ankDor20hipExt15_' model_conf{:} 'Model']);
neu_file =fullfile(DataPath,['TrialNeu_' Misc.side_sel '_hipExt7hipAbd2kne0ankPla40_allModels']);

Misc.IKfile = {[ank_file '.mot'] [hip_file '.mot'] [kne_file '.mot'] [neu_file '.mot']};
Misc.IDfile = {[ank_file '.sto'] [hip_file '.sto'] [kne_file '.sto'] [neu_file '.sto']};

Misc.DofNames_Input={{['ankle_angle_' Misc.side_sel]};{['hip_flexion_' Misc.side_sel]};{['knee_angle_' Misc.side_sel]};...        % for passive forces
{['hip_flexion_' Misc.side_sel] ['hip_adduction_' Misc.side_sel] ['hip_rotation_' Misc.side_sel] ['knee_angle_' Misc.side_sel] ['ankle_angle_' Misc.side_sel]}};   % for neutral

time = [0 10; 0 10; 0 10; 0 10];

Misc.model_path = model_path;
Misc.OutPath    = OutPath;
model_conf ={'rajagopal'};               Misc.model_conf=model_conf; 

Misc.semimem_adjustment=1; % change to "2" to have the lTs recomputed
if Misc.semimem_adjustment==1 
    Misc.updParams.lTs.names      ={'semimem'};
    Misc.updParams.lTs.percentages=[0.98]; % i.e., 0.98 its original value
end

if workFlow_calPass==1
    Misc.Param_read={'selected'};
    Misc.Mode_opti={'passive'}; % passive force-length curve optimization
    Misc.OutName = ['GenericPass_Results'];
    
    Misc.myParams=labelParams(Results_cal.kT.calibratedMRS,Results_cal.params.calibratedMRS,param_label);
    
    Misc.add_Passconstraints=1;
    [Results_PassCal,DatStore_PassCal,Misc_PassCal] = solveMuscleRedundancy_calibrationPassive([],time,[],Misc);
end

%% Simulate walking with prior optimized lMo, lTs, kT, and passive force-length
Misc.IKfile = {fullfile(currentFolder,ExampleFolder,'IK.mot')};
Misc.IDfile = {fullfile(currentFolder,ExampleFolder,'ID.sto')};
Misc.DofNames_Input={['hip_flexion_' Misc.side_sel] ['hip_adduction_' Misc.side_sel] ['hip_rotation_' Misc.side_sel] ['knee_angle_' Misc.side_sel]  ['ankle_angle_' Misc.side_sel]}; 
time = [0 10];

if workFlow_optPass_Fiber==1
Misc.muscleFiberCal= 0;    
Misc.Param_read    = {'selected'};
Misc.lMo_use_ward  = 0;
Misc.OutName       = 'all';

Misc.myParams=labelParams(Results_PassCal.kT.calibratedMRS,Results_PassCal.params.calibratedMRS,param_label);
[Results_all,DatStore,Misc] = solveMuscleRedundancy_ExoCal(model_path,time,OutPath,Misc);
end
%% Plot all results for walking
state_list={'MActivation' 'TForce'};
axis_list ={[0 1] [0 1500]};
title_list={'muscle activations with various workflows' 'muscle-tendon forces with various workflows'};
yaxis_list={'muscle activations ' 'muscle-tendon forces [N]'};

MuscleNames= DatStore.MuscleNames;
NMuscle    = length(MuscleNames);
getSide    = MuscleNames{1}(end-1:end);
selectedMus_noSide=strrep(MuscleNames,getSide,'');

if ~exist('EMG_data','var') || ~exist('US_data','var')  
    EMG_data = ReadMotFile(Misc.EMGfile{1});
    US_data  = ReadMotFile(Misc.USfile{1});
end

extra_frames=5;

for state_opt=1:length(state_list)
fig=figure(state_opt); set(gcf,'color','w','Visible','on'); clf;   

for mus_sel=1:NMuscle
    muscle_name=MuscleNames{mus_sel};
    state_Normal=Results_gen.(state_list{state_opt}).genericMRS(mus_sel,1+extra_frames:end-extra_frames);
    state_Val   =Results_val.(state_list{state_opt}).genericMRS(mus_sel,1+extra_frames:end-extra_frames);
    state_All   =Results_all.(state_list{state_opt}).genericMRS(mus_sel,1+extra_frames:end-extra_frames);
    
    x_axis      =linspace(0,100,length(state_Normal));
%     x_axis = Results_gen.Time.genericMRS(1+extra_frames:end-extra_frames);
    
    subplot(5,8,mus_sel)
    hold on;
    plot(x_axis,state_Normal,'color','g','lineWidth',3,'lineStyle','-');
    plot(x_axis,state_Val,'color','b','lineWidth',3,'lineStyle','--');
    plot(x_axis,state_All,'color','r','lineWidth',3,'lineStyle','-.');
    ylim(axis_list{state_opt});
    set(gca,'FontSize',10);
    title(muscle_name)
end
sgtitle(title_list(state_opt),'fontSize',20)

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
label_y = ylabel(han,yaxis_list{state_opt},'FontSize',17); label_y.Position(1) = -0.05; label_y.Position(2) = 0.5;
label_x = xlabel(han,'gait cycle [%]'       ,'FontSize',17);   label_x.Position(1) = 0.5; label_x.Position(2) = -0.05;
end