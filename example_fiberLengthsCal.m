addpath(genpath('C:\Users\israe\Dropbox\PhDWork\Simulation\code\MuscleRedundancySolver-master\Casadi_installed'));
addpath(genpath('C:\Users\israe\Dropbox\PhDWork\Simulation\code\Exo&Calibration'));
clc; close all; clear Misc
%% Input information
currentFolder = pwd;

% paths and folders
ExampleFolder= 'DataExample';
DataFolder   = 'DataDigitalized';
OutPath      = fullfile(currentFolder,'Results_FiberLengths');  

% experimental data
Misc.IKfile = {fullfile(currentFolder,ExampleFolder,'IK.mot')};
Misc.IDfile = {fullfile(currentFolder,ExampleFolder,'ID.sto')};
model_path  =  fullfile(currentFolder,ExampleFolder,'Model.osim');
           
Misc.EMGfile= {fullfile(DataFolder,'emgFiles.mot')};
Misc.USfile = {fullfile(DataFolder,'ultrasoundFiles.mot')};

% select the leg's side
Misc.side_sel      = 'r';

% input dofs - select the DOFs you want to include in the optimization
Misc.DofNames_Input={['hip_flexion_' Misc.side_sel] ['hip_adduction_' Misc.side_sel] ['hip_rotation_' Misc.side_sel] ['knee_angle_' Misc.side_sel]  ['ankle_angle_' Misc.side_sel]}; 

% ranges
time = [0 10];

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
workFlow_calibration =1;
workFlow_validation  =1;

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
if workFlow_calibration==1
Misc.muscleFiberCal= 1; 
Misc.Param_read    = {'model'}; 
Misc.lMo_use_ward  = 1;
Misc.OutName       = 'ParmEst';

Misc.Set_kT_ByName ={}; % assign kT=35 for all muscles

[Results_cal,DatStore,Misc] = solveMuscleRedundancy_ExoCal(model_path,time,OutPath,Misc);
end

% simulation with previous calibrated parameters: lMo, lTs, and kT. This
% step does not use digitalized muscle fiber lengths
if workFlow_validation==1
Misc.muscleFiberCal= 0;    
Misc.Param_read    = {'selected'};
Misc.lMo_use_ward  = 0;
Misc.OutName       = 'Validation';

Misc.myParams=labelParams(Results_cal.kT.calibratedMRS,Results_cal.params.calibratedMRS,param_label);

[Results_val,DatStore,Misc] = solveMuscleRedundancy_ExoCal(model_path,time,OutPath,Misc); %solveMuscleRedundancy_USDigitalized_upd1
end
%% Quick plotting

MuscleNames= DatStore.MuscleNames;
NMuscle    = length(MuscleNames);
getSide    = MuscleNames{1}(end-1:end);
selectedMus_noSide=strrep(MuscleNames,getSide,'');

AT_value_gen=sum(Results_cal.kT_ATunNormalized.genericMRS);
AT_value_cal=sum(Results_cal.kT_ATunNormalized.calibratedMRS);
disp(['Achilles tendon=' ' generic:' num2str(AT_value_gen,'%4.2f') 'N/mm' ...
                         ' calibrated:' num2str(AT_value_cal,'%4.2f') 'N/mm'])

if ~exist('EMG_data','var') || ~exist('US_data','var')  
    EMG_data = ReadMotFile(Misc.EMGfile{1});
    US_data  = ReadMotFile(Misc.USfile{1});
end

fig=figure(1); set(gcf,'color','w','Visible','on'); clf;    
for mus_sel=1:NMuscle
    muscle_name=MuscleNames{mus_sel};
    state_Normal=Results_gen.lMtildeopt.genericMRS(mus_sel,:);
    state_Cal   =Results_cal.lMtildeopt.genericMRS(mus_sel,:);
    state_Val   =Results_val.lMtildeopt.genericMRS(mus_sel,:);
    x_axis      =linspace(0,100,length(state_Normal));
    
    subplot(5,8,mus_sel)
    hold on;
    plot(x_axis,state_Normal,'color',"#77AC30",'lineWidth',2);
    plot(x_axis,state_Cal,'color',"#0072BD",'lineWidth',2)
    plot(x_axis,state_Val,'color',"#A2142F",'lineWidth',3,'lineStyle',':')
    ylim([0.5 1.5]);
    set(gca,'FontSize',10);
    title(muscle_name)
end

for i=1:length(US_data.names(2:end))
    US_muscle_name=US_data.names(1+i,:);
    ind=find(strcmp(selectedMus_noSide,US_muscle_name));
    x_axis=linspace(0,100,length(US_data.data(:,i+1)));

    subplot(5,8,ind)
    plot(x_axis,US_data.data(:,i+1),'color',"k",'lineWidth',2)
end

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
label_y = ylabel(han,'normalized fiber lengths [lM/lMo]','FontSize',17); label_y.Position(1) = -0.05; label_y.Position(2) = 0.5;
label_x = xlabel(han,'gait cycle [%]'       ,'FontSize',17);   label_x.Position(1) = 0.5; label_x.Position(2) = -0.05;