function [US_data,nUSdata,ind_US,ind_US_AT,ind_USnone,USDigitalizedInterp]=fiberCalibrationSetup(DatStore,Misc,time_opt)
% obtain the indexes and interpolate the fibers

% get ultrasound digitalized data and indexes
All_MuscleNames= DatStore(1).MuscleNames;
getSide        = Misc.side_sel;     

US_file  = ReadMotFile(Misc.USfile{1});
US_data  = US_file.data;
USNames  = US_file.names(2:end);          % lMo and lTs
USNames  = strcat(USNames,['_' getSide]); % update with leg side
nUSdata  = length(USNames);

ind_US=zeros(nUSdata,1);
for i=1:nUSdata
    ind_US(i) = find(strcmp(All_MuscleNames,USNames(i)));
end

USNames_AT={'soleus' 'gasmed' 'gaslat'};
USNames_AT=strcat(USNames_AT,['_' getSide]);
nUSATdata =length(USNames_AT);

ind_US_AT    =zeros(nUSATdata,1);
ind_US_AT_var=zeros(nUSATdata,1);
for i=1:nUSATdata
    ind_US_AT(i)    = find(strcmp(All_MuscleNames,USNames_AT(i)));
    ind_US_AT_var(i)= find(strcmp(USNames        ,USNames_AT(i)));
end

USnone     =setdiff(All_MuscleNames,USNames);
nUSdata_none=length(USnone);
ind_USnone=zeros(nUSdata_none,1);
for i=1:nUSdata_none
    ind_USnone(i) = find(strcmp(All_MuscleNames,USnone(i)));
end

% interpolate US digitalized
US_upd_time=linspace(time_opt(1),time_opt(end),length(US_data(:,1)));
US_data(:,1)=US_upd_time;
USDigitalizedInterp = interp1(US_data(:,1),US_data,time_opt')';
end