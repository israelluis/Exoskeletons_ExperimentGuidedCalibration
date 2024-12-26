function [params,lMo,lTs,FMo,alphao,kT,shift]=ReadMuscleParameters_upd(ModelPath,names,Misc)
% input= path to model and cell array with muscle names
% output= params (5xNMuscles) with  row: (1)  IsomForce (2)OptFiberLength
% 			(3) TendonSlackLength (4) alphao (5) MaxFiberVelocity
% read the model
import org.opensim.modeling.*;
nNames = length(names);
params = zeros(8, nNames);
kT     = zeros(1, nNames);
shift  = zeros(1, nNames);

% read the muscle properties
if strcmp(Misc.Param_read,'model')
model = Model(ModelPath);

muscles = model.getMuscles();

for i = 1:nNames
    kT   =Misc.kT;
    shift=Misc.shift;

    muscle = muscles.get(names{i});
    params(1,i) = muscle.getMaxIsometricForce();  		
    params(2,i) = muscle.getOptimalFiberLength(); 	
    params(3,i) = muscle.getTendonSlackLength();	
    params(4,i) = muscle.getPennationAngleAtOptimalFiberLength(); 
    params(5,i) = muscle.getMaxContractionVelocity();

    params(6,i) = 4.0; % kpe 4.0 % Misc.myParams.kpe(i)
    params(7,i) = 1.0; % so  1.0
    params(8,i) = 0.6; % sM  0.6
end

% get the muscle properties from a labeled structure. User input
elseif strcmp(Misc.Param_read,'selected')
for i = 1:length(names)
    kT(i)       =Misc.myParams.kT(i);
    shift(i)    =getShift(kT(i));
    
    params(1,i) =Misc.myParams.FMo(i);
    params(2,i) =Misc.myParams.lMo(i);
    params(3,i) =Misc.myParams.lTs(i);   
    params(4,i) =Misc.myParams.alphao(i);
    params(5,i) =Misc.myParams.vMmax(i);
    params(6,i) =Misc.myParams.kpe(i);
    params(7,i) =Misc.myParams.so(i);
    params(8,i) =Misc.myParams.sM(i);
end
end

% update params selecting param name and muscle name
% update is performed whereas parameters were read or user input
if  isfield(Misc,'updParams') %|| ~isempty(Misc.updParams)
    updParams_fields=fieldnames(Misc.updParams);
    nFields         =length(updParams_fields);
    for i=1:nFields
        current_field=(updParams_fields(i)); % get current field
        param_ind    =find(strcmp(Misc.param_label,current_field)); %get ind param
        
        names_upd =Misc.updParams.(current_field{:}).names;
        % update if all muscle include
        if strcmp(names_upd,'all');names_upd= cellfun(@(x) x(1:end-2), names, 'un', 0); end
        nMuscles  =length(names_upd);
        
        ind=zeros(nMuscles,1);
        for j=1:nMuscles
            ind(j)=find(strcmp(names,[names_upd{j} '_' Misc.side_sel])); %
        end
        
        % apply change
        updParams_change      =fieldnames(Misc.updParams.(current_field{:}));
        updParams_change_param=updParams_change;
        updParams_change_param(strcmp(updParams_change, 'names')) = [];
        if length(updParams_change_param)>1; disp('you cannot modify a parameter based on value and percentage simultaneously'); end
        if strcmp(updParams_change_param,'values') % replace by a given value
            values=Misc.updParams.(current_field{:}).values;
            params(param_ind,ind)=values;
            
        elseif strcmp(updParams_change_param,'percentages') % modify it based on percentage over the original value
            percentages=Misc.updParams.(current_field{:}).percentages;
            params(param_ind,ind)=params(param_ind,ind).*percentages;
        end
    end
end

% create additional variables with the same information
FMo=params(1,:);
lMo=params(2,:);
lTs=params(3,:);
alphao=params(4,:);

end