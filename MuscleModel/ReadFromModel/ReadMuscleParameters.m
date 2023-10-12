function [params,lMo,lTs,FMo,alphao]=ReadMuscleParameters(ModelPath,names,Misc)
% input= path to model and cell array with muscle names
% output= params (5xNMuscles) with  row: (1)  IsomForce (2)OptFiberLength
% 			(3) TendonSlackLength (4) alphao (5) MaxFiberVelocity

% read the model
import org.opensim.modeling.*;
model = Model(ModelPath);

% read the muscle properties
nNames = length(names);
params = zeros(8, nNames);
muscles = model.getMuscles();


for i = 1:nNames
    muscle = muscles.get(names{i});
    params(1,i) = muscle.getMaxIsometricForce();  		
    params(2,i) = muscle.getOptimalFiberLength(); 	
    params(3,i) = muscle.getTendonSlackLength();	
    params(4,i) = muscle.getPennationAngleAtOptimalFiberLength(); 
    params(5,i) = muscle.getMaxContractionVelocity();

    params(6,i) = 4.0; % kpe 4.0 % Misc.myParams.kpe(i)
    params(7,i) = 1.0; % s0  1.0
    params(8,i) = 0.6; % sM  0.6
end


% create additional variables with the same information
FMo=params(1,:);
lMo=params(2,:);
lTs=params(3,:);
alphao=params(4,:);
end