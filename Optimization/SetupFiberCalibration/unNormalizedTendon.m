function [kT_abs] = unNormalizedTendon(kT_sel,FMo_sel,lTs_sel)
import casadi.*
elo=[2 10]/1000;       % 2 and 10 mm
nMuscle=length(kT_sel);
FT_elo = MX(nMuscle,2); %=zeros(nMuscle,2);
for i=1:nMuscle
lTtilde_elo=(elo./lTs_sel(i))+1;                                   % normalized length
shift      =getShift(kT_sel(i));
fse        =(exp(kT_sel(i).*(lTtilde_elo - 0.995)))/5-0.25+shift;  % normalized force
FT_elo(i,:)=fse.*FMo_sel(i);                                       % un-normalized force
end
elo_mili=elo*1000;                                              % convert to mm
kT_abs   =(FT_elo(:,end)-FT_elo(:,1))./(elo_mili(:,end)-elo_mili(:,1)); % un-normalized stiffness [N/mm]
end