function [myParams] = labelParams(array_kT,array_params,param_label)
    [row,~]=size(array_params);
    for i=1:row
        myParams.([param_label{i}])=array_params(i,:);
    end
    myParams.kT=array_kT;
end