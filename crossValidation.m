
function [cputime,weight,lambda]=crossValidation(model,sample,T,R,P,sigma,x,iterMax,proportion,epsilon,durations,Groups,whichgroup,weights,params)

if (strcmp(model,'noreg'))
    disp('Model without regressors chosen'); 
    [cputime,weight,lambda]=crossValidationNoReg(sample,params.neighbors,params.type,params.distance,T,R,P,sigma,x,iterMax,proportion,epsilon,durations,Groups,whichgroup,weights);    
elseif (strcmp(model,'reg'))
    disp('Model with regressors chosen');  
    [cputime,weight,lambda]=crossValidationReg(sample,params.indexBeta,x,P,params.D,T,R,params.regressor,params.nbLandTypes,proportion,params.sizex,iterMax,sigma,epsilon,params.nbObs,params.isHolidays,params.obsbefore,Groups,durations,whichgroup,weights)
else
    disp('Invalid model name');
end


