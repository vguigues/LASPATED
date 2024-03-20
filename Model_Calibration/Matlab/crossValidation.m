                                            
function [cputime,weight,lambda]=crossValidation(model,sample,T,R,P,sigma,x,iterMax,proportion,epsilon,durations,Groups,whichgroup,weights,params)

if (strcmp(model,'noreg'))
    disp('Model without regressors chosen'); 
    [cputime,weight,lambda]=crossValidationNoReg(sample,params.neighbors,params.type,params.distance,T,R,P,sigma,x,iterMax,proportion,epsilon,durations,Groups,whichgroup,weights, params.delta, params.uppx, params.lowx);    
else
    disp('Invalid model name');
end


