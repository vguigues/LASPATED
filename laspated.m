

function [lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,weight,sigma,iterMax,epsilon,lambda0,params)

if (strcmp(model,'noreg'))
    disp('Model without regressors chosen');
    [lambda,obj]=projectedGradientArmijoFeasible(nbObservations,nbArrivals,params.neighbors,params.type,params.distance,T,R,C,sigma,lambda0,iterMax,params.alpha,epsilon,durations,Groups,whichgroup,weight);
elseif (strcmp(model,'reg'))
    disp('Model with regressors chosen');    
    [lambda,obj]=projectedGradientArmijoFeasible2(nbObservations,nbArrivals,lambda0,C,params.D,T,R,weight,params.regressor,params.nbLandTypes,iterMax,sigma,epsilon,params.indexBeta,Groups,durations,params.sizex,whichgroup);
else
    disp('Invalid model name');
end



