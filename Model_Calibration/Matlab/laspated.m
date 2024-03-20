        
function [lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,sigma,iterMax,epsilon,lambda0,params)
        
function [lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,sigma,iterMax,epsilon,lambda0,params)

if (strcmp(model,'noreg'))
    disp('Model without regressors chosen');
    % projectedGradientArmijoFeasiblev2(nbObservations,nbArrivals,neighbors,type,distance,T,R,C,sigma,x,iterMax,alpha,epsilon,durations,Groups,whichgroup,weights,delta,uppx,lowx);
    [lambda,obj]=projectedGradientArmijoFeasiblev2(nbObservations,nbArrivals,params.neighbors,params.type,params.distance,T,R,C,sigma,lambda0,iterMax,params.alpha,epsilon,durations,Groups,whichgroup,params.weight,params.delta,params.uppx,params.lowx);
elseif (strcmp(model,'reg'))
    disp('Model with regressors chosen');
    [lambda,obj]=projectedGradientArmijoFeasible2(nbObservations,nbArrivals,lambda0,C,params.D,T,R,params.regressor,params.nbRegressors,iterMax,sigma,epsilon,params.indexBeta,durations,params.sizex);
else
    disp('Invalid model name');
end



