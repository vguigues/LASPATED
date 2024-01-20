

function [x,fVal]=projectedGradientArmijoFeasible2(nbObservations,nbArrivals,x,C,D,T,R,regressor,nbRegressors,iterMax,sigma,epsilon,indexBeta,durations,sizex,delta)

k=1;
fVal=[];
bparam=2;
betak=bparam;
lbounds=zeros(sizex,1);
upper=10^(50);
lower=-10^(50)

while ((k<=iterMax)&&((upper-lower)>delta))
      [fold]=oracleObjectiveModel2(nbObservations,nbArrivals,x,C,D,T,R,regressor,nbRegressors,durations,indexBeta);
      [gradient]=oracleGradientModel2(nbObservations,nbArrivals,x,C,D,T,R,regressor,nbRegressors,sizex,durations,indexBeta);
      [z]=projectionRegressor(regressor,indexBeta,nbRegressors,C,D,T,R,epsilon,sizex,x-betak*gradient,lbounds);
      bool=1;
      j=0;
      rhs=gradient'*(x-z);
      while (bool==1)
        zAux=x+(1/(2^j))*(z-x);
        [f]=oracleObjectiveModel2(nbObservations,nbArrivals,zAux,C,D,T,R,regressor,nbRegressors,durations,indexBeta);        
        if (f<=fold-(sigma/2^j)*rhs)
            bool=0;
        else
            j=j+1;
        end
      end
      fVal=[fVal,f];
      x=zAux;
      k=k+1;
      upper=min(upper,f);
      
      betak=bparam/(2^j);
end
