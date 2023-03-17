

function [x,fVal]=projectedGradientArmijoFeasible2(nbObservations,nbArrivals,x,C,D,T,R,weight,regressor,nbLandTypes,iterMax,sigma,epsilon,indexBeta,Groups,durations,sizex,whichgroup)

k=1;
fVal=[];
bparam=2;
betak=bparam;
lbounds=zeros(sizex,1);

while (k<=iterMax)
      [fold]=oracleObjectiveModel2(nbObservations,nbArrivals,x,C,D,T,R,regressor,nbLandTypes,Groups,weight,durations,indexBeta);
      [gradient]=oracleGradientModel2(nbObservations,nbArrivals,x,C,D,T,R,regressor,nbLandTypes,sizex,Groups,whichgroup,weight,durations,indexBeta);
      [z]=projectionRegressor(regressor,indexBeta,nbLandTypes,C,D,T,R,epsilon,sizex,x-betak*gradient,lbounds);
      bool=1;
      j=0;
      rhs=gradient'*(x-z);
      while (bool==1)
        zAux=x+(1/(2^j))*(z-x);
        [f]=oracleObjectiveModel2(nbObservations,nbArrivals,zAux,C,D,T,R,regressor,nbLandTypes,Groups,weight,durations,indexBeta);        
        if (f<=fold-(sigma/2^j)*rhs)
            bool=0;
        else
            j=j+1;
        end
      end
      fVal=[fVal,f];
      x=zAux;
      k=k+1;
      betak=bparam/(2^j);
end
