            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sigma: parameter 0<sigma<1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%betaBar: a positive parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x: starting point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,fVal]=projectedGradientArmijoBoundary2(nbObservations,indexBeta,x,C,D,T,R,regressor,nbLandTypes,sizex,iterMax,sigma,betaBar,epsilon,Groups,durations,sizex,whichgroup)

k=1;
fVal=[];

while (k<=iterMax)
      k
      [fold]=oracleObjectiveModel2(nbObservations,nbCalls,x,C,D,T,R,regressor,nbLandTypes,Groups,weight,durations,indexBeta);
      [gradient]=oracleGradientModel2(nbObservations,nbCalls,x,C,D,T,R,regressor,nbLandTypes,sizex,Groups,whichgroup,weight,durations,indexBeta);                                
      bool=1;
      j=0;
      while (bool==1)
       [z]=projectionRegressor(regressor,indexBeta,nbLandTypes,C,D,T,R,epsilon,sizex,x-(betaBar/(2^j))*gradient);  
       [f]=oracleObjectiveModel2(nbObservations,nbCalls,z,C,D,T,R,regressor,nbLandTypes,Groups,weight,durations,indexBeta);
       rhs=gradient'*(x-z);
       if (f<=fold-sigma*rhs)
            bool=0;
       else
            j=j+1;
       end
      end
      k=k+1;
      fVal=[fVal,f];
      x=z;
end
