
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sigma: parameter 0<sigma<1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%betaBar: a positive parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x: starting point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,fVal]=projectedGradientArmijoBoundaryQuick(nbObservations,nbCalls,neighbors,type,distance,R,sigma,betaBar,x,iterMax,alpha,epsilon,t,g,p)

k=1;
fVal=[];

while (k<=iterMax)
      [fold]=oracleObjectiveModel1Quick(nbObservations,nbCalls,neighbors,type,distance,x,R,alpha,t,g,p);
      [gradient]=oracleGradientModel1Quick(nbObservations,nbCalls,neighbors,type,distance,x,R,alpha,t,g,p);
      bool=1;
      j=0;
      while (bool==1)
        z=max(x-(betaBar/(2^j))*gradient,epsilon); 
        [f]=oracleObjectiveModel1Quick(nbObservations,nbCalls,neighbors,type,distance,z,R,alpha,t,g,p);
        if (f<=fold-sigma*gradient'*(x-z))
            bool=0;
        else
            j=j+1;
        end
      end
      k=k+1;
      fVal=[fVal,f];
      x=z;
end
