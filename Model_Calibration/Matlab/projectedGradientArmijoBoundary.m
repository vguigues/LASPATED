
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sigma: parameter 0<sigma<1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%betaBar: a positive parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x: starting point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,fVal]=projectedGradientArmijoBoundary(nbObservations,nbCalls,neighbors,type,distance,T,R,P,sigma,betaBar,x,iterMax,alpha,epsilon,durations,Groups,whichgroup,weight)

k=1;
fVal=[];

while (k<=iterMax)
      k
      [fold]=oracleObjectiveModel1(nbObservations,nbCalls,neighbors,type,distance,x,T,R,P,alpha,durations,Groups,whichgroup,weight); 
      [gradient]=oracleGradientModel1(nbObservations,nbCalls,neighbors,type,distance,x,T,R,P,alpha,durations,Groups,whichgroup,weight);
      bool=1;
      j=0;
      while (bool==1)
        z=max(x-(betaBar/(2^j))*gradient,epsilon); 
        [f]=oracleObjectiveModel1(nbObservations,nbCalls,neighbors,type,distance,z,T,R,P,alpha,durations,Groups,whichgroup,weight); 
        rhs=sum(sum(sum(gradient.*(x-z))));
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
