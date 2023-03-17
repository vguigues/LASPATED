
                                                  
function [x,fVal]=projectedGradientArmijoFeasibleQuick(nbObservations,nbCalls,neighbors,type,distance,R,sigma,x,iterMax,alpha,epsilon,t,g,p)

k=1;
fVal=[];
bparam=2;
betak=bparam;

while (k<=iterMax)
      %k
      [fold]=oracleObjectiveModel1Quick(nbObservations,nbCalls,neighbors,type,distance,x,R,alpha,t,g,p);
      [gradient]=oracleGradientModel1Quick(nbObservations,nbCalls,neighbors,type,distance,x,R,alpha,t,g,p);
      z=max(x-betak*gradient,epsilon);
      bool=1;
      j=0;
      rhs=gradient'*(x-z);
      while (bool==1)
        zAux=x+(1/(2^j))*(z-x); 
        [f]=oracleObjectiveModel1Quick(nbObservations,nbCalls,neighbors,type,distance,zAux,R,alpha,t,g,p);
        if (f<=fold-(sigma/2^j)*rhs)
            bool=0;
        else
            j=j+1;
        end
      end
      fVal=[fVal,f];
      x=x+(1/2^j)*(z-x);
      k=k+1;
      betak=bparam/(2^j);
end
