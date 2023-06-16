
                                                  
function [x,fVal]=projectedGradientArmijoFeasible(nbObservations,nbArrivals,neighbors,type,distance,T,R,C,sigma,x,iterMax,alpha,epsilon,durations,Groups,whichgroup,weights)

k=1;
fVal=[];
bparam=2;
betak=bparam;

while (k<=iterMax)
      [fold]=oracleObjectiveModel1(nbObservations,nbArrivals,neighbors,type,distance,x,T,R,C,alpha,durations,Groups,whichgroup,weights);
      %fold
      [gradient]=oracleGradientModel1(nbObservations,nbArrivals,neighbors,type,distance,x,T,R,C,alpha,durations,Groups,whichgroup,weights);
      z=max(x-betak*gradient,epsilon);
      bool=1;
      j=0;
      rhs=sum(sum(sum(gradient.*(x-z))));
      while (bool==1)
        %j
        zAux=x+(1/(2^j))*(z-x);
        [f]=oracleObjectiveModel1(nbObservations,nbArrivals,neighbors,type,distance,zAux,T,R,C,alpha,durations,Groups,whichgroup,weights); 
        %f
        if (f<=fold-(sigma/2^j)*rhs)
            bool=0;
        else
            j=j+1;
        end
      end
      fVal=[fVal,f];
      x=zAux;
      betak=bparam/(2^j);
      k=k+1;
      %pause
end
