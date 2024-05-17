function [x,fVal]=projectedGradientArmijoFeasiblev2old(nbObservations,nbArrivals,neighbors,type,distance,T,R,C,sigma,x,iterMax,alpha,epsilon,durations,Groups,whichgroup,weights,delta,uppx,lowx)

k=1;
fVal=[];
bparam=2;
betak=bparam;
upper=10^(10);
lower=-10^(10);

while ((k<=iterMax)&&((upper-lower)>delta))
      [fold]=oracleObjectiveModel1v2(nbObservations,nbArrivals,neighbors,type,distance,x,T,R,C,alpha,durations,Groups,whichgroup,weights);
      %fold
      [gradient]=oracleGradientModel1v2(nbObservations,nbArrivals,neighbors,type,distance,x,T,R,C,alpha,durations,Groups,whichgroup,weights);
      z=max(x-betak*gradient,epsilon);
      bool=1;
      j=0;
      rhs=sum(sum(sum(gradient.*(x-z))));
      while (bool==1)
        %j
        zAux=x+(1/(2^j))*(z-x);
        [f]=oracleObjectiveModel1v2(nbObservations,nbArrivals,neighbors,type,distance,zAux,T,R,C,alpha,durations,Groups,whichgroup,weights); 
        %f
        if (f<=fold-(sigma/2^j)*rhs)
            bool=0;
        else
            j=j+1;
        end
      end
      fVal=[fVal,f];
      x=zAux;
      [gradient]=oracleGradientModel1v2(nbObservations,nbArrivals,neighbors,type,distance,x,T,R,C,alpha,durations,Groups,whichgroup,weights);
      upper=min(upper,f);
      lower=f;
      change=0;
      for t=1:T
        for i=1:R
            for p=1:C
                lower=lower-gradient(p,i,t)*x(p,i,t);
                if (gradient(p,i,t)>0)
                    change=change+gradient(p,i,t)*lowx(p,i,t);
                else
                    change=change+gradient(p,i,t)*uppx(p,i,t);
                end
            end
        end
      end
      lower=lower+change;
      %betak=bparam/(2^j);

      k=k+1;
      %pause
end
