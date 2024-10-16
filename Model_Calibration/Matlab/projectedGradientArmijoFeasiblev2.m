
       
function [x,fVal]=projectedGradientArmijoFeasiblev2(nbObservations,nbArrivals,neighbors,type,distance,T,R,C,sigma,x,iterMax,alpha,epsilon,durations,Groups,whichgroup,weights)

k=1;
fVal=[];
bparam=2;
betak=bparam;
upper=10^(10);
lower=-10^(10);
[fold]=oracleObjectiveModel1v2(nbObservations,nbArrivals,neighbors,type,distance,x,T,R,C,alpha,durations,Groups,whichgroup,weights);
[gradient]=oracleGradientModel1v2(nbObservations,nbArrivals,neighbors,type,distance,x,T,R,C,alpha,durations,Groups,whichgroup,weights);

%while ((k<=iterMax)&&((upper-lower)>delta))
while (k<=iterMax)
%fold
    z=max(x-betak*gradient,epsilon);
    rhs=sum(sum(sum(gradient.*(x-z))));
    [f]=oracleObjectiveModel1v2(nbObservations,nbArrivals,neighbors,type,distance,z,T,R,C,alpha,durations,Groups,whichgroup,weights);
    if (f>fold-sigma*rhs)
        bool=1;
        j=1;
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
        betak=betak/2^(j-1);
    else
        x=z;
        betak=2*betak;
    end
    fold=f;
    [gradient]=oracleGradientModel1v2(nbObservations,nbArrivals,neighbors,type,distance,x,T,R,C,alpha,durations,Groups,whichgroup,weights);
    upper=min(upper,f);
    lower=f;
%     change=0;
%     for t=1:T
%         for i=1:R
%             for p=1:C
%                 lower=lower-gradient(p,i,t)*x(p,i,t);
%                 if (gradient(p,i,t)>0)
%                     change=change+gradient(p,i,t)*lowx(p,i,t);
%                 else
%                     change=change+gradient(p,i,t)*uppx(p,i,t);
%                 end
%             end
%         end
%     end
%     lower=lower+change;
    k=k+1;
end
