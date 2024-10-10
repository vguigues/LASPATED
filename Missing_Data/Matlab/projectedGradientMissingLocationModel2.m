

function [x,fVal]=projectedGradientMissingLocationModel2(nbObservations,nbCalls,nb_missing_calls,neighbors,G,T,R,C,sigma,x,probs,S,iterMax,alpha,epsilon,durations,sample_missing_calls,sample_calls)

k=1;
fVal=[];
bparam=2;
betak=bparam;

[fold]=oracleObjectiveMissingModel2(nbObservations,nbCalls,nb_missing_calls,neighbors,x,T,R,C,G,alpha,durations,S,probs,sample_missing_calls,sample_calls);
[gradientl]=oracleGradientMissingModel2(nbObservations,nbCalls,nb_missing_calls,neighbors,x,G,T,R,C,alpha,durations,S,probs,sample_missing_calls,sample_calls);

while (k<=iterMax)
    k
    z=max(x-betak*gradientl,epsilon);
    rhs=sum(sum(sum(sum(gradientl.*(x-z)))));
    [fL]=oracleObjectiveMissingModel2(nbObservations,nbCalls,nb_missing_calls,neighbors,z,T,R,C,G,alpha,durations,S,probs,sample_missing_calls,sample_calls);
    if (fL>fold-sigma*rhs)
        bool=1;
        j=1;
        while (bool==1)
            zAux=x+(1/(2^j))*(z-x);
            [fL]=oracleObjectiveMissingModel2(nbObservations,nbCalls,nb_missing_calls,neighbors,zAux,T,R,C,G,alpha,durations,S,probs,sample_missing_calls,sample_calls);            
            if (fL<=fold-(sigma/2^j)*rhs)
                bool=0;
            else
                j=j+1;
            end
        end
        x=zAux;
        betak=bparam/2^(j-1);
    else
        x=z;
        j=1;
        betak=2*betak;
    end
    fVal=[fVal,fL];
    fold=fL;
    [gradientl]=oracleGradientMissingModel2(nbObservations,nbCalls,nb_missing_calls,neighbors,x,G,T,R,C,alpha,durations,S,probs,sample_missing_calls,sample_calls);    
    k=k+1;
end
