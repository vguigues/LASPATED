

function [fL]=oracleObjectiveMissingModel2(nbObservations,nbCalls,nb_missing_calls,neighbors,lambda,T,R,P,G,alpha,durations,S,probs,sample_missing_calls,sample_calls)

fL=0;
for g=1:G
    for t=1:T
        for i=1:R
            for p=1:P 
                currentLambda=lambda(t,g,i,p);
                fL=fL+nbObservations(t,g,i,p)*currentLambda*durations(g,t);
            end
        end
   end
end


for g=1:G
    g
    for t=1:T
        for p=1:P
            for n=1:nbObservations(t,g,1,p)
                Sample=mnrnd(sample_missing_calls(t,g,p,n),probs,S);
                u=0;
                for s=1:S
                    aux=1;
                    for k=1:R
                        iaux=exp(-durations(g,t)*lambda(t,g,k,p));
                        iaux=iaux*((durations(g,t)*lambda(t,g,k,p))^(Sample(s,k)+sample_calls(t,g,k,p,n)));
                        denom=factorial(Sample(s,k)+sample_calls(t,g,k,p,n));
                        iaux=iaux/denom;
                        aux=aux*iaux;
                    end
                    u=u+aux;
                end
                fL=fL-log(u/S);
            end
        end
    end
end



