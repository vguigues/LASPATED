

function [gradientl]=oracleGradientMissingModel2(nbObservations,nbCalls,nb_missing_calls,neighbors,lambda,G,T,R,P,alpha,durations,S,probs,sample_missing_calls,sample_calls)

gradientl=zeros(T,G,R,P);
for g=1:G
    for t=1:T
        for i=1:R
            for p=1:P
                gradientl(t,g,i,p)=durations(g,t)*nbObservations(t,g,i,p);
                for n=1:nbObservations(t,g,1,p)
                    %Compute u(c,t,n)
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
                    u=u/S;
                    %Compute u'(c,t,i,n)
                    up=0;
                    for s=1:S
                        faux=exp(-lambda(t,g,i,p)*durations(g,t))*(lambda(t,g,i,p)^(Sample(s,i)+sample_calls(t,g,i,p,n)))*(-durations(g,t)*lambda(t,g,i,p)+Sample(s,i)+sample_calls(t,g,i,p,n));
                        iaux1=1;
                        for k=1:R
                            iaux1=iaux1*((durations(g,t))^(Sample(s,k)+sample_calls(t,g,k,p,n)));
                            denom=factorial(Sample(s,k)+sample_calls(t,g,k,p,n));
                            iaux1=iaux1/denom;
                        end
                        iaux2=1;
                        for k=1:R
                            if (k~=i)
                                iaux2=iaux2*(exp(-durations(g,t)*lambda(t,g,k,p)))*(lambda(t,g,k,p))^(Sample(s,k)+sample_calls(t,g,k,p,n));
                            end
                        end 
                        up=up+faux*iaux1*iaux2;
                    end
                    up=up/S;
                    gradientl(t,g,i,p)=gradientl(t,g,i,p)-(up/u);
                end
            end
        end
    end
end
