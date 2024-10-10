

function [fL,fp]=oracleObjectiveMissing(nbObservations,nbCalls,nb_missing_calls,neighbors,lambda,prob,T,R,P,G,alpha,durations,Groups,whichgroup,weight)

fL=0;
for g=1:G
    for t=1:T
        for i=1:R
            for p=1:P 
                currentLambda=lambda(t,g,i,p);
                fL=fL+nbObservations(t,g,i,p)*currentLambda*durations(g,t)-nbCalls(t,g,i,p)*log(currentLambda);
                for j=1:length(neighbors{1,i})
                    fL=fL+(0.5*alpha)*nbObservations(t,g,i,p)*nbObservations(t,g,neighbors{1,i}(j),p)*((lambda(t,g,i,p)-lambda(t,g,neighbors{1,i}(j),p))^2);
                end
            end
        end
   end
end

fp=0;
for g=1:G
    for t=1:T
        for p=1:P
            somme=0;
            somme2=0;
            for i=1:R
                somme=somme+lambda(t,g,i,p);
                somme2=somme2+nbCalls(t,g,i,p);
            end
            fL=fL-nb_missing_calls(t,g,p)*log(somme);
            fp=fp-nb_missing_calls(t,g,p)*log(prob(p,g,t))-somme2*log(1-prob(p,g,t));
        end
    end
end

for i=1:R
    for p=1:P
        for grindex=1:length(Groups)
            for j=1:length(Groups{1,grindex})
                struct=Groups{1,grindex}(j);
                thisg=struct.g;
                thist=struct.t;
                for itp=1:length(Groups{1,grindex})
                    structtp=Groups{1,grindex}(itp);
                    tpg=structtp.g;
                    tpt=structtp.t;
                    fL=fL+(0.5*weight(grindex))*nbObservations(thist,thisg,i,p)*nbObservations(tpt,tpg,i,p)*((lambda(thist,thisg,i,p)-lambda(tpt,tpg,i,p))^2);
                end
            end
        end
    end
end


for p=1:P
    for grindex=1:length(Groups)
        for j=1:length(Groups{1,grindex})
            struct=Groups{1,grindex}(j);
            thisg=struct.g;
            thist=struct.t;
            for itp=1:length(Groups{1,grindex})
                structtp=Groups{1,grindex}(itp);
                tpg=structtp.g;
                tpt=structtp.t;
                fp=fp+(0.5*weight(grindex))*nbObservations(thist,thisg,1,p)*nbObservations(tpt,tpg,1,p)*((prob(p,thisg,thist)-prob(p,tpg,tpt))^2);
            end
        end
    end
end





