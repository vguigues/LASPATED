

function [gradientl,gradientp]=oracleGradientMissing(nbObservations,nbCalls,nb_missing_calls,neighbors,lambda,prob,G,T,R,P,alpha,durations,Groups,whichgroup,weight)

gradientl=zeros(T,G,R,P);
gradientp=zeros(P,G,T);

somme=zeros(P,G,T);
for g=1:G
    for t=1:T
        for p=1:P
            thissum=0;
            for r=1:R
                thissum=thissum+lambda(t,g,r,p);
            end
            somme(p,g,t)=thissum;
        end
    end
end

for g=1:G
    for t=1:T
        for p=1:P
            gradComponent=nbObservations(t,g,1,p)*durations(g,t)-nb_missing_calls(t,g,p)/somme(p,g,t);
            for r=1:R
                currentLambda=lambda(t,g,r,p);
                gradComponent1=gradComponent-(nbCalls(t,g,r,p)/currentLambda);
                for j=1:length(neighbors{1,r})
                    gradComponent1=gradComponent1+2*alpha*nbObservations(t,g,r,p)*nbObservations(t,g,neighbors{1,r}(j),p)*(lambda(t,g,r,p)-lambda(t,g,neighbors{1,r}(j),p));
                end
                gradientl(t,g,r,p)=gradComponent1;
            end
            
        end
    end
end

for g=1:G
    for t=1:T
        for i=1:R
            for p=1:P
                for j=1:length(Groups{1,whichgroup(t,g)})
                    str=Groups{1,whichgroup(t,g)}(j);
                    gp=str.g;
                    tp=str.t;
                    gradientl(t,g,i,p)=gradientl(t,g,i,p)+2*weight(whichgroup(t,g))*nbObservations(t,g,1,p)*nbObservations(tp,gp,1,p)*(lambda(t,g,i,p)-lambda(tp,gp,i,p));
                end
            end
        end
    end
end


for g=1:G
    for t=1:T
            for p=1:P
                somme=0;
                for r=1:R
                    somme=somme+nbCalls(t,g,r,p);
                end
                gradientp(p,g,t)=-nb_missing_calls(t,g,p)/prob(p,g,t)+((somme)/(1-prob(p,g,t)));              
            end
        end
 end

for g=1:G
     for t=1:T
         for p=1:P
             for j=1:length(Groups{1,whichgroup(t,g)})
                 str=Groups{1,whichgroup(t,g)}(j);
                 gp=str.g;
                 tp=str.t;
                 gradientp(p,g,t)=gradientp(p,g,t)+2*weight(whichgroup(t,g))*nbObservations(t,g,1,p)*nbObservations(tp,gp,1,p)*(prob(p,g,t)-prob(p,gp,tp));
             end
         end
     end
 end
