
function [lambdas_missing,lambdas_nomissing,probs,prob]=lambda_missing_location(T,G,R,P,nbCalls,nbObservations,nb_missing_calls,durations)


lambdas_missing=zeros(T,G,R,P);
lambdas_nomissing=zeros(T,G,R,P);

%Calcul de prob
without_location=0;
with_location=0;
for t=1:T
    for g=1:G
        for p=1:P
            without_location=without_location+nb_missing_calls(t,g,p);
        end
    end
end

for t=1:T
    for g=1:G
        for r=1:R
            for p=1:P
                with_location=with_location+nbCalls(t,g,r,p);
            end
        end
    end
end

prob=without_location/(without_location+with_location);

%Calcul de probs

without_location=zeros(T,G,P);
with_location=zeros(T,G,P);
probs=zeros(T,G,P);
for t=1:T
    for g=1:G
        for p=1:P
            without_location(t,g,p)=without_location(t,g,p)+nb_missing_calls(t,g,p);
        end
    end
end

for t=1:T
    for g=1:G
        for r=1:R
            for p=1:P
                with_location(t,g,p)=with_location(t,g,p)+nbCalls(t,g,r,p);
            end
        end
    end
end

for t=1:T
    for g=1:G
        for p=1:P
            probs(t,g,p)=without_location(t,g,p)/(without_location(t,g,p)+with_location(t,g,p));
        end
    end
end


for t=1:T
    for g=1:G
        for r=1:R
            for p=1:P
                lambdas_nomissing(t,g,r,p)=nbCalls(t,g,r,p)/(nbObservations(t,g,r,p)*durations(t,g));
                lambdas_missing(t,g,r,p)=nbCalls(t,g,r,p)/((1-probs(t,g,p))*nbObservations(t,g,r,p)*durations(t,g));
            end
        end
    end
end
