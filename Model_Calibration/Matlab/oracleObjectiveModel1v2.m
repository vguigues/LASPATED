

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computes the value of the objective for calibrating the model for 
%the process of calls. Model 1 of the paper ... is used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nbObservations(p,i,t) is the number of time period t for
%which we have observations for region i priority p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nbArrivals(p,i,t) is the number of calls for time period t region i
%priority p.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%neighbors{1,i} is the list of neighbors of region i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%type(i) is the type of region i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%distance(i,j) is the distance between centroids of regions i and j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lambda: the point at which we compute the objective.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T: number of time periods during the day.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R: number of regions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%P: number of priorities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%f: value of the objective function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f]=oracleObjectiveModel1v2(nbObservations,nbArrivals,neighbors,type,distance,lambda,T,R,P,alpha,durations,Groups,whichgroup,weight)

f=0;
for t=1:T
        for i=1:R
            for p=1:P 
                currentLambda=lambda(p,i,t);
                f=f+nbObservations(p,i,t)*currentLambda*durations(t)-nbArrivals(p,i,t)*log(currentLambda*durations(t));
                for j=1:length(neighbors{1,i})
                    %Check if the type type(i) of i is the type type(neighbors{1,i}(j)) of neighbors{1,i}(j)
                    if (type(i)==type(neighbors{1,i}(j)))
                        %Add the regularized term with regions i and neighbors{1,i}(j)
                        f=f+(0.5*alpha)*nbObservations(p,i,t)*nbObservations(p,j,t)*((lambda(p,i,t)-lambda(p,neighbors{1,i}(j),t))^2)/(distance(i,neighbors{1,i}(j)));
                    end
                end
            end
        end
end

for i=1:R
    for p=1:P
        for grindex=1:length(Groups)
            for j=1:length(Groups{1,grindex})
                t=Groups{1,grindex}(j);
                for itp=1:length(Groups{1,grindex})
                    tp=Groups{1,grindex}(itp);
                    if (tp~=t)
                        f=f+(0.5*weight(grindex))*nbObservations(p,i,t)*nbObservations(p,i,tp)*((lambda(p,i,t)-lambda(p,i,tp))^2);
                    end
                end
            end
        end
    end
end



