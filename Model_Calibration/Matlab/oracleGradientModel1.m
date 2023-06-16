
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computes the gradient of the objective for calibrating the model for 
%the process of calls. Model 1 of the paper ... is used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nbObservations(p,i,t) is the number of time periods (t,g) for
%which we have observations for region i priority p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nbArrivals(p,i,t) is the number of calls for time period (t,g) region i
%priority p.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%neighbors{1,i} is the list of neighbors of region i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%type(i) is the type of region i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%distance(i,j) is the distance between centroids of regions i and j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%indexLambda(t,g,i,p) is the index of variable lambda(t,g,i,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lambda: the point at which we compute the objective.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T: number of time periods during the day.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%G: number of groups.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R: number of regions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%P: number of priorities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%gradient: gradient of the objective function at lambda.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gradient]=oracleGradientModel1(nbObservations,nbArrivals,neighbors,type,distance,lambda,T,R,P,alpha,durations,Groups,whichgroup,weight)

gradient=zeros(P,R,T);

for t=1:T
        for i=1:R
            for p=1:P
                currentLambda=lambda(p,i,t);
                gradComponent=nbObservations(p,i,t)*durations(t)-(nbArrivals(p,i,t)/currentLambda);
                for j=1:length(neighbors{1,i})
                    %Check if the type type(i) of i is the type type(neighbors{1,i}(j)) of neighbors{1,i}(j)
                    if (type(i)==type(neighbors{1,i}(j)))
                        %Add the regularized term with regions i and neighbors{1,i}(j)
                        gradComponent=gradComponent+2*alpha*(lambda(p,i,t)-lambda(p,neighbors{1,i}(j),t))/(distance(i,neighbors{1,i}(j)));
                    end
                end
                gradient(p,i,t)=gradComponent;
            end
        end
end

for t=1:T
        for i=1:R
            for p=1:P
                for j=1:length(Groups{1,whichgroup(t)})
                    tp=Groups{1,whichgroup(t)}(j);
                    if (tp~=t)
                        gradient(p,i,t)=gradient(p,i,t)+2*weight(whichgroup(t))*(lambda(p,i,t)-lambda(p,i,tp));
                    end
                end
            end
        end
end

