              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computes the value of the objective for calibrating the model for 
%the process of calls. Model 2 of the paper ... is used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%f: value of the objective function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f]=oracleObjectiveModel2(nbObservations,nbArrivals,x,C,D,T,R,regressor,nbRegressors,Groups,weight,durations,indexBeta)

f=0;
for t=1:T
    for d=1:D
        for i=1:R
            for c=1:C
                rates=0;
                for j=1:nbRegressors
                    rates=rates+x(indexBeta(c,d,t,j))*regressor(j,i);
                end
                if (rates==0)
                    rates=0.0000001;
                end
                f=f+nbObservations(c,d,t,i)*rates-nbArrivals(c,d,t,i)*log(rates);
            end
        end
    end
end
 
if (length(Groups)>0)
    for c=1:C
        for m=1:length(Groups)
            for index1=1:length(Groups{1,m})
                d1=Groups{1,m}(index1).d;
                t1=Groups{1,m}(index1).t;
                for index2=1:length(Groups{1,m})
                    if (index2~=index1)
                        d2=Groups{1,m}(index2).d;
                        t2=Groups{1,m}(index2).t;
                        for j=1:nbRegressors
                            f=f+(weight(m)/2)*((x(indexBeta(c,d1,t1,j))/durations(t1))-(x(indexBeta(c,d2,t2,j))/durations(t2)))^2;
                        end
                    end
                end
            end
        end
    end
end
