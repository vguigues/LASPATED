              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computes the value of the objective for calibrating the model for 
%the process of calls. Model 2 of the paper ... is used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%f: value of the objective function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f]=oracleObjectiveModel2(nbObservations,nbArrivals,x,C,D,T,R,regressor,nbRegressors,durations,indexBeta)

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
 
