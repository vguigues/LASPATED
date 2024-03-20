             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computes the gradient of the objective for calibrating the model for 
%the process of calls. Model 1 of the paper ... is used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%gradient: gradient of the objective function at lambda.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gradient]=oracleGradientModel2(nbObservations,nbArrivals,x,C,D,T,R,regressor,nbRegressors,sizex,durations,indexBeta)

gradient=zeros(sizex,1);
for c=1:C
    for d=1:D
        for t=1:T
            for i=1:R
                rates=0;
                for j=1:nbRegressors
                    rates=rates+x(indexBeta(c,d,t,j))*regressor(j,i);
                end
                if (rates==0)
                    rates=0.0000001;
                end
                for j=1:nbRegressors
                    gradient(indexBeta(c,d,t,j))=gradient(indexBeta(c,d,t,j))+nbObservations(c,d,t,i)*regressor(j,i)-nbArrivals(c,d,t,i)*regressor(j,i)/rates;
                end
            end
        end
    end
end




