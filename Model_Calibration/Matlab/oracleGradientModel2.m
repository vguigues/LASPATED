             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computes the gradient of the objective for calibrating the model for 
%the process of calls. Model 1 of the paper ... is used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%gradient: gradient of the objective function at lambda.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gradient]=oracleGradientModel2(nbObservations,nbArrivals,x,C,D,T,R,regressor,nbRegressors,sizex,Groups,whichgroup,weight,durations,indexBeta)

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

if (length(Groups)>0)
    for c=1:C
        for d=1:D
            for t=1:T
                for m=1:length(Groups{1,whichgroup(d,t)})
                    d1=Groups{1,whichgroup(d,t)}(m).d;
                    t1=Groups{1,whichgroup(d,t)}(m).t;
                    if ((d1~=d)||(t1~=t))
                        for j=1:nbRegressors
                            gradient(indexBeta(c,d,t,j))=gradient(indexBeta(c,d,t,j))+(2*weight(whichgroup(d,t))/durations(t))*((gradient(indexBeta(c,d,t,j))/durations(t))-(gradient(indexBeta(c,d1,t1,j))/durations(t1)));
                        end
                    end
                end
            end
        end
    end
end



