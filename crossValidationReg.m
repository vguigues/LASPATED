
%proportion: proportion of data used for validaiton
%(1-proportion): proportion of data used for calibration


function [cputime,weight,x]=crossValidationReg(sample,indexBeta,x,C,D,T,R,regressor,nbLandTypes,proportion,sizex,iterMax,sigma,epsilon,nbObs,isHolidays,obsbefore,Groups,durations,whichgroup,weights)

tic
maxLikelihood=inf;
nbInBlock=floor(nbObs*proportion);

for indexweights=1:length(weights)
    indexweights
    
    likeLihood=0;
    for indexCross=1:floor((1/proportion))
        indexCross
        %Building the calibration data: t=(indexCross-1)*nbInBlock+1:indexCross*nbInBlock 
        nbObservationsCurrent=zeros(C,D,T,R);
        nbCallsCurrent=zeros(C,D,T,R);
        
        for index=(indexCross-1)*nbInBlock+1:indexCross*nbInBlock
            for c=1:C
                for d=1:D
                    for t=1:T
                        for i=1:R
                            day=mod(index,7);
                            if (day==0)
                                day=7;
                            end
                            if (isHolidays{1,index}.b)
                                day=7+isHolidays{1,index}.k;
                            end
                            nbObservationsCurrent(c,day,t,i)=nbObservationsCurrent(c,day,t,i)+1;
                            nbCallsCurrent(c,day,t,i)=nbCallsCurrent(c,day,t,i)+sample{c,day,t,i}(obsbefore(index));      
                        end
                    end
                end
            end
        end    
        beta=epsilon*ones(sizex,1);

        [x,fVal]=projectedGradientArmijoFeasible2(nbObservationsCurrent,nbCallsCurrent,beta,C,D,T,R,weights(indexweights),regressor,nbLandTypes,iterMax,sigma,epsilon,indexBeta,Groups,durations,sizex,whichgroup);

        %Computing the likelihood on the remaining data
        nbObservationsRemaining=zeros(C,D,T,R);
        nbCallsRemaining=zeros(C,D,T,R);        
        for index=1:(indexCross-1)*nbInBlock
            for c=1:C
                for d=1:D
                    for t=1:T
                        for i=1:R
                            day=mod(index,7);
                            if (day==0)
                                day=7;
                            end
                            if (isHolidays{1,index}.b)
                                day=7+isHolidays{1,index}.k;
                            end
                            nbObservationsRemaining(c,day,t,i)=nbObservationsRemaining(c,day,t,i)+1;
                            nbCallsRemaining(c,day,t,i)=nbCallsRemaining(c,day,t,i)+sample{c,day,t,i}(obsbefore(index));      
                        end
                    end
                end
            end
        end
        for index=(indexCross*nbInBlock+1):nbObs
            for c=1:C
                for d=1:D
                    for t=1:T
                        for i=1:R
                            day=mod(index,7);
                            if (day==0)
                                day=7;
                            end
                            if (isHolidays{1,index}.b)
                                day=7+isHolidays{1,index}.k;
                            end
                            nbObservationsRemaining(c,day,t,i)=nbObservationsRemaining(c,day,t,i)+1;
                            nbCallsRemaining(c,day,t,i)=nbCallsRemaining(c,day,t,i)+sample{c,day,t,i}(obsbefore(index));      
                        end
                    end
                end
            end
        end        
        f=0;
        for t=1:T
            for d=1:D
                for i=1:R
                    for c=1:C
                        rates=0;
                        for j=1:(1+nbLandTypes)
                            rates=rates+x(indexBeta(c,d,t,j))*regressor(j,i);
                        end
                        f=f+nbObservationsRemaining(c,d,t,i)*rates-nbCallsRemaining(c,d,t,i)*log(rates);
                    end
                end
            end
        end
        likeLihood=likeLihood+f;
    end
    likeLihood
    if (likeLihood<maxLikelihood)
        maxLikelihood=likeLihood;
        bestweight=weights(indexweights);
    end
end

weight=bestweight;
nbObservations=zeros(C,D,T,R);
nbCalls=zeros(C,D,T,R);
      
for index=1:nbObs
    for t=1:T
        for d=1:D
            for i=1:R
                for c=1:C
                    day=mod(index,7);
                    if (day==0)
                        day=7;
                    end
                    if (isHolidays{1,index}.b)
                        day=7+isHolidays{1,index}.k;
                    end
                    nbObservations(c,day,t,i)=nbObservations(c,day,t,i)+1;
                    nbCalls(c,day,t,i)=nbCalls(c,day,t,i)+sample{c,day,t,i}(obsbefore(index));      
                end
            end
        end
    end
end
[x,fVal]=projectedGradientArmijoFeasible2(nbObservations,nbCalls,x,C,D,T,R,weight,regressor,nbLandTypes,iterMax,sigma,epsilon,indexBeta,Groups,durations,sizex,whichgroup);
cputime=toc;

        
        