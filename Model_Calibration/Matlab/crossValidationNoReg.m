
%proportion: proportion of data used for validation
%(1-proportion): proportion of data used for calibration

function [cputime,weight,lambda]=crossValidationNoReg(sample,neighbors,type,distance,T,R,P,sigma,x,iterMax,proportion,epsilon,durations,Groups,whichgroup,weights, delta, uppx, lowx)

disp('Begin cross validation');
tic
nbObservationsTotal=size(sample,4);
alphas=weights;

nbGroups=length(Groups);
minLoss=inf;
cputime=0;

nbInBlock=floor(nbObservationsTotal*proportion);
 
for indexAlpha=1:length(alphas)
    indexAlpha
    likeLihood=0;
    for indexCross=1:floor((1/proportion))
        indexCross
        %Building the calibration data: j=(indexCross-1)*nbInBlock+1:indexCross*nbInBlock 
        nbObservationsCurrent=zeros(P,R,T);
        nbCallsCurrent=zeros(P,R,T);
        for index=(indexCross-1)*nbInBlock+1:indexCross*nbInBlock
            for t=1:T
                    for i=1:R
                        for p=1:P
                                nbObservationsCurrent(p,i,t)=nbObservationsCurrent(p,i,t)+1;
                                nbCallsCurrent(p,i,t)=nbCallsCurrent(p,i,t)+sample(t,i,p,index);
                        end
                    end
            end
        end
        x=epsilon*ones(P,R,T);
        [lambda,fVal]=projectedGradientArmijoFeasiblev2(nbObservationsCurrent,nbCallsCurrent,neighbors,type,distance,T,R,P,sigma,x,iterMax,alphas(indexAlpha),epsilon,durations,Groups,whichgroup,weights(indexAlpha)*ones(1,nbGroups), delta, uppx, lowx);

        %Computing the likelihood on the remaining data
        nbCallsRemaining=zeros(P,R,T);
        for index=1:(indexCross-1)*nbInBlock
            for t=1:T
                    for i=1:R
                        for p=1:P
                            nbCallsRemaining(p,i,t)=nbCallsRemaining(p,i,t)+sample(t,i,p,index);                            
                        end
                    end
            end
        end
        for index=(indexCross*nbInBlock+1):nbObservationsTotal
            for t=1:T
                    for i=1:R
                        for p=1:P
                                nbCallsRemaining(p,i,t)=nbCallsRemaining(p,i,t)+sample(t,i,p,index);
                        end
                    end
            end
        end
        
        f=0;
        for t=1:T
                for i=1:R
                    for p=1:P
                        currentLambda=lambda(p,i,t);
                        f=f+(nbObservationsTotal-nbInBlock)*currentLambda*durations(t)-nbCallsRemaining(p,i,t)*log(currentLambda);
                    end
                end
        end
        likeLihood=likeLihood+f;
    end
    likeLihood=likeLihood/(floor((1/proportion)));
    if (likeLihood<minLoss)
        minLoss=likeLihood;
        bestAlpha=alphas(indexAlpha);
        bestWeight=weights(indexAlpha);
    end
end

alpha=bestAlpha;
weight=bestWeight;
nbObservationsCurrent=zeros(P,R,T);
nbCallsCurrent=zeros(P,R,T);
for t=1:T
        for i=1:R
            for p=1:P
                for index=1:nbObservationsTotal
                    nbObservationsCurrent(p,i,t)=nbObservationsCurrent(p,i,t)+1;
                    nbCallsCurrent(p,i,t)=nbCallsCurrent(p,i,t)+sample(t,i,p,index);
                end
            end
        end
end
[lambda,fVal]=projectedGradientArmijoFeasiblev2(nbObservationsCurrent,nbCallsCurrent,neighbors,type,distance,T,R,P,sigma,x,iterMax,bestAlpha,epsilon,durations,Groups,whichgroup,bestWeight*ones(1,nbGroups),delta, uppx, lowx);
cputime=toc;

           