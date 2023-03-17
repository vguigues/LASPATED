
%proportion: proportion of data used for validation
%(1-proportion): proportion of data used for calibration

%sample(t,g,i,p,j)

function [cputime,weight,lambda]=crossValidationNoReg(sample,neighbors,type,distance,T,R,P,sigma,x,iterMax,proportion,epsilon,durations,Groups,whichgroup,weights)

tic
%alphas=zeros(1,8);
nbObservationsTotal=size(sample,4);
%weights=[0,0.01,0.05,0.01,0.1,1,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400];
%alphas=[0,0.01,0.05,0.01,0.1,1,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400];

%weights=[0,0.1,1,10,50,100,150,200];
alphas=weights;

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
        [lambda,fVal]=projectedGradientArmijoFeasible(nbObservationsCurrent,nbCallsCurrent,neighbors,type,distance,T,R,P,sigma,x,iterMax,alphas(indexAlpha),epsilon,durations,Groups,whichgroup,weights(indexAlpha));

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
[lambda,fVal]=projectedGradientArmijoFeasible(nbObservationsCurrent,nbCallsCurrent,neighbors,type,distance,T,R,P,sigma,x,iterMax,alphas(indexAlpha),epsilon,durations,Groups,whichgroup,weights(indexAlpha));
cputime=toc;

           