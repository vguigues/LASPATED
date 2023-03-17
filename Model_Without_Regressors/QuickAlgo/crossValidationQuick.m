
%proportion: proportion of data used for validation
%(1-proportion): proportion of data used for calibration

%sampleCalls(t,g,i,p,j)

function [cputime,alpha,lambda]=crossValidationQuick(nbObervationsG,sampleCalls,neighbors,type,distance,T,G,R,P,sigma,betaBar,x,iterMax,proportion,epsilon)

%alphas=[0.001,0.01,0.05,0.1,0.5,1,2,5,10,50,100,1000];
alphas=[0.001,0.5,1];
maxLikelihood=inf;
nbObs=size(sampleCalls,5);
nbInBlock=floor(nbObs*proportion);
tic

for indexAlpha=1:length(alphas)
    indexAlpha
    likeLihood=0;
    for indexCross=1:floor((1/proportion))
        indexCross
        %Building the calibration data: t=(indexCross-1)*nbInBlock+1:indexCross*nbInBlock 
        nbObservationsCurrent=zeros(T,G,R,P);
        nbCallsCurrent=zeros(T,G,R,P);
        for index=(indexCross-1)*nbInBlock+1:indexCross*nbInBlock
            for t=1:T
                for g=1:G
                    for i=1:R
                        for p=1:P
                            if (index<=nbObervationsG(g))
                                nbObservationsCurrent(t,g,i,p)=nbObservationsCurrent(t,g,i,p)+1;
                                nbCallsCurrent(t,g,i,p)=nbCallsCurrent(t,g,i,p)+sampleCalls(t,g,i,p,index);
                            end
                        end
                    end
                end
            end
        end
        lambda=zeros(T,G,R,P);
        for t=1:T
            %t
            for g=1:G
                for p=1:P
                    [lambdaAux,fVal]=projectedGradientArmijoFeasibleQuick(nbObservationsCurrent,nbCallsCurrent,neighbors,type,distance,R,sigma,x,iterMax,alphas(indexAlpha),epsilon,t,g,p);
                    for i=1:R
                        lambda(t,g,i,p)=lambdaAux(i);
                    end
                end
            end
        end
        %Computing the likelihood on the remaining data
        nbCallsRemaining=zeros(T,G,R,P);
        for index=1:(indexCross-1)*nbInBlock
            for t=1:T
                for g=1:G
                    for i=1:R
                        for p=1:P
                            nbCallsRemaining(t,g,i,p)=nbCallsRemaining(t,g,i,p)+sampleCalls(t,g,i,p,index);
                        end
                    end
                end
            end
        end
        for index=(indexCross*nbInBlock+1):nbObs
            for t=1:T
                for g=1:G
                    for i=1:R
                        for p=1:P
                            if (index<=nbObervationsG(g))
                                nbCallsRemaining(t,g,i,p)=nbCallsRemaining(t,g,i,p)+sampleCalls(t,g,i,p,index);
                            end
                        end
                    end
                end
            end
        end
        f=0;
        for t=1:T
            for g=1:G
                for i=1:R
                    for p=1:P
                        currentLambda=lambda(t,g,i,p);
                        f=f+(nbObervationsG(g)-nbObservationsCurrent(t,g,i,p))*currentLambda-nbCallsRemaining(t,g,i,p)*log(currentLambda);
                    end
                end
            end
        end
        likeLihood=likeLihood+f;
    end
    if (likeLihood<maxLikelihood)
        maxLikelihood=likeLihood;
        bestAlpha=alphas(indexAlpha);
    end
end

alpha=bestAlpha;
nbObservationsCurrent=nbObs*ones(T,G,R,P);
nbCallsCurrent=zeros(T,G,R,P);
for index=1:nbObs
    for t=1:T
        for g=1:G
            for i=1:R
                for p=1:P
                    for index=1:nbObervationsG(g)
                        nbObservationsCurrent(t,g,i,p)=nbObservationsCurrent(t,g,i,p)+1;
                        nbCallsCurrent(t,g,i,p)=nbCallsCurrent(t,g,i,p)+sampleCalls(t,g,i,p,index);
                    end
                end
            end
        end
    end
end
cputime=toc;

lambda=zeros(T,G,R,P);
for t=1:T
    for g=1:G
        for p=1:P
            x=epsilon*ones(R,1);
            [lambdaAux,fVal]=projectedGradientArmijoBoundary(nbObservationsCurrent,nbCallsCurrent,neighbors,type,distance,R,sigma,betaBar,x,iterMax,alphas(indexAlpha),epsilon,t,g,p);
            for i=1:R
                lambda(t,g,i,p)=lambdaAux(i);
            end
        end
    end
end