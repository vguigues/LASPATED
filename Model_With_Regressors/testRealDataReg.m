
%Hex7_km
%In line 25, r=24, nan in regressors changed to content of regressors for
%line 24 
%In line 282, r=281, nan in regressors changed to content of regressors for
%line 281 

cd 'C:\gurobi901\win64\matlab\'
gurobi_setup
addpath 'C:\Program Files\Mosek\9.0'
addpath 'C:\Program Files\Mosek\9.0\toolbox\r2015a'
addpath 'C:\Users\vince\Dropbox\Softwares\Heuristic_Ambulance_Dispatch\Call_Model_Calibration';
addpath 'C:\Users\vince\Dropbox\Softwares\Heuristic_Ambulance_Dispatch\Call_Model_Calibration\Model_With_Regressors';

%discretization_dir='C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\New_Data\Hex7_km';
discretization_dir='C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\New_Data\Rect10x10_km';
discretization_dir='C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\New_Data\Districts_km';

%[T,G,R,P,nbLandTypes,nbObservationsG,sample_calls,nbCalls,nbObservations,estimated,type,regressors,neighbors,distance]=read_calls_hol(discretization_dir);
[T,G,R,P,nbLandTypes,nbObservationsG,sample_calls,nbCalls,nbObservations,estimated,type,regressors,neighbors,distance]=read_calls_reg(discretization_dir);
epsilon=10^(-5);
[t1,t2]=size(regressors);
for i=1:t2
    if (sum(regressors(:,i))==0)
        regressors(4,i)=epsilon;
    end
end

C=P;
D=G;
%indexBeta=zeros(C,D,T,nbLandTypes+1);
indexBeta=zeros(C,D,T,1+nbLandTypes);

compteur=1;
for c=1:P
    for d=1:G
        for t=1:T
%           for j=1:(nbLandTypes+1)
            for j=1:(1+nbLandTypes)
                indexBeta(c,d,t,j)=compteur;       
                compteur=compteur+1;
            end
        end
    end
end
%sizex=P*G*T*(nbLandTypes+1);
sizex=P*G*T*(1+nbLandTypes);

Groups=cell(1,G*T);
D=G;
compteur=1;
whichgroup=zeros(D,T);
for g=1:G
    for t=1:T
        aux.d=g;
        aux.t=t;
        whichgroup(g,t)=compteur;
        Groups{1,compteur}=[aux];
    end
end

%Output of the read data function
%sample_calls(T,G,R,P,max(nbObservationsG));
%nbCalls(T,G,R,P);
%nbObservations(T,G,R,P);
    
sample_callsA=sample_calls;
nbCallsA=nbCalls;
nbObservationsA=nbObservations;
C=P;
nbObservations=zeros(P,G,T,R);
nbCalls=zeros(P,G,T,R);
sampleCalls=cell(P,G,T,R);
G=D;
for r=1:R
    for p=1:P
        for g=1:G
            for t=1:T
                nbObservations(p,g,t,r)=nbObservationsA(t,g,r,p);
                %size(nbObservations)
                %init+t
                nbCalls(p,g,t,r)=nbCallsA(t,g,r,p);
                for j=1:min(nbObservationsG)
                    sampleCalls{p,g,t,r}=[sampleCalls{p,g,t,r},sample_callsA(t,g,r,p,j)];
                end
            end
        end
    end
end

epsilon=10^(-);
beta=0.5*ones(sizex,1);
sigma=0.5;
betaBar=1;
iterMax=15;
weight=0;
lbounds=zeros(sizex,1);
durations=0.5*ones(1,T);
regressor=regressors;
%First projected gradient
[x,fVal]=projectedGradientArmijoFeasible2(nbObservations,nbCalls,beta,C,D,T,R,weight,regressor,nbLandTypes,iterMax,sigma,epsilon,indexBeta,Groups,durations,sizex,whichgroup,lbounds);

weeklambda=zeros(T*G,1);
for t=1:T
    for g=1:G
        init=(g-1)*T;
        for r=1:R
            for p=1:P
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rateest=rateest+x(indexBeta(p,g,t,j))*regressors(j,r);
                end
                weeklambda(init+t)=weeklambda(init+t)+rateest;
            end
        end
    end
end
hold on;
plot(weeklambda./0.5*ones(1,T*G),':m');




weeklambdaprioritites=zeros(T*G,P);
for t=1:T
    for g=1:G
        init=(g-1)*T;
        for r=1:R
            for p=1:P
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rateest=rateest+x(indexBeta(p,g,t,j))*regressor(j,r);
                end
                weeklambdaprioritites(init+t,p)=weeklambdaprioritites(init+t,p)+rateest;
            end
        end
    end
end

plot(weeklambdaprioritites(:,1));
hold on;
plot(weeklambdaprioritites(:,2));
hold on;
plot(weeklambdaprioritites(:,3));
legend(['Priority 1';'Priority 2';'Priority 3']);

%Write sum of lambda, summing in t,p.
heatlambdas=zeros(1,R);
for r=1:R
    for g=1:G
        for t=1:T
            for p=1:P
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rateest=rateest+x(indexBeta(p,g,t,j))*regressor(j,r);
                end
                heatlambdas(r)=heatlambdas(r)+rateest;
            end
        end
    end
end
dlmwrite('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\heatLRegRect10x10.txt',heatlambdas);

heatlambdaspriorities=zeros(P,R);
for p=1:P
    for r=1:R
        for g=1:G
            for t=1:T
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rateest=rateest+x(indexBeta(p,g,t,j))*regressor(j,r);
                end
                heatlambdaspriorities(p,r)=heatlambdaspriorities(p,r)+rateest;
            end
        end
    end
end
dlmwrite('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\heatLPRegRect10x10.txt',heatlambdaspriorities);

