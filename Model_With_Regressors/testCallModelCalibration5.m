
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

discretization_dir='C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\New_Data\Rect10x10_km';
[T,G,R,P,nbLandTypes,nbObservationsG,sample_calls,nbCalls,nbObservations,estimated,type,regressors,neighbors,distance]=read_calls_hol(discretization_dir);
C=P;
D=G;
indexBeta=zeros(C,D,T,nbLandTypes+1);

compteur=1;
for c=1:P
    for d=1:G
        for t=1:T
            for j=1:(nbLandTypes+1)
                indexBeta(c,d,t,j)=compteur;       
                compteur=compteur+1;
            end
        end
    end
end
sizex=P*G*T*(nbLandTypes+1);

Groups=cell(1,G*T);
D=G;
compteur=1;
whichgroup=zeros(D,t);
for g=1:G
    for t=1:T
        aux.d=g;
        aux.t=t;
        whichgroup(g,t)=compteur;
        Groups{1,compteur}=[aux];
    end
end

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
                for j=1:max(nbObservationsG)
                    sampleCalls{p,g,t,r}=[sampleCalls{p,g,t,r},sample_callsA(t,g,r,p,j)];
                end
            end
        end
    end
end

epsilon=10^(-5);
beta=0.5*ones(sizex,1);
sigma=0.5;
betaBar=1;
iterMax=30;
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
                    rateest=rateest+x(indexBeta(p,g,t,j))*regressor(j,r);
                end
                weeklambda(init+t)=weeklambda(init+t)+rateest;
            end
        end
    end
end
plot(weeklambda);

weeklambdaprioritites=zeros(T*G,P);
for t=1:T
    for g=1:G
        init=(g-1)*T;
        for r=1:R
            for p=1:P
                rat
                zeest=0;
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatial discretization for test problem in the plane:
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  * 
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  *
%  *    *   *  *  *  *  *  *  *  * 
% 1,2  2,2  *  *  *  *  *  *  *  * 
% 1,1  2,1 3,1 *  *  *  *  *  *  * 
% (Xi,Yi) below for region numbered i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmax=10;
ymax=10;
Nx=10;
Ny=10;
R=Nx*Ny;
C=1;
T=4;
D=7;
Groups=cell(1,2);
whichgroup=zeros(D,T);
for d=1:7
    whichgroup(d,1)=1;
end
for d=1:7
    whichgroup(d,2)=2;
end
for d=1:7
    whichgroup(d,3)=1;
end
for d=1:7
    whichgroup(d,4)=2;
end
Groups{1,1}=[];
for i=1:7
    aux.d=i;
    aux.t=1;
    Groups{1,1}=[Groups{1,1};aux];
end
Groups{1,2}=[];
for i=1:7
    aux.d=i;
    aux.t=2;
    Groups{1,2}=[Groups{1,2};aux];
end
for i=1:7
    aux.d=i;
    aux.t=3;
    Groups{1,1}=[Groups{1,1};aux];
end

for i=1:7
    aux.d=i;
    aux.t=4;
    Groups{1,2}=[Groups{1,2};aux];
end

%nbWeeks=52*10;
nbWeeks=8;
nbYears=floor(nbWeeks/52);
durations=6*ones(1,T);

nbLandTypes=2;
TheoreticalBetas=zeros(C,D,T,nbLandTypes+1);
regressor=zeros(1+nbLandTypes,R);
compteur=1;

for d=1:7
    TheoreticalBetas(1,d,2,1)=0.05;
    TheoreticalBetas(1,d,4,1)=0.05;
                 
    TheoreticalBetas(1,d,1,2)=6;
    TheoreticalBetas(1,d,2,2)=18;
    TheoreticalBetas(1,d,3,2)=6;
    TheoreticalBetas(1,d,4,2)=18;
    
    TheoreticalBetas(1,d,1,3)=3;
    TheoreticalBetas(1,d,2,3)=6;
    TheoreticalBetas(1,d,3,3)=3;
    TheoreticalBetas(1,d,4,3)=6;
end


for i=1:(Ny/2)
    for j=1:(Nx/2)
        type(compteur)=1;
        regressor(1,compteur)=100*rand;
        regressor(2,compteur)=0.5;
        regressor(3,compteur)=0.25;
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        type(compteur)=2;
        regressor(1,compteur)=50*rand;
        regressor(2,compteur)=0.25;
        regressor(3,compteur)=0.5;
        compteur=compteur+1;
    end
end

for i=(Ny/2)+1:Ny
    for j=1:(Nx/2)
        type(compteur)=2;
        regressor(1,compteur)=50*rand;
        regressor(2,compteur)=0.25;
        regressor(3,compteur)=0.5;
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        type(compteur)=1;
        regressor(1,compteur)=100*rand;
        regressor(2,compteur)=0.5;
        regressor(3,compteur)=0.25;
        compteur=compteur+1;
    end
end

sampleCalls=cell(C,D,T,R);
for c=1:C
    for d=1:D
        for t=1:T
            for i=1:R
                sampleCalls{c,d,t,i}=[];     
            end
        end
    end
end
nbObservations=zeros(C,D,T,R);
nbCalls=zeros(C,D,T,R);

C=1;

indexBeta=zeros(C,D,T,nbLandTypes+1);
compteur=1;
for c=1:C
    for d=1:D
        for t=1:T
            for j=1:(nbLandTypes+1)
                indexBeta(c,d,t,j)=compteur;       
                compteur=compteur+1;
            end
        end
    end
end
sizex=C*D*T*(nbLandTypes+1);

ratesemp=zeros(C,D,T,R);

for index=1:nbWeeks*7
    day=mod(index,7);
    if (day==0)
        day=7;
    end
    for c=1:C
        for t=1:T
            for i=1:R
                rate=0;
                for j=1:(nbLandTypes+1)
                    rate=rate+TheoreticalBetas(c,day,t,j)*regressor(j,i);
                end
                thisnbcall=poissrnd(rate,1,1);
                delta=randn(1,1);
                thisnbcall=max(floor(thisnbcall+delta,0);
                sampleCalls{c,day,t,i}=[sampleCalls{c,day,t,i},thisnbcall];
                
                nbObservations(c,day,t,i)=nbObservations(c,day,t,i)+1;
                nbCalls(c,day,t,i)=nbCalls(c,day,t,i)+thisnbcall;
                ratesemp(c,day,t,i)=ratesemp(c,day,t,i)+thisnbcall;
            end
        end
    end
end


for c=1:C
    for d=1:D
        for t=1:T
            for i=1:R
                ratesemp(c,d,t,i)=ratesemp(c,d,t,i)/nbObservations(c,d,t,i);
            end
        end
    end
end

for i=1:Nx
    absc(i)=xmax/(2*Nx)+(xmax/Nx)*(i-1);
end
for i=1:Ny
    ord(i)=ymax/(2*Ny)+(ymax/Ny)*(i-1);
end

% for i=1:Nx
%     for j=1:Ny
%         plot(absc(i),ord(j),'rd');
%         hold on;
%     end
% end

for i=1:R
    Xi=mod(i,Nx);
    if (Xi==0)
       Yi=i/Nx;
    else
       Yi=1+(i-Xi)/Nx;
    end
    if (Xi==0)
        Xi=Nx;
    end
    neighbors{1,i}=[];
    %Possible neighbors are (Xi,Yi-1),(Xi,Yi+1),
    %(Xi-1,Yi),(Xi+1,Yi),
    %(Xi-1,Yi+1),(Xi-1,Yi-1),
    %(Xi+1,Yi+1), and (Xi+1,Yi-1),
    if (Yi-1>=1)
        neighbors{1,i}=[neighbors{1,i},(Yi-2)*Nx+Xi];     
    end
    if (Yi+1<=Ny)
        neighbors{1,i}=[neighbors{1,i},Yi*Nx+Xi];     
    end
    if (Xi-1>=1)
        neighbors{1,i}=[neighbors{1,i},(Yi-1)*Nx+Xi-1];     
    end
    if (Xi+1<=Nx)
        neighbors{1,i}=[neighbors{1,i},(Yi-1)*Nx+Xi+1];     
    end
    if ((Yi+1<=Ny)&&(Xi-1>=1))
        neighbors{1,i}=[neighbors{1,i},Yi*Nx+Xi-1];     
    end
    if ((Yi-1>=1)&&(Xi-1)>=1)
        neighbors{1,i}=[neighbors{1,i},(Yi-2)*Nx+Xi-1];     
    end
    if ((Yi+1<=Ny)&&(Xi+1)<=Nx)
        neighbors{1,i}=[neighbors{1,i},Yi*Nx+Xi+1];     
    end
    if ((Yi-1>=1)&&(Xi+1)<=Nx)
        neighbors{1,i}=[neighbors{1,i},(Yi-2)*Nx+Xi+1];     
    end
    for j=1:R
        Xj=mod(j,Nx);
        if (Xj==0)
          Yj=j/Nx;
        else
          Yj=1+((j-Xj)/Nx);
        end
        if (Xj==0)
            Xj=Nx;
        end
        distance(i,j)=norm([absc(Xi)-absc(Xj);ord(Yi)-ord(Yj)]);
    end
end

% obsdays=zeros(D,1);
% obsbefore=zeros(nbWeeks*7);
% for index=1:nbWeeks*7
%     day=mod(index,7);
%     if (day==0)
%         day=7;
%     end
%     if (isHolidays{1,index}.b)
%         day=7+isHolidays{1,index}.k;
%     end
%     obsdays(day)=obsdays(day)+1;
%     obsbefore(index)=obsdays(day);  
% end

epsilon=10^(-8);
%beta=2*epsilon*ones(sizex,1);
beta=0.5*ones(sizex,1);
sigma=0.5;
betaBar=1;
iterMax=30;
weight=10;
lbounds=zeros(sizex,1);

%First projected gradient
[x,fVal]=projectedGradientArmijoFeasible2(nbObservations,nbCalls,beta,C,D,T,R,weight,regressor,nbLandTypes,iterMax,sigma,epsilon,indexBeta,Groups,durations,sizex,whichgroup,lbounds);

meanrateBeta=0;
rates=[];
ratesest=[];
rateemp=[];
meanrateEmp=0;
%cp=0;
for d=1:D
    for t=1:T
        for i=1:R
            for c=1:C
                rate=0;
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rate=rate+TheoreticalBetas(c,d,t,j)*regressor(j,i);
                    rateest=rateest+x(indexBeta(c,d,t,j))*regressor(j,i);
                    
                    %if (TheoreticalBetas(c,d,t,j)>0)
                    %   meanrateBeta=meanrateBeta+100*(abs(x(indexBeta(c,d,t,j))-TheoreticalBetas(c,d,t,j))/TheoreticalBetas(c,d,t,j));
                    %   cp=cp+1;
                    %end
                end
                rate
                rateest
                rates=[rates,rate];
                ratesest=[ratesest,rateest];    
                meanrateBeta=meanrateBeta+100*(abs(rate-rateest)/rate);
                meanrateEmp=meanrateEmp+100*(abs(rate-ratesemp(c,d,t,i))/rate);
            end
        end
    end
end
%meanrateBeta=meanrateBeta/cp
meanrateBeta=meanrateBeta/(D*T*R*C)
meanrateEmp=meanrateEmp/(D*T*R*C)


ratesT1=[];
ratesE1=[];
for d=1:D
    for t=1:T
        rateTaux=0;
        rateEaux=0;
        for i=1:R
            for c=1:C
                raT=0;
                raE=0;
                for j=1:(1+nbLandTypes)
                    raT=raT+10*TheoreticalBetas(c,d,t,j)*regressor(j,i);
                    raE=raE+x(indexBeta(c,d,t,j))*regressor(j,i);
                end
                rateTaux=rateTaux+raT;
                rateEaux=rateEaux+raE;
            end
        end
        ratesT1=[ratesT1;rateTaux/(R*C)];
        ratesE1=[ratesE1;rateEaux/(R*C)];
    end
end
plot(ratesT1,'-k');
hold on;
plot(ratesE1,'--r');
hold on;
plot(ratesE0,'-.b');
legend(['Theoretical  ';'Covariates   ';'No covariates'])
mean(100*abs((ratesT1-ratesE1)./ratesT1));

 TheoreticalBetas(1,d,2,1)=0.05;
    TheoreticalBetas(1,d,4,1)=0.05;
                 
    TheoreticalBetas(1,d,1,2)=6;
    TheoreticalBetas(1,d,2,2)=18;
    TheoreticalBetas(1,d,3,2)=6;
    TheoreticalBetas(1,d,4,2)=18;
    
    TheoreticalBetas(1,d,1,3)=3;
    TheoreticalBetas(1,d,2,3)=6;
    TheoreticalBetas(1,d,3,3)=3;
    TheoreticalBetas(1,d,4,3)=6;

Relerr=[];
for d=1:D
    for t=1:T
        for j=1:(1+nbLandTypes)
            for c=1:C
                Relerr=[Relerr,abs(TheoreticalBetas(c,d,t,j)-x(indexBeta(c,d,t,j)))];
            end
        end
    end
end
Relerr=Relerr/(C*D*T*(1+nbLandTypes)); 
sample_callsA=sampleCalls;
nbCallsA=nbCalls;
nbObservationsA=nbObservations;
nbObservations=zeros(D*T,R,P);
nbCalls=zeros(D*T,R,P);
sampleCalls=zeros(D*T,R,P,max(nbObservationsG));

for r=1:R
    r
    for p=1:P
        p
        for g=1:D
            g
            init=(g-1)*T;
            for t=1:T
                t
                nbObservations(init+t,r,p)=nbObservationsA(p,g,t,r);
                %size(nbObservations)
                %init+t
                nbCalls(init+t,r,p)=nbCallsA(p,g,t,r);
                for j=1:nbObservationsG(g)
                    sampleCalls(init+t,r,p,j)=sample_callsA{p,g,t,r}(j);
                end
            end
        end
    end
end

T=D*T;
Groups=cell(1,T);
for t=1:T
    Groups{1,t}=[t];
    whichgroup(t)=t;
end



for ind=1:6
alpha=10^(ind-1)
epsilon=0.001;
xnr=ones(T,R,P);
sigma=0.5;
betaBar=1;
iterMax=30;
durations=6*ones(1,T);
weight=0;
proportion=0.2;

%Call model calibration without regressors

[xnr,fVal]=projectedGradientArmijoFeasible(nbObservations,nbCalls,neighbors,type,distance,T,R,P,sigma,xnr,iterMax,alpha,epsilon,durations,Groups,whichgroup,weight)
ratesE0=[];
for t=1:T
    rateEaux=0;
    for i=1:R
        for c=1:C
            rateEaux=rateEaux+xnr(t,i,c)*6;
        end
    end
    ratesE0=[ratesE0;rateEaux/(R*C)];
end
errs(ind)=mean(100*abs((ratesT1-ratesE0)./ratesT1));
end
[cputime,alpha,weight,lambda]=crossValidation(8,sampleCalls,neighbors,type,distance,T,R,P,sigma,betaBar,x,iterMax,proportion,epsilon,durations,Groups,whichgroup);

ratesE0=[];
for t=1:T
    rateEaux=0;
    for i=1:R
        for c=1:C
            rateEaux=rateEaux+lambda(t,i,c)*6;
        end
    end
    ratesE0=[ratesE0;rateEaux/(R*C)];
end
mean(100*abs((ratesT1-ratesE0)./ratesT1));
