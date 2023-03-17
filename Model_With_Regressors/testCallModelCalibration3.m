
cd 'C:\gurobi901\win64\matlab\'
gurobi_setup
addpath 'C:\Program Files\Mosek\9.0'
addpath 'C:\Program Files\Mosek\9.0\toolbox\r2015a'
addpath 'C:\Users\vince\Dropbox\Softwares\Heuristic_Ambulance_Dispatch\Call_Model_Calibration\Model_With_Regressors';

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
D=15;
Groups=cell(1,4);
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
for d=8:15
    whichgroup(d,1)=3;
end
for d=8:15
    whichgroup(d,2)=4;
end
for d=8:15
    whichgroup(d,3)=3;
end
for d=8:15
    whichgroup(d,4)=4;
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
Groups{1,3}=[];
for i=8:15
    aux.d=i;
    aux.t=1;
    Groups{1,3}=[Groups{1,3};aux];
end
Groups{1,4}=[];
for i=8:15
    aux.d=i;
    aux.t=2;
    Groups{1,4}=[Groups{1,4};aux];
end
for i=8:15
    aux.d=i;
    aux.t=3;
    Groups{1,3}=[Groups{1,3};aux];
end
for i=8:15
    aux.d=i;
    aux.t=4;
    Groups{1,4}=[Groups{1,4};aux];
end
nbWeeks=52*4;
nbYears=floor(nbWeeks/52);
durations=6*ones(1,T);
nbHolidaysYear=8;
isHolidays=cell(1,nbWeeks*7);
for j=1:nbWeeks*7
    isHolidays{1,j}.b=0;
    isHolidays{1,j}.k=0;
end
daysH=[1:1:8];
indexH=[8:15];


for year=1:nbYears
    for k=1:length(daysH)
        isHolidays{1,((year-1)*52)*7+daysH(k)}.b=1;
        isHolidays{1,((year-1)*52)*7+daysH(k)}.k=k;
    end
end

for k=1:length(daysH)
    if (nbYears*52*7+daysH(k)<=nbWeeks*7)
          isHolidays{1,nbYears*52*7+daysH(k)}.b=1;
          isHolidays{1,nbYears*52*7+daysH(k)}.k=k;
   end
end

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

for d=8:15
    TheoreticalBetas(1,d,2,1)=0.2;
    TheoreticalBetas(1,d,4,1)=0.2;
    
    TheoreticalBetas(1,d,1,2)=24;
    TheoreticalBetas(1,d,2,2)=72;
    TheoreticalBetas(1,d,3,2)=24;
    TheoreticalBetas(1,d,4,2)=72;
    
    TheoreticalBetas(1,d,1,3)=12;
    TheoreticalBetas(1,d,2,3)=24;
    TheoreticalBetas(1,d,3,3)=12;
    TheoreticalBetas(1,d,4,3)=24;
end

for i=1:(Ny/2)
    for j=1:(Nx/2)
        type(compteur)=1;
        regressor(1,compteur)=50+50*rand;
        %100*rand;
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
        regressor(1,compteur)=50+50*rand;
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

c=1;

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
for index=1:nbWeeks*7
    day=mod(index,7);
    if (day==0)
        day=7;
    end
    if (isHolidays{1,index}.b)
        day=7+isHolidays{1,index}.k;
    end
        
    for c=1:C
        for t=1:T
            for i=1:R
                rate=0;
                for j=1:(nbLandTypes+1)
                    rate=rate+10*TheoreticalBetas(c,day,t,j)*regressor(j,i);
                end
                thisnbcall=poissrnd(rate,1,1);
                thisnbcall=thisnbcall+randn;
                sampleCalls{c,day,t,i}=[sampleCalls{c,day,t,i},thisnbcall];
                nbObservations(c,day,t,i)=nbObservations(c,day,t,i)+1;
                nbCalls(c,day,t,i)=nbCalls(c,day,t,i)+thisnbcall;
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

epsilon=10^(-6);
beta=2*epsilon*ones(sizex,1);
sigma=0.5;
betaBar=1;
iterMax=30;
weights=[0];
lbounds=zeros(sizex,1);

%First projected gradient
errs=[];
for k=1:length(weights)
weight=weights(k);
[x,fVal]=projectedGradientArmijoFeasible2(nbObservations,nbCalls,beta,C,D,T,R,weight,regressor,nbLandTypes,iterMax,sigma,epsilon,indexBeta,Groups,durations,sizex,whichgroup,lbounds);

meanrateBeta=0;
rates=[];
ratesest=[];
%cp=0;
for d=1:D
    for t=1:T
        rate=0;
        rateest=0;        
        for i=1:R
            for c=1:C
                for j=1:(1+nbLandTypes)
                    rate=rate+TheoreticalBetas(c,d,t,j)*regressor(j,i);
                    rateest=rateest+x(indexBeta(c,d,t,j))*regressor(j,i);
                    %if (TheoreticalBetas(c,d,t,j)>0)
                    %   meanrateBeta=meanrateBeta+100*(abs(x(indexBeta(c,d,t,j))-TheoreticalBetas(c,d,t,j))/TheoreticalBetas(c,d,t,j));
                    %   cp=cp+1;
                    %end
                end
            end
        end
        rates=[rates,rate/(C*R)];
        ratesest=[ratesest,rateest/(C*R)];    
        %meanrateBeta=meanrateBeta+100*(abs(rate-rateest)/rate);
    end
end
mean(100*abs(rates-ratesest)./rates)
end

%meanrateBeta=meanrateBeta/cp
meanrateBeta=meanrateBeta/(7*T*R*C)
ms(k)=meanrateBeta;
end
plot(rates);
hold on;
plot(ratesest)

[x,fVal]=projectedGradientArmijoBoundary2(nbObservations,indexBeta,beta,C,D,T,R,regressor,nbLandTypes,sizex,iterMax,sigma,betaBar,epsilon,Groups,durations,sizex,whichgroup);
nbObs=nbWeeks*7;
proportion=0.2;
[weight,x]=crossValidation(sampleCalls,indexBeta,x,C,D,T,R,regressor,nbLandTypes,proportion,betaBar,sizex,iterMax,sigma,epsilon,nbObs);

p=1;
g=1;
meanrateBeta=0;
rates=[];
ratesest=[];
for index=1:nbWeeks*nbYears*7
    for g=1:G
        for t=1:T
            for i=1:R
                for p=1:P
                    rate=0;
                    rateest=0;
                    for j=1:nbLandTypes
                        rate=rate+TheoreticalBetas(t,g,p,j)*regressor(j,i);
                        rateest=rateest+x(indexBeta(t,g,p,j))*regressor(j,i);
                    end
                    %rate
                    %rateest
                    if (index<=nbObservationsG(g))
                        if (isHolidays{1,index}.b)
                            rate=rate+deltaTeoricos(t,isHolidays{1,index}.k,i,p);
                            rateest=rateest+x(indexDelta(t,isHolidays{1,index}.k,i,p));
                            rates=[rates,deltaTeoricos(t,isHolidays{1,index}.k,i,p)];
                            ratesest=[ratesest,x(indexDelta(t,isHolidays{1,index}.k,i,p))];
                        end
                        
                    end
                    meanrateBeta=meanrateBeta+100*(abs(rate-rateest)/rate);     
                end
            end
        end
    end
end
mean(abs((rates-ratesest))./rates)

meanrateBeta=meanrateBeta/(nbWeeks*nbYears*7*R*T);
g=1;
p=1;
j=2;
[TheoreticalBetas(1,g,p,j),TheoreticalBetas(2,g,p,j),TheoreticalBetas(3,g,p,j),TheoreticalBetas(4,g,p,j)]
j=2;
[TheoreticalBetas(1,g,p,j),TheoreticalBetas(2,g,p,j),TheoreticalBetas(3,g,p,j),TheoreticalBetas(4,g,p,j)]

j=1;
[x(indexBeta(1,g,p,j)),x(indexBeta(2,g,p,j)),x(indexBeta(3,g,p,j)),x(indexBeta(4,g,p,j))]
j=2;
[x(indexBeta(1,g,p,j)),x(indexBeta(2,g,p,j)),x(indexBeta(3,g,p,j)),x(indexBeta(4,g,p,j))]

p=1;
g=1;
j=1;
x=[[0:0.1:1],[1:0.1:2],[2:0.1:3],[3:0.1:4]];
y=[TheoreticalBetas(1,g,p,j)*ones(1,11),TheoreticalBetas(2,g,p,j)*ones(1,11),TheoreticalBetas(3,g,p,j)*ones(1,11),TheoreticalBetas(4,g,p,j)*ones(1,11)];
plot(x,y,'k-');

y01=[x001(indexBeta(1,g,p,j))*ones(1,11),x001(indexBeta(2,g,p,j))*ones(1,11),x001(indexBeta(3,g,p,j))*ones(1,11),x001(indexBeta(4,g,p,j))*ones(1,11)];
hold on;
plot(x,y01,'b-.');
y1000=[x1000(indexBeta(1,g,p,j))*ones(1,11),x1000(indexBeta(2,g,p,j))*ones(1,11),x1000(indexBeta(3,g,p,j))*ones(1,11),x1000(indexBeta(4,g,p,j))*ones(1,11)];
hold on;
plot(x,y1000,'r--');

y108=[x108(indexBeta(1,g,p,j))*ones(1,11),x108(indexBeta(2,g,p,j))*ones(1,11),x108(indexBeta(3,g,p,j))*ones(1,11),x108(indexBeta(4,g,p,j))*ones(1,11)];
hold on;
plot(x,y108,'r--');

ylim([0.05 0.25])

xTick = 1:4;
set(gca,'xtick',xTick);
yTick = get(gca,'ytick');
set(gca,'xticklabel',[]);
%xTickLabel = {{'first';'label'},'second',{'third';'label'},'long fourth tick label ','fifth'};
xTickLabel = {'1','2','3','4'};
for k = 1:length(xTick)
    if (k==1)
        xpos=xTick(1)/2;
    else
        xpos=xTick(k-1)+0.5*(xTick(k)-xTick(k-1));
    end
    text(xpos,yTick(1)-0.05*(yTick(end)-yTick(1)),xTickLabel{k},'HorizontalAlignment','center')
end


x=[[0:0.1:1],[1:0.1:2],[2:0.1:3],[3:0.1:4]];
y108=[x108(indexBeta(1,g,p,j))*ones(1,11),x108(indexBeta(2,g,p,j))*ones(1,11),x108(indexBeta(3,g,p,j))*ones(1,11),x108(indexBeta(4,g,p,j))*ones(1,11)];
hold on;
plot(x,y108,'k:');

p=1;
g=1;
i=1;
x=[];
y=[];
for k=1:10
    x=[x,k];
    y=[y,deltaTeoricos(1,k,i,p)];
end
plot(x,y,'k-');

y01=[];
for k=1:10
    y01=[y01,x001(indexDelta(1,k,i,p))];
end
hold on;
plot(x,y01,'b-.');

y1000=[];
for k=1:10
    y1000=[y1000,x1000(indexDelta(1,k,i,p))];
end
hold on;
plot(x,y1000,'r--');

y01=[x001(indexBeta(1,g,p,j))*ones(1,11),x001(indexBeta(2,g,p,j))*ones(1,11),x001(indexBeta(3,g,p,j))*ones(1,11),x001(indexBeta(4,g,p,j))*ones(1,11)];
hold on;
plot(x,y01,'b-.');
y1000=[x1000(indexBeta(1,g,p,j))*ones(1,11),x1000(indexBeta(2,g,p,j))*ones(1,11),x1000(indexBeta(3,g,p,j))*ones(1,11),x1000(indexBeta(4,g,p,j))*ones(1,11)];
hold on;
plot(x,y1000,'r--');

p=1;
g=1;
i=3;


differenceL=[]; 
for index=1:nbWeeks*nbYears
    index
    for g=1:G
        %Day number is (index-1)*G+g isHolidays((index-1)*G+g)
        for t=1:T
            for i=1:R
                for p=1:P
                    rate=0;
                    rateest=0;
                    for j=1:nbLandTypes
                        rate=rate+TheoreticalBetas(t,g,p,j)*regressor(j,i);
                        rateest=rateest+x(indexBeta(t,g,p,j))*regressor(j,i);
                    end
                    rate
                    rateest
                    if (index<=nbObservationsG(g))
                        if (isHolidays{1,(index-1)*G+g}.b)
                            rate=rate+deltaTeoricos(t,isHolidays{1,(index-1)*G+g}.k,i,p);
                            rateest=rateest+x(indexDelta(t,isHolidays{1,(index-1)*G+g}.k,i,p));
                        end
                    end
                    differenceL=[differenceL;abs(rate-rateest)];     
                end
            end
        end
    end
end

%Second projected gradient
[x,fVal]=projectedGradientArmijoFeasible2(nbObservations,nbObservationsHolidays,nbCallsHolidays,nbCallsNoHolidays,neighbors,type,distance,indexBeta,indexDelta,x,T,G,R,P,alpha,regressor,isHolidays,nbLandTypes,nbHolidaysYear,sizex,iterMax,sigma,epsilon);

%Cross validation
proportion=0.2;
[alpha,x]=crossValidation(sampleCalls,neighbors,type,distance,indexBeta,indexDelta,x,T,G,R,P,regressor,isHolidays,nbLandTypes,nbHolidaysYear,proportion,betaBar,sizex,iterMax,sigma,epsilon);


 
