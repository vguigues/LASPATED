             
%file:///C:/gurobi901/win64/docs/refman/matlab_api_overview.html#matlab:MATLAB
%File generator_regressor.cpp in C++
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
Nx=2;
Ny=2;
R=Nx*Ny;
P=1;
nbWeeks=2;
nbYears=1;

T=4*7;
Groups=cell(1,4);
whichgroup=zeros(T,1);
for t=1:T
    whichgroup(t)=mod(t-1,4)+1;
end
Groups{1,1}=[1];
Groups{1,2}=[2];
Groups{1,3}=[3];
Groups{1,4}=[4];
for i=1:4
    for j=1:6
        Groups{1,i}=[Groups{1,i},Groups{1,i}(j)+4];
    end
end
P=1;
nbObervationsTotal=nbWeeks*nbYears;
duration=0.5*ones(1,T);
nbHolidaysYear=10;
holidays=cell(1,nbHolidaysYear);
isHolidays=cell(1,nbWeeks*nbYears);
for j=1:nbWeeks*nbYears
    isHolidays{1,j}.b=0;
    isHolidays{1,j}.k=0;
end

daysH=[12;42;72;93;122;138;243;287;312;350];

for year=1:nbYears
    for k=1:length(daysH)
        isHolidays{1,((year-1)*nbWeeks)*7+daysH(k)}.b=1;
        isHolidays{1,((year-1)*nbWeeks)*7+daysH(k)}.k=k;
    end
end


% holidays{1,1}.day=[];
% holidays{1,1}.week=[];
% for year=1:nbYears
%     day=mod(year,7);
%     if (day==0)
%         day=7;
%     end 
%     day=14;
%     holidays{1,1}.day=[holidays{1,1}.day,day];
%     holidays{1,1}.week=[holidays{1,1}.week,1+(year-1)*nbWeeks];
%     isHolidays{1,((year-1)*nbWeeks)*G+day}.b=1;
%     isHolidays{1,((year-1)*nbWeeks)*G+day}.k=1;
%     
% end
% 
% holidays{1,2}.day=[];
% holidays{1,2}.week=[];
% for year=1:nbYears
%     day=mod(3+year,7);
%     if (day==0)
%         day=7;
%     end 
%     holidays{1,2}.day=[holidays{1,2}.day,day];
%     holidays{1,2}.week=[holidays{1,2}.week,12+(year-1)*nbWeeks];
%     isHolidays{1,(11+(year-1)*nbWeeks)*G+day}.b=1;
%     isHolidays{1,(11+(year-1)*nbWeeks)*G+day}.k=2;
% end
% 
% holidays{1,3}.day=[];
% holidays{1,3}.week=[];
% for year=1:nbYears
%     day=mod(5+year,7);
%     if (day==0)
%         day=7;
%     end 
%     holidays{1,3}.day=[holidays{1,3}.day,day];
%     holidays{1,3}.week=[holidays{1,3}.week,24+(year-1)*nbWeeks];
%     isHolidays{1,(23+(year-1)*nbWeeks)*G+day}.b=1;
%     isHolidays{1,(23+(year-1)*nbWeeks)*G+day}.k=3;
% end
% 
% holidays{1,4}.day=[];
% holidays{1,4}.week=[];
% for year=1:nbYears
%     day=mod(5+year,7);
%     if (day==0)
%         day=7;
%     end 
%     holidays{1,4}.day=[holidays{1,4}.day,day];
%     holidays{1,4}.week=[holidays{1,4}.week,31+(year-1)*nbWeeks];
%     isHolidays{1,(30+(year-1)*nbWeeks)*G+day}.b=1;
%     isHolidays{1,(30+(year-1)*nbWeeks)*G+day}.k=4;
% end
% 
% holidays{1,5}.day=[];
% holidays{1,5}.week=[];
% for year=1:nbYears
%     day=mod(6+year,7);
%     if (day==0)
%         day=7;
%     end 
%     holidays{1,5}.day=[holidays{1,5}.day,day];
%     holidays{1,5}.week=[holidays{1,5}.week,41+(year-1)*nbWeeks];
%     isHolidays{1,(40+(year-1)*nbWeeks)*G+day}.b=1;
%     isHolidays{1,(40+(year-1)*nbWeeks)*G+day}.k=5;
% end
% 
% holidays{1,6}.day=[];
% holidays{1,6}.week=[];
% for year=1:nbYears
%     day=mod(2+year,7);
%     if (day==0)
%         day=7;
%     end 
%     holidays{1,6}.day=[holidays{1,6}.day,day];
%     holidays{1,6}.week=[holidays{1,6}.week,51+(year-1)*nbWeeks];
%     isHolidays{1,(50+(year-1)*nbWeeks)*G+day}.b=1;
%     isHolidays{1,(50+(year-1)*nbWeeks)*G+day}.k=6;
% end

nbLandTypes=2;
betaTeoricos=zeros(T,G,P,nbLandTypes);
for g=1:G
    for p=1:P
%         betaTeoricos(1,g,p,1)=(1/25)*30;
%         betaTeoricos(2,g,p,1)=(1/25)*2;
%         betaTeoricos(3,g,p,1)=(1/25)*30;
%         betaTeoricos(4,g,p,1)=(1/25)*2;
% 
%         betaTeoricos(1,g,p,2)=(1/25)*4;
%         betaTeoricos(2,g,p,2)=(1/25)*2;
%         betaTeoricos(3,g,p,2)=(1/25)*4;
%         betaTeoricos(4,g,p,2)=(1/25)*2;
%         
        
        
    end
end

deltaTeoricos=zeros(T,nbHolidaysYear,R,P);
for t=1:T
    for k=1:nbHolidaysYear
        for i=1:R
            for p=1:P
                if ((i==1)||(i==4))
                if ((k==1)||(k==2))
                    deltaTeoricos(t,k,i,p)=1;
                else
                    deltaTeoricos(t,k,i,p)=10;
                end
                else
                if ((k==1)||(k==2))
                    deltaTeoricos(t,k,i,p)=1;
                else
                    deltaTeoricos(t,k,i,p)=20;
                end
                    
                end
               %deltaTeoricos(t,k,i,p)=rand*duration(t);
            end
        end
    end
end

sampleCalls=zeros(T,G,R,P,7*nbYears*nbWeeks);
nbObservations=zeros(T,G,R,P);
nbObservationsHolidays=zeros(T,nbHolidaysYear,R,P);
nbCallsHolidays=zeros(T,nbHolidaysYear,G,R,P);
nbCallsNoHolidays=zeros(T,G,R,P);
%regressor=10+90*rand(nbLandTypes,R);
regressor=zeros(nbLandTypes,R);
regressor(:,1)=[23;2];
regressor(:,4)=[23;2];
regressor(:,3)=[2;23];
regressor(:,2)=[2;23];
indexWeeks=[1:nbYears*nbWeeks];

g=1;
p=1;

for index=1:nbWeeks*nbYears*7
        %Day number is (index-1)*G+g isHolidays((index-1)*G+g)
        for t=1:T
            for i=1:R
                    rate=0;
                    for j=1:nbLandTypes
                        rate=rate+betaTeoricos(t,g,p,j)*regressor(j,i);
                    end
                    
                    if (index<=nbObservationsG(g))
                        nbObservations(t,g,i,p)=nbObservations(t,g,i,p)+1;
                        if (isHolidays{1,index}.b)
                            sampleCalls(t,g,i,p,index)=poissrnd(rate+deltaTeoricos(t,isHolidays{1,index}.k,i,p),1,1);
                            nbObservationsHolidays(t,isHolidays{1,index}.k,i,p)=nbObservationsHolidays(t,isHolidays{1,index}.k,i,p)+1;
                            nbCallsHolidays(t,isHolidays{1,index}.k,g,i,p)=nbCallsHolidays(t,isHolidays{1,index}.k,g,i,p)+sampleCalls(t,g,i,p,index);
                        else
                            sampleCalls(t,g,i,p,index)=poissrnd(rate,1,1);
                            nbCallsNoHolidays(t,g,i,p)=nbCallsNoHolidays(t,g,i,p)+sampleCalls(t,g,i,p,index);
                        end
                    end
            end
        end
  end

indexBeta=zeros(T,G,P,nbLandTypes);
indexDelta=zeros(T,nbHolidaysYear,R,P);
nbVars=1;
for t=1:T
    for g=1:G
            for p=1:P
                for j=1:nbLandTypes
                    indexBeta(t,g,p,j)=nbVars;
                    nbVars=nbVars+1;
                end
            end         
    end
end

for t=1:T
    for k=1:nbHolidaysYear
        for i=1:R
            for p=1:P
                indexDelta(t,k,i,p)=nbVars;
                nbVars=nbVars+1;
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
type(1)=1;
type(4)=1;
type(2)=2;
type(3)=2;

for i=1:R
    %type(i)=1+floor(4*rand);  
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

epsilon=0.01;
delta=epsilon*ones(T*nbHolidaysYear*R*P,1);
beta=(10^(-40))*ones(T*G*P*nbLandTypes,1);
x=[beta;delta];
%delta=epsilon*ones(T,nbHolidaysYear,R,P);
%beta=epsilon*ones(T,G,P,nbLandTypes);

sigma=0.5;
betaBar=1;
iterMax=300;
alpha=0;
sizex=T*nbHolidaysYear*R*P+T*G*P*nbLandTypes;

%First projected gradient
[x,fVal]=projectedGradientArmijoFeasible2(nbObservations,nbObservationsHolidays,nbCallsHolidays,nbCallsNoHolidays,neighbors,type,distance,indexBeta,indexDelta,x,T,G,R,P,alpha,regressor,isHolidays,nbLandTypes,nbHolidaysYear,sizex,iterMax,sigma,epsilon);
%[x,fVal]=projectedGradientArmijoBoundary2(nbObservations,nbObservationsHolidays,nbCallsHolidays,nbCallsNoHolidays,neighbors,type,distance,indexBeta,indexDelta,x,T,G,R,P,alpha,regressor,isHolidays,nbLandTypes,nbHolidaysYear,sizex,iterMax,sigma,betaBar,epsilon);
x001=x;

alpha=1000;
epsilon=0.01;
delta=epsilon*ones(T*nbHolidaysYear*R*P,1);
beta=epsilon*ones(T*G*P*nbLandTypes,1);
x=[beta;delta];
[x1000,fVal]=projectedGradientArmijoFeasible2(nbObservations,nbObservationsHolidays,nbCallsHolidays,nbCallsNoHolidays,neighbors,type,distance,indexBeta,indexDelta,x,T,G,R,P,alpha,regressor,isHolidays,nbLandTypes,nbHolidaysYear,sizex,iterMax,sigma,epsilon);

alpha=10^(10);
epsilon=0.01;
delta=epsilon*ones(T*nbHolidaysYear*R*P,1);
beta=epsilon*ones(T*G*P*nbLandTypes,1);
x=[beta;delta];
[x108,fVal]=projectedGradientArmijoFeasible2(nbObservations,nbObservationsHolidays,nbCallsHolidays,nbCallsNoHolidays,neighbors,type,distance,indexBeta,indexDelta,x,T,G,R,P,alpha,regressor,isHolidays,nbLandTypes,nbHolidaysYear,sizex,iterMax,sigma,epsilon);


%[x1000,fVal]=projectedGradientArmijoBoundary2(nbObservations,nbObservationsHolidays,nbCallsHolidays,nbCallsNoHolidays,neighbors,type,distance,indexBeta,indexDelta,x,T,G,R,P,alpha,regressor,isHolidays,nbLandTypes,nbHolidaysYear,sizex,iterMax,sigma,betaBar,epsilon);

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
                        rate=rate+betaTeoricos(t,g,p,j)*regressor(j,i);
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
[betaTeoricos(1,g,p,j),betaTeoricos(2,g,p,j),betaTeoricos(3,g,p,j),betaTeoricos(4,g,p,j)]
j=2;
[betaTeoricos(1,g,p,j),betaTeoricos(2,g,p,j),betaTeoricos(3,g,p,j),betaTeoricos(4,g,p,j)]

j=1;
[x(indexBeta(1,g,p,j)),x(indexBeta(2,g,p,j)),x(indexBeta(3,g,p,j)),x(indexBeta(4,g,p,j))]
j=2;
[x(indexBeta(1,g,p,j)),x(indexBeta(2,g,p,j)),x(indexBeta(3,g,p,j)),x(indexBeta(4,g,p,j))]

p=1;
g=1;
j=1;
x=[[0:0.1:1],[1:0.1:2],[2:0.1:3],[3:0.1:4]];
y=[betaTeoricos(1,g,p,j)*ones(1,11),betaTeoricos(2,g,p,j)*ones(1,11),betaTeoricos(3,g,p,j)*ones(1,11),betaTeoricos(4,g,p,j)*ones(1,11)];
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
                        rate=rate+betaTeoricos(t,g,p,j)*regressor(j,i);
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


 
