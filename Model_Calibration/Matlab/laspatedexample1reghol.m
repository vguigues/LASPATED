
function [weight,lambda,obj]=laspatedex1reg()

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
nbWeeks=20;
nbYears=1;
durations=ones(1,T);
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
    TheoreticalBetas(1,d,2,1)=0.1;
    TheoreticalBetas(1,d,4,1)=0.1;
    
    TheoreticalBetas(1,d,1,2)=12;
    TheoreticalBetas(1,d,2,2)=36;
    TheoreticalBetas(1,d,3,2)=12;
    TheoreticalBetas(1,d,4,2)=36;
    
    TheoreticalBetas(1,d,1,3)=6;
    TheoreticalBetas(1,d,2,3)=12;
    TheoreticalBetas(1,d,3,3)=6;
    TheoreticalBetas(1,d,4,3)=12;
end

bs=rand(1,20);
rs=rand(1,20);

for i=1:(Ny/2)
    for j=1:(Nx/2)
        type(compteur)=1;
        indx=mod(compteur,10);
        if (indx==0)
            indx=10;
            indy=compteur/10;
        else
            indy=1+floor(compteur/10);
        end
        regressor(1,compteur)=computex1(indx,indy,1,bs,rs);
        regressor(2,compteur)=0.5;
        regressor(3,compteur)=0.25;
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        type(compteur)=2;
        indx=mod(compteur,10);
        if (indx==0)
            indx=10;
            indy=compteur/10;
        else
            indy=1+floor(compteur/10);
        end
        regressor(1,compteur)=computex1(indx,indy,2,bs,rs);
        regressor(2,compteur)=0.25;
        regressor(3,compteur)=0.5;
        compteur=compteur+1;
    end
end

for i=(Ny/2)+1:Ny
    for j=1:(Nx/2)
        type(compteur)=2;
        indx=mod(compteur,10);
        if (indx==0)
            indx=10;
            indy=compteur/10;
        else
            indy=1+floor(compteur/10);
        end
        regressor(1,compteur)=computex1(indx,indy,2,bs,rs);
        regressor(2,compteur)=0.25;
        regressor(3,compteur)=0.5;
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        type(compteur)=1;
        indx=mod(compteur,10);
        if (indx==0)
            indx=10;
            indy=compteur/10;
        else
            indy=1+floor(compteur/10);
        end
        regressor(1,compteur)=computex1(indx,indy,1,bs,rs);
        regressor(2,compteur)=0.5;
        regressor(3,compteur)=0.25;
        compteur=compteur+1;
    end
end

sample=cell(C,D,T,R);
for c=1:C
    for d=1:D
        for t=1:T
            for i=1:R
                sample{c,d,t,i}=[];     
            end
        end
    end
end
nbObservations=zeros(C,D,T,R);
nbArrivals=zeros(C,D,T,R);
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
                    rate=rate+TheoreticalBetas(c,day,t,j)*regressor(j,i);
                end
                thisnbcall=poissrnd(rate,1,1);
                sample{c,day,t,i}=[sample{c,day,t,i},thisnbcall];
                nbObservations(c,day,t,i)=nbObservations(c,day,t,i)+1;
                nbArrivals(c,day,t,i)=nbArrivals(c,day,t,i)+thisnbcall;
            end
        end
    end
end
epsilon=10^(-6);
beta=2*epsilon*ones(sizex,1);
sigma=0.5;
iterMax=30;
Groups=[];
whichgroup=[];
nbGroups=length(Groups);
params.indexBeta=indexBeta;
params.D=D;
params.regressor=regressor;
params.nbRegressors=1+nbLandTypes;
params.sizex=sizex;
params.nbObs=nbWeeks*7;
lbounds=zeros(sizex,1);
model='reg';
[lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,sigma,iterMax,epsilon,beta,params);


obsdays=zeros(D,1);
obsbefore=zeros(nbWeeks*7);
 for index=1:nbWeeks*7
     day=mod(index,7);
     if (day==0)
         day=7;
     end
     if (isHolidays{1,index}.b)
         day=7+isHolidays{1,index}.k;
     end
     obsdays(day)=obsdays(day)+1;
     obsbefore(index)=obsdays(day);  
 end

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
 
differenceL2=[];
for c=1:C
    for d=1:D
        for t=1:T
            for i=1:R
                ratee=0;
                rate=0;
                for j=1:(nbLandTypes+1)
                    ratee=ratee+lambda(indexBeta(c,d,t,j))*regressor(j,i);
                    rate=rate+TheoreticalBetas(c,d,t,j)*regressor(j,i);
                end
                differenceL2=[differenceL2;abs(rate-ratee)/rate];
            end
        end
    end
end
errcov=mean(differenceL2);

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
        %distance(i,j)=norm([absc(Xi)-absc(Xj);ord(Yi)-ord(Yj)]);
        %Setting the distance to one is equivalent to not considering
        %distances in penalizations. Rather we use the region type filled
        %above.
        distance(i,j)=1;
    end
end

T=60;
epsilon=0.001;
lambda0=epsilon*ones(C,R,T);
sigma=0.5;
betaBar=1;
iterMax=30;
model='noreg';
params.neighbors=neighbors;
params.type=type;
params.distance=distance;
params.alpha=1;

 
Groups{1,1}=[];
Groups{1,2}=[];
for t=1:14
    Groups{1,1}=[Groups{1,1},2*t-1];
    Groups{1,2}=[Groups{1,2},2*t];
end

Groups{1,3}=[];
Groups{1,4}=[];
for t=15:30
    Groups{1,3}=[Groups{1,3},2*t-1];
    Groups{1,4}=[Groups{1,4},2*t];
end

weights=[0,0.01,0.05,0.1,1,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400];
alphas=[0,0.01,0.05,0.1,1,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400];

 for i=1:100
     compteur=1;
     for d=1:D
         for t=1:4
             nbArrivals2(1,i,compteur)=nbArrivals(1,d,t,i);
             compteur=compteur+1;
         end
     end
 end
 nbObervationsTotal=nbWeeks*nbYears;
 nbArrivals=nbArrivals2;
 nbObservations=nbObervationsTotal*ones(1,100,15*4);
 
 nbGroups=4;
 whichgroup=zeros(T,1);
 for t=1:28
     whichgroup(t)=mod(t-1,2)+1;
 end
 t=29;
 for d=8:15
     whichgroup(t)=3;
     whichgroup(t+1)=4;
     whichgroup(t+2)=3;
     whichgroup(t+3)=4;
     t=t+4;
 end
 
durations=ones(1,60);
err=[];
for it=1:length(weights)
    it
    params.alpha=0;
    params.weight=weights(it)*ones(1,nbGroups);
    [lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,sigma,iterMax,epsilon,lambda0,params);
    differenceL2=[];
    for t=1:T
        for i=1:R
            for p=1:C
                delta=mod(t-1,4)+1;
                d=ceil(t/4);
                rate=0;
                for j=1:(nbLandTypes+1)
                    rate=rate+TheoreticalBetas(c,d,delta,j)*regressor(j,i);
                end
                differenceL2=[differenceL2;abs(rate-lambda(p,i,t))/rate];
            end
        end
    end
    err=[err,mean(differenceL2)];
end
errnosp=err;

err=[];
for it=1:length(weights)
    it
    params.alpha=alphas(it);
    params.weight=weights(it)*ones(1,nbGroups);
    [lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,sigma,iterMax,epsilon,lambda0,params);
    differenceL2=[];
    for t=1:T
        for i=1:R
            for p=1:C
                delta=mod(t-1,4)+1;
                d=ceil(t/4);
                rate=0;
                for j=1:(nbLandTypes+1)
                    rate=rate+TheoreticalBetas(c,d,delta,j)*regressor(j,i);
                end
                differenceL2=[differenceL2;abs(rate-lambda(p,i,t))/rate];
            end
        end
    end
    err=[err,mean(differenceL2)];
end
errscol=err

for i=1:(Ny/2)
    for j=1:(Nx/2)
        type(compteur)=1;
        indx=mod(compteur,10);
        if (indx==0)
            indx=10;
            indy=compteur/10;
        else
            indy=1+floor(compteur/10);
        end
        regressor(1,compteur)=computex1(indx,indy,1,bs,rs);
        regressor(2,compteur)=0.5;
        regressor(3,compteur)=0.25;
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        type(compteur)=1;
        indx=mod(compteur,10);
        if (indx==0)
            indx=10;
            indy=compteur/10;
        else
            indy=1+floor(compteur/10);
        end
        regressor(1,compteur)=computex1(indx,indy,2,bs,rs);
        regressor(2,compteur)=0.25;
        regressor(3,compteur)=0.5;
        compteur=compteur+1;
    end
end

for i=(Ny/2)+1:Ny
    for j=1:(Nx/2)
        type(compteur)=1;
        indx=mod(compteur,10);
        if (indx==0)
            indx=10;
            indy=compteur/10;
        else
            indy=1+floor(compteur/10);
        end
        regressor(1,compteur)=computex1(indx,indy,2,bs,rs);
        regressor(2,compteur)=0.25;
        regressor(3,compteur)=0.5;
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        type(compteur)=1;
        indx=mod(compteur,10);
        if (indx==0)
            indx=10;
            indy=compteur/10;
        else
            indy=1+floor(compteur/10);
        end
        regressor(1,compteur)=computex1(indx,indy,1,bs,rs);
        regressor(2,compteur)=0.5;
        regressor(3,compteur)=0.25;
        compteur=compteur+1;
    end
end

err=[];
for it=1:length(weights)
    it
    params.alpha=alphas(it);
    params.weight=weights(it)*ones(1,nbGroups);
    [lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,sigma,iterMax,epsilon,lambda0,params);
    differenceL2=[];
    for t=1:T
        for i=1:R
            for p=1:C
                delta=mod(t-1,4)+1;
                d=ceil(t/4);
                rate=0;
                for j=1:(nbLandTypes+1)
                    rate=rate+TheoreticalBetas(c,d,delta,j)*regressor(j,i);
                end
                differenceL2=[differenceL2;abs(rate-lambda(p,i,t))/rate];
            end
        end
    end
    err=[err,mean(differenceL2)];
end
errsnocol=err

[min(errnosp) min(errscol) min(errsnocol) errnosp(1)]
