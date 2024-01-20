
function [err,err1,err2,err3]=laspatedexample1noreg(nbperiods,nbWeeks)

xmax=10;
ymax=10;
Nx=10;
Ny=10;
R=Nx*Ny;
T=nbperiods;
Groups=cell(1,4);
nbGroups=4;
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

C=1;
nbYears=1;
nbObervationsTotal=nbWeeks*nbYears;
maxObs=nbObervationsTotal;
durations=ones(1,T);
nbVars=1;
sample=zeros(T,R,C,nbObervationsTotal);
compteur=1;

for k=1:length(Groups{1,1})
            Theoreticallambdar(Groups{1,1}(k),1)=0.5;
end
for k=1:length(Groups{1,2})
            Theoreticallambdar(Groups{1,2}(k),1)=0.1;
end
for k=1:length(Groups{1,3})
            Theoreticallambdar(Groups{1,3}(k),1)=0.5;
end
for k=1:length(Groups{1,4})
            Theoreticallambdar(Groups{1,4}(k),1)=0.1;
end

for k=1:length(Groups{1,1})
            Theoreticallambdab(Groups{1,1}(k),1)=0.1;
end
for k=1:length(Groups{1,2})
            Theoreticallambdab(Groups{1,2}(k),1)=0.5;
end
for k=1:length(Groups{1,3})
            Theoreticallambdab(Groups{1,3}(k),1)=0.1;
end
for k=1:length(Groups{1,4})
            Theoreticallambdab(Groups{1,4}(k),1)=0.5;
end


for i=1:(Ny/2)
    for j=1:(Nx/2)
        type(compteur)=1;
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        type(compteur)=2;
        compteur=compteur+1;
    end
end

for i=(Ny/2)+1:Ny
    for j=1:(Nx/2)
        type(compteur)=2;
        compteur=compteur+1;
    end
    for j=1:(Nx/2)
        type(compteur)=1;
        compteur=compteur+1;
    end
end

rmap=[1,2,3,4,5,11,12,13,14,15,21,22,23,24,25,31,32,33,34,35,41,42,43,44,45,56,57,58,59,60,66,67,68,69,70,76,77,78,79,80,86,87,88,89,90,96,97,98,99,100];
bmap=[6,7,8,9,10,16,17,18,19,20,26,27,28,29,30,36,37,38,39,40,46,47,48,49,50,51,52,53,54,55,61,62,63,64,65,71,72,73,74,75,81,82,83,84,85,91,92,93,94,95];

nbArrivals=zeros(C,R,T);
nbObservations=nbObervationsTotal*ones(C,R,T);
% for t=1:T
%         for p=1:C
%             thissample=poissrnd(Theoreticallambdar(t,p)*durations(t)*50,1,nbObervationsTotal);
%             for i=1:nbObervationsTotal
%                 samplesites=1+floor(50*rand(1,thissample(i)));
%                 for j=1:thissample(i)
%                     if (samplesites(j)==51)
%                         samplesites(j)=50;
%                     end
%                     sample(t,rmap(samplesites(j)),p,i)=sample(t,rmap(samplesites(j)),p,i)+1;
%                     nbArrivals(p,rmap(samplesites(j)),t)=nbArrivals(p,rmap(samplesites(j)),t)+1;
%                 end
%             end
%             thissample=poissrnd(Theoreticallambdab(t,p)*durations(t)*50,1,nbObervationsTotal);
%             for i=1:nbObervationsTotal
%                 samplesites=1+floor(50*rand(1,thissample(i)));
%                 for j=1:thissample(i)
%                     if (samplesites(j)==51)
%                         samplesites(j)=50;
%                     end
%                     sample(t,bmap(samplesites(j)),p,i)=sample(t,bmap(samplesites(j)),p,i)+1;
%                     nbArrivals(p,bmap(samplesites(j)),t)=nbArrivals(p,bmap(samplesites(j)),t)+1;
%                 end
%             end
% 
%         end
% end

nbArrivals=zeros(C,R,T);
nbObservations=nbObervationsTotal*ones(C,R,T);
for t=1:T
        for p=1:C
            for r=1:100
                if (type(r)==1)
                    thissample=poissrnd(Theoreticallambdar(t,p)*durations(t),1,nbObervationsTotal);
                else
                    thissample=poissrnd(Theoreticallambdab(t,p)*durations(t),1,nbObervationsTotal);
                end
                for j=1:nbObervationsTotal
                    nbArrivals(p,r,t)=nbArrivals(p,r,t)+thissample(j);
                    sample(t,r,p,j)=thissample(j);
                end
            end
        end
end

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

epsilon=0.001;
lambda0=epsilon*ones(C,R,T);
sigma=0.5;
betaBar=1;
iterMax=30;
model='noreg';
params.neighbors=neighbors;
params.type=type;
params.distance=distance;
%weights=[0,0.01,0.05,0.1,1,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400];
%alphas=[0,0.01,0.05,0.1,1,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400];


weights=[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,30,40,50,60,70,80,90,100];
alphas=[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,30,40,50,60,70,80,90,100];

%weights=[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4];
%alphas=[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4];

%weights=[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,1,2,3,4,5,6,7,8,9,10];
%alphas=[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,1,2,3,4,5,6,7,8,9,10];

weights=[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.4]/10;
alphas=[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.4]/10;

params.delta=1;
params.uppx=1000*ones(C,R,T);
params.lowx=epsilon*ones(C,R,T);

err=[];
%for it=1:10
for it=1:length(weights)
    it
    params.alpha=alphas(it);
    params.weight=weights(it)*ones(1,nbGroups);

    [lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,sigma,iterMax,epsilon,lambda0,params);
    differenceL2=[];
    for t=1:T
            for i=1:R
                for p=1:C
                    if (type(i)==1)
                        differenceL2=[differenceL2;abs(Theoreticallambdar(t,p)-lambda(p,i,t))/Theoreticallambdar(t,p)];
                    else
                        differenceL2=[differenceL2;abs(Theoreticallambdab(t,p)-lambda(p,i,t))/Theoreticallambdab(t,p)];
                    end
                end
            end
    end
    err=[err,mean(differenceL2)]; 
end

err1=[];
for it=1:length(weights)
%for it=1:10
    it
    params.alpha=0;
    params.weight=weights(it)*ones(1,nbGroups);

    [lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,sigma,iterMax,epsilon,lambda0,params);
    differenceL2=[];
    for t=1:T
            for i=1:R
                for p=1:C
                    if (type(i)==1)
                        differenceL2=[differenceL2;abs(Theoreticallambdar(t,p)-lambda(p,i,t))/Theoreticallambdar(t,p)];
                    else
                        differenceL2=[differenceL2;abs(Theoreticallambdab(t,p)-lambda(p,i,t))/Theoreticallambdab(t,p)];
                    end
                end
            end
    end
    err1=[err1,mean(differenceL2)]; 
end


%%%%%%%%%%%%%%%%%%%%%%
%2 time groups 
%%%%%%%%%%%%%%%%%%%%%%

Groups=cell(1,2);
nbGroups=2;
whichgroup=zeros(T,1);
for t=1:T
    whichgroup(t)=mod(t-1,2)+1;
end
Groups{1,1}=[1];
Groups{1,2}=[2];
for i=1:2
    for j=1:13
        Groups{1,i}=[Groups{1,i},Groups{1,i}(j)+2];
    end
end

err2=[];
for it=1:length(weights)
%for it=1:10
    it
    params.alpha=alphas(it);
    params.weight=weights(it)*ones(1,nbGroups);
    [lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,sigma,iterMax,epsilon,lambda0,params);
    differenceL2=[];
    for t=1:T
            for i=1:R
                for p=1:C
                    if (type(i)==1)
                        differenceL2=[differenceL2;abs(Theoreticallambdar(t,p)-lambda(p,i,t))/Theoreticallambdar(t,p)];
                    else
                        differenceL2=[differenceL2;abs(Theoreticallambdab(t,p)-lambda(p,i,t))/Theoreticallambdab(t,p)];
                    end
                end
            end
    end
    err2=[err2,mean(differenceL2)]; 
end

err3=[];
for it=1:length(weights)
%for it=1:10
    it
    params.alpha=0;
    params.weight=weights(it)*ones(1,nbGroups);
    [lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,sigma,iterMax,epsilon,lambda0,params);
    differenceL2=[];
    for t=1:T
            for i=1:R
                for p=1:C
                    if (type(i)==1)
                        differenceL2=[differenceL2;abs(Theoreticallambdar(t,p)-lambda(p,i,t))/Theoreticallambdar(t,p)];
                    else
                        differenceL2=[differenceL2;abs(Theoreticallambdab(t,p)-lambda(p,i,t))/Theoreticallambdab(t,p)];
                    end
                end
            end
    end
    err3=[err3,mean(differenceL2)]; 
end

differenceL2=[];
    for t=1:T
            for i=1:R
                for p=1:C
                    if (type(i)==1)
                        differenceL2=[differenceL2;abs(Theoreticallambdar(t,p)-estimated(t,i,p))/Theoreticallambdar(t,p)];
                    else
                        differenceL2=[differenceL2;abs(Theoreticallambdab(t,p)-estimated(t,i,p))/Theoreticallambdab(t,p)];
                    end
                end
            end
    end
mean(differenceL2) 

weights=[0.01,0.03,0.04,0.05,0.06,0.1];
alphas=[0.01,0.03,0.04,0.05,0.06,0.1];

weights=[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,1,2,3,4,5,6,7,8,9,10];
alphas=[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,1,2,3,4,5,6,7,8,9,10];

weights=[0.02,0.03,0.04,0.05];
alphas=[0.02,0.03,0.04,0.05];

alphas=[0.03,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.4];


%lw=length(weights);
%weights=weights(1:5);
proportion=0.2;
[cputime,weightcv,lambdacv]=crossValidation(model,sample,T,R,C,sigma,lambda0,iterMax,proportion,epsilon,durations,Groups,whichgroup,weights,params);

errcv=[];
params.alpha=weightcv;
params.weight=weightcv*ones(1,nbGroups);
[lambdacv,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,sigma,iterMax,epsilon,lambda0,params);
differenceL2=[];
for t=1:T
    for i=1:R
        for p=1:C
            if (type(i)==1)
                differenceL2=[differenceL2;abs(Theoreticallambdar(t,p)-lambdacv(p,i,t))/Theoreticallambdar(t,p)];
            else
                differenceL2=[differenceL2;abs(Theoreticallambdab(t,p)-lambdacv(p,i,t))/Theoreticallambdab(t,p)];
            end
        end
    end
end
errcv=mean(differenceL2);

weights=[0,0.01,0.05,0.1,1,5,10,20,30,40,50,60,70,80,90,100];
alphas=[0,0.01,0.05,0.1,1,5,10,20,30,40,50,60,70,80,90,100];

weights=[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,1,2,3,4,5,6,7,8,9,10];
alphas=[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,1,2,3,4,5,6,7,8,9,10];


weights=[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.4];
alphas=[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.4];


endw=length(weights);

plot(weights(1:endw),err(1:endw),'k-','LineWidth',1.5);
hold on;
plot(weights(1:endw),err1(1:endw),'-.r','LineWidth',1.5);
hold on;
plot(weights(1:endw),err2(1:endw),':b','LineWidth',1.5);
hold on;
plot(weights(1:endw),err3(1:endw),'--k','LineWidth',1.5);
hold on;
plot(weights(1:endw),errcv*ones(1,endw),'-r','LineWidth',1.5);

legend(['4  groups, neighbors   ';'4  groups, no neighbors';'2  groups, neighbors   ';'2  groups, no neighbors';'Cross validation       ']);

legend(['4  groups, neighbors   ';'4  groups, no neighbors';'2  groups, neighbors   ';'2  groups, no neighbors']);

[a,b]=min(err);
hold on;
plot(weights(b),a,'kd','LineWidth',1.5);
[a1,b1]=min(err1);
hold on;
plot(weights(b1),a1,'-.rd','LineWidth',1.5);
hold on;
[a2,b2]=min(err2);
plot(weights(b2),a2,':bd','LineWidth',1.5);
[a3,b3]=min(err3);
hold on;
plot(weights(b3),a3,'--kd','LineWidth',1.5);

legend(['4  groups, neighbors   ';'4  groups, no neighbors';'2  groups, neighbors   ';'2  groups, no neighbors']);
%legend(['4  groups, neighbors   ';'4  groups, no neighbors';'2  groups, neighbors   ';'2  groups, no neighbors';'Cross validation       ']);

for t=1:T
    for i=1:R
        for p=1:C
            estimated(t,i,p)=nbArrivals(p,i,t)/(nbObservations(p,i,t)*durations(t));
        end
    end
end

differenceL2=[];
for t=1:T
    for i=1:R
        for p=1:C
            if (type(i)==1)
                differenceL2=[differenceL2;abs(Theoreticallambdar(t,p)-estimated(t,i,p))/Theoreticallambdar(t,p)];
            else
                differenceL2=[differenceL2;abs(Theoreticallambdab(t,p)-estimated(t,i,p))/Theoreticallambdab(t,p)];
            end
        end
    end
end
mean(differenceL2)

aux=[];
for t=1:T
    for i=1:R
        if (type(i)==1)
            aux=[aux;estimated(t,i,1)-Theoreticallambdar(t,1)];
        elseif (type(i)==2)
            aux=[aux;estimated(t,i,1)-Theoreticallambdab(t,1)];
        end
    end
end
mean(aux);

figure();

Groups=cell(1,2);
nbGroups=2;
whichgroup=zeros(T,1);
for t=1:T
    whichgroup(t)=mod(t-1,2)+1;
end
Groups{1,1}=[1];
Groups{1,2}=[2];
for i=1:2
    for j=1:13
        Groups{1,i}=[Groups{1,i},Groups{1,i}(j)+2];
    end
end

params.alpha=alphas(b2);
params.weight=weights(b2)*ones(1,nbGroups);
[lambda,obj]=laspated(model,nbObservations,nbArrivals,T,R,C,durations,Groups,whichgroup,sigma,iterMax,epsilon,lambda0,params);

figure()
xT=[];
xr=[];
xemp=[];
yT=[];
yr=[];
yemp=[];
for t=1:T
    xT=[xT,[t-1:0.1:t-1+0.9]];
    xr=[xr,[t-1:0.1:t-1+0.9]];
    xemp=[xemp,[t-1:0.1:t-1+0.9]];
    yT=[yT,Theoreticallambdar(t,1)*ones(1,10)];
    yr=[yr,lambda(1,1,t)*ones(1,10)];
    yemp=[yemp,estimated(t,1,1)*ones(1,10)];      
end
plot(xT,yT,'-k');
hold on;
plot(xr,yr,'--b');
hold on;
plot(xemp,yemp,'-.r');
ylim([-0.1 1.1*max(yemp)])
xlabel('Time periods');
ylabel('Intensities');
legend(['True intensity       ';'Regularized estimator';'Empirical estimator  '])


figure()

xT=[];
xr=[];
xemp=[];
yT=[];
yr=[];
yemp=[];

for t=1:T
    xT=[xT,[t-1:0.1:t-1+0.9]];
    xr=[xr,[t-1:0.1:t-1+0.9]];
    xemp=[xemp,[t-1:0.1:t-1+0.9]];

    yT=[yT,Theoreticallambdab(t,1)*ones(1,10)];
    yr=[yr,lambda(1,6,t)*ones(1,10)];
    yemp=[yemp,estimated(t,6,1)*ones(1,10)];      
end
plot(xT,yT,'-k');
hold on;
plot(xr,yr,'--b');
hold on;
plot(xemp,yemp,'-.r');
ylim([-0.1 1.1*max(yemp)])
xlabel('Time periods');
ylabel('Intensities');
legend(['True intensity       ';'Regularized estimator';'Empirical estimator  '])

