cd 'C:\gurobi901\win64\matlab\'
gurobi_setup
addpath 'C:\Program Files\Mosek\9.0'
addpath 'C:\Program Files\Mosek\9.0\toolbox\r2015a'
addpath 'C:\Users\vince\Dropbox\Softwares\Heuristic_Ambulance_Dispatch\Call_Model_Calibration\Model_Without_Regressors';
addpath 'C:\Users\vince\Dropbox\Softwares\Heuristic_Ambulance_Dispatch\Call_Model_Calibration';

%Rect 10x10 alpha=1000
%Hex 7 alpha=10000
%Hex 8 alpha=10000

%Test with real data
discretization_dir='C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\New_Data\Hex7_km';
%[T,G,R,P,nbLandTypes,nbObservationsG,sample_calls,nbCalls,nbObservations,estimated,type,regressors,neighbors,distance]=read_calls(discretization_dir);
[T,G,R,P,nbLandTypes,nbObservationsG,nbCalls,estimated,type,regressors,neighbors,distance]=read_calls(discretization_dir);


sample_callsA=sample_calls;
nbCallsA=nbCalls;
nbObservationsA=nbObservations;
C=P;
nbObservations=zeros(7*T,R,P);
nbCalls=zeros(7*T,R,P);
sampleCalls=zeros(7*T,R,P,104);
%sampleCalls=zeros(7*T,R,P,min(nbObservationsG));
nbObervationsTotal=min(nbObservationsG)*7;
for r=1:R
    for p=1:P
        for g=1:G
            init=(g-1)*T;
            for t=1:T
                nbObservations(init+t,r,p)=nbObservationsA(t,g,r,p);
                %size(nbObservations)
                %init+t
                nbCalls(init+t,r,p)=nbCallsA(t,g,r,p);
                for j=1:104
                    sampleCalls(init+t,r,p,j)=sample_callsA(t,g,r,p,j);
                    %sampleCalls(init+t,r,p,j)=sample_callsA{p,g,t,r}(j);
                end
            end
        end
    end
end

T=7*T;
Groups=cell(1,T);
for t=1:T
    Groups{1,t}=[t];
    whichgroup(t)=t;
end
alpha=1000;
epsilon=0.001;
x=ones(T,R,P);
sigma=0.5;
betaBar=1;
iterMax=20;
weight=0.1;
durations=0.5*ones(1,T);

x=ones(T,R,P);
nbObservationsTotal=size(sampleCalls,4);
proportion=0.2;
%[cputime,alpha,weight,lambda]=crossValidation(nbObservationsTotal,sampleCalls,neighbors,type,distance,T,R,P,sigma,betaBar,x,iterMax,proportion,epsilon,durations,Groups,whichgroup);


%Call model calibration without regressors

%[xnr,fVal]=projectedGradientArmijoBoundary(nbObservations,nbCalls,neighbors,type,distance,T,R,P,sigma,betaBar,x,iterMax,alpha,epsilon,durations,Groups,whichgroup,weight)
alpha=1000;
[x,fVal]= projectedGradientArmijoFeasible(nbObservations,nbCalls,neighbors,type,distance,T,R,P,sigma,x,iterMax,alpha,epsilon,durations,Groups,whichgroup,1000)

weeklambda=zeros(T,1);
for t=1:T
        for r=1:R
            for p=1:P
                weeklambda(t)=weeklambda(t)+x(t,r,p);
                %weeklambda(t)=weeklambda(t)+estimated(t,g,r,p);
            end
        end
end
plot(weeklambda);

estTarget=zeros(48*7,1);
for t=1:48
    for g=1:7
        for r=1:R
            for p=1:3
                estTarget((g-1)*48+t,1)=estTarget((g-1)*48+t,1)+estimated(t,g,r,p);
            end
        end
    end
end

estTarget=zeros(48*7,1);
for t=1:48
    for g=1:7
        for r=1:R
            for p=1:3
                estTarget((g-1)*48+t,1)=estTarget((g-1)*48+t,1)+estimated(t,g,r,p);
            end
        end
    end
end


LHD((g-1)*48+t,(p-1)*160+r)

estimated=zeros(T,G,R,P);
for t=1:T
    for g=1:7
        for r=1:R
            for p=1:3
                estimated(t,g,r,p)=LH7((g-1)*48+t,(p-1)*297+r);
            end
        end
    end
end


heatmap=zeros(1,R);
for t=1:48
    for g=1:7
        for r=1:R
            for p=1:3
                %rateest=0;
                %for j=1:(1+nbLandTypes)
                %    rateest=rateest+LR10s((g-1)*48+t,(p-1)*4916+r);
                %end
                %heatmap(p,r)=heatmap(p,r)+LR10((g-1)*48+t,(p-1)*76+r);  
                heatmap(1,r)=heatmap(1,r)+estimated(t,g,r,p);  
            end
        end
    end
end
dlmwrite('C:\Users\vince\Desktop\Figures_Calls\heatEmpR10.txt',heatmap);

%x(indexBeta(c,d,t,j))
%x(t,i,p)

errR10=inf;
lambda100R10=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda100_Rect10.txt');
thiserr=0;
for t=1:48
    for g=1:7
        thisest=0;
        for r=1:76
            for p=1:3
                thisest=thisest+lambda100R10((g-1)*48+t,(p-1)*76+r);
            end
        end
        thiserr=thiserr+abs(2*estTarget((g-1)*48+t)-thisest);
    end
end
if (thiserr<errR10)
    pen=100;
    errR10=thiserr;
    LR10=lambda100R10;
end

lambda1000R10=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda1000_Rect10.txt');
thiserr=0;
for t=1:48
    for g=1:7
        thisest=0;
        for r=1:76
            for p=1:3
                thisest=thisest+lambda1000R10((g-1)*48+t,(p-1)*76+r);
            end
        end
        thiserr=thiserr+abs(2*estTarget((g-1)*48+t)-thisest);
    end
end
if (thiserr<errR10)
    pen=1000;
    errR10=thiserr;
    LR10=lambda1000R10;
end


lambda10R10=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda10_Rect10.txt');
thiserr=0;
for t=1:48
    for g=1:7
        thisest=0;
        for r=1:76
            for p=1:3
                thisest=thisest+lambda10R10((g-1)*48+t,(p-1)*76+r);
            end
        end
        thiserr=thiserr+abs(2*estTarget((g-1)*48+t)-thisest);
    end
end
if (thiserr<errR10)
    pen=10;
    errR10=thiserr;
    LR10=lambda10R10;
end

lambda1R10=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda1_Rect10.txt');
thiserr=0;
for t=1:48
    for g=1:7
        thisest=0;
        for r=1:76
            for p=1:3
                thisest=thisest+lambda1R10((g-1)*48+t,(p-1)*76+r);
            end
        end
        thiserr=thiserr+abs(2*estTarget((g-1)*48+t)-thisest);
    end
end
if (thiserr<errR10)
    pen=1;
    errR10=thiserr;
    LR10=lambda1R10;
end

LR10s=zeros(48*7,3);
for t=1:48
    for g=1:7
        for r=1:76
            for p=1:3
                %thisest=thisest+lambda1000R10((g-1)*48+t,(p-1)*76+r);
                %thisest=thisest+lambda1R10((g-1)*48+t,(p-1)*76+r);
                LR10s((g-1)*48+t,p)=LR10s((g-1)*48+t,p)+LR10((g-1)*48+t,(p-1)*76+r);
            end
        end
    end
end


Emps=zeros(48*7,3);
for t=1:48
    for g=1:7
        for r=1:4916
            for p=1:3
                %thisest=thisest+lambda1000R10((g-1)*48+t,(p-1)*76+r);
                %thisest=thisest+lambda1R10((g-1)*48+t,(p-1)*76+r);
                Emps((g-1)*48+t,p)=Emps((g-1)*48+t,p)+estimated(t,g,r,p);
            end
        end
    end
end


type = zeros(R,1);
%regressors = zeros(nbLandTypes, R);
neighbors=cell(1,R);
distance = zeros(R,R);

neighbors_file = fopen(neighbors_file_path);
continueBool=1;

while continueBool
    str=fgetl(neighbors_file);
    [line]=strsplit (str);
    if (length(str) == 3 & str=='END')
        continueBool=0;
    else
        i=str2num(line{1,1})+1;
        coordinates_centers(i).lat=str2num(line{1,2});
        coordinates_centers(i).long=str2num(line{1,3});
%         land_type=str2num(line{1,4})+1;
%         type(i)=land_type;
%         neighbors{1,i}=[];
%         for l=1:nbLandTypes
%             regressors(l,i)=str2num(line{1,l+4});
%         end
%         for l=5+nbLandTypes:2:length(line)
%             j=str2num(line{1,l})+1;
%             d=str2num(line{1,l+1});
%             neighbors{1,i}=[neighbors{1,i},j];
%             distance(i,j)=d;
%         end
    end
end

fclose(neighbors_file);

estimated=zeros(T,7,160,3);
for t=1:48
    for g=1:7
        init=(g-1)*48;
        for r=1:160
            for p=1:3
                estimated(t,g,r,p)=LHD((g-1)*48+t,(p-1)*160+r);
            end
        end
    end
end


LCov=zeros(48*7,3);
estimated=zeros(T,7,76,3);
for t=1:48
    for g=1:7
        init=(g-1)*48;
        for r=1:76
            for p=1:3
                rateest=0;
                for j=1:(1+nbLandTypes)
                    %rateest=rateest+LCH7(indexBeta(p,g,t,j))*regressors(j,r);
                    rateest=rateest+lambda1CovR10(indexBeta(p,g,t,j))*regressors(j,r);
                end
                estimated(t,g,r,p)=rateest;
            end
        end
    end
end

LR100Ps=zeros(48*7,3);
for t=1:48
    for g=1:7
        thisest=0;
        for r=1:4916
            for p=1:3
                %thisest=thisest+lambda1000R10((g-1)*48+t,(p-1)*76+r);
                %thisest=thisest+lambda1R10((g-1)*48+t,(p-1)*76+r);
                LR100Ps((g-1)*48+t,p)=LR100Ps((g-1)*48+t,p)+LR100((g-1)*48+t,(p-1)*4916+r);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%

heatlambdas=zeros(1,76);
for r=1:76
    for t=1:48
        for g=1:7
            for p=1:3
                heatlambdas(r)=heatlambdas(r)+LR10((g-1)*48+t,(p-1)*76+r);
            end
        end
    end
end
dlmwrite('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\figures\Heatmaps\heatRect10x10.txt',heatlambdas);


heatlambdaspriorities=zeros(3,76);
for r=1:76
    for t=1:48
        for g=1:7
            for p=1:3
                heatlambdaspriorities(p,r)=heatlambdaspriorities(p,r)+LR10((g-1)*48+t,(p-1)*76+r);
            end
        end
    end
end
dlmwrite('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\figures\Heatmaps\heatPRect10x10.txt',heatlambdaspriorities);




estTargetPs=zeros(48*7,1);
for t=1:48
    for g=1:7
        for r=1:76
            for p=1:3
                estTargetPs((g-1)*48+t)=estTargetPs((g-1)*48+t)+estimated(t,g,r,p);
            end
        end
    end
end

weeklambdaprioritites=zeros(T*G,P);
for t=1:T
    for g=1:G
        init=(g-1)*T;
        for r=1:R
            for p=1:P
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rateest=rateest+x(indexBeta(p,g,t,j))*regressors(j,r);
                end
                weeklambdaprioritites(init+t,p)=weeklambdaprioritites(init+t,p)+rateest;
            end
        end
    end
end

p=3
plot(2*estTargetPs(:,p),'-k');
hold on;
plot(2*LR10Ps(:,p),'--r');
hold on;
plot(2*LR100Ps(:,p),'-.b');
hold on;
plot(2*weeklambdaprioritites(:,p));
hold on;





%%%%%%

errR100=inf;
lambda1R100=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda1_Rect100.txt');
thiserr=0;

for t=1:48
    for g=1:7
        thisest=0;
        for r=1:4916
            for p=1:3
                thisest=thisest+lambda1R100((g-1)*48+t,(p-1)*4916+r);
            end
        end
        thiserr=thiserr+abs(2*estTarget((g-1)*48+t)-thisest);
    end
end
if (thiserr<errR100)
    pen=1;
    errR100=thiserr;
    LR100=lambda1R100;
end

lambda1000R100=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda1000_Rect100.txt');
thiserr=0;

for t=1:48
    for g=1:7
        thisest=0;
        for r=1:4916
            for p=1:3
                thisest=thisest+lambda1000R100((g-1)*48+t,(p-1)*4916+r);
            end
        end
        thiserr=thiserr+abs(2*estTarget((g-1)*48+t)-thisest);
    end
end
if (thiserr<errR100)
    pen=1000;
    errR100=thiserr;
    LR100=lambda1000R100;
end

LR100P=zeros(48*7,3);
for t=1:48
    for g=1:7
        thisest=0;
        for r=1:4916
            for p=1:3
                thisest=thisest+LR100((g-1)*48+t,(p-1)*4916+r);
                LR100P((g-1)*48+t,p)=LR100P((g-1)*48+t,p)+LR100((g-1)*48+t,(p-1)*4916+r);
                %thisest=thisest+lambda1R10((g-1)*48+t,(p-1)*76+r);
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errH7=inf;
lambda1Hex7=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda1_hex7.txt');
thiserr=0;
for t=1:48
    for g=1:7
        thisest=0;
        for r=1:297
            for p=1:3
                thisest=thisest+lambda1Hex7((g-1)*48+t,(p-1)*297+r);
            end
        end
        thiserr=thiserr+abs(2*estTarget((g-1)*48+t)-thisest);
    end
end
if (thiserr<errH7)
    pen=1;
    errH7=thiserr;
    LH7=lambda1Hex7;
end

lambda104Hex7=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda10000_hex7.txt');
thiserr=0;
for t=1:48
    for g=1:7
        thisest=0;
        for r=1:297
            for p=1:3
                thisest=thisest+lambda104Hex7((g-1)*48+t,(p-1)*297+r);
            end
        end
        thiserr=thiserr+abs(2*estTarget((g-1)*48+t)-thisest);
    end
end
if (thiserr<errH7)
    pen=104;
    errH7=thiserr;
    LH7=lambda104Hex7;
end

LH7P=zeros(48*7,1);

for t=1:48
    for g=1:7
        thisest=0;
        for r=1:297
            for p=1:3
                thisest=thisest+LH7((g-1)*48+t,(p-1)*297+r);
            end
        end
        LH7P((g-1)*48+t)=thisest;
    end
end

heatmap=zeros(1,R);
for t=1:48
    for g=1:7
        for r=1:R
            for p=1:3
                %rateest=0;
                %for j=1:(1+nbLandTypes)
                %    rateest=rateest+LR100((g-1)*48+t,(p-1)*4916+r);
                %end
                heatmap(r)=heatmap(r)+LH7((g-1)*48+t,(p-1)*R+r);  
            end
        end
    end
end
dlmwrite('C:\Users\vince\Desktop\Figures_Calls\heatH7.txt',heatmap);


LH8P=zeros(48*7,1);
for t=1:48
    for g=1:7
        thisest=0;
        for r=1:297
            for p=1:3
                thisest=thisest+LH8((g-1)*48+t,(p-1)*297+r);
            end
        end
        LH8P((g-1)*48+t)=thisest;
    end
end

LCov=zeros(48*7,3);
for t=1:48
    for g=1:7
        init=(g-1)*48;
        thiserr=0
        for r=1:76
            for p=1:3
                rateest=0;
                for j=1:(1+nbLandTypes)
                    %rateest=rateest+LCH7(indexBeta(p,g,t,j))*regressors(j,r);
                    rateest=rateest+lambda1CovR10(indexBeta(p,g,t,j))*regressors(j,r);
                end
                thiserr=thiserr+rateest;
            end
        end
        LCov((g-1)*48+t)=thiserr;
    end
end


LH7Ps=zeros(48*7,3);
for t=1:48
    for g=1:7
        for r=1:297
            for p=1:3
                LH7Ps((g-1)*48+t,p)=LH7Ps((g-1)*48+t,p)+LH7((g-1)*48+t,(p-1)*297+r);
            end
        end
    end
end

errH8=inf;
lambda1Hex8=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda1_hex8.txt');
thiserr=0;
for t=1:48
    for g=1:7
        thisest=0;
        for r=1:1811
            for p=1:3
                thisest=thisest+lambda1Hex8((g-1)*48+t,(p-1)*1811+r);
            end
        end
        thiserr=thiserr+abs(2*estTarget((g-1)*48+t)-thisest);
    end
end
if (thiserr<errH8)
    pen=1;
    errH8=thiserr;
    LH8=lambda1Hex8;
end

lambda104Hex8=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda10000_hex8.txt');
thiserr=0;
for t=1:48
    for g=1:7
        thisest=0;
        for r=1:1811
            for p=1:3
                thisest=thisest+lambda104Hex8((g-1)*48+t,(p-1)*1811+r);
            end
        end
        thiserr=thiserr+abs(2*estTarget((g-1)*48+t)-thisest);
    end
end
if (thiserr<errH8)
    pen=104;
    errH8=thiserr;
    LH8=lambda104Hex8;
end


LH8Ps=zeros(48*7,1);
for t=1:48
    for g=1:7
        thisest=0;
        for r=1:1811
            for p=1:3
                %thisest=thisest+lambda104Hex8((g-1)*48+t,(p-1)*1811+r);
                thisest=thisest+LH8((g-1)*48+t,(p-1)*1811+r);
            end
        end
        LH8P((g-1)*48+t)=thisest;
    end
end


lambda1Hex8=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda1_hex8.txt');
LH8Ps=zeros(48*7,3);
for t=1:48
    for g=1:7
        for r=1:1811
            for p=1:3
                LH8Ps((g-1)*48+t,p)=LH8Ps((g-1)*48+t,p)+LH8((g-1)*48+t,(p-1)*1811+r);
                %thisest=thisest+LH8((g-1)*48+t,(p-1)*1811+r);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%

errD=inf;
lambda1D=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda1_District.txt');
thiserr=0;
for t=1:48
    for g=1:7
        thisest=0;
        for r=1:160
            for p=1:3
                thisest=thisest+lambda1D((g-1)*48+t,(p-1)*160+r);
            end
        end
        thiserr=thiserr+abs(2*estTarget((g-1)*48+t)-thisest);
    end
end
if (thiserr<errD)
    pen=1;
    errD=thiserr;
    LHD=lambda1D;
end

lambda103D=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda1000_District.txt');
thiserr=0;
for t=1:48
    for g=1:7
        thisest=0;
        for r=1:160
            for p=1:3
                thiest=thisest+lambda103D((g-1)*48+t,(p-1)*160+r);
            end
        end
        thiserr=thiserr+abs(2*estTarget((g-1)*48+t)-thisest);
    end
end
if (thiserr<errD)
    pen=103;
    errD=thiserr;
    LHD=lambda103D;
end

LDP=zeros(48*7,1);
for t=1:48
    for g=1:7
        thisest=0;
        for r=1:160
            for p=1:3
                thisest=thisest+LHD((g-1)*48+t,(p-1)*160+r);
            end
        end
        LDP((g-1)*48+t)=thisest;
    end
end

LDPs=zeros(48*7,3);
for t=1:48
    for g=1:7
        for r=1:160
            for p=1:3
                LDPs((g-1)*48+t,p)=LDPs((g-1)*48+t,p)+LHD((g-1)*48+t,(p-1)*160+r);        
            end
        end
    end
end


plot(estTarget,'-k');
hold on;
plot(LR10P,'--r');
hold on;
plot(LR100P,'-.b');
hold on;
plot(LH7P,'-.k');
hold on:
plot(LH8P,'-.k');
hold on:
plot(LDisP,'-.k');

%Cov R10
lambda1CovR10=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda1_CovRect10.txt');

%Cov H
errCov=inf
lambda1000_CovDist=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda1000_CovDist.txt');
thiserr=0

for t=1:48
    for g=1:7
        init=(g-1)*48;
        thisest=0;
        for r=1:76
            for p=1:3
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rateest=rateest+lambda1000_CovDist(indexBeta(p,g,t,j))*regressors(j,r);
                end
                thisest=thisest+rateest;
            end
        end
        thiserr=thiserr+abs(estTarget((g-1)*48+t)-thisest);
    end
end
if (thiserr<errCov)
    pen=1000;
    errCov=thiserr;
    LCov=lambda1000_CovDist;
end

lambda100_CovDist=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda100_CovDist.txt');
thiserr=0
for t=1:48
    for g=1:7
        init=(g-1)*48;
        thisest=0
        for r=1:160
            for p=1:3
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rateest=rateest+lambda100_CovDist(indexBeta(p,g,t,j))*regressors(j,r);
                end
                thisest=thisest+rateest;
            end
        end
        thiserr=thiserr+abs(estTarget((g-1)*48+t)-thisest);
    end
end
if (thiserr<errCov)
    pen=100;
    errCov=thiserr;
    LCov=lambda100_CovDist;
end

lambda10_CovDist=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda10_CovDist.txt');

thiserr=0
for t=1:48
    for g=1:7
        init=(g-1)*48;
        thisest=0
        for r=1:160
            for p=1:3
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rateest=rateest+lambda10_CovDist(indexBeta(p,g,t,j))*regressors(j,r);
                end
                thisest=thisest+rateest;
            end
        end
        thiserr=thiserr+abs(estTarget((g-1)*48+t)-thisest);
    end
end
if (thiserr<errCov)
    pen=10;
    errCov=thiserr;
    LCov=lambda10_CovDist;
end

lambda1_CovDist=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda1_CovDist.txt');
thiserr=0
for t=1:48
    for g=1:7
        init=(g-1)*48;
        thisest=0
        for r=1:160
            for p=1:3
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rateest=rateest+lambda1_CovDist(indexBeta(p,g,t,j))*regressors(j,r);
                end
                thisest=thisest+rateest;
            end
        end
        thiserr=thiserr+abs(estTarget((g-1)*48+t)-thisest);
    end
end
if (thiserr<errCov)
    pen=1;
    errCov=thiserr;
    LCov=lambda1_CovDist;
end

for t=1:T
    for r=1:R
        for g=1:7
            for p=1:3
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rateest=rateest+LCov(indexBeta(p,g,t,j))*regressors(j,r);
                end
                estimated(t,g,r,p)=rateest;    
            end
        end
    end
end

LCov


LCovs=zeros(48*7,3);
for t=1:48
    for g=1:7
        for r=1:160
            for p=1:3
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rateest=rateest+LCov(indexBeta(p,g,t,j))*regressors(j,r);
                end
                LCovs((g-1)*48+t,p)=LCovs((g-1)*48+t,p)+rateest;
            end
        end
    end
end




errCov=inf;
lambda1000=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda1000_CovHex7.txt');
LCH7=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda1000_CovHex7.txt');
thiserr=0;
for t=1:48
    for g=1:7
        init=(g-1)*48;
        thisest=0;
        for r=1:297
            for p=1:3
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rateest=rateest+lambda1000(indexBeta(p,g,t,j))*regressors(j,r);
                end
                thisest=thisest+rateest;
            end
        end
        thiserr=thiserr+abs(estTarget((g-1)*48+t)-thisest);
    end
end
if (thiserr<errCov)
    pen=1000;
    errCov=thiserr;
    LCH7=lambda1000;
end

lambda10=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda10_CovHex7.txt');
LCH7=dlmread('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\lambda10_CovHex7.txt');
thiserr=0;
for t=1:48
    for g=1:7
        init=(g-1)*48;
        thisest=0;
        for r=1:297
            for p=1:3
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rateest=rateest+lambda10(indexBeta(p,g,t,j))*regressors(j,r);
                end
                thisest=thisest+rateest;
            end
        end
        thiserr=thiserr+abs(estTarget((g-1)*48+t)-thisest);
    end
end
if (thiserr<errCov)
    pen=10;
    errCov=thiserr;
    LCH7=lambda10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


LH7s=zeros(48*7,1);
for t=1:48
    for g=1:7
        thisest=0;
        for r=1:297
            for p=1:3
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rateest=rateest+2*LCH7(indexBeta(p,g,t,j))*regressors(j,r);
                end
                thisest=thisest+rateest; 
            end
        end
        LH7s((g-1)*48+t)=thisest;
    end
end



%Cov District

legend(['Empirical';'Rect 10  ';'Rect 100 ';'Hex 7    '])


weeklambdaprioritites=zeros(T,P);
for t=1:T
    for r=1:R
        for p=1:P
            weeklambdaprioritites(t,p)=weeklambdaprioritites(t,p)+x(t,r,p);  
        end
    end
end
plot(weeklambdaprioritites(:,1));
hold on;
plot(weeklambdaprioritites(:,2));
hold on;
plot(weeklambdaprioritites(:,3));
legend(['Priority 1';'Priority 2';'Priority 3']);

x=ones(T,R,P);
nbObservationsTotal=size(sampleCalls,4);
proportion=0.2;
[cputime,alpha,weight,lambda]=crossValidation(nbObservationsTotal,sampleCalls,neighbors,type,distance,T,R,P,sigma,betaBar,x,iterMax,proportion,epsilon,durations,Groups,whichgroup);

%Write sum of lambda, summing in t,p.


heatlambdas=zeros(1,76);
for r=1:76
    for t=1:48
        for g=1:7
            for p=1:3
                heatlambdas(r)=heatlambdas(r)+LR10((g-1)*48+t,(p-1)*76+r);
            end
        end
    end
end
dlmwrite('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\figures\Heatmaps\heatRect10x10.txt',heatlambdas);


heatlambdaspriorities=zeros(3,4916);
for r=1:4916
    for t=1:48
        for g=1:7
            for p=1:3
                heatlambdaspriorities(p,r)=heatlambdaspriorities(p,r)+LR100((g-1)*48+t,(p-1)*4916+r);
            end
        end
    end
end
dlmwrite('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\figures\Heatmaps\heatPRect100x100.txt',heatlambdaspriorities);


heatlambdas=zeros(1,4916);
for r=1:4916
    for t=1:48
        for g=1:7
            for p=1:3
                heatlambdas(r)=heatlambdas(r)+LR100((g-1)*48+t,(p-1)*4916+r);
            end
        end
    end
end
dlmwrite('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\figures\Heatmaps\heatRect100x100.txt',heatlambdas);


LR100Ps=zeros(48*7,3);
for t=1:48
    for g=1:7
        for r=1:4916
            for p=1:3
                %thisest=thisest+lambda1000R10((g-1)*48+t,(p-1)*76+r);
                %thisest=thisest+lambda1R10((g-1)*48+t,(p-1)*76+r);
                LR100Ps((g-1)*48+t,p)=LR100Ps((g-1)*48+t,p)+LR100((g-1)*48+t,(p-1)*4916+r);
            end
        end
    end
end


heatlambdaspriorities=zeros(3,4916);
for r=1:4916
    for t=1:48
        for g=1:7
            for p=1:3
                heatlambdaspriorities(p,r)=heatlambdaspriorities(p,r)+LR100((g-1)*48+t,(p-1)*4916+r);
            end
        end
    end
end
dlmwrite('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\figures\Heatmaps\heatPRect100x100.txt',heatlambdaspriorities);

heatlambdas=zeros(1,4916);
for t=1:48
    for g=1:7
        for r=1:4916
            for p=1:3
                heatlambdas(r)=heatlambdas(r)+estimated(t,g,r,p);
            end
        end
    end
end
dlmwrite('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\figures\Heatmaps\heatEmpRect100x100.txt',heatlambdas);

heatlambdaspriorities=zeros(3,4916);
for t=1:48
    for g=1:7
        for r=1:4916
            for p=1:3
                heatlambdaspriorities(p,r)=heatlambdaspriorities(p,r)+estimated(t,g,r,p);
            end
        end
    end
end
dlmwrite('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\figures\Heatmaps\heatPEmpRect100x100.txt',heatlambdaspriorities);



estTarget=zeros(48*7,3);
for t=1:48
    for g=1:7
        for r=1:R
            for p=1:3
                estTarget((g-1)*48+t,p)=estTarget((g-1)*48+t,p)+estimated(t,g,r,p);
            end
        end
    end
end

heatlambdasprioritites=zeros(1,R);
for t=1:48
    for g=1:7
        for r=1:160
            for p=1:3
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rateest=rateest+2*LCov(indexBeta(p,g,t,j))*regressors(j,r);
                end
                heatlambdasprioritites(r)=heatlambdasprioritites(r)+rateest;
            end
        end
    end
end
dlmwrite('C:\Users\vince\Desktop\Figures_Calls\heatCovDistrict.txt',heatlambdasprioritites);

heatlambdasprioritites=zeros(1,R);
for t=1:48
    for g=1:7
        for r=1:160
            for p=1:3
                heatlambdasprioritites(r)=heatlambdasprioritites(r)+LHD((g-1)*48+t,(p-1)*160+r);        
            end
        end
    end
end
dlmwrite('C:\Users\vince\Desktop\Figures_Calls\heatDistrict.txt',heatlambdasprioritites);


heatlambdasprioritites=zeros(3,R);
for t=1:48
    for g=1:7
        for r=1:R
            for p=1:3
                heatlambdasprioritites(p,r)=heatlambdasprioritites(p,r)+2*estimated(t,g,r,p);        
            end
        end
    end
end
dlmwrite('C:\Users\vince\Desktop\Figures_Calls\heatPEmpH8.txt',heatlambdasprioritites);


estimated(t,g,r,p)


heatlambdasprioritites=zeros(1,76);
for t=1:48
    for g=1:7
        init=(g-1)*T;
        for r=1:76
            for p=1:3
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rateest=rateest+lambda1CovR10(indexBeta(p,g,t,j))*regressors(j,r);
                end
                heatlambdasprioritites(1,r)=heatlambdasprioritites(1,r)+rateest;
            end
        end
    end
end
dlmwrite('C:\Users\vince\Desktop\Figures_Calls\heatCovRect10x10.txt',heatlambdasprioritites);



weeklambdaprioritites=zeros(T*G,P);
for t=1:T
    for g=1:G
        init=(g-1)*T;
        for r=1:R
            for p=1:P
                rateest=0;
                for j=1:(1+nbLandTypes)
                    rateest=rateest+x(indexBeta(p,g,t,j))*regressors(j,r);
                end
                weeklambdaprioritites(init+t,p)=weeklambdaprioritites(init+t,p)+rateest;
            end
        end
    end
end


heatlambdas=zeros(1,76);
for r=1:76
    for t=1:48
        for g=1:7
            for p=1:3
                heatlambdas(r)=heatlambdas(r)+LR10((g-1)*48+t,(p-1)*76+r);
            end
        end
    end
end
dlmwrite('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\figures\Heatmaps\heatRect10x10.txt',heatlambdas);


heatlambdaspriorities=zeros(3,76);
for r=1:76
    for t=1:48
        for g=1:7
            for p=1:3
                heatlambdaspriorities(p,r)=heatlambdaspriorities(p,r)+LR10((g-1)*48+t,(p-1)*76+r);
            end
        end
    end
end
dlmwrite('C:\Users\vince\Dropbox\Articles_Math\Anton_Vincent_Yao\Call_Models\figures\Heatmaps\heatPRect10x10.txt',heatlambdaspriorities);

