addpath 'C:\Users\vince\Dropbox\Softwares\LASPATED\Model_Calibration\Matlab'
cd 'C:\gurobi901\win64\matlab\'
gurobi_setup

p1=[];
p2=[];
p3=[];
for g=1:G
    for t=1:T
        p1=[p1;probs(t,g,1)];
        p2=[p2;probs(t,g,2)];
        p3=[p3;probs(t,g,3)];
    end
end

plot([1:T*G],p1,'-k');
hold on;
plot([1:T*G],p2,'--r');
hold on;
plot([1:T*G],p3,'-.b');
hold on;
plot([1:T*G],prob*ones(1,T*G),'-k');
legend(['High priority calls        ';'Intermediate priority calls';'Low priority calls         ';'p                          ']);
ylim([0 0.7])
xlabel('Time');
ylabel('Probability of not reporting the location');


for g=1:G
    for t=1:T
        lambda=0;
        for r=1:R
            for p=1:P
                lambda=lambda+lambdas_missing(t,g,r,p);
            end
        end
        %All p
        meanmissing((g-1)*T+t)=lambda;
        q5missing((g-1)*T+t)=poissinv(0.05,lambda);
        q95missing((g-1)*T+t)=poissinv(0.95,lambda);
        %p1
        lambda=0;
        for r=1:R
            lambda=lambda+lambdas_missing(t,g,r,1);
        end
        meanmissing1((g-1)*T+t)=lambda;
        q5missing1((g-1)*T+t)=poissinv(0.05,lambda);
        q95missing1((g-1)*T+t)=poissinv(0.95,lambda);
        %p2
        lambda=0;
        for r=1:R
            lambda=lambda+lambdas_missing(t,g,r,2);
        end
        meanmissing2((g-1)*T+t)=lambda;
        q5missing2((g-1)*T+t)=poissinv(0.05,lambda);
        q95missing2((g-1)*T+t)=poissinv(0.95,lambda);
        %p3
        lambda=0;
        for r=1:R
            lambda=lambda+lambdas_missing(t,g,r,3);
        end
        meanmissing3((g-1)*T+t)=lambda;
        q5missing3((g-1)*T+t)=poissinv(0.05,lambda);
        q95missing3((g-1)*T+t)=poissinv(0.95,lambda);
        
        
        lambda=0;
        for r=1:R
            for p=1:P
                lambda=lambda+lambdas_nomissing(t,g,r,p);
            end
        end
        %All p
        meannomissing((g-1)*T+t)=lambda;
        q5nomissing((g-1)*T+t)=poissinv(0.05,lambda);
        q95nomissing((g-1)*T+t)=poissinv(0.95,lambda);
        %p1
        lambda=0;
        for r=1:R
            lambda=lambda+lambdas_nomissing(t,g,r,1);
        end
        meannomissing1((g-1)*T+t)=lambda;
        q5nomissing1((g-1)*T+t)=poissinv(0.05,lambda);
        q95nomissing1((g-1)*T+t)=poissinv(0.95,lambda);
        %p2
        lambda=0;
        for r=1:R
            lambda=lambda+lambdas_nomissing(t,g,r,2);
        end
        meannomissing2((g-1)*T+t)=lambda;
        q5nomissing2((g-1)*T+t)=poissinv(0.05,lambda);
        q95nomissing2((g-1)*T+t)=poissinv(0.95,lambda);
        %p3
        lambda=0;
        for r=1:R
            lambda=lambda+lambdas_nomissing(t,g,r,3);
        end
        meannomissing3((g-1)*T+t)=lambda;
        q5nomissing3((g-1)*T+t)=poissinv(0.05,lambda);
        q95nomissing3((g-1)*T+t)=poissinv(0.95,lambda);
    end
end

plot([1:T*G],meanmissing3,'-k');
hold on;
plot([1:T*G],q5missing3,'--k');
hold on;
plot([1:T*G],q95missing3,'--k');
hold on;
plot([1:T*G],meannomissing3,'-r');
hold on;
plot([1:T*G],q5nomissing3,'--r');
hold on;
plot([1:T*G],q95nomissing3,'--r');
legend(['Mean  - Corrected  ';'q0.05 - Corrected  ';'q0.95 - Corrected  ';'Mean  - Uncorrected';'q0.05 - Uncorrected';'q0.95 - Uncorrected']);
ylim([0 55])

lambda=0;
for g=1:G
    for t=1:T
        for r=1:R
            for p=1:P
                lambda=lambda+lambdas_missing(t,g,r,p);
            end
        end
    end
end
lambdan=0;
for g=1:G
    for t=1:T
        for r=1:R
            for p=1:P
                lambdan=lambdan+lambdas_nomissing(t,g,r,p);
            end
        end
    end
end
x=[4000:20:8000];
y=poisspdf(x,lambda);
bar(x,y,0.2,'r');
x=[4000:20:8000];
y=poisspdf(x,lambdan);
hold on;
bar(x+0.5,y,0.2,'k');
xlabel('Number of calls');
ylabel('Probabilities');


prior=1;
lambda=0;
for g=1:G
    for t=1:T
        for r=1:R
            lambda=lambda+lambdas_missing(t,g,r,prior);
        end
    end
end
lambdan=0;
for g=1:G
    for t=1:T
        for r=1:R
            lambdan=lambdan+lambdas_nomissing(t,g,r,prior);
        end
    end
end
x=[1800:5:3200];
y=poisspdf(x,lambda);
bar(x,y,0.2,'r');
x=[1800:5:3200];
y=poisspdf(x,lambdan);
hold on;
bar(x+0.5,y,0.2,'k');
xlabel('Number of calls');
ylabel('Probabilities');
legend(['Corrected  ';'Uncorrected']);

lambda=0;
prior=2;
for g=1:G
    for t=1:T
        for r=1:R
            lambda=lambda+lambdas_missing(t,g,r,prior);
        end
    end
end
lambdan=0;
for g=1:G
    for t=1:T
        for r=1:R
            lambdan=lambdan+lambdas_nomissing(t,g,r,prior);
        end
    end
end
x=[500:5:1500];
y=poisspdf(x,lambda);
bar(x,y,0.2,'r');
x=[500:5:1500];
y=poisspdf(x,lambdan);
hold on;
bar(x+0.5,y,0.2,'k');
xlabel('Number of calls');
ylabel('Probabilities');

prior=3;
lambda=0;
for g=1:G
    for t=1:T
        for r=1:R
            lambda=lambda+lambdas_missing(t,g,r,prior);
        end
    end
end
lambdan=0;
for g=1:G
    for t=1:T
        for r=1:R
            lambdan=lambdan+lambdas_nomissing(t,g,r,prior);
        end
    end
end
x=[1000:10:3000];
y=poisspdf(x,lambda);
bar(x,y,0.2,'r');
x=[1000:10:3000];
y=poisspdf(x,lambdan);
hold on;
bar(x+0.5,y,0.2,'k');
xlabel('Number of calls');
ylabel('Probabilities');

%Lambda missing, nomissing rect 10x10

fileorig = fopen('C:\Users\vince\Dropbox\Softwares\LASPATED\Model_Calibration\Matlab\lambdar10.txt','w');
filem = fopen('C:\Users\vince\Dropbox\Softwares\LASPATED\Model_Calibration\Matlab\lambdamr10.txt','w');

for g=1:G
    for t=1:T
        for r=1:R 
            for p=1:P    
               fprintf(fileorig,'%d %d %d %d %f\n',g,t,r,p,lambdas_nomissing(t,g,r,p));
               fprintf(filem,'%d %d %d %d %f\n',g,t,r,p,lambdas_missing(t,g,r,p));
            end
        end
    end
end
fclose(fileorig);
fclose(filem);

fileorig = fopen('C:\Users\vince\Dropbox\Softwares\LASPATED\Model_Calibration\Matlab\lambdar100.txt','w');
filem = fopen('C:\Users\vince\Dropbox\Softwares\LASPATED\Model_Calibration\Matlab\lambdamr100.txt','w');

for g=1:G
    for t=1:T
        for r=1:R 
            for p=1:P    
               fprintf(fileorig,'%d %d %d %d %f\n',g,t,r,p,lambdas_nomissing(t,g,r,p));
               fprintf(filem,'%d %d %d %d %f\n',g,t,r,p,lambdas_missing(t,g,r,p));
            end
        end
    end
end
fclose(fileorig);
fclose(filem);

fileorig = fopen('C:\Users\vince\Dropbox\Softwares\LASPATED\Model_Calibration\Matlab\lambdah8.txt','w');
filem = fopen('C:\Users\vince\Dropbox\Softwares\LASPATED\Model_Calibration\Matlab\lambdamh8.txt','w');

for g=1:G
    for t=1:T
        for r=1:R 
            for p=1:P    
               fprintf(fileorig,'%d %d %d %d %f\n',g,t,r,p,lambdas_nomissing(t,g,r,p));
               fprintf(filem,'%d %d %d %d %f\n',g,t,r,p,lambdas_missing(t,g,r,p));
            end
        end
    end
end
fclose(fileorig);
fclose(filem);

fileorig = fopen('C:\Users\vince\Dropbox\Softwares\LASPATED\Model_Calibration\Matlab\lambdadis.txt','w');
filem = fopen('C:\Users\vince\Dropbox\Softwares\LASPATED\Model_Calibration\Matlab\lambdamdis.txt','w');

for g=1:G
    for t=1:T
        for r=1:R 
            for p=1:P    
               fprintf(fileorig,'%d %d %d %d %f\n',g,t,r,p,lambdas_nomissing(t,g,r,p));
               fprintf(filem,'%d %d %d %d %f\n',g,t,r,p,lambdas_missing(t,g,r,p));
            end
        end
    end
end
fclose(fileorig);
fclose(filem);

S=10;
C=P;
sigma=0.5;
iterMax=10;
epsilon=0.0001;
x=epsilon*ones(T,G,R,P);
prob=epsilon*ones(P,G,T);
alpha=1;
durations=0.5*ones(G,T);
nbGroups=8;
Groups=cell(1,8);
whichgroup=zeros(T,G);
weights=ones(1,nbGroups);
for i=1:8
    Groups{1,i}=[];
end

for g=1:5
    for t=13:20
        str.g=g;
        str.t=t;
        whichgroup(t,g)=1;
        Groups{1,1}=[Groups{1,1},str];
    end
end


for g=1:5
    for t=21:36
        str.g=g;
        str.t=t;
        whichgroup(t,g)=2;
        Groups{1,2}=[Groups{1,2},str];
    end
end

for g=1:5
    for t=37:44
        str.g=g;
        str.t=t;
        whichgroup(t,g)=3;
        Groups{1,3}=[Groups{1,3},str];
    end
end

for g=1:4
    for t=45:48
        str.g=g;
        str.t=t;
        whichgroup(t,g)=4;
        Groups{1,4}=[Groups{1,4},str];
    end
end

for t=45:48
      str.g=7;
      str.t=t;
      whichgroup(t,7)=4;
      Groups{1,4}=[Groups{1,4},str];
end

for g=1:5
    for t=1:12
        str.g=g;
        str.t=t;
        whichgroup(t,g)=4;
        Groups{1,4}=[Groups{1,4},str];
    end
end

for t=45:48
    str.g=5;
    str.t=t;
    whichgroup(t,5)=5;
    Groups{1,5}=[Groups{1,5},str];
end

for t=45:48
    str.g=6;
    str.t=t;
    whichgroup(t,6)=5;
    Groups{1,5}=[Groups{1,5},str];
end

for t=1:12
    str.g=6;
    str.t=t;
    whichgroup(t,6)=5;
    Groups{1,5}=[Groups{1,5},str];
end

for t=1:12
    str.g=7;
    str.t=t;
    whichgroup(t,7)=5;
    Groups{1,5}=[Groups{1,5},str];
end

for t=13:20
    str.g=6;
    str.t=t;
    whichgroup(t,6)=6;
    Groups{1,6}=[Groups{1,6},str];
end

for t=13:20
    str.g=7;
    str.t=t;
    whichgroup(t,7)=6;
    Groups{1,6}=[Groups{1,6},str];
end


for t=21:36
    str.g=6;
    str.t=t;
    whichgroup(t,6)=7;
    Groups{1,7}=[Groups{1,7},str];
end

for t=21:36
    str.g=7;
    str.t=t;
    whichgroup(t,7)=7;
    Groups{1,7}=[Groups{1,7},str];
end

for t=37:44
    str.g=6;
    str.t=t;
    whichgroup(t,6)=8;
    Groups{1,8}=[Groups{1,8},str];
end

for t=37:44
    str.g=7;
    str.t=t;
    whichgroup(t,7)=8;
    Groups{1,8}=[Groups{1,8},str];
end
probs=(1/R)*ones(1,R);

[x,fVal]=projectedGradientMissingLocation(nbObservations,nbCalls,nb_missing_calls,neighbors,G,T,R,C,sigma,x,prob,iterMax,alpha,epsilon,durations,Groups,whichgroup,weights);
[x,fVal]=projectedGradientMissingLocationModel2(nbObservations,nbCalls,nb_missing_calls,neighbors,G,T,R,C,sigma,x,probs,S,iterMax,alpha,epsilon,durations,sample_missing_calls,sample_calls);

% lambda=epsilon*ones(T,G,R,P);
% delta=0.0001;
% [fL]=oracleObjectiveMissingModel2(nbObservations,nbCalls,nb_missing_calls,neighbors,lambda,T,R,P,G,alpha,durations,S,probs,sample_missing_calls,sample_calls);
% [gradientl]=oracleGradientMissingModel2(nbObservations,nbCalls,nb_missing_calls,neighbors,lambda,G,T,R,P,alpha,durations,S,probs,sample_missing_calls,sample_calls);
% gradapprox=zeros(T,G,R,P);
% lambdaaux=lambda;
% for t=1:T
%     g=1;
%     r=1;
%     p=1;
% %    for g=1:G
% %        for r=1:R
% %            for p=1:P
%                 size(lambda)
%                 size(lambdaaux)
%                 lambdaaux(t,g,r,p)=lambda(t,g,r,p)+delta;
%                 [fL]=oracleObjectiveMissingModel2(   nbObservations,nbCalls,nb_missing_calls,neighbors,lambda,   T,R,P,G,alpha,durations,S,probs,sample_missing_calls,sample_calls)
%                 [fLaux]=oracleObjectiveMissingModel2(nbObservations,nbCalls,nb_missing_calls,neighbors,lambdaaux,T,R,P,G,alpha,durations,S,probs,sample_missing_calls,sample_calls);
%                 gradapprox(t,g,r,p)=(-fL+fLaux)/delta;      
%                 lambdaaux(t,g,r,p)=lambda(t,g,r,p);
% %            end
% %        end
% %    end
% end
% 
% err=[];
% for t=1:T
%     err=[err;gradapprox(t,g,r,p)-gradientl(t,g,r,p)];
% end
% plot(err);



