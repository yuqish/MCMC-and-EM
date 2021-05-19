clc;
clear;

y = readtable('coal-mine.csv');
tau = table2array(y);

%parameters
hyper=0.1;
rho=0.05;

N=10000;

plotind=1;

for d=[2,3,4,5] %no.of break points=d-1
tinit=[];
delta=round((1963-1851)/d);
temp=1851;
for i=1:d-1
    temp = temp+delta;
    tinit = [tinit; temp*ones(1,N)];
end
tvec=[1851*ones(1,N);tinit;1963*ones(1,N)];
lambdavec=zeros(d,N);
theta=zeros(1,N);
theta(1)=gamrnd(2,hyper);
lambdavec(:,1) = gamrnd(2,theta(1),d,1);
acc=0;
for i=1:N-1
    %draw lambda,theta with Gibbs
    for j=1:d
        n=sum(numel(find(tau>=tvec(j,i) & tau<tvec(j+1,i))));
        lambdavec(j,i+1)=gamrnd(n+2,1/(tvec(j+1,i)-tvec(j,i)+theta(i)));
    end
    theta(i+1)=gamrnd(2*d+2,1/(sum(lambdavec(:,i+1))+hyper));
    % Draw t with RW-MH
    for j=2:d
        Li=getL(lambdavec(:,i+1),tvec(:,i),tau);
        R=rho*(tvec(j+1,i)-tvec(j-1,i));
        tprop=tvec(j,i)+floor((2*R+1)*rand)-R;
        tvecTemp=tvec(:,i);
        tvecTemp(j)=tprop;
        Lp=getL(lambdavec(:,i+1),tvecTemp,tau);
        alpha=0;
        if (tprop>tvec(j-1,i+1) && tprop<tvec(j+1,i)) %check boundary condition
            alpha=min(1,Lp/Li);
        end
        if (rand<alpha)
            acc=acc+1;
            tvec(j,i+1)=tprop;
        else
            tvec(j,i+1)=tvec(j,i);
        end
    end
end

figure(1);
acc_rate=acc/(N*(d-1))
subplot(2,4,plotind);
plot(1:N,tvec);
xlabel('iteration');
ylabel('t');
title(['Number of Breakpoints = ',num2str(d-1)])

subplot(2,4,plotind+4);
histogram(tvec(2:end-1,:));
xlabel('t');
ylabel('count');
title(['Number of Breakpoints = ',num2str(d-1)])

figure(2);
subplot(2,4,plotind);
plot(1:N,theta);
xlabel('iteration');
ylabel('\theta');
title(['Number of Breakpoints = ',num2str(d-1)])

subplot(2,4,plotind+4);
histogram(theta);
xlabel('\theta');
ylabel('count');
title(['Number of Breakpoints = ',num2str(d-1)])

figure(3);
subplot(2,4,plotind);
plot(1:N,lambdavec);
xlabel('iteration');
ylabel('\lambda');
title(['Number of Breakpoints = ',num2str(d-1)])

subplot(2,4,plotind+4);
histogram(lambdavec);
xlabel('\lambda');
ylabel('count');
title(['Number of Breakpoints = ',num2str(d-1)])
plotind=plotind+1;
end

function L=getL(lambdavec,tvec,tau)
d=length(tvec)-1;
nvec=zeros(d,1);
for i = 1:d
    nvec(i) = sum(tau >= tvec(i) & tau < tvec(i+1));
end
L=prod(lambdavec.^nvec.*diff(tvec).*exp(-lambdavec.*diff(tvec)));
end