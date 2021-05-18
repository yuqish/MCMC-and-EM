clc;
clear;
close all;

y = readtable('coal-mine.csv');
tau = table2array(y);
d=4;
hyper=1;
N=1000;
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
rho=1;
for i=1:N-1
    % Draw t with RW-MH
    for j=2:d
        Li=getL(lambdavec(j,i),tvec(j,i),tvec(j+1,i),tau);
        R=rho*(tvec(j+1,i)-tvec(j-1,i));
        tprop=tvec(j,i)+floor((2*R+1)*rand)-R;
        Lp=getL(lambdavec(j,i),tprop,tvec(j+1,i),tau);
        
        if (rand<Lp/Li)
            tvec(j,i+1)=tprop;
        else
            tvec(j,i+1)=tvec(j,i);
        end
    end
    
    %draw lambda,theta with Gibbs
    for j=1:d
        n=sum(numel(find(tau>=tvec(j,i+1) & tau<tvec(j+1,i+1))));
        lambdavec(j,i+1)=gamrnd(n+2,1/(tvec(j+1,i+1)-tvec(j,i+1)+theta(i)));
    end
    theta(i+1)=gamrnd(2*d+2,1/(sum(lambdavec(:,i+1))+hyper));
end


function L=getL(lambda,t,t1,tau)
n=sum(numel(find(tau>=t & tau<t1)));
L=lambda^n*(t1-t)*exp(-lambda*(t1-t));
end