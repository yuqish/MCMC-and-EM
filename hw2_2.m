% Generating data 
y = csvread('mixture-observations.csv');


% EM algorithm 
g1=@(y) normpdf(y,1,2);
g0=@(y) normpdf(y,0,1);
M = 100;
theta = zeros(1,M);
theta(1) = 0.5;
for m = 1:(M - 1)
    theta(m + 1) = mean(theta(m).*g1(y)./(g0(y).*(1-theta(m))+g1(y).*theta(m))); 
end

% Plotting EM learning curves 

figure
plot(theta);

title('EM learning curve');
xlabel('Iteration number');
ylabel('\theta');
legend('EM estimates');