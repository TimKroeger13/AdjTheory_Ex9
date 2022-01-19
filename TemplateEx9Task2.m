%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%   Exercise 9: Adjustment Calculation - part IV  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 11, 2018
%   Last changes   : January 11, 2022
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;
format longG;

%--------------------------------------------------------------------------
%   Task 2 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   Observations and initial values for the unknowns
%--------------------------------------------------------------------------
%Load data
data = load("Exercise9Task2.txt");

%Plot the data
x = data(:,1);
y = data(:,2);

plot(x,y,'-b')
hold on;
scatter(x,y,'b','filled')
Buffer = 5; % in pecrent
xlim([min(x)-(max(x)- min(x)) * Buffer / 100, ...
    max(x)+(max(x)- min(x)) * Buffer / 100]);
ylim([-min(y)-(max(y)- min(y)) * Buffer / 100, ...
    max(y)+(max(y)- min(y)) * Buffer / 100]);
stem(x,y,'b')
hold off;
xlabel('X values');
ylabel('Y values');
title('Task 2: Unknown values as a plot');
legend_cell = {'','Data points'};
legend(legend_cell,'Location','northwest')

%functional model:
%a*exp(b*x) + c;

%Error-free values
%x

%Vector of observations
L = y;

%Number of observations
no_n = length(y);

%Initial values for the unknowns  
a = 1.7637;
b = log(2);
c = 0;

%Vector of initial values for the unknowns
X_0 = [a,b,c]';

%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n - no_u;

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
s_y = 0.15;  %[m]
S_LL = eye(no_n) * s_y^2;

%Theoretical standard deviation
sigma_0 = 1;     %a priori

%Cofactor matrix of the observations
Q_LL = (1/sigma_0^2) * S_LL;

%Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
epsilon = 10^-8;
delta = 10^-12;
max_x_hat = inf;
Check2 = inf;

%Number of iterations
iteration = 0;

while max_x_hat > epsilon || Check2 > delta
    
     %Observations as functions of the approximations for the unknowns
     L_0 = a*exp(b*x) + c;
     
     %Vector of reduced observations
     l = L - L_0;
    
     %Design matrix with the elements from the Jacobian matrix J
     A(:,1) = exp(b*x);
     A(:,2) = a*x .* exp(b*x);
     A(:,3) = 1;
    
     %Normal matrix
     N = A' * P * A;
     
     %Vector of right hand side of normal equations
     n =  A' * P * l;
    
     %Inversion of normal matrix / Cofactor matrix of the unknowns
     Q_xx = inv(N);
    
     %Solution of the normal equations
     x_hat = Q_xx * n;
       
     %Update
     X_0 = X_0 + x_hat;

     a = X_0(1);
     b = X_0(2);
     c = X_0(3);
    
     %Check 1
     max_x_hat = max(abs(x_hat));
     
     %Vector of residuals
     v = A*x_hat - l;
 
     %Vector of adjusted observations
     L_hat = L+v;
    
     %Objective function
     vTPv = v'*P*v;
    
     %Functional relationships without the observations
     Psy_X_hat = a*exp(b*x) + c;

     %Check 2
     Check2 = max(abs(L_hat - Psy_X_hat));
    
     %Update number of iterations
     iteration = iteration+1;
  
end

if Check2<=delta
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end

%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);

%VC matrix of adjusted unknowns
S_XX_hat = s_0^2 * Q_xx;

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));

%Cofactor matrix of adjusted observations
Q_LL_hat = A * Q_xx * A';

%VC matrix of adjusted observations
S_LL_hat = s_0^2 * Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));

%Cofactor matrix of the residuals
Q_vv = Q_LL - Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2 * Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));

%--------------------------------------------------------------------------
%  Plot
%--------------------------------------------------------------------------


%Residual Plot

Fu = @(x) a*exp(b*x) + c;
fplot(Fu,'red')
xlabel('X-values');
ylabel('Y-values');
title('Task 2: exponational groth plot');
Buffer = 5; % in pecrent
xlim([min(x)-(max(x)- min(x)) * Buffer / 100, ...
    max(x)+(max(x)- min(x)) * Buffer / 100]);
ylim([-min(L_0)-(max(L_0)- min(L_0)) * Buffer / 100, ...
    max(L_0)+(max(L_0)- min(L_0)) * Buffer / 100]);
hold on;
plot(x,Psy_X_hat-v,'bo','linewidth',1);
for i = 1:no_n

    plot([x(i),x(i)],[Psy_X_hat(i),Psy_X_hat(i)-v(i)],'blue')

end
legend_cell = {'Functional model','Data points'};
legend(legend_cell,'Location','northwest')
hold off;







bar(sort(v),'blue')

Buffer = 5; % in pecrent
xlim([0 + 0 * Buffer / 100, ...
    length(v) + length(v) * Buffer / 100]);
ylim([min(sort(v)) + min(sort(v)) * Buffer / 100, ...
    max(sort(v)) + max(sort(v)) * Buffer / 100]);
title('Task 2: Resudial Plot');
ylabel('Difference to the observation');
hold off;



%Tables
table(L,L_hat,s_L_hat,v,s_v,'VariableNames',...
{'L'...
'L_hat',...
's_L_hat',...
'v',...
's_v'})

table([1.7637,log(2),0]',X_0,s_X,'VariableNames',...
{'X',...
'X_hat',...
's_X_hat'},...
'RowNames',{'a','b','c'})








