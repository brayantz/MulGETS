% fit the equation of y=a+bX^n, n is the oder of X

function [ahat]=funct(x,y,or)

mdl = @(a,x)(a(1)+abs(a(2))*x.^or); % function 

a0 = [1;1]; % intial parameter values 

[ahat] = nlinfit(x,y,mdl,a0); % nonlinear least-squares regression.
end
