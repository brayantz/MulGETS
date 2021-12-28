function [p00,p10,total] = transition(occ)
%
% adapted by Francois Brissette from transition.m written by Robert Leconte
%
% This function calculates transition probabilities p00 and p10 from a matrix of occurence occ of 
% any dimension as long as rows are years and columns days. Takes into
% account NaN as long as they are positive or negtive values exceeding 2 or
% smaller than -2
%
% total is just a verification parameter.  If there are no NaN total
% should be equal to : (j-1)*i 
%

%
[i,j]=size(occ);
occ1=occ(:,1:j-1);
occ2=occ(:,2:j);

%
% subtract occ1-occ2:   -1 correspond to n01  1 correspond to n10
%
diffocc=occ1-occ2;

[k]=find(diffocc == 1);
n01=length(k);
[k]=find(diffocc == -1);
n10=length(k);
%
% scalar multiplication of occ1 and occ2.  1 correspond to n11
%
prodocc=occ1.*occ2;
[k]=find(prodocc == 1);
n11=length(k);
%

[k]=find(prodocc == 0 & diffocc == 0);  % finds n00
n00=length(k);

p10=n10/(n10+n11);
p00=n00/(n00+n01);
total=n01+n11+n00+n10;
end


