function [mat,K,occ] = spatial_iterate_occ(C,nstations,random,...
    transitions,n,Tolerance,maxiter)
%
% this function iterates from a starting correlation matrix, and modifies it
% so that it becomes the correlation matrix of random numbers giving
% an occurence correlation matrix equal to the original
%
%

kiter=0.1;   % iteration parameter in calculation of new estimate of matrix 'mat'
mat=C;  % start with observed correlation matrix as best estimate of random numbers correlation matrix
occ = genocc(mat, transitions, random, n, nstations);
K=corrcoef(occ);
K= forcecorr(K);

val=max(max(abs(K-C)));
%count=0;
ii=0;
while val>Tolerance && ii<maxiter % change value of random number matrix, only if computations are to continue
    ii=ii+1;         % 'live' counter fo follow progress on screen
    
%     for i=1:nstations
%         for j=1:nstations
%             if i==j
%                 mat(i,j)=1;
%             else
%                 mat(i,j)=mat(i,j)+kiter*(C(i,j)-K(i,j));
%             end
%         end
%     end
    
    mat=mat+kiter*(C-K);
    mat(1:(nstations+1):end) = 1; % make sure diagonal 1.
    mat(mat>1)=0.999; % avoid solutions with correlation coefficient greater than 1
    
    occ = genocc(mat, transitions, random, n, nstations);
    K=corrcoef(occ);
    K = forcecorr(K);
    
    val=max(max(abs(K-C)));
    
end
%mat=mat./sqrt(diag(mat)*diag(mat)');   % have to really think whether or not this line is needed.

end


%%
function occ = genocc(mat, transitions, random, n, nstations)
% p00 and p10 are vector, seed is tf (wet or dry).

% insure that matrix is positive semi-definite (all eigen values >= 0).  if
% not, cholesky deconposition will not work
[V,D] = eig(mat);
i=find(D<0);      % find eigenvalues smaller than 0
if length(i)>=1
    %count=count+1;
    D(i)=0.000001;    % they should be only a little bit below zero
    Q=(V*D)/V;   % reconstruct matrix, which becomes a covariance matrix diagonal slighly larger than 1
    % we keep on iterating with the covariance matrix and
    % only noramlize at the end, allows for better
    % accuracy
    %mat=Q./sqrt(diag(Q)*diag(Q)');   % transform covariance matrix in correlation matrix that is positive definite
    mat=Q;
end

U=chol(mat);
L=U';   % lower diagonal matrix  L*U=C

values=L*random;
values=values';

corr_random=values;

mean_random=mean(corr_random);
std_random=std(corr_random);
unit_vector=ones(n,1);
mean_matrix=unit_vector*mean_random;
std_matrix=unit_vector*std_random;
norm_random=(corr_random-mean_matrix)./(std_matrix);   % this matrix contains the N(0,1) serially independant but spatially correlated random numbers
unif_random = 0.5*erfc(-norm_random./sqrt(2)); % unifom random

%
% at this point, for each station,:a series of length 'n' of precip
% occurence will be generated p00 p10
%
%     occ2=zeros(n,nstations);


%     for i=1:nstations
%         for j=2:n
%             if occ(j-1,i)==0
%                 if unif_random(j,i)>transitions(1,i)
%                     occ(j,i)=1;
%                 end
%             else
%                 if unif_random(j,i)>transitions(2,i)
% 	                occ(j,i)=1;
%                 end
%             end
%         end
%     end


p00 = repmat(transitions(1,:), n, 1);
p10 = repmat(transitions(2,:), n, 1);
tf00 = unif_random > p00; % wet (1) or dry (0).
tf10 = unif_random > p10;

tf = tf00 | tf10;
pcond = false(n, nstations);
pcond(2:end,:) = tf(1:end-1,:);
occ = tf00;
occ(pcond) = tf10(pcond);
end






