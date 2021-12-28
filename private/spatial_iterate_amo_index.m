function [mat,K,am] = spatial_iterate_amo_index(C,nstations,random,n,...
    threshold,occ,phat,dist,Tolerance,maxiter)
%
% this function iterates from a starting correlation matrix, and modifies it
% so that it becomes the correlation matrix of random numbers giving
% an amount correlation matrix equal to the original
%
%

kiter=0.1;   % iteration parameter in calculation of new estimate of matrix 'mat'
mat=C;  % start with observed correlation matrix as best estimate of random numbers correlation matrix
am = genam(mat, phat, occ, dist, random, n, nstations, threshold);
K=corrcoef(am);
K=forcecorr(K);

val=max(max(abs(K-C)));
%count=0;
ii=0;
while val>Tolerance && ii<maxiter
    ii=ii+1;         % 'live' counter fo follow progress on screen
    
    %     for i=1:nstations
    %         for j=1:n
    %             if occurences(j,i)==1
    %                 if dist==1 % exponential distribution
    %                     prec(j,i)=-log(1-randU01(j,i))/phat.alpha(j,i)+ threshold;
    %                 elseif dist==2 % gamma distribution
    %                     prec(j,i)=gaminv(randU01(j,i),phat.alpha(j,i),phat.beta(j,i))+ threshold;
    %                 end
    %             end
    %         end
    %     end
    mat=mat+kiter*(C-K);
    mat(1:(nstations+1):end) = 1; % make sure diagonal 1.
    mat(mat>1)=0.999; % avoid solutions with correlation coefficient greater than 1
    
    am = genam(mat, phat, occ, dist, random, n, nstations, threshold);
    K=corrcoef(am);
    K=forcecorr(K);
    
    val=max(max(abs(K-C)));
end

end


function am = genam(mat, phat, occ, dist, random, n, nstations, threshold)

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
% at this point, for each station we will generate a series of length
% 'n' of precipitation using a matrix of occurrrences already
% determined
%

am = zeros(n, nstations);
if dist == 1
    for i = 1:nstations
        tf = occ(:,i) == 1;
        if any(tf)
            m = -log(1 - unif_random(tf, i))./phat.par1(tf, i);
            m(m < 0) = 0;
            am(tf, i) = m + threshold;
        end
    end
elseif dist == 2
    for i = 1:nstations
        tf = occ(:,i) == 1;
        if any(tf)
            m = gaminv(unif_random(tf, i), phat.par1(tf, i), phat.par2(tf, i));
            m(m < 0) = 0;
            am(tf, i) = m + threshold;
        end
    end
end

end













