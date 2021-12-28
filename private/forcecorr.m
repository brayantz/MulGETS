function m = forcecorr(m)
% Square matrix of correlations.

m(isnan(m)) = 0.001; % Dry months, force small value.
m(m <= 0) = 0.001; % No sense negative corr.
m(1:(size(m, 1)+1):end) = 1; % make sure diagonal 1.

end
