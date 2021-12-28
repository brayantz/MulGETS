function generation = baseMulgets(observation, years_sim, pivotYear, dist, graph)
% *************************************************************************
% Multisite weather generator of the Ecole de technologie superieure
%                               (MulGETS)
% *************************************************************************
%
% MulGETS is Matlab-based multisite weather generator for generating
% spatially correlated daily precipitation, maximum and minimum
% temperatures (Tmax and Tmin) series. The algorithm originates from the
% Wilks approach proposed in 1998. In 2007, Francois Brissette et al.
% presented an algorithm for efficient generation of multisite
% precipitation data following the Wilks approach. Afterward, A minor tune
% was conducted in 2013 after a further evaluation over various climates. A
% component for generating multisite temperature was also added.
%
% The Matlab code was programmed by Dr. Jie Chen, followed a work version
% provided by Prof. Francois Brissette, at the Ecole de technologie
% superieure, University of Quebec
% Contact email: jie.chen@etsmtl.ca (Jie Chen)
%
% ****************************
% Input data
% ****************************
% The input data consists of daily precipitation, Tmax and Tmin for
% multisite, meteorological data shall be separated by stations with a matlab
% structure named "observation", the data of each station shall be provided
% with the order of year, month, day, Tmax, Tmin and precipitation
% Missing data should be assigned as NaN.
%
% ****************************
% Output data
% ****************************
% The output consists of daily precipitation, Tmax and Tmin, the generated
% meteorological time series is separated by stations with a matlab structure
% named "generation", the order of each station is year, month, day, Tmax,
% Tmin and precipitation
%
% References:
% (1) Wilks, D. S., 1998. Multisite generalization of a daily stochastic
% precipitation generation model. J. Hydrol. 210, 178-191.
% (2) Brissette, F.P., Khalili, M., Leconte, R., 2007. Efficient stochastic
% generation of multi-site sythetic precipitation data. J. Hydrol. 345,
% 121-133.
% (3) Chen, J., Brissette, F. P.,Leconte, R., Caron, A., 2012. A versatile
% weather generator for daily precipitation and temperature. Transactions
% of the ASABE. 55(3): 895-906.
% (4) Chen, J., Brissette, F. P.,Zhang, X.-C., 2014. A multi-site stochastic
% weather generator for daily precipitation and temperature. Transaction of
% the ASABE (Accepted)


% declare several important variables
% load(filenamein)
sitename = fields(observation);
nstations=numel(sitename);
months={'Jan' ;'Feb'; 'Mar'; 'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'};
length_month=[31 28 31 30 31 30 31 31 30 31 30 31]; % days of each month
begin_month=[1 32 60 91 121 152 182 213 244 274 305 335]; % the first julin day of each month
Tolerance=0.001;% convergence threshold
maxiter=200; % convergence threshold
% graph = 1;
threshold = 0.1;

%% basic processing of inputs
% find the start year and the end year for each station
begin_year=zeros(1,nstations);
end_year=begin_year;
for i=1:nstations
    begin_year(1,i)=eval(['observation.S',num2str(i),'(1,1);']); % begin year of the dataset
    end_year(1,i)=eval(['observation.S',num2str(i),'(end,1);']); % end year of the dataset
end
%
% define the stations from S to Sn
for i=1:nstations
    stationname{i,1}=['S',num2str(i)];
end

% ************************************************************************
%         START MULTISITE PRECIPITATION GENERATION (SIX STEPS)
% *************************************************************************

%% step 1: determination of weather generator parameters: p00, p10 (eq.3)
% and precip distribution function on a monthly scale
%
[trans_prob, Sam, Socc]=trans_proba(observation,stationname,nstations,threshold,...
    begin_month,length_month,dist,months);

%% step 2: computation of correlation matrices of precipatation occurrence
% and amounts (eq.1)
[corr_occ,corr_amounts]=corr_precip(Sam, Socc, stationname,nstations,threshold,...
    begin_month,length_month,months,begin_year,end_year);

%% step 3: generate spatially correlated precipitation occurrence
% the step 3 includes two sub-steps:
% (1) automatic determination of correlation matrix of random number using
% algorithm (eq.7).Diagonalization and repalcement of negative eigenvalues
% if matrix is non-positive definite (eqs.5 & 6)
% (2) generate precipitation occurrence using first-order Markov chain with
% Cholesky factorization (eq.2)
[corr_occ_rand,corr_occ_gen,occurrences_gen]=multisite_occ_generation(corr_occ,...
    trans_prob,nstations,length_month,years_sim,months,Tolerance,maxiter,graph);%, folderOut);

%% step 4: establish link between occurrence index (eq.10) and average
% precipitation amounts for each station and construct the multi-exponential
% or multi-gamma distribution for each station (eq.11)
[phat_multiexp,site_bins,occ_index,SamSeason,SoccSeason,corrOcc, occurences_genSeason]=multi_distribution_par(Sam, Socc,...
    stationname,occurrences_gen,nstations,threshold,dist,graph);%,folderOut);

%% step 5: generate precipitation amounts based on the occurrence index
[seasonal_cor,corr_amounts_rand_index,corr_amounts_gen_index,...
    precip_gen_index]=multisite_occ_index(SamSeason,SoccSeason,corrOcc,phat_multiexp,...
    corr_amounts,stationname,occurences_genSeason,nstations,threshold,years_sim,...
    Tolerance,maxiter,site_bins,begin_year,end_year,dist,graph);%, folderOut);

% *************************************************************************
%                   END MULTISITE PRECIPITATION GENERATION
% *************************************************************************

generation=station_precip_only(precip_gen_index,years_sim,pivotYear,nstations,length_month);

end



