function TTpr_gen = mulgets(TTpr_obs, dist, nYears, iniYear, graph)
% Multisite weather generator of the Ecole de technologie superieure - MulGETS
% (Precipitation only).
% 
% Input arguments
%   TTpr_obs    - Observed precipitation as timetable.
%   Dist        - Distribution type: 1 (Multi-exponential) or 2 (Multi-gamma).
%   nYears      - Number of years to generate.
%   pivotYear   - Initial year.
%   graph       - Display graphs.
%
% Output arguments
%   TTpr_gen    - Generated precipitation as timetable.


%

[nTimes, nSites] = size(TTpr_obs);
sitesId = strcat({'S'}, arrayfun(@(x) {num2str(x)}, 1:nSites));
dv = datevec(TTpr_obs.Properties.RowTimes);
dv(:, 4:6) = [];

Obs = [];
nn = nan(nTimes, 2);
for i = 1:nSites
    Obs.(sitesId{i}) = [dv, nn, TTpr_obs{:,i}];
end

%
Gen = modifiedMulgets(Obs, nYears, iniYear, dist, graph);
sitesName = TTpr_obs.Properties.VariableNames;
TTpr_gen = array2timetable(Gen.pr, 'RowTimes', Gen.dt, 'VariableNames', sitesName);

end
