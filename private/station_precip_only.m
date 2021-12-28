% separate the generated precipitation and temperature by station

function generation=station_separate_precip(data,years_sim,pivotYear,nstations,length_month)

% combine seasonal precip
startDate = datetime([pivotYear, 1, 1]);
endYear = datetime([pivotYear + years_sim - 1, 12, 31]);
dt = (startDate:days:endYear)';
dv = datevec(dt);

nDays= size(dt, 1); % the number of days to be generated
pr=zeros(nDays,nstations);

%
tf = dv(:,2) == 12; % dec
pr(tf,:) = data(1).pre(1:years_sim*length_month(12),:);

tf = dv(:,2) == 1; % jan
pr(tf,:) = data(1).pre(years_sim*length_month(12) + 1:years_sim*length_month(12)+years_sim*length_month(1),:);

tf = dv(:,2) == 2 & dv(:,3) ~= 29; % feb
pr(tf,:) = data(1).pre(years_sim*length_month(12) + years_sim*length_month(1) + 1:end,:);

%
tf = dv(:,2) == 3 ; % mar
pr(tf,:) = data(2).pre(1:years_sim*length_month(3),:);

tf = dv(:,2) == 4; % apr
pr(tf,:) = data(2).pre(years_sim*length_month(3) + 1:years_sim*length_month(3) + years_sim*length_month(4),:);

tf = dv(:,2) == 5; % may
pr(tf,:) = data(2).pre(years_sim*length_month(3) + years_sim*length_month(4) + 1:end,:);

%
tf = dv(:,2) == 6; % jun
pr(tf,:) = data(3).pre(1:years_sim*length_month(6),:);

tf = dv(:,2) == 7; % jul
pr(tf,:) = data(3).pre(years_sim*length_month(6) + 1:years_sim*length_month(6) + years_sim*length_month(7),:);

tf = dv(:,2) == 8; % aug
pr(tf,:) = data(3).pre(years_sim*length_month(6) + years_sim*length_month(7) + 1:end,:);

%
tf = dv(:,2) == 9; % sep
pr(tf,:) = data(4).pre(1:years_sim*length_month(9),:);

tf = dv(:,2) == 10; % oct
pr(tf,:) = data(4).pre(years_sim*length_month(9) + 1:years_sim*length_month(9) + years_sim*length_month(10),:);

tf = dv(:,2) == 11; %nor
pr(tf,:) = data(4).pre(years_sim*length_month(9) + years_sim*length_month(10) + 1:end,:);

generation.pr = pr;
generation.dv = dv;
generation.dt = dt;

% %% put precip and temperature in a file for each station
% nn = nan(size(series, 1), 1);
% for i=1:nstations
%     station = ['S' num2str(i)];
%     S=zeros(nDays,3);
%     S=addymd(S);
%     S(:,4)=nn; % generated Tmax
%     S(:,5)=nn; % generated Tmin
%     S(:,6)=series(:,i); % generated precip
%     generation.(station)=S;
end
