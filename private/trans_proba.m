%
% calculate monthly p00 and p10 parameters for precip occurrence, and the
% mean and stdev of the exponential and gamma distribution for precip amounts

function [trans_prob, Sam, Socc]=trans_proba(observation,stationname,nstations,...
    threshold,begin_month,length_month,dist,months)


Sp = [];
% extract precip for each station
for i=1:nstations
    m = observation.(stationname{i}); 
    [dat]=feb29_treat(m); % treat Feb 29th, the model does not account for Feb 29th
    pre=dat(:,6);
    pre=reshape(pre,365,[]);
    pre=pre';
    Sp.(stationname{i}) = pre;   
end

% preliminary treatment create matrices 'shortname'occ for occurences, with
% NaN as missing data using the threshold specified above

Sam = [];
Socc = [];
for i=1:nstations
    am = Sp.(stationname{i});    
    occ = nan(size(am));    
    occ(am > threshold) = 1;
    occ(am <= threshold) = 0;
    am(am <= threshold) = 0;
    Sam.(stationname{i}) = am;
    Socc.(stationname{i}) = occ;
end

%
% calculate transition matrices and exponential or gamma distribution
%
trans_prob = struct('prob', cell(12, 1), 'month', cell(12, 1));
%phat_value = [];
for i=1:nstations
    for imonth=1:12
        % occurence
        occ = Socc.(stationname{i});
        %am = Sam.(stationname{i});
        occmonth=occ(:,begin_month(imonth):begin_month(imonth)+length_month(imonth)-1);
        %ammonth=am(:,begin_month(imonth):begin_month(imonth)+length_month(imonth)-1);

        [p00,p10,~] = transition(occmonth);

        trans_prob(imonth).prob(1,i) = p00;
        trans_prob(imonth).prob(2,i) = p10;
        trans_prob(imonth).month = months{imonth};

        % exponential or gamma distribution parameters
%         qq=find(ammonth > 0);
%         if dist==1
%             phat=1/mean(ammonth(qq)-threshold);
%             phat_value(imonth).phat(1,i)=phat;
%         elseif dist==2
%             phat=gamfit(ammonth(qq)-threshold);
%             phat_value(imonth).phat(1:2,i)=phat; 
%         end
%         phat_value(imonth).month=months{imonth};
       
    end
    
end

end






