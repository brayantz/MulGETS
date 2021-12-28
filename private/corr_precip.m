%
% calculate the spatial correlation matrices for precip occurences and
% amounts for each month  

function [corr_occ,corr_amounts]=corr_precip(Sam, Socc, stationname,nstations,...
    threshold,begin_month,length_month,months,begin_year,end_year)

% % extract precip for each station
% for i=1:nstations
%     S=eval(['observation.S',num2str(i),';']); 
%     [dat]=feb29_treat(S); % treat Feb 29th, the model does not account for Feb 29th
%     pre=dat(:,6);
%     pre=reshape(pre,365,[]);
%     pre=pre';
%     eval([char(stationname(i)) '=pre;']);   
% end
% 
% 
% % preliminary treatment create matrices 'stationname'occ for occurences, with
% % NaN as missing data and 'stationname'amounts fo0r amounts with NaN as
% % missing data, all using a threshold sepcified above
% 
% for i=1:nstations
%     eval(['mat=' char(stationname(i)) ';' ]);    
%     amounts=mat;
%     occ=nan(size(mat));
%     j=find(amounts>threshold);
%     occ(j)=1;
%     
%     j=find(amounts<=threshold);
%     occ(j)=0;
%     amounts(j)=0;
%     eval([char(stationname(i)) 'amounts=amounts;' ]);
%     eval([char(stationname(i)) 'occ=occ;' ]);    
% end

%
% calculate monthly correlation matrices (same dates for all stations)
%
corr_occ = struct('cor', cell(12, 1), 'month', cell(12, 1));
corr_amounts = struct('cor', cell(12, 1), 'month', cell(12, 1));

for imonth = 1:12   
    occ_stn = cell(1, nstations);
    am_stn = cell(1, nstations);
    for istn = 1:nstations
        occ = Socc.(stationname{istn});
        occ = occ(:,begin_month(imonth):begin_month(imonth)+length_month(imonth)-1);
        occ = occ(:);
        occ_stn{istn} = occ;
        
        am = Sam.(stationname{istn});
        am = am(:,begin_month(imonth):begin_month(imonth)+length_month(imonth)-1);
        am = am(:);
        am_stn{istn} = am;
    end
    occ_stn = horzcat(occ_stn{:});
    ro = corrcoef(occ_stn, 'Rows', 'pairwise');
    ro = forcecorr(ro);
            
    am_stn = horzcat(am_stn{:});
    ra = corrcoef(am_stn, 'Rows', 'pairwise');
    ra = forcecorr(ra);
    
    corr_occ(imonth).cor = ro;
    corr_amounts(imonth).cor = ra;
    
    corr_occ(imonth).month = months{imonth};
    corr_amounts(imonth).month = months{imonth};
end


% for i=1:nstations
%     for j=1:nstations
%         eval(['station_i_occ=' char(stationname(i)) 'occ;' ]);
%         eval(['station_i_amounts=' char(stationname(i)) 'amounts;' ]);
%         eval(['station_j_occ=' char(stationname(j)) 'occ;' ]);
%         eval(['station_j_amounts=' char(stationname(j)) 'amounts;' ]);        
%         
%         % find intersection of year period between stations       
%         yearstart=max([begin_year(i) begin_year(j)]);
%         yearend=min([end_year(i) end_year(j)]);
%         years_i=[begin_year(i):end_year(i)];
%         years_j=[begin_year(j):end_year(j)];
%         ifind=find(years_i>=yearstart & years_i<=yearend);
%         jfind=find(years_j>=yearstart & years_j<=yearend);
%         nyears=length(ifind);
%         station_i_occ=station_i_occ(ifind,:);
%         station_i_amounts=station_i_amounts(ifind,:);
%         station_j_occ=station_j_occ(jfind,:);
%         station_j_amounts=station_j_amounts(jfind,:);
%         
%         % calculate monthly correlations for both occurences and amounts
%         
%         for imonth=1:12
%            % start with occurence
%            if nyears~=0   % if there is an intersection between station data
%                occ_i=station_i_occ(:,begin_month(imonth):begin_month(imonth)+length_month(imonth)-1);
%                occ_j=station_j_occ(:,begin_month(imonth):begin_month(imonth)+length_month(imonth)-1);
%                occ_i=reshape(occ_i',length_month(imonth)*nyears,1);
%                occ_j=reshape(occ_j',length_month(imonth)*nyears,1);
%                
%                ij=find(occ_i>=0 & occ_j>=0 );  % find all elements where there is no missing values in both
%                corr=corrcoef(occ_i(ij),occ_j(ij));
%                corr=corr(1,2);
%            else
%                corr=-999;
%            end
%            eval(['corr_occ_' char(months(imonth)) '(i,j)=corr;']);
%            corr_occ(imonth).cor(i,j)=corr;
%            corr_occ(imonth).month=char(months(imonth));
%            
%            % do amounts
%            if nyears~=0   % if there is an intersection between station data
%                amounts_i=station_i_amounts(:,begin_month(imonth):begin_month(imonth)+length_month(imonth)-1);
%                amounts_j=station_j_amounts(:,begin_month(imonth):begin_month(imonth)+length_month(imonth)-1);
%                amounts_i=reshape(amounts_i',length_month(imonth)*nyears,1);
%                amounts_j=reshape(amounts_j',length_month(imonth)*nyears,1);
%                
%                ij=find(amounts_i>=0 & amounts_j>=0 );  % find all elements where there is precip data in both
%                corr=corrcoef(amounts_i(ij),amounts_j(ij));               
%                corr=corr(1,2);
%            else
%                corr=-999;
%            end
%            corr_amounts(imonth).cor(i,j)=corr;
%            corr_amounts(imonth).month=char(months(imonth));
%         end
%         
%     end
% end

end
