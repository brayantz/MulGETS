%
% generate precipitation amounts based on the occurrence index

function [seasonal_cor,corr_amounts_rand_index,corr_amounts_gen_index,...
    precip_gen_index] = multisite_occ_index(SamSeason,SoccSeason,corrOcc,phat_multiexp,corr_amounts,...
    stationname,occurences_genSeason,nstations,threshold,years_sim,Tolerance,...
    maxiter,site_bins,begin_year,end_year,dist,graph)%,folderOut)

% calculate seasonal correlation for precip amount
seasons={'DJF' ;'MAM'; 'JJA'; 'SON'};
nseasons = numel(seasons);
seasonal_cor = struct('cor', cell(nseasons, 1), 'season', cell(nseasons, 1));
for i = 1:nseasons
    am = SamSeason.(seasons{i});
    corr = corrcoef(am, 'Rows', 'pairwise');
    corr = forcecorr(corr);
    
    seasonal_cor(i).cor=corr;
    seasonal_cor(i).season= seasons{i};
end


corr_amounts_rand_index = struct('cor', cell(nseasons, 1), 'season', cell(nseasons, 1));
corr_amounts_gen_index = struct('cor', cell(nseasons, 1), 'season', cell(nseasons, 1));
precip_gen_index = struct('pre', cell(nseasons, 1),  'season', cell(nseasons, 1));
for ijk=1:nseasons
        
    C=seasonal_cor(ijk).cor;
    
    % obtain correlation matrix of generated precip occurrence.
    tempo = corrOcc.(seasons{ijk});
    
    % calculate moran for each occurence
    occ = occurences_genSeason.(seasons{ijk});
    nStepsGen = size(occ, 1);
    
    moran=nan(nStepsGen,nstations);
    for i=1:nstations
        tempo2 = tempo(i,:);
        tempo2(i) = [];
        for j=1:nStepsGen
            thisOcc = occ(j,:);
            hasPrecip = thisOcc(i) == 1;
            if hasPrecip
                thisOcc(i) = [];
                moran(j,i)=thisOcc*tempo2'/sum(tempo2);
               % to avoid rewriting below
            end
        end
    end
    
    
    % load precip parameters as a function of moran
    
%     for i=1:nstations
%         alpha(i,:).alpha=phat_multiexp(i).(seasons{ijk})(1,:);
%         if dist==2 % gamma distibution
%             beta(i,:).beta=phat_multiexp(i).(seasons{ijk})(2,:);
%         end
%     end
    
    % calculate right away the values of lambda based on moran to avoid doing
    % it every step of the way during the iterations
    phatGen = [];
    if dist==1 % exponential distribution
        par1 = nan(nStepsGen,nstations); 
        for i = 1:nstations
            par1Bin = phat_multiexp(i).(seasons{ijk})(1,:).';
            bins = site_bins(i).(seasons{ijk});
            [~, ~, loc] = histcounts(moran(:,i), bins);
            hasLoc = occ(:,i);
            loc(~hasLoc) = [];
            par1(hasLoc, i) = par1Bin(loc);% Will return nans if moran is nan
        end
        phatGen.par1 = par1;
                
%         for i=1:length_season(ijk)*years_sim
%             for j=1:nstations
%                 bins=eval(['site_bins(j).bins',char(seasons(ijk))]);
%                 for ii=1:size(bins,2)-1;
%                     if moran(i,j)>=bins(ii) & moran(i,j)<bins(ii+1)
%                         t1=ii;
%                         phat.alpha(i,j)=alpha(j).alpha(t1);
%                     end
%                 end
%             end
%         end
    elseif dist==2 % gamma distribution
        par1 = nan(nStepsGen,nstations);
        par2 = nan(nStepsGen,nstations);
        for i = 1:nstations
            par1Bin = phat_multiexp(i).(seasons{ijk})(1,:).';
            par2Bin = phat_multiexp(i).(seasons{ijk})(2,:).';
            bins = site_bins(i).(seasons{ijk});
            [~, ~, loc] = histcounts(moran(:,i), bins); % find moran in ranges (bins).
            hasLoc = occ(:,i);
            loc(~hasLoc) = [];
            par1(hasLoc, i) = par1Bin(loc);
            par2(hasLoc, i) = par2Bin(loc);
        end
        phatGen.par1 = par1;
        phatGen.par2 = par2;
        
%         for i=1:length_season(ijk)*years_sim
%             for j=1:nstations
%                 bins=eval(['site_bins(j).bins',char(seasons(ijk))]);
%                 for ii=1:size(bins,2)-1;
%                     if moran(i,j)>=bins(ii) & moran(i,j)<bins(ii+1)
%                         t1=ii;
%                         phat.alpha(i,j)=alpha(j).alpha(t1);
%                         phat.beta(i,j)=beta(j).beta(t1);
%                     end
%                 end
%             end
%         end
    end
    
    % generate random number
    random=zeros(nstations,nStepsGen);
    for i=1:nstations
        random(i,:)=randn(1,nStepsGen);
    end
    
    
    % automatic determination of correlation matrix of random number using
    % algorithm of eq.7
    
%     tttt=cputime;
    [M,K,precip]= spatial_iterate_amo_index(C,nstations,random,...
        nStepsGen,threshold,occ,phatGen,dist,Tolerance,maxiter);
%     time_in_second=cputime-tttt;
    
    % correlation of random numbers needed
    corr_amounts_rand_index(ijk).cor=M;
    corr_amounts_rand_index(ijk).season=seasons{ijk};
    % resulting correlation of amount
    corr_amounts_gen_index(ijk).cor=K;
    corr_amounts_gen_index(ijk).season=seasons{ijk};
    % generated precip amount
    precip_gen_index(ijk).pre=precip;
    precip_gen_index(ijk).season=seasons{ijk};
end

%
% produce graphics
%
if graph==1
    % correlation of generated precip amount
    plotheight = 48;
    plotwidth = 45;
    subplotsx = 2;
    subplotsy = 2;
    leftedge = 3;
    rightedge = 2;
    topedge = 2;
    bottomedge = 4;
    spacex = 3;
    spacey = 7;
    sub_pos = subplot_pos(plotwidth, plotheight, leftedge, rightedge, bottomedge,...
        topedge, subplotsx, subplotsy, spacex, spacey);
    sub_pos = fliplr(sub_pos);
    
    figure
    set(gcf, 'Position', [445, 222, 602, 385])
    set(gcf, 'Color', 'w')
    
    for ii=1:nseasons
        tf = ~logical(tril(true(size(seasonal_cor(ii).cor))));
        C=seasonal_cor(ii).cor;
        CC = C(tf);
        K=corr_amounts_gen_index(ii).cor;
        KK=K(tf);
        
        axes('Position', sub_pos{ii},...
            'XTick', [0, 0.5, 1],...
            'YTick', [0, 0.5, 1],...
            'XLim', [0, 1],...
            'YLim', [0, 1],...
            'Box', 'on',...
            'FontSize', 8);
        hold on
        
        plot(CC,KK,'o');
        plot([0, 1],[0, 1], 'Color', [0.9290, 0.6940, 0.1250])
        
        if ii == 3
            xlabel('\bf Observed correlation of amounts', 'HorizontalAlignment', 'left');
        end
       
        if ii == 1
            ylabel('\bf Generated correlation of amounts', 'HorizontalAlignment', 'right');
        end
        
        title(seasons{ii}, 'VerticalAlignment', 'baseline')

    end
%     figure_size = get(gcf,'position');
%     set(gcf,'PaperPosition', figure_size/100);
%     fileName = fullfile(folderOut, 'Corr-Amount.jpg');
%     print(gcf,'-djpeg', fileName, '-r500');
    
    %     % needed correlation for generating precip amount
    %     figure
    %     for ii=1:4
    %         C=seasonal_cor(1,ii).cor;
    %         CC=reshape(triu(C,1),1,[]);
    %         i=find(CC~=0);
    %         CC=CC(i);
    %
    %         M=corr_amounts_rand_index(1,ii).cor;
    %         MM=reshape(triu(M,1),1,[]);
    %         i=find(MM~=0);
    %         MM=MM(i);
    %
    %         subplot(2,2,ii)
    %         plot(CC,MM,'o',[0 1],[0 1]);
    %         xlabel('observed correlation');
    %         ylabel('random numbers correlation');
    %         axis([0 1 0 1])
    %         axis square
    %         set(gca,'FontSize',12)
    %         set(gcf,'Color',[1 1 1])
    %     end
end
end
