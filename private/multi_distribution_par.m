%
% establish link between occurrence index (eq.10) and average precip
% amounts (and stdev when using gamma distribution) for each station and
% construct the multi-exponential (or multi-gamma) distribution for
% each station (eq.11)
%

function [phat_multiexp,site_bins,occ_index,SamSeason,SoccSeason,corrOcc, occurences_genSeason]=...
    multi_distribution_par(Sam, Socc, stationname,...
    occurrences_gen,nstations,threshold,dist,graph)%, folderOut)

%
seasons={'DJF' ;'MAM'; 'JJA'; 'SON';};
nSeasons = numel(seasons);
begin_season=[1 91 183 275];
length_season=[90 92 92 91];

% generate season precip occrrence (the spatial correlation is induced)
occurences_genSeason.DJF=[occurrences_gen(12).occ; occurrences_gen(1).occ; occurrences_gen(2).occ];
occurences_genSeason.MAM=[occurrences_gen(3).occ; occurrences_gen(4).occ; occurrences_gen(5).occ];
occurences_genSeason.JJA=[occurrences_gen(6).occ; occurrences_gen(7).occ; occurrences_gen(8).occ];
occurences_genSeason.SON=[occurrences_gen(9).occ; occurrences_gen(10).occ; occurrences_gen(11).occ];

% Weigths Corr
corrDJF = corrcoef(occurences_genSeason.DJF);
corrOcc.DJF=forcecorr(corrDJF);

corrMAM = corrcoef(occurences_genSeason.MAM);
corrOcc.MAM=forcecorr(corrMAM);

corrJJA = corrcoef(occurences_genSeason.JJA);
corrOcc.JJA=forcecorr(corrJJA);

corrSON=corrcoef(occurences_genSeason.SON);
corrOcc.SON=forcecorr(corrSON);


% [n, stations] by season arrays of amount and occurences.
SamSeason = [];
SoccSeason = [];
for i = 1:nSeasons
    cam = cell(1, nstations);
    cocc = cell(1, nstations);
    for j = 1:nstations
        m = Sam.(stationname{j});
        o = Socc.(stationname{j});
        
        m = [m(:,335:365), m(:,1:334)];
        o = [o(:,335:365), m(:,1:334)];
        
        m = m(:,begin_season(i):begin_season(i)+length_season(i)-1);
        o = o(:,begin_season(i):begin_season(i)+length_season(i)-1);
        
        m = m(:);
        o = o(:);
        cam{j} = m;
        cocc{j} = o;
    end
    
    cam = horzcat(cam{:});
    cocc = horzcat(cocc{:});
      
    SamSeason.(seasons{i}) = cam;
    SoccSeason.(seasons{i}) = cocc;
end


% Determine ranges of bins for k index for each station and season with obs
% data.
phat_multiexp = struct('DJF', cell(nstations, 1), 'MAM', cell(nstations, 1),...
    'JJA', cell(nstations, 1), 'SON', cell(nstations, 1));
site_bins = struct('DJF', cell(nstations, 1), 'MAM', cell(nstations, 1),...
    'JJA', cell(nstations, 1), 'SON', cell(nstations, 1));

for ijk=1:nstations % loop for each station
    for qq=1:4   % loop for each season
        
        % Compute k index for all days with precip (morancum).
        corSeason = corrOcc.(seasons{qq});
        corSeason = corSeason(ijk,:);
        corSeason(ijk) = [];
        
        precipdata = SamSeason.(seasons{qq});
        occdata = SoccSeason.(seasons{qq});
        isOcc = occdata == 1;
        precipdata(isOcc) = precipdata(isOcc) - threshold; % for fited function (mean and std), recover threshold when applied inverse func.
               
        c = size(precipdata, 1);
        moran=nan(c,1);
        for i=1:c   % for each day
            if precipdata(i,ijk)>0 % Outside loop!
                thisPrecip = precipdata(i,:);
                hasPrecip = thisPrecip > 0;
                thisOcc = zeros(1, nstations);
                thisOcc(hasPrecip) = 1;
                thisOcc(ijk) = [];
                moran(i)=thisOcc*corSeason'/sum(corSeason);
            end
        end
        isDry=isnan(moran);
        hasEnoughEvents = sum(~isDry) > 10;

        if hasEnoughEvents
            raincum=precipdata(~isDry,ijk); % For dry months may be empty
            morancum=moran(~isDry);
        else % Not enough ocurrences, will generate smalls ammount of pr.
            upperLim = 5;%generate amounts between th and 5 mm;
            nEvents = 15;
            raincum = (threshold - upperLim)*rand(nEvents, 1) + threshold; 
            morancum = rand(nEvents, 1);%rand k values. 
        end
        
        % Find the ranges of k bins considering a min of 50 elements.
        nRanges = 10;
        nMin = 30;
        nMinRanges = 3;
        bins = 0:1/nRanges:1;
        bins(end) = 1.001; %
        nElem = histcounts(morancum, bins);
        while any(nElem < nMin) && nRanges > nMinRanges %Minimun 3 ranges.
            bins = 0:1/nRanges:1;
            bins(end) = 1.001; % Due histcount behaviour.
            nElem = histcounts(morancum, bins);
            nRanges = nRanges - 1;
        end
        site_bins(ijk).(seasons{qq})=bins; % bins for each season
        
        
        % Compute parameters for each range.
        qqq=length(bins);
        bincenter=(bins(2:qqq)+bins(1:qqq-1))*0.5;
        
        if dist==1 % exponential distribution
            binMean=zeros(size(bincenter));
            for i=1:qqq-1
                j=morancum>=bins(i) & morancum<bins(i+1);
                binMean(i)=mean(raincum(j)); % calculate the mean precip for each bin
            end
            Bm.(seasons{qq}) = binMean;
        elseif dist==2 % gamma distribution
            binMean=zeros(size(bincenter));
            binStd=zeros(size(bincenter));
            for i=1:qqq-1
                j=morancum>=bins(i) & morancum<bins(i+1);
                binMean(i)=mean(raincum(j));
                binStd(i)=std(raincum(j));
            end
            Bm.(seasons{qq}) = binMean;
            Bs.(seasons{qq}) = binStd;
            
        end
        Bcent.(seasons{qq})=bincenter;
        
        
        % best-line fit and normalize data so that mean of multi-exponential
        % distribution fits that of the observed one
        
        % using a function of y=a+bX^n to fit the relationship between
        % occurrence index and the coresponging mean precip, n is integer
        % value between 2 and 15, which is optimized by minimizing the
        % relative difference between the real and fucntion-generated mean
        % precip
        orderExp=1:20;
        dif_tot=zeros(size(orderExp));
        for or=orderExp
            [P]=funct(bincenter,binMean,or); % estimate the parameters of a and b
            yMean=P(1)+P(2)*bincenter.^or; % apply the equation
            dif=abs((yMean-binMean)./binMean); % calculate the relative difference
            dif_tot(1,or)=sum(dif);
        end
        [~, indMin] = min(dif_tot);
        or_min=orderExp(indMin); % find the best order
        [P]=funct(bincenter,binMean,or_min); % estimate the parameters of a and b again
        yMean=P(1)+P(2)*bincenter.^or_min; % apply the equation again
        Ym.(seasons{qq})= yMean;
        
        if dist==2 % gamma distribution
            orderExp=1:20;
            dif_tot=zeros(size(orderExp));
            for or=orderExp
                [P]=funct(bincenter,binStd,or); % estimate the parameters of a and b
                yStd=P(1)+P(2)*bincenter.^or; % apply the equation
                dif=abs((yStd-binStd)./binStd); % calculate the relative difference
                dif_tot(1,or)=sum(dif);
            end
            [~, indMin] = min(dif_tot);
            or_min=orderExp(indMin); % find the best order
            [P]=funct(bincenter,binStd,or_min); % estimate the parameters of a and b again
            yStd=P(1)+P(2)*bincenter.^or_min; % apply the equation again
            Ys.(seasons{qq}) = yStd;
        end
        
        
        % Adjust means returned from curve fitted(Y)
        % use nbins from simulated occurences instead of observed ones,
        % this will insure precipitation means that fit that of observed
        
        if hasEnoughEvents
            occ = occurences_genSeason.(seasons{qq});
            hasOcc = occ(:,ijk)==1;
            occ = occ(hasOcc,:);
            occ(:,ijk) = []; % delete current stn.
            
            genLength = size(occ, 1);
            moran = nan(genLength, 1);
            for g=1:genLength
                moran(g)=(occ(g,:)*corSeason')/sum(corSeason);
            end
            
            nbins = nan(1, qqq-1);
            for i=1:qqq-1
                j=find(moran>=bins(i) & moran<bins(i+1));
                nbins(i)=length(j);
            end
            ntotal=sum(nbins);
            pbins=nbins/ntotal;
            
            tf=precipdata(:,ijk)>0; % Global mean of season
            raindata=precipdata(tf,ijk);
            correction_factor=mean(raindata)/sum(pbins.*yMean);  % make sure the precipitation mean stays the same
        else
            correction_factor = 1;
        end

        mean_value=yMean*correction_factor; % bias correction for fitted mean precip
        
        
        % transfer mean and stdev to exponential and gamma distribution
        % parameters
        phat_value = nan(1, qqq-1);
        if dist==1 % exponential distribution
            phat_value= 1./mean_value; % Lambda
        elseif dist==2
            phat_value(1,:)=mean_value.^2./yStd; % a
            phat_value(2,:)=yStd/mean_value; % b
        end
        phat_multiexp(ijk).(seasons{qq})=phat_value;
        
        %         if dist==2 % gamma distibution
        %             Pstd.(seasons{qq})=yStd; % From fited curve
        %             Pstd2.(seasons{qq})=Pstdev2; % From obs data in bins
        %         end
        
    end
    
    
    % Save occ_index?
    occ_index.(stationname{ijk}).Bm=Bm; % bin mean
    occ_index.(stationname{ijk}).Ym=Ym; % bin mean adjusted
    occ_index.(stationname{ijk}).Bcent=Bcent; % bin center
    
    if dist==2
        occ_index.(stationname{ijk}).Ys=Ys; 
        occ_index.(stationname{ijk}).Bs=Bs; % bin std
    end
    %
    % end of main computing, beginning of plotting routine
    %
    if graph==1

        plotheight = 48;
        plotwidth = 60;
        leftedge = 3;
        rightedge = 2;
        topedge = 2;
        bottomedge = 4;
        spacex = 4;
        spacey = 7;
        if dist==1 % exponential distibution (only plot mean precip)
            subplotsx = 2;
            subplotsy = 2;
        elseif dist==2 % gamma distribution (plot stdev also)
            subplotsx = 4;
            subplotsy = 2;
        end
        sub_pos = subplot_pos(plotwidth, plotheight, leftedge, rightedge, bottomedge,...
            topedge, subplotsx, subplotsy, spacex, spacey);
        sub_pos = fliplr(sub_pos);
    
        maxvalue=max([max(Ym.DJF),max(Bm.DJF),max(Ym.MAM),...
            max(Bm.MAM),max(Ym.JJA),max(Bm.JJA),...
            max(Ym.SON),max(Bm.SON)]);
        
        figure
        set(gcf, 'Position', [445, 222, 602, 385])
        set(gcf, 'Color', 'w')
        
        % DJF
        axes('Position', sub_pos{1},...
            'XTick', [0 0.25 0.5 0.75 1],...
            'XLim', [0, 1],...
            'YLim', [0, maxvalue+1],...
            'Box', 'on',...
            'FontSize', 8);
        hold on
        title('DJF', 'VerticalAlignment', 'baseline')
        plot(Bcent.DJF,Ym.DJF,Bcent.DJF,Bm.DJF,'o')
        ylabel('\bf Mean precip(mm)')
        xlabel('\bf Ocurrence index')
        grid on
        
        % MAM
        axes('Position', sub_pos{2},...
            'XTick', [0 0.25 0.5 0.75 1],...
            'XLim', [0, 1],...
            'YLim', [0, maxvalue+1],...
            'Box', 'on',...
            'FontSize', 8);
        hold on
        title('MAM','VerticalAlignment', 'baseline')
        plot(Bcent.MAM,Ym.MAM,Bcent.MAM,Bm.MAM,'o')
        ylabel('\bf Mean precip(mm)')
        xlabel('\bf Ocurrence index')
        grid on

        % JJA
        axes('Position', sub_pos{3},...
            'XTick', [0 0.25 0.5 0.75 1],...
            'XLim', [0, 1],...
            'YLim', [0, maxvalue+1],...
            'Box', 'on',...
            'FontSize', 8);
        hold on
        title('JJA','VerticalAlignment', 'baseline')
        plot(Bcent.JJA,Ym.JJA,Bcent.JJA,Bm.JJA,'o')
        ylabel('\bf Mean precip(mm)')
        xlabel('\bf Ocurrence index')
        grid on
        
        % SON
        axes('Position', sub_pos{4},...
            'XTick', [0 0.25 0.5 0.75 1],...
            'XLim', [0, 1],...
            'YLim', [0, maxvalue+1],...
            'Box', 'on',...
            'FontSize', 8);
        hold on
        title('SON', 'VerticalAlignment', 'baseline')
        plot(Bcent.SON,Ym.SON,Bcent.SON,Bm.SON,'o')
        ylabel('\bf Mean precip(mm)')
        xlabel('\bf Ocurrence index')
        grid on
        
        if dist==2
            maxvalue=max([max(Ys.DJF),max(Bs.DJF),max(Ys.MAM),...
                max(Bs.MAM),max(Ys.JJA),max(Bs.JJA),...
                max(Ys.SON),max(Bs.SON)]);
            
            % DJF
            axes('Position', sub_pos{5},...
                'XTick', [0 0.25 0.5 0.75 1],...
                'XLim', [0, 1],...
                'YLim', [0, maxvalue+1],...
                'Box', 'on',...
                'FontSize', 8);
            hold on
            title('DJF', 'VerticalAlignment', 'baseline')
            plot(Bcent.DJF,Ys.DJF,Bcent.DJF,Bs.DJF,'o')
            ylabel('\bf Stdev of precip(mm)')
            xlabel('\bf Ocurrence index')
            grid on
                        
            % MAM
            axes('Position', sub_pos{6},...
                'XTick', [0 0.25 0.5 0.75 1],...
                'XLim', [0, 1],...
                'YLim', [0, maxvalue+1],...
                'Box', 'on',...
                'FontSize', 8);
            hold on
            title('MAM','VerticalAlignment', 'baseline')
            plot(Bcent.MAM,Ys.MAM,Bcent.MAM,Bs.MAM,'o')
            ylabel('\bf Stdev of precip(mm)')
            xlabel('\bf Ocurrence index')
            grid on
            
            % JJA
            axes('Position', sub_pos{7},...
                'XTick', [0 0.25 0.5 0.75 1],...
                'XLim', [0, 1],...
                'YLim', [0, maxvalue+1],...
                'Box', 'on',...
                'FontSize', 8);
            hold on
            title('JJA','VerticalAlignment', 'baseline')
            plot(Bcent.JJA,Ys.JJA,Bcent.JJA,Bs.JJA,'o')
            ylabel('\bf Stdev of precip(mm)')
            xlabel('\bf Ocurrence index')
            grid on
            
            % SON
            axes('Position', sub_pos{8},...
                'XTick', [0 0.25 0.5 0.75 1],...
                'XLim', [0, 1],...
                'YLim', [0, maxvalue+1],...
                'Box', 'on',...
                'FontSize', 8);
            hold on
            title('SON', 'VerticalAlignment', 'baseline')
            plot(Bcent.SON,Ys.SON,Bcent.SON,Bs.SON,'o')
            ylabel('\bf Stdev of precip(mm)')
            xlabel('\bf Ocurrence index')
            grid on
        end
%         figure_size = get(gcf,'position');
%         set(gcf,'PaperPosition', figure_size/100);
%         fileName = fullfile(folderOut, sprintf('Fit-K-Index_%s.jpg',stationname{ijk}));
%         print(gcf,'-djpeg', fileName, '-r500');
    end
    
    
%     eval(['occ_index.S',num2str(ijk),'.PmeanDJF=PmeanDJF;']);
%     eval(['occ_index.S',num2str(ijk),'.PmeanMAM=PmeanMAM;']);
%     eval(['occ_index.S',num2str(ijk),'.PmeanJJA=PmeanJJA;']);
%     eval(['occ_index.S',num2str(ijk),'.PmeanSON=PmeanSON;']);
%     
%     eval(['occ_index.S',num2str(ijk),'.meanfitDJF=meanfitDJF;']);
%     eval(['occ_index.S',num2str(ijk),'.meanfitMAM=meanfitMAM;']);
%     eval(['occ_index.S',num2str(ijk),'.meanfitJJA=meanfitJJA;']);
%     eval(['occ_index.S',num2str(ijk),'.meanfitSON=meanfitSON;']);
%     
%     eval(['occ_index.S',num2str(ijk),'.bincenterDJF=bincenterDJF;']);
%     eval(['occ_index.S',num2str(ijk),'.bincenterMAM=bincenterMAM;']);
%     eval(['occ_index.S',num2str(ijk),'.bincenterJJA=bincenterJJA;']);
%     eval(['occ_index.S',num2str(ijk),'.bincenterSON=bincenterSON;']);
%     
%     if dist==2
%         eval(['occ_index.S',num2str(ijk),'.PstdevDJF=PstdevDJF;']);
%         eval(['occ_index.S',num2str(ijk),'.PstdevMAM=PstdevMAM;']);
%         eval(['occ_index.S',num2str(ijk),'.PstdevJJA=PstdevJJA;']);
%         eval(['occ_index.S',num2str(ijk),'.PstdevSON=PstdevSON;']);
%         
%         eval(['occ_index.S',num2str(ijk),'.Pstdev2DJF=Pstdev2DJF;']);
%         eval(['occ_index.S',num2str(ijk),'.Pstdev2MAM=Pstdev2MAM;']);
%         eval(['occ_index.S',num2str(ijk),'.Pstdev2JJA=Pstdev2JJA;']);
%         eval(['occ_index.S',num2str(ijk),'.Pstdev2SON=Pstdev2SON;']);
%     end
    
end
end



