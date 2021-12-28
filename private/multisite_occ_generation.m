% generate spatially correlated precip occurrence with two steps: 
% (1) automatic determination of correlation matrix of random number using
% algorithm (eq.7).Diagonalization and repalcement of negative eigenvalues
% if matrix is non-positive definite (eqs.5 & 6)
% (2) generate precipitation occurrence using first-order Markov chain with
% Cholesky factorization (eq.2)

function [corr_occ_rand,corr_occ_gen,occurrences_gen]=multisite_occ_generation(...
    corr_occ,trans_prob,nstations,length_month,years_sim,months,Tolerance,...
    maxiter,graph)%, folderOut)

n=length_month*years_sim;       % number of days in the iterative process

% obtain correlation matrix of occurences
corr_occ_rand = struct('cor', cell(12,1), 'month', cell(12, 1));
corr_occ_gen = struct('cor', cell(12,1), 'month', cell(12, 1));
occurrences_gen = struct('occ', cell(12,1), 'month', cell(12, 1));

for q=1:12    % do for each month
    C=corr_occ(q).cor;
	transitions=trans_prob(q).prob;    
    
	% since random numbers generated have a normal distribution, each p00 and
	% p10 has to be recalculated according to a normal number	
	
% 	transitions_normal=zeros(size(transitions));
%     for i=1:nstations
%         transitions_normal(1,i)=-sqrt(2)*erfcinv(2*transitions(1,i));
%         transitions_normal(2,i)=-sqrt(2)*erfcinv(2*transitions(2,i));
%     end
    
    % produce independently normally ditributed random number
	random=zeros(nstations,n(q));
	for i=1:nstations
        random(i,:)=randn(1,n(q));
	end
	
	% automatic determination of correlation matrix of random number using
    % algorithm of eq.7
	
	%tttt=cputime;
	[M,K,occurrences]= spatial_iterate_occ(C,nstations,random,...
        transitions,n(q),Tolerance,maxiter);
	%time_in_second=cputime-tttt;
	
	%'correlation of random numbers needed'
%     eval(['corr_occ_rand_' char(months(q)) '=M;' ]);
% 	%'resulting correlation of occurence'
%     eval(['corr_occ_gen_' char(months(q)) '=K;' ]);
%     eval(['occurrences_' char(months(q)) '=occurrences;' ]);
      
	corr_occ_rand(q).cor=M;
    corr_occ_rand(q).month=months{q};  
    
    corr_occ_gen(q).cor=K;
    corr_occ_gen(q).month=months{q};   
    
    occurrences_gen(q).occ=occurrences;
    occurrences_gen(q).month=months{q};
end

%
% produce graphics
%
if graph==1
    % correlation of generated precip occurrence
    plotheight = 48;
    plotwidth = 45;
    subplotsx = 4;
    subplotsy = 3;
    leftedge = 3;
    rightedge = 2;
    topedge = 2;
    bottomedge = 4;
    spacex = 3;
    spacey = 4;
    sub_pos = subplot_pos(plotwidth, plotheight, leftedge, rightedge, bottomedge,...
        topedge, subplotsx, subplotsy, spacex, spacey);
    sub_pos = fliplr(sub_pos);
    
    figure
    set(gcf, 'Position', [445, 222, 602, 385])
    set(gcf, 'Color', 'w')
    for ii=1:12
        
        tf = ~logical(tril(true(size(corr_occ(ii).cor))));
        C=corr_occ(ii).cor;
        CC = C(tf);
        K=corr_occ_gen(ii).cor;
        KK= K(tf);
        
        axes('Position', sub_pos{ii},...
            'XTick', [0, 0.5, 1],...
            'YTick', [0, 0.5, 1],...
            'XLim', [0, 1],...
            'YLim', [0, 1],...[0.9290 0.6940 0.1250]
            'Box', 'on',...
            'FontSize', 8);
        hold on
        
        plot(CC,KK,'o');
        plot([0, 1],[0, 1], 'Color', [0.9290, 0.6940, 0.1250])
        
        if ii == 5
            ylabel('\bf Generated correlation of amounts');
        end
        
        if ii == 10
            xlabel('\bf Observed correlation of amounts', 'HorizontalAlignment', 'left');
        end
      
        
        title(months{ii}, 'VerticalAlignment', 'baseline')
        %set(gca,'FontSize',12)
    end
%     figure_size = get(gcf,'position');
%     set(gcf,'PaperPosition', figure_size/100);
%     fileName = fullfile(folderOut, 'Corr-Ocurrences.jpg');
%     print(gcf,'-djpeg', fileName, '-r500');


%     % needed correlation for generating precip occurrence
%     figure
%     for ii=1:12
%         C=corr_occ(1,ii).cor;
%         CC=reshape(triu(C,1),1,[]);
%         i=find(CC~=0);
%         CC=CC(i);
% 
%         M=corr_occ_rand(1,ii).cor;
%         MM=reshape(triu(M,1),1,[]);
%         i=find(MM~=0);
%         MM=MM(i);
% 
%         subplot(3,4,ii)
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



