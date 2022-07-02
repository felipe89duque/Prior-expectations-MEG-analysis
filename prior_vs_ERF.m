clear variables; clc, close all;
% grandavg_timelock_planar_prior = load('grandaveragetimelock_planar_prior.mat');

load('grandavg_timelock_planar_prior_distributed.mat');


%% Parameters

% Pick 'faces' or 'houses'
stimulus = 'faces';

% time limits/latency for cluster analysis
latency = [-0.2, 0.7];

% number of plots for topoplot figure (number of bins to devide time for
% clusters visualisation
num_plots = 1;
gif = 0;


%% 
if strcmp(stimulus,'faces')
   prior_1 = prior_1_faces;
   prior_2 = prior_2_faces;
   prior_3 = prior_3_faces;
   
elseif strcmp(stimulus,'houses')
   prior_1 = prior_1_houses;
   prior_2 = prior_2_houses;
   prior_3 = prior_3_houses;
end

for pp=1:24
    prior_3{pp} = rmfield(prior_2{pp},'ntrials');
    prior_2{pp} = rmfield(prior_2{pp},'ntrials');
    prior_1{pp} = rmfield(prior_1{pp},'ntrials');
end
for pp=1:24
    prior_3{pp} = rmfield(prior_2{pp},'sampleinfo');
    prior_2{pp} = rmfield(prior_2{pp},'sampleinfo');
    prior_1{pp} = rmfield(prior_1{pp},'sampleinfo');
end
for pp=1:24
    prior_3{pp} = rmfield(prior_2{pp},'fsample');
    prior_2{pp} = rmfield(prior_2{pp},'fsample');
    prior_1{pp} = rmfield(prior_1{pp},'fsample');
end

%% calculate grand average for each condition

% make a plot
% Then take the difference of the averages using ft_math
cfg           = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
% raweffectFICvsFC = ft_math(cfg, grandavg_high_coher, grandavg_low_coher);
% 
high_minus_low = cell(size(prior_1));
for i = 1:length(prior_1)
    high_minus_low{i} = ft_math(cfg, prior_3{i}, prior_1{i});
%     high_minus_low_houses{i} = ft_math(cfg, high_coher_houses{i}, low_coher_houses{i});
end

% calculate grand average for each condition
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
grandavg_high_minus_low = ft_timelockgrandaverage(cfg,high_minus_low{:});

%% Permutation test based on cluster statistics

cfg = [];
cfg.method      = 'triangulation';
cfg.layout      = 'CTF270Nott.lay';                % specify layout of channels
cfg.feedback    = 'no';                              % show a neighbour plot
neighbours      = ft_prepare_neighbours(cfg, prior_1{1}); % define neighbouring channels

cfg = [];
cfg.channel     = 'MEG';
cfg.neighbours  = neighbours; % defined as above
cfg.latency     = latency;
cfg.avgovertime = 'no';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.computecritval = 'yes';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.01;
cfg.tail = 0;
cfg.correctm    = 'cluster';
cfg.numrandomization = 1000;%'all';  

% difference between low and high coherence
Nsub = length(prior_1);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number


stat = ft_timelockstatistics(cfg,prior_3{:},prior_1{:});
%%
% % Make a vector of all p-values associated with the clusters from ft_timelockstatistics.
pos_cluster_pvals = [stat.posclusters(:).prob];
% 
% % Then, find which clusters are deemed interesting to visualize, here we use a cutoff criterion based on the
% % cluster-associated p-value, and take a 5% two-sided cutoff (i.e. 0.025 for the positive and negative clusters,
% % respectively
pos_clust = find(pos_cluster_pvals < cfg.alpha);
pos       = ismember(stat.posclusterslabelmat, pos_clust);
% 
% cluster1 = pos | neg;
% cluster2 = stat.prob <= cfg.alpha;
% 
% % and now for the negative clusters...
neg_cluster_pvals = [stat.negclusters(:).prob];
neg_clust         = find(neg_cluster_pvals < cfg.alpha);
neg               = ismember(stat.negclusterslabelmat, neg_clust);

mask_pos = stat.stat.*pos;
mask_neg = stat.stat.*neg;

factors = unique(factor(num_plots));
rows = max(factors(end),num_plots/factors(end));
cols = num_plots/rows;
time_divisions = linspace(cfg.latency(1), cfg.latency(2), num_plots +1);
cont=1;


for i = time_divisions
    if i ~= cfg.latency(2) 
        time_index = find(stat.time>=i,1);
        t_steps = floor(length(stat.time)/num_plots);
        time = time_index:time_index+t_steps-1;
        
         pos_masked = squeeze(nanmean(mask_pos(:,time),2));
         pos_masked(find(pos_masked == 0)) = NaN;
         
         neg_masked = squeeze(nanmean(mask_neg(:,time),2));
         neg_masked(find(neg_masked == 0)) = NaN;
        if ~gif
            figure(1); %pos
            subplot(cols,rows,cont);
            topoplot(pos_masked,'layout','CTF270Nott.lay','maplimits',[-2.576,2.576],'style','straight');
            %colormap(flipud(pink))
            title(['time = [',num2str(stat.time(time(1))),', ',num2str(stat.time(time(end))),']']);

            figure(2); %neg
            subplot(cols,rows,cont);
            topoplot(neg_masked,'layout','CTF270Nott.lay','maplimits',[-2.576,2.576],'style','straight');
            %colormap(pink)
            title(['time = [',num2str(stat.time(time(1))),', ',num2str(stat.time(time(end))),']']);
        else
           figure(3); topoplot(squeeze(nanmean(mask_pos(:,time),2)),'layout','CTF270Nott.lay','maplimits',[-2.576,2.576],'style','straight');
           %colormap(flipud(pink))
           title(['time: ',newline,num2str(stat.time(time(1)))]);
           saveas(gcf,strcat('gif_pos_clust_',stimulus,'_',num2str(cont),'_prior.png'));
        end
        
        cont = cont +1;
    end
%     figure;topoplot(squeeze(nanmean(mask_neg(:,time),2)),'layout','CTF270Nott.lay');
end

pos_stat = stat.stat.*mask_pos;
neg_stat = stat.stat.*mask_neg;
% figure;plot(stat.time,nanmean(pos_stat,1));title('mean t-statistic over positive cluster');
% figure;plot(stat.time,nanmean(neg_stat,1));title('mean t-statistic over negative cluster');

%%
% number of rectangles to plot, showing the clusters with statistical
% significance
step_pos_stat = nanmean(pos_stat,1) > 0;
pos_clust_lims = diff(step_pos_stat);
pos_clust_start = stat.time(find(pos_clust_lims > 0));
pos_clust_end = stat.time(find(pos_clust_lims < 0));

if isempty(pos_clust_start) & (pos_clust~=0)
    pos_clust_start = latency(1);
end
if isempty(pos_clust_end) & (pos_clust~=0)
    pos_clust_end = latency(2);
end

% check the order of the clusters starts and ends to make them into pairs
diff_start_end = length(pos_clust_start) - length(pos_clust_end);

if abs(diff_start_end) > 1
    error('There is an inequivalent number of clusters starting and ending');
elseif diff_start_end == 1
    pos_clust_end = [pos_clust_end, latency(2)];
elseif diff_start_end == -1
    pos_clust_start = [ latency(1), pos_clust_start];
else % same length
    if ~isempty(pos_clust_start)
        if pos_clust_end(1) < pos_clust_start(1)
            pos_clust_start = [latency(1), pos_clust_start];
        end
        if pos_clust_start(end) > pos_clust_end(end)
            pos_clust_end = [pos_clust_end, latency(2)];
        end
    end
end
num_pos_clusts = length(pos_clust_start);

% negative clusters
step_neg_stat = nanmean(neg_stat,1) > 0;
neg_clust_lims = diff(step_neg_stat);
neg_clust_start = stat.time(find(neg_clust_lims > 0));
neg_clust_end = stat.time(find(neg_clust_lims < 0));

if isempty(neg_clust_start) & (neg_clust~=0)
    neg_clust_start = latency(1);
end
if isempty(neg_clust_end) & (neg_clust~=0)
    neg_clust_end = latency(2);
end

% check the order of the clusters starts and ends to make them into pairs
diff_start_end = length(neg_clust_start) - length(neg_clust_end);

if abs(diff_start_end) > 1
    error('There is an inequivalent number of clusters starting and ending');
elseif diff_start_end == 1
    neg_clust_end = [neg_clust_end, latency(2)];
elseif diff_start_end == -1
    neg_clust_start = [ latency(1), neg_clust_start];
else % same length
    if ~isempty(neg_clust_start)
        if neg_clust_end(1) < neg_clust_start(1)
            neg_clust_start = [latency(1), neg_clust_start];
        end
        if neg_clust_start(end) > neg_clust_end(end)
            neg_clust_end = [neg_clust_end, latency(2)];
        end
    end
end
num_neg_clusts = length(neg_clust_start);

 %%   
    
        


grandavg_colapsed_sensors = nanmean(grandavg_high_minus_low.avg,1);
y_lims = [min(grandavg_colapsed_sensors), max(grandavg_colapsed_sensors)];
figure;hold on;
% make sure of how many rectangles there should be (number of clusters)
for i = 1:length(pos_clust_start)
    rectangle('Position',[pos_clust_start(i) y_lims(1) (pos_clust_end(i)-pos_clust_start(i)) diff(y_lims)],'FaceColor',[0.7 0.7 0.7]);
end
for i = 1:length(neg_clust_start)
    rectangle('Position',[neg_clust_start(i) y_lims(1) (neg_clust_end(i)-neg_clust_start(i)) diff(y_lims)],'FaceColor',[0.9 0.7 0.7]);
end
plot(grandavg_high_minus_low.time,grandavg_colapsed_sensors);
title(['High prior minus low prior ',stimulus])
% min_x = stat.time(find(nanmean(pos_stat,1) > 0,1));
% max_x = stat.time(find(nanmean(pos_stat,1) > 0,1,'last'));
% plot([min_x,min_x],y_lims,'k--');
% plot([max_x,max_x],y_lims,'k--');
set(gca,'FontSize',16);set(gcf,'Color','w');



% cfg = [];
% %cfg.style     = 'blank';
% cfg.layout = 'CTF270Nott.lay';
% % cfg.xlim = [-0.2:0.1:0.6];
% cfg.highlight = 'on';
% cfg.maskparameter = mask_pos;
% cfg.highlightchannel = find(stat.mask);
% cfg.comment   = 'no';
% figure; ft_topoplotER(cfg, grandavg_high_minus_low_faces)
% title('Nonparametric: significant with cluster-based multiple comparison correction')
% 
% %% 
% figure;plot([1,2],[mean(grandavg_low_coher.avg,'all'),mean(grandavg_high_coher.avg,'all')]);
