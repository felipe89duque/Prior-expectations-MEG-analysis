clear variables; clc, close all;
% grandavg_timelock_planar_prior = load('grandaveragetimelock_planar_prior.mat');

load('GAtimelock_planar_stimcoherv2_distributed.mat');



%% Parameters

% Pick 'faces' or 'houses'
stimulus = 'houses';

% time limits/latency for cluster analysis
latency = [-0.2, 0.7];

% number of plots for topoplot figure (number of bins to devide time for
% clusters visualisation
num_plots = 8;
gif = 0;


%% 
if strcmp(stimulus,'faces')
   high_coher = high_coher_faces;
   low_coher = low_coher_faces;
   
elseif strcmp(stimulus,'houses')
   high_coher = high_coher_houses;
   low_coher = low_coher_houses;
end

for pp=1:24
    low_coher{pp} = rmfield(low_coher{pp},'ntrials');
    high_coher{pp} = rmfield(high_coher{pp},'ntrials');
end
for pp=1:24
    low_coher{pp} = rmfield(low_coher{pp},'sampleinfo');
    high_coher{pp} = rmfield(high_coher{pp},'sampleinfo');
end
for pp=1:24
    low_coher{pp} = rmfield(low_coher{pp},'fsample');
    high_coher{pp} = rmfield(high_coher{pp},'fsample');
end

%% calculate grand average for each condition

% make a plot
% Then take the difference of the averages using ft_math
cfg           = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
% raweffectFICvsFC = ft_math(cfg, grandavg_high_coher, grandavg_low_coher);
% 
high_minus_low = cell(size(high_coher));
for i = 1:length(high_coher)
    high_minus_low{i} = ft_math(cfg, high_coher{i}, low_coher{i});
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
neighbours      = ft_prepare_neighbours(cfg, high_coher{1}); % define neighbouring channels

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
Nsub = length(high_coher);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number


stat = ft_timelockstatistics(cfg,high_coher{:},low_coher{:});
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
mask_neg(find(neg == 0)) = NaN;

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
         
        if ~gif
            figure(1); %pos
            subplot(cols,rows,cont);
           
            try
                topoplot(pos_masked,'layout','CTF270Nott.lay','maplimits',[-2.576,2.576],'style','straight');
            catch
                %topoplot(zeros(size(mask_pos,1),1),'layout','CTF270Nott.lay');
            end
            %colormap(default);
            title(['time = [',num2str(stat.time(time(1))),', ',num2str(stat.time(time(end))),']']);

%             figure(2); %neg
%             subplot(cols,rows,cont);
%             try
%                 topoplot(squeeze(nanmean(mask_neg(:,time),2)),'layout','CTF270Nott.lay');%,'maplimits',[2.576/t_steps,2.576]);
%             
%             catch
%             end
%             %colormap(flipud(pink))
%             title(['time = [',num2str(stat.time(time(1))),', ',num2str(stat.time(time(end))),']']);
        else
           figure(3); topoplot(pos_masked,'layout','CTF270Nott.lay','maplimits',[-2.576,2.576],'style','straight');
           %colormap(flipud(pink))
           title(['time:',newline,num2str(stat.time(time(1)))]);
           saveas(gcf,strcat('gif_pos_clust_',stimulus,'_',num2str(cont),'.png'));
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

if isempty(pos_clust_start) && (pos_clust~=0)
    pos_clust_start = latency(1);
end
if isempty(pos_clust_end) && (pos_clust~=0)
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

step_neg_stat = nanmean(neg_stat,1) > 0;
 %%   plot ERFs
        

% Grand average difference
grandavg_colapsed_sensors = nanmean(grandavg_high_minus_low.avg,1);
y_lims = [min(grandavg_colapsed_sensors), max(grandavg_colapsed_sensors)];
figure;hold on;

for i = 1:length(pos_clust_start)
    rectangle('Position',[pos_clust_start(i) y_lims(1) (pos_clust_end(i)-pos_clust_start(i)) diff(y_lims)],'FaceColor',[0.7 0.7 0.7]);
end

plot(grandavg_high_minus_low.time,grandavg_colapsed_sensors);
ylim tight;
title('High coherence minus low coherence')
set(gca,'FontSize',16);set(gcf,'Color','w');

% plot each subject
figure;hold on;
for i = 1:length(high_coher)
    plot(high_minus_low{i}.time,nanmean(high_minus_low{i}.avg,1));
end
ylim tight;
title('High coherence minus low coherence')
set(gca,'FontSize',16);set(gcf,'Color','w');

% Plot ERF for each sensor
% calculate grand average for each condition
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
grandavg_high_coher = ft_timelockgrandaverage(cfg,high_coher{:});
grandavg_low_coher = ft_timelockgrandaverage(cfg,low_coher{:});

cfg = [];
cfg.layout = 'CTF270Nott.lay';
cfg.layout = ft_prepare_layout(cfg);
figure; ft_multiplotER(cfg,grandavg_high_coher,grandavg_low_coher);
figure; ft_multiplotER(cfg,grandavg_high_minus_low);

% select the individual subject data from the time points and calculate the mean
timesel = find(grandavg_high_coher.time >= latency(1) & grandavg_high_coher.time <= latency(2));
for isub = 1:length(high_coher)
    avg_high_coher(isub) = mean(high_coher{isub}.avg(:,timesel),'all');
    avg_low_coher(isub)  = mean(low_coher{isub}.avg(:,timesel),'all');
end

% plot to see the effect in each subject
M = [avg_high_coher(:) avg_low_coher(:)];
figure; plot(M', 'o-'); xlim([0.5 2.5])
legend({'subj1', 'subj2', 'subj3', 'subj4', 'subj5', 'subj6', ...
        'subj7', 'subj8', 'subj9', 'subj10'}, 'location', 'EastOutside');
% Same but with the grand average among subjects
mean_high_coher = mean(grandavg_high_coher.avg(:,timesel),'all');
mean_low_coher  = mean(grandavg_low_coher.avg(:,timesel),'all');
grand_M = [mean_high_coher mean_low_coher];
figure; plot(grand_M', 'o-'); xlim([0.5 2.5])

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
