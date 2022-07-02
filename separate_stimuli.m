clear;clc;

load('GAtimelock_planar_stimcoherv2.mat');

[num_subjects, num_conditions, num_sensors, time] = size(timelockprestim.avg);

prior_1_faces = repmat({timelockprestim},1,num_subjects);
prior_2_faces = repmat({timelockprestim},1,num_subjects);
prior_3_faces = repmat({timelockprestim},1,num_subjects);

prior_1_houses = repmat({timelockprestim},1,num_subjects);
prior_2_houses = repmat({timelockprestim},1,num_subjects);
prior_3_houses = repmat({timelockprestim},1,num_subjects);

for i = 1:num_subjects
   prior_1_faces{i}.avg = squeeze(timelockprestim.avg(i,1,:,:));
   prior_1_faces{i}.ntrials = timelockprestim.ntrials(i,1);
   
   prior_2_faces{i}.avg = squeeze(timelockprestim.avg(i,2,:,:));
   prior_2_faces{i}.ntrials = timelockprestim.ntrials(i,2);

   prior_3_faces{i}.avg = squeeze(timelockprestim.avg(i,3,:,:)); 
   prior_3_faces{i}.ntrials = timelockprestim.ntrials(i,3);
   
   prior_1_houses{i}.avg = squeeze(timelockprestim.avg(i,4,:,:)); 
   prior_1_houses{i}.ntrials = timelockprestim.ntrials(i,4);
   
   prior_2_houses{i}.avg = squeeze(timelockprestim.avg(i,5,:,:)); 
   prior_2_houses{i}.ntrials = timelockprestim.ntrials(i,5);
   
   prior_3_houses{i}.avg = squeeze(timelockprestim.avg(i,6,:,:)); 
   prior_3_houses{i}.ntrials = timelockprestim.ntrials(i,6);
end

save('GAtimelock_planar_stimcoherv2_distributed.mat',...
    'low_coher_faces','mid_coher_faces','high_coher_faces',...
    'low_coher_houses','mid_coher_houses','high_coher_houses');

%%
load('grandaveragetimelock_planar_prior.mat');

[num_subjects, num_conditions, num_sensors, time] = size(timelockprecuebl.avg);

prior_1_faces = repmat({timelockprestim},1,num_subjects);
prior_2_faces = repmat({timelockprestim},1,num_subjects);
prior_3_faces = repmat({timelockprestim},1,num_subjects);

prior_1_houses = repmat({timelockprestim},1,num_subjects);
prior_2_houses = repmat({timelockprestim},1,num_subjects);
prior_3_houses = repmat({timelockprestim},1,num_subjects);

for i = 1:num_subjects
   prior_1_faces{i}.avg = squeeze(timelockprestim.avg(i,1,:,:));
   prior_1_faces{i}.ntrials = timelockprestim.ntrials(i,1);
   
   prior_2_faces{i}.avg = squeeze(timelockprestim.avg(i,2,:,:));
   prior_2_faces{i}.ntrials = timelockprestim.ntrials(i,2);

   prior_3_faces{i}.avg = squeeze(timelockprestim.avg(i,3,:,:)); 
   prior_3_faces{i}.ntrials = timelockprestim.ntrials(i,3);
   
   prior_1_houses{i}.avg = squeeze(timelockprestim.avg(i,4,:,:)); 
   prior_1_houses{i}.ntrials = timelockprestim.ntrials(i,4);
   
   prior_2_houses{i}.avg = squeeze(timelockprestim.avg(i,5,:,:)); 
   prior_2_houses{i}.ntrials = timelockprestim.ntrials(i,5);
   
   prior_3_houses{i}.avg = squeeze(timelockprestim.avg(i,6,:,:)); 
   prior_3_houses{i}.ntrials = timelockprestim.ntrials(i,6);
end

save('grandavg_timelock_planar_prior_distributed.mat',...
    'prior_1_faces','prior_2_faces','prior_3_faces',...
    'prior_1_houses','prior_2_houses','prior_3_houses');
