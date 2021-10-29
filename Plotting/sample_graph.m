clear 
addpath ~/GitHub/BeckmannPhytoplankton/Functions
bioinf_fcns  = bioinf_functions; 

cd ~/GitHub/BeckmannPhytoplankton/output
load('workspace.mat','eco_params','env_forcing')


% load('pruned_tree.mat');
load('last_year.mat');


sample_index = find(phylogeny_table.TLN);

sample_index = randsample(sample_index,1e5);

srcnode  = phylogeny_table.srcnode(sample_index);
snknode  = phylogeny_table.snknode(sample_index);
gentime  = phylogeny_table.divtime(snknode)-phylogeny_table.divtime(srcnode);
 divtime = phylogeny_table.divtime(sample_index);
deadtime = phylogeny_table.deadtime(sample_index);
   t_opt = phylogeny_table.t_opt(sample_index);
 rgbgene = phylogeny_table.rgb_genome(sample_index,:);
 

% normalise rgb gene (centered on origin)
maxclr  = max(abs(rgbgene));
minclr  = -maxclr;
plotclr = (rgbgene-minclr)./(maxclr-minclr);


scatter(divtime,t_opt,15,plotclr,'filled','MarkerFacealpha',1)





