function final_year_phylogeny(location)

    disp(' ')
    disp('Extract sub-tree for max 1e6 agents from last year of simulation')
    
    cd(location)
    load('workspace.mat','env_forcing','eco_params')

    %% load times of death
    ds2 = datastore('time_of_dying_*.csv');
    ds2.SelectedFormats       = {'%u32','%f'};
    time_of_dying = tall(ds2);
    disp(' ')
    disp('Loading times of death...')
    deaths_table = gather(time_of_dying);
    deadtime = table2array(gather(deaths_table(:,2)));
    
    % get index of agents still alive in final year
    i_lastyr = find(deadtime>env_forcing.tmax-env_forcing.daysperyear);
    
    %% sample agents still alive in final year
    nsample = min(1e6,numel(i_lastyr)); % maximum 1e6 terminal leaf nodes
    isample = sort(randsample(numel(i_lastyr),nsample)); % sample without replacement
    i_lastyr = i_lastyr(isample);

    %% assign metadata to tall array (not yet loaded into memory)
    ds1 = datastore('phylogeny*.csv');
    for i=[1 2 5]
        ds1.SelectedFormats{i}       = '%u32';
    end
    if eco_params.bioinf_verbose
        for i=(6+eco_params.nrgb):size(ds1.VariableNames,2)
            ds1.SelectedFormats{i}       = '%u64';
        end
    end

    phylogeny_tall = tall(ds1);
    % then load srcnodes for entire graph
    disp(' ')
    disp('Loading srcnode...')
    % get srcnode index
    srcnode = table2array(gather(phylogeny_tall(:,1)));
    %% use srcnodes to backtrack descent of each extant agent to source node

    % create cell array with cell for each sampled agent (terminal leaf node)
    nodeid = cell(nsample,1);
    for i=1:nsample
        nodeid{i}=i_lastyr(i); % index of terminal leaf node
        % iterate back through generations
        while nodeid{i}(end)~=1 % until root is reached
            nodeid{i}(end+1)=srcnode(nodeid{i}(end));
        end
    end
    % find unique nodes among all paths back from sampled terminal nodes
    [nodes] = unique(horzcat(nodeid{:}));
    % this is indexes all the nodes in the phylogenetic tree

    clear nodeid srcnode
    %% extract all meta data for sampled graph (base tree of sampled nodes)
    disp(' ')
    disp('Loading metadata for phylogenetic tree')
    
    % read in sampled tree from tall array
    phylogeny_table = gather(phylogeny_tall(nodes,:));
   
    % merge rgb columns into simgle variable
    if eco_params.bioinf_verbose
        phylogeny_table = mergevars(phylogeny_table,[(6+eco_params.nrgb):size(phylogeny_table,2)]);
        table_headings = {'srcnode' 'snknode' 'divtime' 't_opt' 'otu' 'rgb_genome' 'binary_genome'};
    else
        table_headings = {'srcnode' 'snknode' 'divtime' 't_opt' 'otu' 'rgb_genome'};
    end
    
    % merge rgb columns into simgle variable
    phylogeny_table = mergevars(phylogeny_table,5+[1:eco_params.nrgb]);
    
    % assign names to variables
    phylogeny_table.Properties.VariableNames = table_headings;
    
    %% reindex src and snknode IDs to make snk nodes consectutive
    srcnode = phylogeny_table.srcnode;
    snknode = phylogeny_table.snknode;
    % make snk node indices consecutive
    [~,snknode2] = sort(snknode);
    % reindex src nodes to match snk nodes
    ii = sparse(size(snknode2,2),1); % create blank sparse matrix
    ii(snknode,1) = snknode2;        % map original snk index to new snk index
    srcnode2 = full(ii(srcnode));    % convert src nodes using snknode mapping
    phylogeny_table.srcnode = srcnode2;
    phylogeny_table.snknode = snknode2;
    
    %% add times of death  to table
    deadtime = deadtime(nodes);
    phylogeny_table = addvars( phylogeny_table , deadtime,'After', 3 );
    
    
    %% identify sampled nodes
    phylogeny_table.TLN = zeros(size(srcnode2));
    phylogeny_table.TLN(ii(i_lastyr))=true;
    
    %% write out pruned tree data to new .mat file
    save('last_year.mat','phylogeny_table','-v7.3')
    
    
    %%
end












% end