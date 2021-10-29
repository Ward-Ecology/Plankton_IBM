function process_phylogeny(location)

    disp(' ')
    disp('Extract sub-tree for max 1e6 agents alive at end of simulation')

    cd(location)
    load('workspace.mat','phylogeny','genome','eco_params')

    %%
    % %% extract and sample indices of extant agents and their matching bioinformatic data
    phylogeny_index = phylogeny.phylogeny_index;

    % % identify extant agents in phylogeny index (i.e. exclude empty elements)
    iii = find(phylogeny_index); % find nonzero rows of phylogeny_index

    % % extract phylogeny index and matching rgb and genome data
    phylogeny_index2 = phylogeny_index( iii );
    genome2          = genome(iii,:);
    
    clear phylogeny
    %% sample terminal leaf nodes
    nsample = min(1e6,size(phylogeny_index2,1)); % maximum 1e6 terminal leaf nodes
    isample = sort(randsample(numel(phylogeny_index2),nsample)); % sample without replacement
    % isample gives location of extant agents in phylogeny_index2 and matching rgb and genome darrays

    % extract sample of extant agents and corresponding rgb and genome data
    phylogeny_terminal = phylogeny_index2(isample);
%     genome_terminal =          genome2(isample,:);

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
    srcnode = table2array(gather(phylogeny_tall(:,1)));

    %% use srcnodes to backtrack descent of each extant agent to source node

    % create cell array with cell for each sampled agent (terminal leaf node)
    nodeid = cell(nsample,1);
    for i=1:nsample
        nodeid{i}=phylogeny_terminal(i); % index of terminal leaf node
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
   
    % merge rgb columns into single variable
    if eco_params.bioinf_verbose
        phylogeny_table = mergevars(phylogeny_table,[(6+eco_params.nrgb):size(phylogeny_table,2)]);
        table_headings = {'srcnode' 'snknode' 'divtime' 't_opt' 'otu' 'rgb_genome' 'binary_genome'};
    else
        table_headings = {'srcnode' 'snknode' 'divtime' 't_opt' 'otu' 'rgb_genome'};
    end
    % merge rgb columns into simgle variable
    phylogeny_table = mergevars(phylogeny_table,[6:(5+eco_params.nrgb)]);
    
    % assign names to variables
    phylogeny_table.Properties.VariableNames = table_headings;
    
    %%
%     % add binary genome of terminal nodes as additional variable
%     % (padded with zeros for non-terminal nodes)
%     phylogeny_table.binary_genome = uint64(zeros(numel(nodes),eco_params.ngenome));
%     % process and add large binary genome array column by column
%     % (too large to fit full padded intermediate array in memory)
%     for i=1:eco_params.ngenome
%         binary_genome = uint64(zeros(max(phylogeny_index2),1));
%         % load bioinformatic data into array using locations specified by phylogeny_index2
%         binary_genome(phylogeny_index2)    = genome2(:,i);
%         % extract sub tree nodes to save in table
%         phylogeny_table.binary_genome(:,i) = binary_genome(nodes);
%     end
%     
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
    
    
    
    %% load times of death and add to table
    ds2 = datastore('time_of_dying_*.csv');
    ds2.SelectedFormats       = {'%u32','%f'};
    time_of_dying = tall(ds2);
    disp(' ')
    disp('Loading times of death...')
    deaths_table = gather(time_of_dying(nodes,:));
    deadtime = table2array(gather(deaths_table(:,2)));
    phylogeny_table = addvars( phylogeny_table , deadtime,'After', 3 );
    
    
    %% identify terminal leaf nodes
    phylogeny_table.TLN = zeros(size(srcnode2));
    phylogeny_table.TLN(ii(phylogeny_terminal))=true;
    
    %% write out pruned tree data to new .mat file
    save('pruned_tree.mat','phylogeny_table','-v7.3')

    %%
end








