%% Time-step IBM
function [yout y tout eco_params bioinformatics phylogeny] = timestep_IBM(y,env_forcing,eco_params,bioinformatics,bioinf_fcns,phylogeny)



    % create blank arays for export
    % calculate average biomass value within each bin for easy comparison to MCM
    yout = zeros(env_forcing.tmax,size(eco_params.T_bins,2)+1);
    % outout time vector
    tout = 1:env_forcing.tmax;
    
    t1=tic;
    for t=0:eco_params.delta_t:env_forcing.tmax
        % evaluate equations at time=t
        [y eco_params bioinformatics phylogeny] = MUSE_IBM_eqs(t,y,env_forcing,eco_params,bioinformatics,bioinf_fcns,phylogeny);
            
        % print simulation information at set intervals
        if rem(t,1)==0 & t>0 
            if rem(t,30)==0
                time1=toc(t1);
            end
            
            % get P cell biomasses
            P   = y(2:end-1);

            % get index of extant agents lying in range of thermal bins used by MCM
            Topt=eco_params.T_optimum;
            iextant_inbin = P>0 & eco_params.mintemp<=Topt & eco_params.maxtemp>=Topt;

            % bin phenotypes into intervals used by MCM thermal optima
            [a,~,b]=histcounts(Topt(iextant_inbin),eco_params.T_bins);
            % get total biomass in each phenotypic bin
            c = accumarray(b,P(iextant_inbin),size(a'));

            % collate state variables for output
            yout(t,:)=[y(1);c;y(end)];

            % display computation info
%             clc
            if rem(t,30)==0
                disp(['Y:' num2str(ceil(t/360),'%03i')...
                    ', M:' num2str(ceil((rem(t-1,360)+1)/30),'%02i')...
                    ', D:' num2str(rem(t-1,30)+1,'%02i') ...
                    ', Prog. time ' num2str(time1,'%3.2f') ' sec' ...
                    ])
            end
        end
        
        
        % write phylogenetic data to file every 30 days
        if eco_params.use_bioinf & eco_params.bioinf_verbose
            if t==env_forcing.tmax || rem(t,30)==0 & t>0
                disp(['Writing phylogenetic data to {phylogeny.csv}...'])
                t2=tic;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % write phylogenetic data to output .csv file
                phylodata  = phylogeny.phylodata;
                ii = 1:find(phylogeny.phylogeny_index>0,1,'last'); % remove trailing zeros
                phylogeny_index = phylogeny.phylogeny_index(ii); % index of extant cells
                phylogeny = rmfield(phylogeny,'phylodata');
                szphy=size(phylodata); % record size of preallocated array
                
                phylodata = phylodata( find(phylodata(:,1)), : ) ; % trim empty rows from phylogeny array
                fname1=['output/phylogeny_' num2str(floor((t-1)/env_forcing.daysperyear)+1,'%04i') '.csv'];
                [fid,~] = fopen(fname1, 'a'); % append yearly file (or create new if not existing)
                fprintf(fid,eco_params.fmt,phylodata'); % print phylogeny array
                fclose(fid); 
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % cache times at which each agent ceases to exist (dies or divides into two new agents)
                [indx_of_dying,isort] = sort(phylogeny.time_of_dying(:,1));    % index of unwritten dead cells
                time_of_dying         =      phylogeny.time_of_dying(isort,2); % time of death or division
                                
                % remove zeros
                ii = find(indx_of_dying);                 
                indx_of_dying = indx_of_dying(ii);
                time_of_dying = time_of_dying(ii);
                
                % index of cells to cache (all those up to first extant cell)
                first_extant = min(phylogeny_index(phylogeny_index>0));
                cache_dead = find(indx_of_dying<first_extant);
                % index of cells not to cache (all those after first extant cell)
                pass_dead  = find(indx_of_dying>first_extant);
                
                % unless at end of run, in which case write all extant cells
                if t==env_forcing.tmax
                    % combine extant and recently died
                    indx_of_dying     = [indx_of_dying;phylogeny_index];
                    time_of_dying     = [time_of_dying;zeros(size(phylogeny_index))+t];
                    % re-sort
                    [~,ii] = sort(indx_of_dying);                 
                    indx_of_dying = indx_of_dying(ii);
                    time_of_dying = time_of_dying(ii);
                    % cache all and pass none
                    cache_dead    = find(indx_of_dying);
                    pass_dead     = [];
                end
                % write times of death and division
                if numel(cache_dead)>0
                    fname2=['output/time_of_dying_' num2str(floor((t-1)/env_forcing.daysperyear)+1,'%04i') '.csv'];
                    [fid,~] = fopen(fname2, 'a'); % append yearly file (or create new if not existing)
                    fprintf(fid,'%i,%f\n',[indx_of_dying(cache_dead) time_of_dying(cache_dead)].'); % print phylogeny array
                    fclose(fid);
                    % reset times of death
                    phylogeny.time_of_dying = phylogeny.time_of_dying.*0;
                    % pass unwritten to blank array
                    phylogeny.time_of_dying(1:numel(pass_dead),:) = [indx_of_dying(pass_dead) time_of_dying(pass_dead)];
                end
                
                % update n written
                eco_params.nwritten(3) = numel(pass_dead);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % reset phylodata and delete written death indices
                nrows=size(phylodata,1);
                % note n written and reset phylogeny to all zeros
                eco_params.nwritten(2) = eco_params.nwritten(1);
                % clear surplus variables and save workspace at end of each
                % year
                clear phylodata binarydata P Topt b iextant*
                if t==env_forcing.tmax || rem(t,env_forcing.daysperyear)==0
                    tic
                    !rm output/workspace.mat
                    genome = bioinformatics.genome;
                    save('output/workspace.mat','phylogeny','genome','eco_params','env_forcing','tout')
                    clear genome
                    toc
                end
                phylogeny.phylodata  = zeros(szphy);
                
                time2=toc(t2);
                disp(['Written ' ...
                      num2str(nrows) ' rows of phylogenetic data to {phylogeny.csv} in ' ...
                      num2str(time2,'%4.3f') ' sec'])
            end              
        end
        
        if rem(t,30)==0 & t>0
            t1=tic;
        end
    end
end