function [ bioinf_fcns ] = bioinf_functions( )
    bioinf_fcns.rgb_mutate     = @rgb_mutate;
    bioinf_fcns.gene_mutate    = @gene_mutate;
    bioinf_fcns.bioinf_inherit = @bioinf_inherit;
    bioinf_fcns.print_genomes  = @print_genomes;
    bioinf_fcns.maxrows        = @maxrows;
    bioinf_fcns.maxcols        = @maxcols;
    
end

% 'biomass' is [n_population x 1] vector
 
% 'rgb' is [n_population x n_rgb] matrix
% n_population is number of phenotypes
% n_rgb is number of orthogonal colour dimensions 
%   (3 for basic colour, use more for more discrimiatory power).
 
 
 

 
%% Apply single-point mutations to all genes in all extant populations
function [rgb] = rgb_mutate(rgb,biomass)
    % function to apply random mutations to each 'colour' of each extant
    % population. 'rgb' is matrix of neutral rgb tracers. 'biomass' is
    % population/community biomass
 
    % get locations of non-zero biomass
    kk=find(reshape(biomass,[],1)); 
 
    nrgb =size(rgb   ,2); % number of rgb tracers
    
    %% RGB TRACER
    if nrgb>0
        % extract tracers for extant populations
        extant=full(rgb(kk,:));
        % add random mutation
        extant = extant+randn(size(extant));
        % place back in full rgb array
        rgb(kk,:)=extant;
    end
    
end
 %%
 function [genome] = gene_mutate(genome,newbio,pmut)

    % get locations of non-zero biomass
    kk=find(reshape(newbio,[],1)); 

    ngenes=size(genome,2); % number of 53 bit genes
    
    %% BINARY GENOME
    if ngenes>0
        % extract genes for extant populations
        extant=full(genome(kk,:));
        
        npopn=numel(kk);      % number of extant populations
        nbit = 53; % uint64 (maximum 53 bits so integers also fit within float32)
        
        
        usePMut = true; %use logarithmic gradient for probability of mutating each gene
        if usePMut
            %%%%%% FROM Test_clock.m
            % identify genes in which to flip single base (1 in each population)
            i_gene_mut = rand([npopn ngenes])<pmut; %random number for each gene in each population
            igene = find(i_gene_mut); % indices of each gene selected - produces 1d vector ordered row by row
            % identify base to flip in those genes 
            ibase=randi(nbit,size(igene)); % bit to flip in each gene selected (normal dist)
            
            %get values of genes to mutate
            extantLin = extant(:);
            uint_values = extantLin(igene);

            % get original value of mutated bases (bits)
            mut=bitget(uint_values,ibase);
            % mutate them
            extantLin(igene) = bitset(uint_values,ibase,1-mut);

            %reconvert back to 2d vector
            extant = reshape(extantLin, [npopn, ngenes]);

            %%%%%%
        else
            ww=ones(1,ngenes); 
            %%%%%% old method (uniform selection)
            igene = randsample(ngenes, npopn, true, ww)'; % identify gene in which to flip single base
            % (each gene can flip >1 base)
            ibase = randi(nbit,npopn,1);    % identify base to flip in that gene

            % get index of genes to mutate within overall extant population
            genind=sub2ind([npopn ngenes],1:npopn,igene);
            
            % extract fp_values of mutant genes of extant populations
            uint_values=extant(genind);

            % get original value of mutated bases (bits)
            mut=bitget(uint_values(:),ibase(:));
            
            % get new genome after flipping a single bit
            extant(genind)=bitset(uint_values(:),ibase(:),1-mut);

            %%%%%%
        end



        
        % place back in full genome array
        genome(kk,:)=extant;
    end
    
end
 %%
function [binary_genome] = print_genomes(decimal_genome)
    
    nbit = 53; % uint64
    
    npopn =size(decimal_genome,1); % sample size
    ngenes=size(decimal_genome,2); % number of 64 bit genes
 
    % get original value of mutated bases (i.e. binary bits)
    genome_strings=dec2bin(reshape(decimal_genome',[],1),nbit);
    
    binary_genome=reshape(genome_strings',nbit*ngenes,npopn)';
    
  
end

%% Select inheritence of rgb genes based on dominant biomass source (perfect clonal or mutant)
% rgb_tracer function is based on mutation as fraction of growth rate
% ggr is gross growth rate, Pmut is mutation matrix
function [jj] = bioinf_inherit(P,mu,Pmut)
    %% Identify maximum source to each population (for each location)
    jmax=maxcols(P,mu,Pmut);

    % jmax is index of maximum product in intermediate dimension (of matrix multiplication)
    % (i.e. population [matrix column] that dominates biomass flux (gross growth and mutation) into each population)
    kk=sub2ind(size(jmax),repmat((1:size(jmax,1))',[1,size(jmax,2)]),jmax); % Convert this to linear index
    % reshape index arrays and subtract index
    % jj=sparse(reshape(jj,1,[])-(1:numel(jj))); % reduced memory version (subtract linear index and make spares)
    jj=reshape(kk,1,[]);
       
end
           
%% Functions to efficiently identify dominant source
% see https://uk.mathworks.com/matlabcentral/answers/375408-find-maximum-intermediate-product-in-matrix-multiplication?s_tid=prof_contriblnk
function K = maxrows(A,B)
    AT = A.'; % gets the product data contiguous in memory
    K = zeros(size(B));
    zi= 1:size(B,1);
    for j=1:size(B,2) % for each population, find dominant source location at each location
        [v,k] = max(AT.*B(:,j),[],1); % Uses implicit array expansion (all influxes to each location)
        k(v==0) = zi(v==0); % Set index to local point where maximum incoming flux == 0
        K(:,j) = k(:);
    end
end

function K = maxcols(P,mu,mutmat)
    P   = P.*rand(size(P));
    muP = (mu.*P)';
    
    B  = muP;
    AT = mutmat.'; % gets the product data contiguous in memory

    K = zeros(size(B));
    zi= 1:size(B,2);
    for j=1:size(B,1) % for each location, find dominant source population for each population
        sources = B(j,:).*AT; 
        dindx = find(speye(size(sources)));
        sources(dindx) = sources(dindx) + P; % add local continuing biomass to diagonal
        [v,k] = max(sources,[],2); % Uses implicit array expansion
        k(v==0) = zi(v==0); % Set index to local point where maximum incoming flux == 0
        K(j,:) = k(:);
    end
end
 
%%
 




