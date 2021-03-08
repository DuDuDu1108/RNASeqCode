
%%%%%%%% Logic for comparing my RNA-seq with other RNA-seq

%%% A: Identify candidate genes (annotated, expressed across species)
    % 1: get the translation between Hg, Ms, Cy, entrez, ensembl, gene_symbol
    % 2: get the annotated orthologues among Hg, Cy, Ms (15220 genes)
    % 3: assign expressoion value to all the overlap genes and PCA them
    % 4: get the expressed orthologues (at least in one cell, log2(RPM+1) > 4
    % in at least one sample)

%%% B: Normalize expression (log2(X+1) -> zscore separately)
    % 1: read in all the datasets (my RNA-seq choose TPM, their choose RPM,
    % because they use SC3 RNA seq only identify mRNA 3' end, similar to TPM)
    % 2: log2(X+1) and then zscore them separately

%%% C: Visual analysis
    % 1: find corresponding normalized expression for all candiate genes
    % 2: PCA analysis (PC1 is referring to species, can check; PC2 is referring
    % to developmental stage? Draw with PC2 and PC3)
    % 3: other plots?


%%%%%%%% Additional information:
% v1_HgCyMs_translation (version 1) ---- all the genes align through
% species (18491 genes)
% v2_HgCyMs_translation (version 2) ---- genes have orthologs in three
% species (15150 genes, some genes couldn't find ensembl id through DAVID)
% v3_HgCyMs_translation (version 3) ---- v2 plus corresponding genes
% expression (order is Hg8, Cy421, Ms108)
% v4_HgCyMs_translation (version 4) ---- fill v3 empty cell with value 0
% v5_HgCyMs_translation (version 5) ---- zscore(expression level)
% separately in species and then combine them together


addpath(genpath('/Users/dududu/CRAZY_SCIENNY/RNA-Seq'));
addpath(genpath('/Users/dududu/CRAZY_SCIENNY/MATLAB'));

%%%%%%%% The part need to run everytime

%%% read in all the data source
    
% Read Cy and Ms RNA-seq files and transfer them into matrix
% (Cy is monkey and Ms is mice).
% From paper "A developmental coordinate of pluripotency among mice,
% monkeys and humans"
very_raw_Cy = readtable('ReadInFiles/From Nakamura 2016/GSE74767_SC3seq_Cy_ProcessedData.txt');
raw_Cy = table2cell(very_raw_Cy(:,:)); Cy_SampleNames = very_raw_Cy.Properties.VariableNames;
very_raw_Ms = readtable('ReadInFiles/From Nakamura 2016/GSE74767_SC3seq_Ms_ProcessedData.txt');
raw_Ms = table2cell(very_raw_Ms(:,:)); Ms_SampleNames = very_raw_Ms.Properties.VariableNames;

% Read my Hg RNA-seq files
%%% add column name in 'My Hg RNA-seq data.xls'
very_raw_Hg = readtable('ReadInFiles/MyCells/My Hg RNA-seq data TPM.xls');
raw_Hg = table2cell(very_raw_Hg(:,:)); Hg_SampleNames = very_raw_Hg.Properties.VariableNames;

load v5-Hg-Ms-Cy-translation.mat v5_HgCyMs_translation
load v4-Hg-Ms-Cy-translation.mat v4_HgCyMs_translation
load Cy_CellNames.mat Cy_CellNames; load Ms_CellNames.mat Ms_CellNames; load Hg_CellNames.mat Hg_CellNames


%% AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

% Identify candidate genes (annotated, expressed across species)


%% 1. Get the translation between Hg, Ms, Cy, entrez, ensembl, gene_symbol


    %%% I have original two are CyToHg and MsToHg from Nakamura paper
        CyToHg = table2cell(readtable('ReadInFiles/From Nakamura 2016/CyToHg.xlsx'));
        MsToHg = table2cell(readtable('ReadInFiles/From Nakamura 2016/MsToHg.xlsx'));
    
    %%% combine Hg_entrez_id from CyToHg and MsToHg, unique and sort
        templist = cell(length(CyToHg)+length(MsToHg),1);
        templist(1:length(CyToHg),1) = CyToHg(:,1);
        templist((length(CyToHg)+1):(length(CyToHg)+length(MsToHg)),1) = MsToHg(:,1);
        templist = cell2mat(templist);
        templist = sort(unique(templist),'ascend');
        templist = num2cell(templist); % (all the human entrez_id appeared)
    
    %%% Put corresponding Hg, Cy, Ms, entrez/gene_symbol/ensembl together
    
        % (column order is: Hg_entrez_id, Hg_gene_symbol, Hg_ensembl_id,
        % Cy_entrez_id, Cy_gene_symbol, Ms_entrez_id, Ms_gene_symbol)

        % by compare text
        v1_HgCyMs_translation = cell(length(templist),7);
        v1_HgCyMs_translation(:,1) = templist(:,1);
        for v1Idx = 1:length(v1_HgCyMs_translation)
            for CyIdx = 1:length(CyToHg)
                Cycompare = strcmp(num2str(v1_HgCyMs_translation{v1Idx,1}), num2str(CyToHg{CyIdx,1}));
                if Cycompare == 1
                    v1_HgCyMs_translation{v1Idx,2} = CyToHg{CyIdx,2}; % Hg_gene_symbol
                    v1_HgCyMs_translation(v1Idx,4:5) = CyToHg(CyIdx,3:4); % Cy_entrez_id, Cy_gene_symbol
                end
            end
            for MsIdx = 1:length(MsToHg)
                Mscompare = strcmp(num2str(v1_HgCyMs_translation{v1Idx,1}), num2str(MsToHg{MsIdx,1}));
                if Mscompare == 1
                    v1_HgCyMs_translation(v1Idx,6:7) = MsToHg(MsIdx,3:4); % Ms_entrez_id, Ms_gene_symbol
                end
            end
        end
        save v1-Hg-Ms-Cy-translation.mat v1_HgCyMs_translation
        writecell(v1_HgCyMs_translation, 'temp_v1_HgCyMs_translation.xls')
        %     % by find index (v1Idx = 2027, Ms has repeated genes, faster)
        %     v1_HgCyMs_translation2 = cell(length(templist),7);
        %     v1_HgCyMs_translation2(:,1) = templist(:,1);
        %     for v1Idx = 1:length(v1_HgCyMs_translation2)
        %         CyIdx = find([CyToHg{:,1}] == v1_HgCyMs_translation2{v1Idx,1});
        %         if ~isempty(CyIdx)
        %             v1_HgCyMs_translation2{v1Idx,2} = CyToHg{CyIdx,2}; % Hg_gene_symbol
        %             v1_HgCyMs_translation2(v1Idx,4:5) = CyToHg(CyIdx,3:4); % Cy_entrez_id, Cy_gene_symbol
        %         end
        %         MsIdx = find([MsToHg{:,1}] == v1_HgCyMs_translation2{v1Idx,1});
        %         if ~isempty(MsIdx)
        %             v1_HgCyMs_translation2(v1Idx,6:7) = MsToHg(MsIdx,3:4); % Cy_entrez_id, Cy_gene_symbol
        %         end
        %     end

    %%% insert one line at the beginning label what each column is in 'v1_HgCyMs_translation.xls'
    
    
    
%% 2. Get the annotated orthologues among Hg, Cy, Ms (15220 genes in paper, I only got 15150 genes)

    %%% use DAVID webtool translate Hg_entrez_id to ensembl_id (will lose
    %%% some)
    
    
    %%% fill emsembl id to v1_HgCyMs_translation
        EntrezToEnsembl = table2cell(readtable('ReadInFiles/HgEntrezToEnsembl.xlsx'));
        for v1Idx = 1:length(v1_HgCyMs_translation)
            EmsemblIdx = find([EntrezToEnsembl{:,1}] == v1_HgCyMs_translation{v1Idx,1});
            if ~isempty(EmsemblIdx)
               v1_HgCyMs_translation{v1Idx,3} = EntrezToEnsembl{EmsemblIdx,2}; % Hg_ensembl_id
            end
        end
        writecell(v1_HgCyMs_translation, 'HgCyMs_translation_All.xls')
    %%% delete 'v1_HgCyMs_translation.xls'
    %%% insert one line at the beginning label what each column is in 'All_HgCyMs_translation.xls'
    
    %%% get rid of genes that not across all three species
        v2_HgCyMs_translation = cell(length(templist),7);
        v2Idx = 1;
        for v1Idx = 1:length(v1_HgCyMs_translation)
            EmptyIdx = cellfun(@isempty,v1_HgCyMs_translation(v1Idx,:));
            allNonEmpty = all(~EmptyIdx);
            if allNonEmpty
                v2_HgCyMs_translation(v2Idx,:) = v1_HgCyMs_translation(v1Idx,:);
                v2Idx = v2Idx + 1;
            end
        end
        writecell(v2_HgCyMs_translation, 'HgCyMs_translation_Overlap.xls')
    %%% insert one line at the beginning label what each column is in 'Overlap_HgCyMs_translation.xls'
    
    
    
%% 3. Assign expressoion value to all the overlap genes and PCA them


    %%% read in all the data source
    
        % Read Cy and Ms RNA-seq files and transfer them into matrix
        % (Cy is monkey and Ms is mice).
        % From paper "A developmental coordinate of pluripotency among mice,
        % monkeys and humans"
        very_raw_Cy = readtable('ReadInFiles/From Nakamura 2016/GSE74767_SC3seq_Cy_ProcessedData.txt');
        raw_Cy = table2cell(very_raw_Cy(:,:)); Cy_SampleNames = very_raw_Cy.Properties.VariableNames;
        very_raw_Ms = readtable('ReadInFiles/From Nakamura 2016/GSE74767_SC3seq_Ms_ProcessedData.txt');
        raw_Ms = table2cell(very_raw_Ms(:,:)); Ms_SampleNames = very_raw_Ms.Properties.VariableNames;

        % Read my Hg RNA-seq files
        %%% add column name in 'My Hg RNA-seq data.xls'
        very_raw_Hg = readtable('ReadInFiles/MyCells/My Hg RNA-seq data.xls');
        raw_Hg = table2cell(very_raw_Hg(:,:)); Hg_SampleNames = very_raw_Hg.Properties.VariableNames;
        
        
    %%% find the corresponding expression level
    
        %%% make a form for gene name (7 columns) and corresponding
        %%% expression level: row order is the same as
        %%% 'v2_HgCyMs_translation'; column order is the same as gene label(7 columns,1:7), Hg(8,8:15),
        %%% Cy(421,16:436), Ms(108,437:544)
    
        %%% read in 'HgCyMs_translation_Overlap.xls'
            v2_HgCyMs_translation = readtable('HgCyMs_translation_Overlap.xls');
        %%% fill in with values
            v3_HgCyMs_translation = cell(height(v2_HgCyMs_translation), size(raw_Hg,2)+size(raw_Cy,2)+size(raw_Ms,2)+2);
            v3_HgCyMs_translation(:,1:7) = table2cell(v2_HgCyMs_translation(:,:));
            for v3Idx = 1:length(v3_HgCyMs_translation)

                % Hg genes (8:15)
                for HgIdx = 1:length(raw_Hg)
                    Hgcompare = strcmp(v3_HgCyMs_translation{v3Idx,3}, raw_Hg{HgIdx,1});
                    if Hgcompare == 1
                        v3_HgCyMs_translation(v3Idx,8:15) = raw_Hg(HgIdx,2:9);
                    end
                end

                % Cy genes (16:436)
                CyIdx = find([raw_Cy{:,1}] == v3_HgCyMs_translation{v3Idx,4});
                if ~isempty(CyIdx)
                    v3_HgCyMs_translation(v3Idx,16:436) = raw_Cy(CyIdx,3:423);
                end

                % Ms genes (437:544)
                MsIdx = find([raw_Ms{:,1}] == v3_HgCyMs_translation{v3Idx,6});
                if ~isempty(MsIdx)
                    v3_HgCyMs_translation(v3Idx,437:544) = raw_Ms(MsIdx,3:110);
                end


            end
            save v3-Hg-Ms-Cy-translation.mat v3_HgCyMs_translation % (next time just need to load)
            
            
        %%% find the corresponding Zscore(expression level) seaprately in
        %%% different species
            
            % fill all the empty cell with 0
            v4_HgCyMs_translation = v3_HgCyMs_translation;
            EmptyIdx = cellfun(@isempty, v4_HgCyMs_translation);
            v4_HgCyMs_translation(EmptyIdx) = {0};
            save v4-Hg-Ms-Cy-translation.mat v4_HgCyMs_translation % (next time just need to load)
            
            % Zscore separately in species (column order is the same as gene label(7 columns,1:7), Hg(8,8:15),
            % Cy(421,16:436), Ms(108,437:544))
            v5_HgCyMs_translation = cell(length(v4_HgCyMs_translation), size(raw_Hg,2)+size(raw_Cy,2)+size(raw_Ms,2)+2);
            v5_HgCyMs_translation(:,1:7) = v4_HgCyMs_translation(:,1:7);
            temp_mat = cell2mat(v4_HgCyMs_translation(:,8:end));
            v5_HgCyMs_translation(:,8:15) = num2cell(zscore(log2(temp_mat(:,1:8)+1),[],2)); % Hg
            v5_HgCyMs_translation(:,16:436) = num2cell(zscore(temp_mat(:,9:429),[],2)); % Cy
            v5_HgCyMs_translation(:,437:544) = num2cell(zscore(temp_mat(:,430:537),[],2)); % Ms
            save v5-Hg-Ms-Cy-translation.mat v5_HgCyMs_translation % (next time just need to load)
            
          
%% (extra) Find each sample name corresponding cell type (for PCA label, run everytime)

        clear CyCtoS MsCtoS HgCtoS Cy_CellNames Ms_CellNames Hg_CellNames
    %%% get the 1-to-1 map for celltype and samplename
        Sampleinfo = table2cell(readtable('ReadInFiles/From Nakamura 2016/SampleTable.xlsx'));
        Cycount = 0; Mscount = 0; Hgcount = 0;
        for ii = 1:length(Sampleinfo)
            compareCy = strcmp(Sampleinfo{ii,4}, 'Cynomolgus Monkey');
            compareMs = strcmp(Sampleinfo{ii,4}, 'Mouse');
            compareHg = strcmp(Sampleinfo{ii,4}, 'Human');
            if compareCy == 1
                Cycount = Cycount + 1;
                CyCtoS{Cycount,1} = Sampleinfo{ii,3}; % 1st column is cell type name
                CyCtoS{Cycount,2} = Sampleinfo{ii,1}; % 3rd column is sample name
            end
            if compareMs == 1
                Mscount = Mscount + 1;
                MsCtoS{Mscount,1} = Sampleinfo{ii,3};
                MsCtoS{Mscount,2} = Sampleinfo{ii,1};
            end
            if compareHg == 1
                Hgcount = Hgcount + 1;
                HgCtoS{Hgcount,1} = Sampleinfo{ii,3};
                HgCtoS{Hgcount,2} = Sampleinfo{ii,1};
            end
        end
    %%% assign each sample column number in RNAseq and put in 3rd column of
    %%% (Specie)CtoS
        for mapid = 1:length(CyCtoS) % Cy
            for rnaseqid = 3:length(Cy_SampleNames)
                compare = strcmp(CyCtoS{mapid,2}, Cy_SampleNames{1,rnaseqid});
                if compare == 1
                    CyCtoS{mapid,3} = rnaseqid;
                end
            end
        end
        for mapid = 1:length(MsCtoS) % Ms
            for rnaseqid = 3:length(Ms_SampleNames)
                compare = strcmp(MsCtoS{mapid,2}, Ms_SampleNames{1,rnaseqid});
                if compare == 1
                    MsCtoS{mapid,3} = rnaseqid;
                end
            end
        end
        
    %%% get cell type names and sample column number from RNAseq
        Cy_CellNames = unique(CyCtoS(:,1)); Ms_CellNames = unique(MsCtoS(:,1)); Hg_CellNames = unique(HgCtoS(:,1));
        for Cyid = 1:length(Cy_CellNames) % Cy
            Cy_CellNames{Cyid,2} = [];
            for mapid = 1:length(CyCtoS)
                compare = strcmp(CyCtoS{mapid,1}, Cy_CellNames{Cyid,1});
                if compare == 1
                    Cy_CellNames{Cyid,2} = [Cy_CellNames{Cyid,2},CyCtoS{mapid,3}];
                end
            end
        end
        for Msid = 1:length(Ms_CellNames) % Ms
            Ms_CellNames{Msid,2} = [];
            for mapid = 1:length(MsCtoS)
                compare = strcmp(MsCtoS{mapid,1}, Ms_CellNames{Msid,1});
                if compare == 1
                    Ms_CellNames{Msid,2} = [Ms_CellNames{Msid,2},MsCtoS{mapid,3}];
                end
            end
        end
        
    %%% save the cell sample relationship
        save Cy_CellNames.mat Cy_CellNames; save Ms_CellNames.mat Ms_CellNames; save Hg_CellNames.mat Hg_CellNames
        
        