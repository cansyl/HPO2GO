% HPO2GO v1.1
%
% Mapping between Human Phenotype Ontology (HPO) and Gene Ontology (GO)
% terms for the prediction of gene/protein - function - phenotype - disease
% associations.
%
% ------------------------------------------------------------------------
% 
% HPO2GO: Prediction of Human Phenotype Ontology Term Associations with
% Cross Ontology Annotation Co-occurrences
%
% Author: Tunca Dogan1,2,3,*
% 
% 1 Cancer Systems Biology Laboratory (CanSyL), Graduate School of
%   Informatics, METU, Ankara, 06800, Turkey
% 2 Department of Health Informatics, Graduate School of Informatics, METU,
%   Ankara, 06800, Turkey
% 3 European Molecular Biology Laboratory, European Bioinformatics
%   Institute (EMBL-EBI), Hinxton, Cambridge, CB10 1SD, UK
%
% ------------------------------------------------------------------------
%
% Copyright (C) 2018 CanSyL
% 
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program. If not, see http://www.gnu.org/licenses/.
%
% ------------------------------------------------------------------------
%
%
%
% ----------------------------------
% Source Code of the Whole Analysis
% ----------------------------------


% Loading HPO annotation data:

[HPO_gene_id,HPO_gene_symbol,HPO_Name,HPO_ID]=textread('HPO_gene_to_phenotype_annotation_01_2017_ALL_SOURCES_ALL_FREQUENCIES.txt', '%s %s %s %s', 'delimiter', '\t','headerlines',1);
HPO_gene_annotation=[HPO_gene_id HPO_gene_symbol HPO_Name HPO_ID];
HPO_gene_annotation=uniqueRowsCA(HPO_gene_annotation);
save HPO2GO_Files/HPO_test_HPO_gene_annotation.mat HPO_gene_annotation -v7
HPO_gene_symbol_unique=unique(HPO_gene_annotation(:,2));
dlmcell('HPO2GO_Files/HPO_gene_symbol_unique.txt',HPO_gene_symbol_unique);

% Loading all GO annoatation data for human proteins with manual evidence code:

[~,GO_UniProt_acc_all,GO_gene_symbol_all,~,GO_term_id_all,GO_term_name_all,GO_term_aspect_all,~,~,~,~,~,~]=textread('GOA_UniProt_human_protein_annotation.tsv', '%s %s %s %s %s %s %s %s %s %s %s %s %s', 'delimiter', '\t','headerlines',1);
GO_annot_manual_human_all=[GO_gene_symbol_all GO_term_id_all];
GO_annot_manual_human_all=uniqueRowsCA(GO_annot_manual_human_all);
save HPO2GO_Files/HPO_test_GO_annotation_all.mat GO_annot_manual_human_all -v7
GO_annot_manual_human_all_addcol=[GO_gene_symbol_all GO_term_id_all GO_term_name_all GO_term_aspect_all];
GO_annot_manual_human_all_addcol=uniqueRowsCA(GO_annot_manual_human_all_addcol);
save HPO2GO_Files/HPO_test_GO_annotation_all_addcol.mat GO_annot_manual_human_all_addcol -v7

% Saving the processed GO annotation file will additional columns as text:

GO_annot_manual_human_all_addcol_txt=GO_annot_manual_human_all_addcol';
fid=fopen('HPO2GO_Files/GO_annot_human_proteins_UniProtGOA_01_2017.txt', 'w');
fprintf(fid, '%s\t%s\t%s\t%s\n', GO_annot_manual_human_all_addcol_txt{:});
fclose(fid);

GO_annot_manual_human_all_uniprotacc=[GO_UniProt_acc_all GO_gene_symbol_all GO_term_id_all];
GO_annot_manual_human_all_uniprotacc=uniqueRowsCA(GO_annot_manual_human_all_uniprotacc);
save HPO2GO_Files/HPO_test_GO_annotation_all_uniprotacc.mat GO_annot_manual_human_all_uniprotacc -v7

% Selecting the GO annotations for the HPO annotated protein list:

[Lia,Locb]=ismember(GO_annot_manual_human_all(:,1),HPO_gene_symbol_unique);
GO_annot_manual_human=GO_annot_manual_human_all(Lia==1,:);
save HPO2GO_Files/HPO_test_GO_annotation.mat GO_annot_manual_human -v7

% Generating the HPO and GO term pairs:

to=0;
HPO_GO_mapping_ind=zeros(2000000,2);
for i=1:length(HPO_gene_symbol_unique);
    disp(['Mapping HPO terms with GO terms with gene #: ', num2str(i), ' / ', num2str(length(HPO_gene_symbol_unique))])
    [Lia,Locb]=ismember(HPO_gene_annotation(:,2),HPO_gene_symbol_unique(i,1));
    [Lia2,Locb2]=ismember(GO_annot_manual_human(:,1),HPO_gene_symbol_unique(i,1));
    if sum(Lia2)>0
        siz=sum(Lia)*sum(Lia2);
        [A,B]=meshgrid(find(Lia==1),find(Lia2==1));
        c=cat(2,A',B');
        d=reshape(c,[],2);
        HPO_GO_mapping_ind(to+1:to+siz,:)=d;
        to=to+siz;
    end
end
HPO_GO_mapping_ind(HPO_GO_mapping_ind(:,1)==0,:)=[];
HPO_GO_mapping_terms=HPO_gene_annotation(HPO_GO_mapping_ind(:,1),4);
HPO_GO_mapping_terms(:,2)=GO_annot_manual_human(HPO_GO_mapping_ind(:,2),2);
save HPO2GO_Files/HPO_GO_mapping_terms.mat HPO_GO_mapping_terms -v7

% Generating CAFA id vs. gene symbol vs. entry name vs. UniProt acc. mapping file:

% (downloading all swissprot human proteins with necessary columns)

[Humprot_UniProt_acc,Humprot_gene_symbol,Humprot_entry_name]=textread('HPO2GO_Files/Protein_acc_name_mapping.tab', '%s %s %s', 'delimiter', '\t','headerlines',1);
Humprot_gene_symbol=strtok(Humprot_gene_symbol);
save HPO2GO_Files/Humprot_mapping_variables.mat Humprot_UniProt_acc Humprot_gene_symbol Humprot_entry_name -v7
[CAFA3_Targets_human_headers,~]=fastaread('CAFA3_targets/Target files/target.9606.fasta');
CAFA3_Targets_human_headers=CAFA3_Targets_human_headers';
CAFA3_Targets_human_idmapping=cellfun(@(x) strsplit(x,' '), CAFA3_Targets_human_headers(:),'uni',0);
CAFA3_Targets_human_idmapping=vertcat(CAFA3_Targets_human_idmapping{:,1});
save HPO2GO_Files/CAFA3_HPO_target_variables.mat CAFA3_Targets_human_headers CAFA3_Targets_human_idmapping -v7
[~,Locb]=ismember(CAFA3_Targets_human_idmapping(:,2),Humprot_entry_name);
Locb(Locb==0,1)=length(Humprot_entry_name)+1;
Humprot_entry_name(end+1,1)=cellstr(' ');
Humprot_UniProt_acc(end+1,1)=cellstr(' ');
Humprot_gene_symbol(end+1,1)=cellstr(' ');
CAFA3_Targets_human_all_mappings=CAFA3_Targets_human_idmapping;
CAFA3_Targets_human_all_mappings(:,3)=Humprot_UniProt_acc(Locb,1);
CAFA3_Targets_human_all_mappings(:,4)=Humprot_gene_symbol(Locb,1);
save HPO2GO_Files/CAFA3_Targets_human_all_mappings.mat CAFA3_Targets_human_all_mappings -v7

% Calculating the co-occurrence similarity between ontology terms:

[HPO_GO_mapping_terms_unique,I,J]=uniqueRowsCA(HPO_GO_mapping_terms);
HPO_GO_mapping_terms_freq=sort(J);
HPO_GO_mapping_terms_pair_hist=hist(HPO_GO_mapping_terms_freq,0.5:1:max(HPO_GO_mapping_terms_freq)+0.5);
HPO_GO_mapping_terms_pair_hist=HPO_GO_mapping_terms_pair_hist(1,1:end-1)';
save HPO2GO_Files/HPO_GO_mapping_terms_unique.mat HPO_GO_mapping_terms_unique -v7
save HPO2GO_Files/HPO_GO_mapping_terms_pair_hist.mat HPO_GO_mapping_terms_pair_hist -v7

[HPO_terms_unique,~,N]=unique(HPO_gene_annotation(:,4));
HPO_terms_annot_hist=hist(N,0.5:1:max(N)+0.5);
HPO_terms_annot_hist=HPO_terms_annot_hist(1,1:end-1)';
save HPO2GO_Files/HPO_terms_unique.mat HPO_terms_unique -v7
save HPO2GO_Files/HPO_terms_annot_hist.mat HPO_terms_annot_hist -v7

[GO_terms_unique,~,N]=unique(GO_annot_manual_human(:,2));
GO_terms_annot_hist=hist(N,0.5:1:max(N)+0.5);
GO_terms_annot_hist=GO_terms_annot_hist(1,1:end-1)';
save HPO2GO_Files/GO_terms_unique.mat GO_terms_unique -v7
save HPO2GO_Files/GO_terms_annot_hist.mat GO_terms_annot_hist -v7

[~,HPO_GO_mapping_unique_array_ind]=ismember(HPO_GO_mapping_terms_unique(:,1),HPO_terms_unique);
[~,HPO_GO_mapping_unique_array_ind(:,2)]=ismember(HPO_GO_mapping_terms_unique(:,2),GO_terms_unique);

HPO_GO_mapping_sematic_sim_high_det=zeros(length(HPO_GO_mapping_terms_unique),5);
for i=1:length(HPO_GO_mapping_terms_unique);
    disp(['Calculating the semantic similarity between the HPO ang GO term pair #: ', num2str(i), ' / ', num2str(length(HPO_GO_mapping_terms_unique))])
    t1=HPO_GO_mapping_terms_pair_hist(i,1);
    t2=HPO_terms_annot_hist(HPO_GO_mapping_unique_array_ind(i,1),1);
    t3=GO_terms_annot_hist(HPO_GO_mapping_unique_array_ind(i,2),1);
    HPO_GO_mapping_sematic_sim_high_det(i,1:5)=[i (2*t1)/(t2+t3) t1 t2 t3]; % columns: indice of mapping, semantic similarity of the mapped terms, # of co-annotated genes, total # of annotation of the mapped HPO term on different genes, total # of annotation of the mapped GO term on different genes
end
save HPO2GO_Files/HPO_GO_mapping_sematic_sim_high_det.mat HPO_GO_mapping_sematic_sim_high_det -v7

HPO_GO_Raw_Mapping_merged=[HPO_GO_mapping_terms_unique num2cell(HPO_GO_mapping_sematic_sim_high_det(:,2:end))];
save HPO2GO_Files/HPO_GO_Raw_Mapping_merged.mat HPO_GO_Raw_Mapping_merged -v7

% Saving the original raw mapping file as text:

HPO_GO_Raw_Mapping_merged_txt=HPO_GO_Raw_Mapping_merged';
fid=fopen('HPO_GO_Raw_Original_Mapping.txt', 'w');
fprintf(fid, '%s\t%s\t%.4f\t%d\t%d\t%d\n', HPO_GO_Raw_Mapping_merged_txt{:});
fclose(fid);

% Extracting the statistics for the original mapping - part of Table 1:

HPO_GO_mapping_sematic_sim=HPO_GO_mapping_sematic_sim_high_det(:,2);
S=1.1;
Mapping_stat=0;
for i=1:10;
    S=S-0.1;
    Mapping_stat(i,1)=S;
    Mapping_stat(i,2)=length(find(HPO_GO_mapping_sematic_sim>=S));
    Mapping_stat(i,3)=length(unique(HPO_GO_mapping_terms_unique(HPO_GO_mapping_sematic_sim>=S,1)));
    Mapping_stat(i,4)=length(unique(HPO_GO_mapping_terms_unique(HPO_GO_mapping_sematic_sim>=S,2)));
    Mapping_stat(i,5)=length(find(ismember(HPO_gene_annotation(:,4),unique(HPO_GO_mapping_terms_unique(HPO_GO_mapping_sematic_sim>=S,1)))==1));
    Mapping_stat(i,6)=length(find(ismember(GO_annot_manual_human(:,2),unique(HPO_GO_mapping_terms_unique(HPO_GO_mapping_sematic_sim>=S,2)))==1));
end
Mapping_table_title={['Threshold'], ['# of unique mappings'], ['# of unique HPO terms'], ['# of unique GO terms'], ['Total # of annot with HPO term'], ['Total # of annot with GO term']};
Mapping_stat_cel=num2cell(Mapping_stat);
Mapping_table=vertcat(Mapping_table_title,Mapping_stat_cel);
save HPO2GO_Files/HPO_sematic_sim_Mapping_stat_table.mat Mapping_table -v7

% Plotting the original co-occurrence similarity distributions for different numbers of co-annotated
% genes - Figure 3 (only when m=0) and part of Figure 1:

n=[1 5 25 50];
m=[0 0.1 0.2];
for i=1:4;
    for j=1:3;
        P=HPO_GO_mapping_sematic_sim_high_det(HPO_GO_mapping_sematic_sim_high_det(:,3)>=(n(1,i)),2);
        X=max(hist(P,0.00005:0.01:1.00005));
        P=P(P>=m(1,j));
        figure;hist(P,0.00005:0.01:1.00005);
        hold on;
        x_lab=xlabel('Co-occurrence similarity bins');
        y_lab=ylabel('# of unique mappings');
        set(x_lab,'FontSize',14);
        set(y_lab,'FontSize',14);
        eval(['title(''Co-occurrence Similarity Distribution (n >= ', num2str(n(1,i)), ')'',''FontSize'',16);'])
        grid on
        axis([0 0.3 0 (X*(11/10))])
        hold off;
    end
end

% Generating the randomized set (for comparison with the original set, to determine the thresholds):

GO_annot_manual_human_rand=randsample(GO_annot_manual_human(:,1),length(GO_annot_manual_human(:,1)));
GO_annot_manual_human_rand(:,2)=randsample(GO_annot_manual_human(:,2),length(GO_annot_manual_human(:,2)));
GO_annot_manual_human_rand=uniqueRowsCA(GO_annot_manual_human_rand);
HPO_gene_annotation_rand=randsample(HPO_gene_annotation(:,2),length(HPO_gene_annotation(:,2)));
HPO_gene_annotation_rand(:,2)=randsample(HPO_gene_annotation(:,4),length(HPO_gene_annotation(:,4)));
HPO_gene_annotation_rand=uniqueRowsCA(HPO_gene_annotation_rand);

to=0;
HPO_GO_mapping_ind_rand=zeros(2000000,2);
for i=1:length(HPO_gene_symbol_unique);
    disp(['Mapping HPO terms with GO terms with gene #: ', num2str(i), ' / ', num2str(length(HPO_gene_symbol_unique))])
    [Lia,Locb]=ismember(HPO_gene_annotation_rand(:,1),HPO_gene_symbol_unique(i,1));
    [Lia2,Locb2]=ismember(GO_annot_manual_human_rand(:,1),HPO_gene_symbol_unique(i,1));
    if sum(Lia2)>0
        siz=sum(Lia)*sum(Lia2);
        [A,B]=meshgrid(find(Lia==1),find(Lia2==1));
        c=cat(2,A',B');
        d=reshape(c,[],2);
        HPO_GO_mapping_ind_rand(to+1:to+siz,:)=d;
        to=to+siz;
    end
end
HPO_GO_mapping_ind_rand(HPO_GO_mapping_ind_rand(:,1)==0,:)=[];
HPO_GO_mapping_terms_rand=HPO_gene_annotation_rand(HPO_GO_mapping_ind_rand(:,1),2);
HPO_GO_mapping_terms_rand(:,2)=GO_annot_manual_human_rand(HPO_GO_mapping_ind_rand(:,2),2);
save HPO2GO_Files/HPO_GO_mapping_terms_rand.mat HPO_GO_mapping_terms_rand -v7

[HPO_GO_mapping_terms_unique_rand,I,J]=uniqueRowsCA(HPO_GO_mapping_terms_rand);
HPO_GO_mapping_terms_freq_rand=sort(J);
HPO_GO_mapping_terms_pair_hist_rand=hist(HPO_GO_mapping_terms_freq_rand,0.5:1:max(HPO_GO_mapping_terms_freq_rand)+0.5);
HPO_GO_mapping_terms_pair_hist_rand=HPO_GO_mapping_terms_pair_hist_rand(1,1:end-1)';
save HPO2GO_Files/HPO_GO_mapping_terms_unique_rand.mat HPO_GO_mapping_terms_unique_rand -v7
save HPO2GO_Files/HPO_GO_mapping_terms_pair_hist_rand.mat HPO_GO_mapping_terms_pair_hist_rand -v7

[HPO_terms_unique_rand,~,N]=unique(HPO_gene_annotation_rand(:,2));
HPO_terms_annot_hist_rand=hist(N,0.5:1:max(N)+0.5);
HPO_terms_annot_hist_rand=HPO_terms_annot_hist_rand(1,1:end-1)';
[GO_terms_unique_rand,~,N]=unique(GO_annot_manual_human_rand(:,2));
GO_terms_annot_hist_rand=hist(N,0.5:1:max(N)+0.5);
GO_terms_annot_hist_rand=GO_terms_annot_hist_rand(1,1:end-1)';
[~,HPO_GO_mapping_unique_array_ind_rand]=ismember(HPO_GO_mapping_terms_unique_rand(:,1),HPO_terms_unique_rand);
[~,HPO_GO_mapping_unique_array_ind_rand(:,2)]=ismember(HPO_GO_mapping_terms_unique_rand(:,2),GO_terms_unique_rand);

HPO_GO_mapping_sematic_sim_high_det_rand=zeros(length(HPO_GO_mapping_terms_unique_rand),5);
for i=1:length(HPO_GO_mapping_terms_unique_rand);
    disp(['Calculating the semantic similarity between the HPO ang GO term pair #: ', num2str(i), ' / ', num2str(length(HPO_GO_mapping_terms_unique_rand))])
    t1=HPO_GO_mapping_terms_pair_hist_rand(i,1);
    t2=HPO_terms_annot_hist_rand(HPO_GO_mapping_unique_array_ind_rand(i,1),1);
    t3=GO_terms_annot_hist_rand(HPO_GO_mapping_unique_array_ind_rand(i,2),1);
    HPO_GO_mapping_sematic_sim_high_det_rand(i,1:5)=[i (2*t1)/(t2+t3) t1 t2 t3];  % columns: indice of mapping, semantic similarity of the mapped terms, no of co-annotation of the mapped terms on different genes, total no of annotation of the mapped HPO term on different genes, total no of annotation of the mapped GO term on different genes
end
save HPO2GO_Files/HPO_GO_mapping_sematic_sim_high_det_rand.mat HPO_GO_mapping_sematic_sim_high_det_rand -v7

HPO_GO_random_mapping_merge=[HPO_GO_mapping_terms_unique_rand num2cell(HPO_GO_mapping_sematic_sim_high_det_rand(:,2:end))];
save HPO2GO_Files/HPO_GO_random_mapping_merge.mat HPO_GO_random_mapping_merge -v7

% Saving the randomized mapping file as text:

HPO_GO_random_mapping_merge_txt=HPO_GO_random_mapping_merge';
fid=fopen('HPO_GO_Random_Mapping.txt', 'w');
fprintf(fid, '%s\t%s\t%.4f\t%d\t%d\t%d\n', HPO_GO_random_mapping_merge_txt{:});
fclose(fid);

% Extracting the statistics for the randomized mapping - part of Table 1:

HPO_GO_mapping_sematic_sim_rand=HPO_GO_mapping_sematic_sim_high_det_rand(:,2);
S=1.1;
Mapping_stat_rand=0;
for i=1:10;
    S=S-0.1;
    Mapping_stat_rand(i,1)=S;
    Mapping_stat_rand(i,2)=length(find(HPO_GO_mapping_sematic_sim_rand>=S));
    Mapping_stat_rand(i,3)=length(unique(HPO_GO_mapping_terms_unique_rand(HPO_GO_mapping_sematic_sim_rand>=S,1)));
    Mapping_stat_rand(i,4)=length(unique(HPO_GO_mapping_terms_unique_rand(HPO_GO_mapping_sematic_sim_rand>=S,2)));
    Mapping_stat_rand(i,5)=length(find(ismember(HPO_gene_annotation_rand(:,2),unique(HPO_GO_mapping_terms_unique_rand(HPO_GO_mapping_sematic_sim_rand>=S,1)))==1));
    Mapping_stat_rand(i,6)=length(find(ismember(GO_annot_manual_human_rand(:,2),unique(HPO_GO_mapping_terms_unique_rand(HPO_GO_mapping_sematic_sim_rand>=S,2)))==1));
end
Mapping_table_title_rand={['Threshold'], ['# of unique mappings'], ['# of unique HPO terms'], ['# of unique GO terms'], ['Total # of annot with HPO term'], ['Total # of annot with GO term']};
Mapping_stat_cel_rand=num2cell(Mapping_stat_rand);
Mapping_table_rand=vertcat(Mapping_table_title_rand,Mapping_stat_cel_rand);
save HPO2GO_Files/HPO_sematic_sim_Mapping_stat_table_rand.mat Mapping_table_rand -v7

% Plotting the randomized co-occurrence similarity distributions for different numbers of co-annotated genes:

n=[1 5 25 75];
for i=1:4;
    P=HPO_GO_mapping_sematic_sim_high_det_rand(HPO_GO_mapping_sematic_sim_high_det_rand(:,3)>=(n(1,i)),2);
    figure;hist(P,0.00005:0.01:1.00005);
    hold on;
    x_lab=xlabel('Co-occurrence similarity bins');
    y_lab=ylabel('# of unique mappings');
    set(x_lab,'FontSize',14);
    set(y_lab,'FontSize',14);
    eval(['title(''Co-occurrence Similarity Distribution of Random Mapping (n >= ', num2str(n(1,i)), ')'',''FontSize'',16);'])
    grid on
    axis([0 0.3 0 (max(hist(P,0.00005:0.01:1.00005))*(11/10))])
    hold off;
end

% Plotting the number of mappings against threshold selections (normal vs. random) - Figure 4:

load HPO2GO_Files/HPO_GO_mapping_sematic_sim_high_det.mat
load HPO2GO_Files/HPO_GO_mapping_sematic_sim_high_det_rand.mat
n=0;
for j=1:5;
    n=n+1;
    S=1.1;
    for i=1:20;
        if S>0.11 % using 0.11 is crucial since matlab has a weird behaviour if 0.1 is used
            S=S-0.1;
        else
            S=S-0.01;
        end
        plo_nor(j,21-i)=length(find(HPO_GO_mapping_sematic_sim_high_det(HPO_GO_mapping_sematic_sim_high_det(:,2)>=S,3)>=n));
        plo_rand(j,21-i)=length(find(HPO_GO_mapping_sematic_sim_high_det_rand(HPO_GO_mapping_sematic_sim_high_det_rand(:,2)>=S,3)>=n));
    end
    thres=[0:0.01:0.09 0.1:0.1:1];
    figure;
    hold on;
    F=plot(thres,log(plo_nor(j,:)));
    set(F,'Color',[0.2 0.2 1],'LineWidth',2)
    R=plot(thres,log(plo_rand(j,:)));
    set(R,'Color',[1 0.2 0.2],'LineWidth',2)
    legendo=legend('Original','Randomized');
    set(legendo,'FontSize',14,'Location','northeast');
    x_lab=xlabel('Co-occurrence Similarity Threshold');
    y_lab=ylabel('Log ( total # of unique mappings )');
    set(x_lab,'FontSize',14);
    set(y_lab,'FontSize',14);
    eval(['title(''Cumulative # of Mappings - Original vs. Randomized Distributions (n >= ', num2str(n), ')'',''FontSize'',16);'])
    grid on
    axis([0 1 (min(log(plo_rand(j,:)))*(4/5)) (max(log(plo_nor(j,:)))*(6/5))])
    hold off;
end

% Significance test (kstest and ttest):

% (to differentiate between normal and random at different similarity thresholds (S) and number of co-occurences on different genes (n) )

sig_kstest=NaN(5,7);
sig_ttest=NaN(5,7);
n=0;
for j=1:5;
    n=n+1;
    S=-0.1;
    for i=1:7;
        S=S+0.1;
        dist1=HPO_GO_mapping_sematic_sim_high_det(HPO_GO_mapping_sematic_sim_high_det(:,2)>=S,:);
        dist1=sort(dist1(dist1(:,3)>=n,2));
        dist2=HPO_GO_mapping_sematic_sim_high_det_rand(HPO_GO_mapping_sematic_sim_high_det_rand(:,2)>=S,:);
        dist2=sort(dist2(dist2(:,3)>=n,2));
        
        dist1_hist=hist(dist1,20);
        dist2_hist=hist(dist2,20);
        dist1_hist(dist1_hist==0)=NaN;  % if not set to NaN zeros in both histograms are treated as if they are the same actual values and the p-value is incorrectly adjusted in these cases
        dist2_hist(dist2_hist==0)=NaN;
        if length(dist1)>10 && length(dist2)>10
            [~,sig_kstest(j,i)]=kstest2(dist1_hist,dist2_hist);
            [~,sig_ttest(j,i)]=ttest2(dist1_hist,dist2_hist,'Vartype','unequal');
        else
            sig_kstest(j,i)=NaN;
            sig_ttest(j,i)=NaN;
        end
    end
end
save HPO2GO_Files/HPO_KStest_result_table.mat sig_kstest -v7
save HPO2GO_Files/HPO_ttest_result_table.mat sig_ttest -v7

% (checking the number of unique HPO and GO terms we used at different thresholds)

n=0;
for j=1:10;
    n=n+1;
    S=1.1;
    for i=1:20
        if S>0.1001 % using 0.11 is crucial since matlab has a weird behaviour if 0.1 is used
            S=S-0.1;
        else
            S=S-0.01;
        end
        HPO_GO_mapping_sematic_sim_high_det_thres=HPO_GO_mapping_sematic_sim_high_det(HPO_GO_mapping_sematic_sim_high_det(:,2)>S,:);
        HPO_GO_mapping_sematic_sim_high_det_thres_membnum=HPO_GO_mapping_sematic_sim_high_det_thres(HPO_GO_mapping_sematic_sim_high_det_thres(:,3)>=n,:);
        HPO_GO_mapping_terms_unique_thres_membnum=HPO_GO_mapping_terms_unique(HPO_GO_mapping_sematic_sim_high_det_thres_membnum(:,1),:);
        num_HPO_terms(n,21-i)=length(unique(HPO_GO_mapping_terms_unique_thres_membnum(:,1)));
        num_GO_terms(n,21-i)=length(unique(HPO_GO_mapping_terms_unique_thres_membnum(:,2)));
    end
end

% (considering sig_Sstest result: n>=2 -observed number of co-annotated genes- and 
% S>=0.1 -co-occurrence similarity thereshold- was observed to be sufficient to assume significance)

% Saving the finalized HPO2GO mapping file:

n=2;
S=0.1;

load HPO2GO_Files/HPO_GO_mapping_sematic_sim_high_det.mat
load HPO2GO_Files/HPO_GO_mapping_terms.mat
load HPO2GO_Files/HPO_GO_mapping_terms_unique.mat
high_det_thres=HPO_GO_mapping_sematic_sim_high_det(HPO_GO_mapping_sematic_sim_high_det(:,2)>=S,:);
high_det_thres=high_det_thres(high_det_thres(:,3)>=n,:);
HPO2GO_Finalized_Mapping=HPO_GO_mapping_terms_unique(high_det_thres(:,1),:);
HPO2GO_Finalized_Mapping(:,3:4)=num2cell(HPO_GO_mapping_sematic_sim_high_det(high_det_thres(:,1),2:3));

% (deleting the mappings to HPO terms that are not belong to phenotypic abnormality sub-ontology)

[HPO_ID_not_abnor,HPO_Name_not_abnor]=textread('HPO2GO_Files/HPO_terms_not_abnormality.txt', '%s %s', 'delimiter', '\t');
[Lia,Locb]=ismember(HPO2GO_Finalized_Mapping(:,1),HPO_ID_not_abnor);
HPO2GO_Finalized_Mapping(Lia==1,:)=[];

save HPO2GO_Files/HPO2GO_Finalized_Mapping.mat HPO2GO_Finalized_Mapping -v7
length(unique(HPO2GO_Finalized_Mapping(:,1)))
length(unique(HPO2GO_Finalized_Mapping(:,2)))

figure;
hold on;
hist(cell2mat(HPO2GO_Finalized_Mapping(:,3)),0:0.01:1)
axis([0 1 0 max(hist(cell2mat(HPO2GO_Finalized_Mapping(:,3)),0:0.01:1))*1.1])
hold off;
figure;
hold on;
hist(HPO_GO_mapping_sematic_sim_high_det(:,2),0:0.01:1)
axis([0 1 0 max(hist(HPO_GO_mapping_sematic_sim_high_det(:,2),0:0.01:1))*1.1])
hold off;
figure;
hold on;
hist(cell2mat(HPO2GO_Finalized_Mapping(:,4)),0.0005:1:max(cell2mat(HPO2GO_Finalized_Mapping(:,4)))+0.0005)
axis([0 30 0 max(hist(cell2mat(HPO2GO_Finalized_Mapping(:,4)),0.0005:1:max(cell2mat(HPO2GO_Finalized_Mapping(:,4)))+0.0005))*1.1])
hold off;
figure;
hold on;
hist(HPO_GO_mapping_sematic_sim_high_det(:,3),0.0005:1:max(HPO_GO_mapping_sematic_sim_high_det(:,3))+0.0005)
axis([0 30 0 max(hist(HPO_GO_mapping_sematic_sim_high_det(:,3),0.0005:1:max(HPO_GO_mapping_sematic_sim_high_det(:,3))+0.0005))*1.1])
hold off;

% Saving the finalized HPO2GO mapping file as text:

HPO2GO_Finalized_Mapping_txt=HPO2GO_Finalized_Mapping';
fid=fopen('HPO2GO_Finalized_Mapping.txt', 'w');
fprintf(fid, '%s\t%s\t%.4f\t%d\n', HPO2GO_Finalized_Mapping_txt{:});
fclose(fid);


% Comparing HPO2GO mappings with the manually generated HPO to GO associations:

[HPO_manual_GO_HPOID,HPO_manual_GO_GOID]=textread('HPO2GO_Files/HPO_manual_GO_associations_03_2018.txt', '%s %s', 'delimiter', '\t');
HPO_manual_GO=[HPO_manual_GO_HPOID HPO_manual_GO_GOID];
length(HPO_manual_GO)
length(unique(HPO_manual_GO_HPOID))
length(unique(HPO_manual_GO_GOID))
HPO_manual_GO=uniqueRowsCA(HPO_manual_GO);
HPO_manual_GO_txt=HPO_manual_GO';
fid=fopen('HPO2GO_Files/HPO_manual_GO_associations_03_2018.txt', 'w');
fprintf(fid, '%s\t%s\n', HPO_manual_GO_txt{:});
fclose(fid);
save HPO2GO_Files/HPO_manual_GO_associations.mat HPO_manual_GO -v7
[Lia,Locb]=ismember(HPO_manual_GO(:,2),unique(HPO_manual_GO(:,2)));
Locb_hist=hist(Locb,0.5:1:(max(Locb)+0.5));
Locb_hist(Locb_hist>10)
HPO_manual_GO_unique=unique(HPO_manual_GO(:,2));
HPO_manual_GO_unique(find(Locb_hist>10),1)

load HPO2GO_Files/HPO2GO_Finalized_Mapping.mat
[Lia,Locb]=ismember(HPO_manual_GO(:,1),HPO2GO_Finalized_Mapping(:,1));
sum(Lia)
length(unique(HPO_manual_GO(Lia==1,1)))
[Lia,Locb]=ismember(HPO_manual_GO(:,2),HPO2GO_Finalized_Mapping(:,2));
sum(Lia)
length(unique(HPO_manual_GO(Lia==1,2)))
HPO_manual_GO_merged=strcat(HPO_manual_GO(:,1),{' '},HPO_manual_GO(:,2));
HPO2GO_Finalized_Mapping_merged=strcat(HPO2GO_Finalized_Mapping(:,1),{' '},HPO2GO_Finalized_Mapping(:,2));
[Lia,Locb]=ismember(HPO_manual_GO_merged,HPO2GO_Finalized_Mapping_merged);
sum(Lia)


% Predicting HPO terms for CAFA3 human target proteins:

load HPO2GO_Files/HPO_GO_mapping_sematic_sim_high_det.mat
load HPO2GO_Files/HPO_GO_mapping_terms_unique.mat
load HPO2GO_Files/CAFA3_Targets_human_all_mappings.mat
load HPO2GO_Files/CAFA3_HPO_target_variables.mat
load HPO2GO_Files/Humprot_mapping_variables.mat
load HPO2GO_Files/HPO_GO_mapping_terms.mat
load HPO2GO_Files/HPO_test_GO_annotation.mat
load HPO2GO_Files/HPO_test_HPO_gene_annotation.mat
load HPO2GO_Files/HPO_test_GO_annotation_all.mat
load HPO2GO_Files/HPO2GO_Finalized_Mapping.mat

mapGO_unique=unique(HPO2GO_Finalized_Mapping(:,2));
[~,Locb0]=ismember(HPO2GO_Finalized_Mapping(:,2),mapGO_unique);
[Lia,Locb]=ismember(GO_annot_manual_human_all(:,2),mapGO_unique);
to=0;
CAFA3_HPO_predictions=cell(2000000,3);
for i=1:length(mapGO_unique);
    disp(['Predicting HPO terms for the mapped GO term #: ', num2str(i), ' / ', num2str(length(mapGO_unique))])
    genesymbol_temp=GO_annot_manual_human_all(Locb==i,1);
    HPO_temp=HPO2GO_Finalized_Mapping(Locb0==i,1);
    score_temp=HPO2GO_Finalized_Mapping(Locb0==i,3);
    if sum(Lia)>0
        siz=length(genesymbol_temp)*length(HPO_temp);
        [A,B]=meshgrid(1:length(genesymbol_temp),1:length(HPO_temp));
        c=cat(2,A',B');
        d=reshape(c,[],2);
        CAFA3_HPO_predictions(to+1:to+siz,1)=genesymbol_temp(d(:,1),1);
        CAFA3_HPO_predictions(to+1:to+siz,2)=HPO_temp(d(:,2),1);
        CAFA3_HPO_predictions(to+1:to+siz,3)=score_temp(d(:,2),1);
        to=to+siz;
    end
end
CAFA3_HPO_predictions(cellfun(@isempty,CAFA3_HPO_predictions(:,1))==1,:)=[];
[CAFA3_HPO_predictions_un,I,C]=uniqueRowsCA(CAFA3_HPO_predictions(:,1:2));
% score_unique=zeros(length(CAFA3_HPO_predictions_un),1); % Alternative 1) Extremely slow
% for i=1:length(CAFA3_HPO_predictions_un)
%     disp(['Calculating the max score for the prediction #: ', num2str(i), ' / ', num2str(length(CAFA3_HPO_predictions_un))])
%     score_unique(i,1)=max(cell2mat(CAFA3_HPO_predictions(C==i,3)));
% end
% score=cell2mat(CAFA3_HPO_predictions(:,3));
% score_unique=arrayfun(@(x) max(score(C==x,1)), 1:max(C)); % Alternative 2) Slow
score=cell2mat(CAFA3_HPO_predictions(:,3));
score_unique=splitapply(@max,score,C);  % Alternative 3) Very fast but requires Matlab 2017b
CAFA3_HPO_predictions_un(:,3)=num2cell(score_unique);
CAFA3_HPO_predictions=CAFA3_HPO_predictions_un;
CAFA3_HPO_predictions(ismember(CAFA3_HPO_predictions(:,1),'2xchrna4-3xchrnb2_human')==1,1)=cellstr('CHRNA4');
CAFA3_HPO_predictions(ismember(CAFA3_HPO_predictions(:,1),'3xchrna4-2xchrnb2_human')==1,1)=cellstr('CHRNA4');

mapsymbol_unique=unique(CAFA3_HPO_predictions(:,1));
[Lia0,Locb0]=ismember(CAFA3_HPO_predictions(:,1),mapsymbol_unique);
[Lia,Locb]=ismember(mapsymbol_unique,CAFA3_Targets_human_all_mappings(:,4));
Locb(Locb==0,1)=length(CAFA3_Targets_human_all_mappings)+1;
CAFA3_Targets_human_all_mappings(end+1,:)=cellstr(' ');
mapCAFAid_unique=CAFA3_Targets_human_all_mappings(Locb,1);
CAFA3_HPO_predictions_CAFAid_genesymbol_HPOid_sco=mapCAFAid_unique(Locb0,1);
CAFA3_HPO_predictions_CAFAid_genesymbol_HPOid_sco(:,2:4)=CAFA3_HPO_predictions;
save HPO2GO_Files/CAFA3_HPO_predictions_CAFAid_genesymbol_HPOid_sco.mat CAFA3_HPO_predictions_CAFAid_genesymbol_HPOid_sco -v7
CAFA3_pred_eval_save=CAFA3_HPO_predictions_CAFAid_genesymbol_HPOid_sco;
CAFA3_pred_eval_save(cellfun(@isempty,CAFA3_pred_eval_save(:,1))==1,:)=[];
CAFA3_pred_eval_save=uniqueRowsCA(CAFA3_pred_eval_save);
save HPO2GO_Files/CAFA3_HPO_predictions_semantic.mat CAFA3_pred_eval_save -v7
CAFA3_HPO_target_predictions=CAFA3_pred_eval_save(:,[1 3 4]);
save HPO2GO_Files/CAFA3_HPO_target_predictions.mat CAFA3_HPO_target_predictions -v7
length(CAFA3_HPO_target_predictions)
length(unique(CAFA3_HPO_target_predictions(:,1)))
length(unique(CAFA3_HPO_target_predictions(:,2)))

% Saving CAFA3 HPO target predictions in a text file:

CAFA3_pred_eval_save_txt=CAFA3_pred_eval_save(:,[1 3 4])';
fid=fopen('CAFA3_HPO_target_predictions.txt', 'w');
fprintf(fid, '%s\t%s\t%.2f\n', CAFA3_pred_eval_save_txt{:});
fclose(fid);



% Generating the finalized HPO2protein predictions:

load HPO2GO_Files/HPO_test_GO_annotation_all_uniprotacc.mat
[Lia,Locb]=ismember(CAFA3_HPO_predictions_CAFAid_genesymbol_HPOid_sco(:,2),GO_annot_manual_human_all_uniprotacc(:,2));
HPO2protein_predictions=GO_annot_manual_human_all_uniprotacc(Locb,1);
HPO2protein_predictions(:,2:3)=CAFA3_HPO_predictions_CAFAid_genesymbol_HPOid_sco(:,3:4);
Lia=strncmp('EBI-',HPO2protein_predictions(:,1),4);
HPO2protein_predictions(Lia==1,:)=[];

% (deleting the predictions with terms that are not belong to phenotypic abnormality sub-ontology)

[HPO_ID_not_abnor,~]=textread('HPO2GO_Files/HPO_terms_not_abnormality.txt', '%s %s', 'delimiter', '\t');
[Lia,~]=ismember(HPO2protein_predictions(:,2),HPO_ID_not_abnor);
HPO2protein_predictions(Lia==1,:)=[];

save HPO2GO_Files/HPO2protein_predictions.mat HPO2protein_predictions -v7
length(HPO2protein_predictions)
length(unique(HPO2protein_predictions(:,1)))
length(unique(HPO2protein_predictions(:,2)))

% Saving HPO2protein predictions in a text file:

HPO2protein_predictions_txt=HPO2protein_predictions';
fid=fopen('HPO2protein_predictions.txt', 'w');
fprintf(fid, '%s\t%s\t%.4f\n', HPO2protein_predictions_txt{:});
fclose(fid);

% Calculating how many of the targets already have disease annotation in UniProt:

list=unique(HPO2protein_predictions(:,1));
dlmcell('hpo_list.txt',list)
% (load the list to UniProt ID mapping and obtain the disease annotation file)
[CAFA3_HPO_UniProt_acc,CAFA3_HPO_gene_name,CAFA3_HPO_disease_annot,CAFA3_HPO_orphanet_annot]=textread('CAFA3_HPO_predicted_target_list_w_disease_annot.tab', '%s %s %s %s', 'delimiter', '\t', 'headerlines', 1, 'bufsize', 10000);
CAFA3_HPO_disease_annot_sort=sort(CAFA3_HPO_disease_annot);
CAFA3_HPO_disease_annot_sort(cellfun(@isempty,CAFA3_HPO_disease_annot_sort(:,1))==1,:)=[];
length(list)-length(CAFA3_HPO_disease_annot_sort);
% (14545 out of 18134 proteins with predictions have no existing disease annotation in UniProt)



% Performance test on CAFA2 targets

% (download CAFA2 files from: https://figshare.com/articles/Supplementary_Data_for_CAFA2/2059944)

[CAFA2_HPO_ground_geneCAFAid,CAFA2_HPO_ground_HPOid]=textread('CCAFA2/CAFA2-master/benchmark/groundtruth/propagated_HPO.txt', '%s %s', 'delimiter', '\t');
CAFA2_HPO_ground_geneCAFAid_unique=unique(CAFA2_HPO_ground_geneCAFAid);

[CAFA2_Targets_human_headers,~]=fastaread('CAFA2/CAFA2_Supplementary_data/data/CAFA2-targets/eukarya/sp_species.9606.tfa');
CAFA2_Targets_human_headers=CAFA2_Targets_human_headers';
CAFA2_Targets_human_idmapping=cellfun(@(x) strsplit(x,' '), CAFA2_Targets_human_headers(:),'uni',0);
CAFA2_Targets_human_idmapping=vertcat(CAFA2_Targets_human_idmapping{:,1});
[Lia,~]=ismember(CAFA2_Targets_human_idmapping(:,1),CAFA2_HPO_ground_geneCAFAid_unique);
save HPO2GO_Files/CAFA2_HPO_target_variables.mat CAFA2_Targets_human_headers CAFA2_Targets_human_idmapping -v7
CAFA2_HPO_ground_proteinname=CAFA2_Targets_human_idmapping(Lia==1,2);

load HPO2GO_Files/Humprot_mapping_variables.mat
[~,Locb]=ismember(CAFA2_Targets_human_idmapping(:,2),Humprot_entry_name);
Locb(Locb==0,1)=length(Humprot_entry_name)+1;
Humprot_entry_name(end+1,1)=cellstr(' ');
Humprot_UniProt_acc(end+1,1)=cellstr(' ');
Humprot_gene_symbol(end+1,1)=cellstr(' ');
CAFA2_Targets_human_all_mappings=CAFA2_Targets_human_idmapping;
CAFA2_Targets_human_all_mappings(:,3)=Humprot_UniProt_acc(Locb,1);
CAFA2_Targets_human_all_mappings(:,4)=Humprot_gene_symbol(Locb,1);
save HPO2GO_Files/CAFA2_Targets_human_all_mappings.mat CAFA2_Targets_human_all_mappings -v7

% (id mapping between CAFA2_HPO_ground_proteinname and associated gene symbols, loaded as: CAFA2_HPO_ground_proteinname_genesymbol)

[CAFA2_HPO_ground_proteinname_genesymbol(:,1),idx]=sort(CAFA2_HPO_ground_proteinname_genesymbol(:,1));
CAFA2_HPO_ground_proteinname_genesymbol(:,2)=CAFA2_HPO_ground_proteinname_genesymbol(idx,2);
isequal(CAFA2_HPO_ground_proteinname_genesymbol(:,1),CAFA2_HPO_ground_proteinname(:,1))
[~,Locb]=ismember(CAFA2_HPO_ground_geneCAFAid,CAFA2_HPO_ground_geneCAFAid_unique);
CAFA2_HPO_ground_genesymbol=CAFA2_HPO_ground_proteinname_genesymbol(Locb,2);
CAFA2_HPO_ground_annot_genesymbol_HPOid=CAFA2_HPO_ground_genesymbol;
CAFA2_HPO_ground_annot_genesymbol_HPOid(:,2)=CAFA2_HPO_ground_HPOid;
CAFA2_HPO_ground_annot_CAFAid_genesymbol_HPOid=CAFA2_HPO_ground_geneCAFAid;
CAFA2_HPO_ground_annot_CAFAid_genesymbol_HPOid(:,2:3)=CAFA2_HPO_ground_annot_genesymbol_HPOid;
save HPO2GO_Files/CAFA2_HPO_ground_variables.mat CAFA2_HPO_ground_annot_CAFAid_genesymbol_HPOid CAFA2_HPO_ground_annot_genesymbol_HPOid CAFA2_HPO_ground_geneCAFAid CAFA2_HPO_ground_geneCAFAid_unique CAFA2_HPO_ground_genesymbol CAFA2_HPO_ground_HPOid CAFA2_HPO_ground_proteinname CAFA2_HPO_ground_proteinname_genesymbol -v7

% (loading CAFA2 trainnig set for HPO)

[CAFA2_HPO_training_UniProtid,CAFA2_HPO_training_HPOid]=textread('CAFA2/CAFA2_Supplementary_data/data/HPO-t0/hpoa.hp', '%s %s', 'delimiter', '\t');
CAFA2_HPO_training_UniProtid_unique=unique(CAFA2_HPO_training_UniProtid);

% (id mapping between CAFA2_HPO_training_UniProtid_unique and associated gene symbols, loaded as: CAFA2_HPO_training_UniProtid_genesymbol)

[CAFA2_HPO_training_UniProtid_genesymbol(:,1),idx]=sort(CAFA2_HPO_training_UniProtid_genesymbol(:,1));
CAFA2_HPO_training_UniProtid_genesymbol(:,2)=CAFA2_HPO_training_UniProtid_genesymbol(idx,2);
[Lia,Locb]=ismember(CAFA2_HPO_training_UniProtid,CAFA2_HPO_training_UniProtid_genesymbol(:,1));
CAFA2_HPO_training_UniProtid_onlywithgenesym=CAFA2_HPO_training_UniProtid_genesymbol(Locb(Locb>0),1);
CAFA2_HPO_training_genesymbol_onlywithgenesym=CAFA2_HPO_training_UniProtid_genesymbol(Locb(Locb>0),2);
CAFA2_HPO_training_HPOid_onlywithgenesym=CAFA2_HPO_training_HPOid(Locb>0,1);
CAFA2_HPO_training_annot_genesymbol_HPOid=CAFA2_HPO_training_genesymbol_onlywithgenesym;
CAFA2_HPO_training_annot_genesymbol_HPOid(:,2)=CAFA2_HPO_training_HPOid_onlywithgenesym;
CAFA2_HPO_training_annot_genesymbol_HPOid_unique=uniqueRowsCA(CAFA2_HPO_training_annot_genesymbol_HPOid);
save HPO2GO_Files/HPO_CAFA2_training_variables.mat CAFA2_HPO_training_UniProtid CAFA2_HPO_training_HPOid CAFA2_HPO_training_UniProtid_unique CAFA2_HPO_training_UniProtid_genesymbol CAFA2_HPO_training_UniProtid_onlywithgenesym CAFA2_HPO_training_genesymbol_onlywithgenesym CAFA2_HPO_training_HPOid_onlywithgenesym CAFA2_HPO_training_annot_genesymbol_HPOid CAFA2_HPO_training_annot_genesymbol_HPOid_unique -v7

CAFA2_HPO_training_annot_genesymbol_unique=unique(CAFA2_HPO_training_annot_genesymbol_HPOid_unique(:,1));

[GO_UniProt_acc_all,GO_gene_symbol_all,GO_term_id_all,GO_term_evid_all]=textread('HPO2GO_Files/GOA_UniProt_human_protein_annotation_2014_01.txt', '%s %s %s %s', 'delimiter', '\t');
ind_IEA=find(strcmp(GO_term_evid_all,'IEA')==1);
GO_annot_manual_human_all=[GO_gene_symbol_all GO_term_id_all];
GO_annot_manual_human_all(ind_IEA,:)=[];
GO_annot_manual_human_all=uniqueRowsCA(GO_annot_manual_human_all);
save HPO2GO_Files/HPO_test_GO_annotation_all.mat GO_annot_manual_human_all -v7
GO_annot_manual_human_all_addcol=[GO_gene_symbol_all GO_term_id_all GO_term_evid_all];
GO_annot_manual_human_all_addcol(ind_IEA,:)=[];
GO_annot_manual_human_all_addcol=uniqueRowsCA(GO_annot_manual_human_all_addcol);
save HPO2GO_Files/HPO_test_GO_annotation_all_addcol.mat GO_annot_manual_human_all_addcol -v7

[Lia,Locb]=ismember(GO_annot_manual_human_all(:,1),CAFA2_HPO_training_annot_genesymbol_unique);
CAFA2_GO_annot_manual_human=GO_annot_manual_human_all(Lia==1,:);
save HPO2GO_Files/CAFA2_HPO_test_GO_annotation.mat CAFA2_GO_annot_manual_human -v7

to=0;
CAFA2_HPO_GO_mapping_ind=zeros(2000000,2);
for i=1:length(CAFA2_HPO_training_annot_genesymbol_unique);
    disp(['Mapping HPO terms with GO terms with gene #: ', num2str(i), ' / ', num2str(length(CAFA2_HPO_training_annot_genesymbol_unique))])
    [Lia,Locb]=ismember(CAFA2_HPO_training_annot_genesymbol_HPOid_unique(:,1),CAFA2_HPO_training_annot_genesymbol_unique(i,1));
    [Lia2,Locb2]=ismember(CAFA2_GO_annot_manual_human(:,1),CAFA2_HPO_training_annot_genesymbol_unique(i,1));
    if sum(Lia2)>0
        siz=sum(Lia)*sum(Lia2);
        [A,B]=meshgrid(find(Lia==1),find(Lia2==1));
        c=cat(2,A',B');
        d=reshape(c,[],2);
        CAFA2_HPO_GO_mapping_ind(to+1:to+siz,:)=d;
        to=to+siz;
    end
end
CAFA2_HPO_GO_mapping_ind(CAFA2_HPO_GO_mapping_ind(:,1)==0,:)=[];
CAFA2_HPO_GO_mapping_terms=CAFA2_HPO_training_annot_genesymbol_HPOid_unique(CAFA2_HPO_GO_mapping_ind(:,1),2);
CAFA2_HPO_GO_mapping_terms(:,2)=CAFA2_GO_annot_manual_human(CAFA2_HPO_GO_mapping_ind(:,2),2);
save HPO2GO_Files/CAFA2_HPO_GO_mapping_terms.mat CAFA2_HPO_GO_mapping_terms -v7

[CAFA2_HPO_GO_mapping_terms_unique,I,J]=uniqueRowsCA(CAFA2_HPO_GO_mapping_terms);
CAFA2_HPO_GO_mapping_terms_freq=sort(J);
CAFA2_HPO_GO_mapping_terms_pair_hist=hist(CAFA2_HPO_GO_mapping_terms_freq,0.5:1:max(CAFA2_HPO_GO_mapping_terms_freq)+0.5);
CAFA2_HPO_GO_mapping_terms_pair_hist=CAFA2_HPO_GO_mapping_terms_pair_hist(1,1:end-1)';
save HPO2GO_Files/CAFA2_HPO_GO_mapping_terms_unique.mat CAFA2_HPO_GO_mapping_terms_unique -v7
save HPO2GO_Files/CAFA2_HPO_GO_mapping_terms_pair_hist.mat CAFA2_HPO_GO_mapping_terms_pair_hist -v7

[CAFA2_HPO_terms_unique,~,N]=unique(CAFA2_HPO_training_annot_genesymbol_HPOid_unique(:,2));
CAFA2_HPO_terms_annot_hist=hist(N,0.5:1:max(N)+0.5);
CAFA2_HPO_terms_annot_hist=CAFA2_HPO_terms_annot_hist(1,1:end-1)';
save HPO2GO_Files/CAFA2_HPO_terms_unique.mat CAFA2_HPO_terms_unique -v7
save HPO2GO_Files/CAFA2_HPO_terms_annot_hist.mat CAFA2_HPO_terms_annot_hist -v7

[CAFA2_GO_terms_unique,~,N]=unique(CAFA2_GO_annot_manual_human(:,2));
CAFA2_GO_terms_annot_hist=hist(N,0.5:1:max(N)+0.5);
CAFA2_GO_terms_annot_hist=CAFA2_GO_terms_annot_hist(1,1:end-1)';
save HPO2GO_Files/CAFA2_GO_terms_unique.mat CAFA2_GO_terms_unique -v7
save HPO2GO_Files/CAFA2_GO_terms_annot_hist.mat CAFA2_GO_terms_annot_hist -v7

[~,CAFA2_HPO_GO_mapping_unique_array_ind]=ismember(CAFA2_HPO_GO_mapping_terms_unique(:,1),CAFA2_HPO_terms_unique);
[~,CAFA2_HPO_GO_mapping_unique_array_ind(:,2)]=ismember(CAFA2_HPO_GO_mapping_terms_unique(:,2),CAFA2_GO_terms_unique);

CAFA2_HPO_GO_mapping_sematic_sim_high_det=zeros(length(CAFA2_HPO_GO_mapping_terms_unique),5);
for i=1:length(CAFA2_HPO_GO_mapping_terms_unique);
    disp(['Calculating the semantic similarity between the HPO ang GO term pair #: ', num2str(i), ' / ', num2str(length(CAFA2_HPO_GO_mapping_terms_unique))])
    t1=CAFA2_HPO_GO_mapping_terms_pair_hist(i,1);
    t2=CAFA2_HPO_terms_annot_hist(CAFA2_HPO_GO_mapping_unique_array_ind(i,1),1);
    t3=CAFA2_GO_terms_annot_hist(CAFA2_HPO_GO_mapping_unique_array_ind(i,2),1);
    CAFA2_HPO_GO_mapping_sematic_sim_high_det(i,1:5)=[i (2*t1)/(t2+t3) t1 t2 t3]; % columns: indice of mapping, semantic similarity of the mapped terms, no of co-annotation of the mapped terms on different genes, total no of annotation of the mapped HPO term on different genes, total no of annotation of the mapped GO term on different genes
end
save HPO2GO_Files/CAFA2_HPO_GO_mapping_sematic_sim_high_det.mat CAFA2_HPO_GO_mapping_sematic_sim_high_det -v7

% (using the parameters: n>=2 and S>=0.1 as previously obtained from randomization test)

n=2;
S=0.1;

high_det_thres=CAFA2_HPO_GO_mapping_sematic_sim_high_det(CAFA2_HPO_GO_mapping_sematic_sim_high_det(:,2)>=S,:);
high_det_thres=high_det_thres(high_det_thres(:,3)>=n,:);
CAFA2_HPO_GO_selected_mappings=CAFA2_HPO_GO_mapping_terms_unique(high_det_thres(:,1),:);
CAFA2_HPO_GO_selected_mappings(:,3)=num2cell(CAFA2_HPO_GO_mapping_sematic_sim_high_det(high_det_thres(:,1),2));
save HPO2GO_Files/CAFA2_HPO_GO_selected_mappings.mat CAFA2_HPO_GO_selected_mappings -v7
length(CAFA2_HPO_GO_selected_mappings)
length(unique(CAFA2_HPO_GO_selected_mappings(:,1)))
length(unique(CAFA2_HPO_GO_selected_mappings(:,2)))

% Predicting HPO terms for CAFA2 human target proteins:

% (the variable containing the ground-truth / benchmark set with respective HPO annotations: CAFA2_HPO_ground_annot_CAFAid_genesymbol_HPOid)

load HPO2GO_Files/CAFA2_HPO_test_GO_annotation_all.mat
load HPO2GO_Files/CAFA2_HPO_GO_selected_mappings.mat
load HPO2GO_Files/CAFA2_HPO_GO_mapping_terms_unique.mat
mapGO_unique=unique(CAFA2_HPO_GO_selected_mappings(:,2));
[~,Locb0]=ismember(CAFA2_HPO_GO_selected_mappings(:,2),mapGO_unique);
[Lia,Locb]=ismember(GO_annot_manual_human_all(:,2),mapGO_unique);
to=0;
CAFA2_HPO_predictions=cell(2000000,3);
for i=1:length(mapGO_unique);
    disp(['Predicting HPO terms for the mapped GO term #: ', num2str(i), ' / ', num2str(length(mapGO_unique))])
    genesymbol_temp=GO_annot_manual_human_all(Locb==i,1);
    HPO_temp=CAFA2_HPO_GO_selected_mappings(Locb0==i,1);
    score_temp=CAFA2_HPO_GO_selected_mappings(Locb0==i,3);
    if sum(Lia)>0
        siz=length(genesymbol_temp)*length(HPO_temp);
        [A,B]=meshgrid(1:length(genesymbol_temp),1:length(HPO_temp));
        c=cat(2,A',B');
        d=reshape(c,[],2);
        CAFA2_HPO_predictions(to+1:to+siz,1)=genesymbol_temp(d(:,1),1);
        CAFA2_HPO_predictions(to+1:to+siz,2)=HPO_temp(d(:,2),1);
        CAFA2_HPO_predictions(to+1:to+siz,3)=score_temp(d(:,2),1);
        to=to+siz;
    end
end
CAFA2_HPO_predictions(cellfun(@isempty,CAFA2_HPO_predictions(:,1))==1,:)=[];
[CAFA2_HPO_predictions_un,I,C]=uniqueRowsCA(CAFA2_HPO_predictions(:,1:2));
% score_unique=zeros(length(CAFA2_HPO_predictions_un),1); % Alternative 1) Extremely slow
% for i=1:length(CAFA2_HPO_predictions_un)
%     disp(['Calculating the max score for the prediction #: ', num2str(i), ' / ', num2str(length(CAFA2_HPO_predictions_un))])
%     score_unique(i,1)=max(cell2mat(CAFA2_HPO_predictions(C==i,3)));
% end
% score=cell2mat(CAFA2_HPO_predictions(:,3));
% score_unique=arrayfun(@(x) max(score(C==x,1)), 1:max(C)); - Alternative 2) Slow
score=cell2mat(CAFA2_HPO_predictions(:,3));
score_unique=splitapply(@max,score,C);  % Alternative 3) Very fast but requires Matlab 2017b
CAFA2_HPO_predictions_un(:,3)=num2cell(score_unique);
CAFA2_HPO_predictions=CAFA2_HPO_predictions_un;
CAFA2_HPO_predictions(ismember(CAFA2_HPO_predictions(:,1),'2xchrna4-3xchrnb2_human')==1,1)=cellstr('CHRNA4');
CAFA2_HPO_predictions(ismember(CAFA2_HPO_predictions(:,1),'3xchrna4-2xchrnb2_human')==1,1)=cellstr('CHRNA4');

load HPO2GO_Files/CAFA2_Targets_human_all_mappings.mat
mapsymbol_unique=unique(CAFA2_HPO_predictions(:,1));
[Lia0,Locb0]=ismember(CAFA2_HPO_predictions(:,1),mapsymbol_unique);
[Lia,Locb]=ismember(mapsymbol_unique,CAFA2_Targets_human_all_mappings(:,4));
Locb(Locb==0,1)=length(CAFA2_Targets_human_all_mappings)+1;
CAFA2_Targets_human_all_mappings(end+1,:)=cellstr(' ');
mapCAFAid_unique=CAFA2_Targets_human_all_mappings(Locb,1);
CAFA2_HPO_predictions_CAFAid_genesymbol_HPOid_sco=mapCAFAid_unique(Locb0,1);
CAFA2_HPO_predictions_CAFAid_genesymbol_HPOid_sco(:,2:4)=CAFA2_HPO_predictions;
save HPO2GO_Files/CAFA2_HPO_predictions_CAFAid_genesymbol_HPOid_sco.mat CAFA2_HPO_predictions_CAFAid_genesymbol_HPOid_sco -v7
CAFA2_pred_eval_save=CAFA2_HPO_predictions_CAFAid_genesymbol_HPOid_sco;
CAFA2_pred_eval_save(cellfun(@isempty,CAFA2_pred_eval_save(:,1))==1,:)=[];

CAFA2_pred_eval_save134=CAFA2_pred_eval_save(:,[1 3 4]);
CAFA2_pred_eval_save_txt=CAFA2_pred_eval_save134';
fid=fopen('HPO2GO_Files/CAFA2_HPO_target_predictions.txt', 'w');
fprintf(fid, '%s\t%s\t%.4f\n', CAFA2_pred_eval_save_txt{:});
fclose(fid);

CAFA2_HPO_target_predictions=CAFA2_pred_eval_save(:,[1 3 4]);
save HPO2GO_Files/CAFA2_HPO_target_predictions.mat CAFA2_HPO_target_predictions -v7



% Performance test on CAFA2 targets using propagated HPO annotations and GO annotation with all evidence codes (HPOprop2GOall)

% Propagating the HPO annotations to the root:

[CAFA2_HPO_training_UniProtid,CAFA2_HPO_training_HPOid]=textread('CAFA2/CAFA2_Supplementary_data/data/HPO-t0/hpoa.hp', '%s %s', 'delimiter', '\t');
CAFA2_HPO_training=[CAFA2_HPO_training_UniProtid CAFA2_HPO_training_HPOid];
CAFA2_HPO_training_unique=uniqueRowsCA(CAFA2_HPO_training);
CAFA2_HPO_training_UniProtid_unique=unique(CAFA2_HPO_training_UniProtid);
CAFA2_HPO_training_HPOid_unique=unique(CAFA2_HPO_training_HPOid);
ont=pfp_ontbuild('CAFA2/CAFA2-master/ontology/hpo_v810-17-SEP-2013.obo');
term_list=cellstr(vertcat(ont.term.id));
[Lia,Locb]=ismember(CAFA2_HPO_training_HPOid,term_list);
[Lia2,Locb2]=ismember(CAFA2_HPO_training_UniProtid,CAFA2_HPO_training_UniProtid_unique);
source_annot_mat=zeros(length(CAFA2_HPO_training_UniProtid_unique),length(term_list));
row_nonz=Locb2(Locb>0);
col_nonz=Locb(Locb>0);
for i=1:length(row_nonz);
    disp(['Ind: ', num2str(i), ' / ', num2str(length(row_nonz))])
    source_annot_mat(row_nonz(i,1),col_nonz(i,1))=1;
end
source_annot_mat_log=logical(source_annot_mat);
A=pfp_annotprop(ont.DAG, source_annot_mat_log);
[A_ind_row,A_ind_col]=find(A==1);
CAFA2_HPO_training_prop=CAFA2_HPO_training_UniProtid_unique(A_ind_row,1);
CAFA2_HPO_training_prop(:,2)=term_list(A_ind_col,1);
CAFA2_HPO_training_prop(strcmp(CAFA2_HPO_training_prop(:,2),'HP:0000001')==1,:)=[];
CAFA2_HPO_training_propagated=vertcat(CAFA2_HPO_training_prop,CAFA2_HPO_training);
CAFA2_HPO_training_propagated=uniqueRowsCA(CAFA2_HPO_training_propagated);
CAFA2_HPO_training_propagated_unique=unique(CAFA2_HPO_training_propagated(:,1));
% (download the uniprot acc, entry name, gene symbol for the proteins in 'CAFA2_HPO_training_propagated_unique' as gene_list.txt)
[CAFA2_HPO_training_propagated_unique_UniProtID,CAFA2_HPO_training_propagated_unique_entry_name,CAFA2_HPO_training_propagated_unique_genesymbol]=textread('HPO2GO_Files/gene_list.txt', '%s %s %s', 'delimiter', '\t', 'headerlines', 1);
[Lia,Locb]=ismember(CAFA2_HPO_training_propagated(:,1),CAFA2_HPO_training_propagated_unique_UniProtID);
CAFA2_HPO_training_propagated=CAFA2_HPO_training_propagated(Lia==1,:);
CAFA2_HPO_training_propagated(:,3)=CAFA2_HPO_training_propagated_unique_genesymbol(Locb(Lia==1,1),1);
save HPO2GO_Files/CAFA2_HPO_training_propagated.mat CAFA2_HPO_training_propagated -v7

CAFA2_HPO_training_propagated_txt=CAFA2_HPO_training_propagated';
fid=fopen('HPO2GO_Files/CAFA2_HPO_training_propagated.txt', 'w');
fprintf(fid, '%s\t%s\t%s\n', CAFA2_HPO_training_propagated_txt{:});
fclose(fid);

% Generating the GO annotation arrays:

load HPO2GO_Files/CAFA2_HPO_training_propagated.mat
load HPO2GO_Files/CAFA2_Targets_human_all_mappings.mat

[GO_UniProt_acc_all,GO_gene_symbol_all,GO_term_id_all,GO_term_evid_all]=textread('HPO2GO_Files/GOA_UniProt_human_protein_annotation_2014_01.txt', '%s %s %s %s', 'delimiter', '\t');
GO_annot_manual_human_all=[GO_gene_symbol_all GO_term_id_all];
GO_annot_manual_human_all=uniqueRowsCA(GO_annot_manual_human_all);
save HPO2GO_Files/HPOprop2GOall_HPO_test_GO_annotation_all.mat GO_annot_manual_human_all -v7
GO_annot_manual_human_all_addcol=[GO_gene_symbol_all GO_term_id_all GO_term_evid_all];
GO_annot_manual_human_all_addcol=uniqueRowsCA(GO_annot_manual_human_all_addcol);
save HPO2GO_Files/HPOprop2GOall_HPO_test_GO_annotation_all_addcol.mat GO_annot_manual_human_all_addcol -v7

benchmark_hpo=pfp_loaditem('CAFA2/CAFA2-master/benchmark/lists/hpo_HUMAN_type1.txt', 'char');
[Lia,Locb]=ismember(CAFA2_Targets_human_all_mappings(:,1),benchmark_hpo);
CAFA2_Targets_human_all_mappings_CAFA2hpobench=CAFA2_Targets_human_all_mappings(Lia==1,:);
save HPO2GO_Files/HPOprop2GOall_CAFA2_Targets_human_all_mappings_CAFA2hpobench.mat CAFA2_Targets_human_all_mappings_CAFA2hpobench -v7
[Lia,Locb]=ismember(GO_annot_manual_human_all(:,1),CAFA2_Targets_human_all_mappings_CAFA2hpobench(:,4));
GO_annot_manual_human_all_CAFA2hpobench=GO_annot_manual_human_all(Lia==1,:);
save HPO2GO_Files/HPOprop2GOall_HPO_test_GO_annotation_all_CAFA2hpobench.mat GO_annot_manual_human_all_CAFA2hpobench -v7

[HPOprop_UniProt_id,HPOprop_ID,HPOprop_gene_symbol]=textread('CAFA2_HPO_training_propagated.txt', '%s %s %s', 'delimiter', '\t');
HPOprop_gene_annotation=[HPOprop_UniProt_id HPOprop_gene_symbol HPOprop_ID];
HPOprop_gene_annotation=uniqueRowsCA(HPOprop_gene_annotation);
save HPO2GO_Files/HPOprop2GOall_HPOprop_test_HPOprop_gene_annotation.mat HPOprop_gene_annotation -v7
HPOprop_gene_symbol_unique=unique(HPOprop_gene_annotation(:,2));


% Selecting the GO annotations for the HPOprop annotated protein list:


[Lia,Locb]=ismember(GO_annot_manual_human_all(:,1),HPOprop_gene_symbol_unique);
GO_annot_manual_human=GO_annot_manual_human_all(Lia==1,:);
save HPO2GO_Files/HPOprop2GOall_HPOprop_test_GO_annotation.mat GO_annot_manual_human -v7

% Generating the HPOprop and GO term pairs:

to=0;
HPOprop_GO_mapping_ind=zeros(10000000,2);
for i=1:length(HPOprop_gene_symbol_unique);
    disp(['Mapping HPOprop terms with GO terms with gene #: ', num2str(i), ' / ', num2str(length(HPOprop_gene_symbol_unique))])
    [Lia,Locb]=ismember(HPOprop_gene_annotation(:,2),HPOprop_gene_symbol_unique(i,1));
    [Lia2,Locb2]=ismember(GO_annot_manual_human(:,1),HPOprop_gene_symbol_unique(i,1));
    if sum(Lia2)>0
        siz=sum(Lia)*sum(Lia2);
        [A,B]=meshgrid(find(Lia==1),find(Lia2==1));
        c=cat(2,A',B');
        d=reshape(c,[],2);
        HPOprop_GO_mapping_ind(to+1:to+siz,:)=d;
        to=to+siz;
    end
end
HPOprop_GO_mapping_ind(HPOprop_GO_mapping_ind(:,1)==0,:)=[];
HPOprop_GO_mapping_terms=HPOprop_gene_annotation(HPOprop_GO_mapping_ind(:,1),3);
HPOprop_GO_mapping_terms(:,2)=GO_annot_manual_human(HPOprop_GO_mapping_ind(:,2),2);
save HPO2GO_Files/HPOprop2GOall_HPOprop_GO_mapping_terms.mat HPOprop_GO_mapping_terms -v7


% Calculating the co-occurrence similarity between ontology terms:

[HPOprop_GO_mapping_terms_unique,I,J]=uniqueRowsCA(HPOprop_GO_mapping_terms);
HPOprop_GO_mapping_terms_freq=sort(J);
HPOprop_GO_mapping_terms_pair_hist=hist(HPOprop_GO_mapping_terms_freq,0.5:1:max(HPOprop_GO_mapping_terms_freq)+0.5);
HPOprop_GO_mapping_terms_pair_hist=HPOprop_GO_mapping_terms_pair_hist(1,1:end-1)';
save HPO2GO_Files/HPOprop2GOall_HPOprop_GO_mapping_terms_unique.mat HPOprop_GO_mapping_terms_unique -v7
save HPO2GO_Files/HPOprop2GOall_HPOprop_GO_mapping_terms_pair_hist.mat HPOprop_GO_mapping_terms_pair_hist -v7

[HPOprop_terms_unique,~,N]=unique(HPOprop_gene_annotation(:,3));
HPOprop_terms_annot_hist=hist(N,0.5:1:max(N)+0.5);
HPOprop_terms_annot_hist=HPOprop_terms_annot_hist(1,1:end-1)';
save HPO2GO_Files/HPOprop2GOall_HPOprop_terms_unique.mat HPOprop_terms_unique -v7
save HPO2GO_Files/HPOprop2GOall_HPOprop_terms_annot_hist.mat HPOprop_terms_annot_hist -v7

[GO_terms_unique,~,N]=unique(GO_annot_manual_human(:,2));
GO_terms_annot_hist=hist(N,0.5:1:max(N)+0.5);
GO_terms_annot_hist=GO_terms_annot_hist(1,1:end-1)';
save HPO2GO_Files/HPOprop2GOall_GO_terms_unique.mat GO_terms_unique -v7
save HPO2GO_Files/HPOprop2GOall_GO_terms_annot_hist.mat GO_terms_annot_hist -v7

[~,HPOprop_GO_mapping_unique_array_ind]=ismember(HPOprop_GO_mapping_terms_unique(:,1),HPOprop_terms_unique);
[~,HPOprop_GO_mapping_unique_array_ind(:,2)]=ismember(HPOprop_GO_mapping_terms_unique(:,2),GO_terms_unique);

HPOprop_GO_mapping_sematic_sim_high_det=zeros(length(HPOprop_GO_mapping_terms_unique),5);
for i=1:length(HPOprop_GO_mapping_terms_unique);
    disp(['Calculating the semantic similarity between the HPOprop ang GO term pair #: ', num2str(i), ' / ', num2str(length(HPOprop_GO_mapping_terms_unique))])
    t1=HPOprop_GO_mapping_terms_pair_hist(i,1);
    t2=HPOprop_terms_annot_hist(HPOprop_GO_mapping_unique_array_ind(i,1),1);
    t3=GO_terms_annot_hist(HPOprop_GO_mapping_unique_array_ind(i,2),1);
    HPOprop_GO_mapping_sematic_sim_high_det(i,1:5)=[i (2*t1)/(t2+t3) t1 t2 t3]; % columns: indice of mapping, semantic similarity of the mapped terms, # of co-annotated genes, total # of annotation of the mapped HPOprop term on different genes, total # of annotation of the mapped GO term on different genes
end
save HPO2GO_Files/HPOprop2GOall_HPOprop_GO_mapping_sematic_sim_high_det.mat HPOprop_GO_mapping_sematic_sim_high_det -v7

HPOprop_GO_Raw_Mapping_merged=[HPOprop_GO_mapping_terms_unique num2cell(HPOprop_GO_mapping_sematic_sim_high_det(:,2:end))];
save HPO2GO_Files/HPOprop2GOall_HPOprop_GO_Raw_Mapping_merged.mat HPOprop_GO_Raw_Mapping_merged -v7
HPOprop_GO_Raw_Mapping_merged_txt=HPOprop_GO_Raw_Mapping_merged';
fid=fopen('HPO2GO_Files/HPOprop2GOall_HPOprop_GO_Raw_Mapping_merged.txt', 'w');
fprintf(fid, '%s\t%s\t%.4f\t%d\t%d\t%d\n', HPOprop_GO_Raw_Mapping_merged_txt{:});
fclose(fid);

% (large-scale performance test with different parameters)

load HPO2GO_Files/HPOprop2GOall_HPOprop_GO_mapping_sematic_sim_high_det.mat
load HPO2GO_Files/HPOprop2GOall_HPOprop_GO_mapping_terms_unique.mat
load HPO2GO_Files/HPO_test_GO_annotation_all_CAFA2hpobench.mat
load CAFA2/CAFA2-master/ontology/HPO.mat
load CAFA2/CAFA2-master/benchmark/groundtruth/hpoa.mat
n_mat=[2 5 10 20 30 40 50 60 70 80 90 100 110 120 130 140 145 150 155 160 165 170 175 180 185 190 195 200 210 220 230 250 250 300 400 500];
S_mat=0:0.01:1;

for j=1:length(n_mat)
    disp(['Analyzing the performance for n: ', num2str(j), ' / ', num2str(length(n_mat))])
    load HPO2GO_Files/CAFA2_Targets_human_all_mappings.mat
    n=n_mat(1,j);
    S=0.05;
    high_det_thres=HPOprop_GO_mapping_sematic_sim_high_det(HPOprop_GO_mapping_sematic_sim_high_det(:,2)>=S,:);
    high_det_thres=high_det_thres(high_det_thres(:,3)>=n,:);
    HPOprop2GOall_Mapping=HPOprop_GO_mapping_terms_unique(high_det_thres(:,1),:);
    HPOprop2GOall_Mapping(:,3:4)=num2cell(HPOprop_GO_mapping_sematic_sim_high_det(high_det_thres(:,1),2:3));  
    mapGO_unique=unique(HPOprop2GOall_Mapping(:,2));
    [~,Locb0]=ismember(HPOprop2GOall_Mapping(:,2),mapGO_unique);
    [Lia,Locb]=ismember(GO_annot_manual_human_all_CAFA2hpobench(:,2),mapGO_unique);
    to=0;
    CAFA2_HPOprop_predictions=cell(15000000,3);
    for i=1:length(mapGO_unique)
        disp(['Predicting HPO terms for the mapped GO term #: ', num2str(i), ' / ', num2str(length(mapGO_unique))])
        genesymbol_temp=GO_annot_manual_human_all_CAFA2hpobench(Locb==i,1);
        HPOprop_temp=HPOprop2GOall_Mapping(Locb0==i,1);
        score_temp=HPOprop2GOall_Mapping(Locb0==i,3);
        if sum(Lia)>0
            siz=length(genesymbol_temp)*length(HPOprop_temp);
            [A,B]=meshgrid(1:length(genesymbol_temp),1:length(HPOprop_temp));
            c=cat(2,A',B');
            d=reshape(c,[],2);
            CAFA2_HPOprop_predictions(to+1:to+siz,1)=genesymbol_temp(d(:,1),1);
            CAFA2_HPOprop_predictions(to+1:to+siz,2)=HPOprop_temp(d(:,2),1);
            CAFA2_HPOprop_predictions(to+1:to+siz,3)=score_temp(d(:,2),1);
            to=to+siz;
        end
    end
    CAFA2_HPOprop_predictions(cellfun(@isempty,CAFA2_HPOprop_predictions(:,1))==1,:)=[];
    [CAFA2_HPOprop_predictions_un,I,C]=uniqueRowsCA(CAFA2_HPOprop_predictions(:,1:2));
    score=cell2mat(CAFA2_HPOprop_predictions(:,3));
    score_unique=splitapply(@max,score,C);  % Alternative 3) Very fast but requires Matlab 2017b
    CAFA2_HPOprop_predictions_un(:,3)=num2cell(score_unique);
    CAFA2_HPOprop_predictions=CAFA2_HPOprop_predictions_un;
    CAFA2_HPOprop_predictions(ismember(CAFA2_HPOprop_predictions(:,1),'2xchrna4-3xchrnb2_human')==1,1)=cellstr('CHRNA4');
    CAFA2_HPOprop_predictions(ismember(CAFA2_HPOprop_predictions(:,1),'3xchrna4-2xchrnb2_human')==1,1)=cellstr('CHRNA4');
    mapsymbol_unique=unique(CAFA2_HPOprop_predictions(:,1));
    [Lia0,Locb0]=ismember(CAFA2_HPOprop_predictions(:,1),mapsymbol_unique);
    [Lia,Locb]=ismember(mapsymbol_unique,CAFA2_Targets_human_all_mappings(:,4));
    Locb(Locb==0,1)=length(CAFA2_Targets_human_all_mappings)+1;
    CAFA2_Targets_human_all_mappings(end+1,:)=cellstr(' ');
    mapCAFAid_unique=CAFA2_Targets_human_all_mappings(Locb,1);
    CAFA2_HPOprop_predictions_CAFAid_genesymbol_HPOpropid_sco=mapCAFAid_unique(Locb0,1);
    CAFA2_HPOprop_predictions_CAFAid_genesymbol_HPOpropid_sco(:,2:4)=CAFA2_HPOprop_predictions;
    CAFA2_pred_eval_save=CAFA2_HPOprop_predictions_CAFAid_genesymbol_HPOpropid_sco;
    CAFA2_pred_eval_save(cellfun(@isempty,CAFA2_pred_eval_save(:,1))==1,:)=[];
    CAFA2_HPOprop_target_predictions=CAFA2_pred_eval_save(:,[1 3 4]);
    ind_greaterone=cell2mat(CAFA2_HPOprop_target_predictions(:,3))>1;
    CAFA2_HPOprop_target_predictions(ind_greaterone,3)=num2cell(1.0000);
    clear CAFA2_Targets_human_all_mappings
    
    CAFA2_HPOprop_target_predictions_txt=CAFA2_HPOprop_target_predictions';
    fid=fopen('HPO2GO_Files/CAFA2_HPOprop_target_predictions_test.txt', 'w');
    fprintf(fid, '%s\t%s\t%.4f\n', CAFA2_HPOprop_target_predictions_txt{:});
    fclose(fid);
    
    pred=cafa_import('HPO2GO_Files/CAFA2_HPOprop_target_predictions_test.txt', HPO, false);
    benchmark_hpo=pfp_loaditem('CAFA2/CAFA2-master/benchmark/lists/hpo_HUMAN_type1.txt', 'char');
    fmax_hpo = pfp_seqmetric(benchmark_hpo, pred, oa, 'fmax');
    wfmax_hpo = pfp_seqmetric(benchmark_hpo, pred, oa, 'wfmax', 'w', eia);
    smin_hpo = pfp_seqmetric(benchmark_hpo, pred, oa, 'smin', 'w', eia);
    nsmin_hpo = pfp_seqmetric(benchmark_hpo, pred, oa, 'nsmin', 'w', eia);
    pr_hpo = pfp_seqmetric(benchmark_hpo, pred, oa, 'pr');
    wpr_hpo = pfp_seqmetric(benchmark_hpo, pred, oa, 'wpr', 'w', eia);
    fmax_hpo_all=(2*pr_hpo(:,1).*pr_hpo(:,2))./(pr_hpo(:,1)+pr_hpo(:,2));
    ev_hpo = cafa_eval_term_auc('HPOprop2GO', benchmark_hpo, pred, oa,'full');
    auc_hpo=mean(ev_hpo.auc(isnan(ev_hpo.auc)==0));
    benchmark=pfp_loaditem('CAFA2/CAFA2-master/benchmark/groundtruth/propagated_hpo.txt', 'char');
    benchmark_fix=benchmark(1:2:(end-1),1);
    benchmark_fix_uni=unique(benchmark_fix);
    fmax_limitedknow_benchmark = pfp_seqmetric(benchmark_fix_uni, pred, oa, 'fmax');
    perf_HPOprop2GOall(j,1)=n_mat(1,j);
    perf_HPOprop2GOall(j,2)=fmax_hpo;
    perf_HPOprop2GOall(j,3)=wfmax_hpo;
    perf_HPOprop2GOall(j,4)=smin_hpo;
    perf_HPOprop2GOall(j,5)=nsmin_hpo;
    perf_HPOprop2GOall(j,6)=auc_hpo;
    indthres=find(fmax_hpo_all==max(fmax_hpo_all));
    perf_HPOprop2GOall(j,7)=S_mat(max(indthres));
end
save HPO2GO_Files/perf_HPOprop2GOall.mat perf_HPOprop2GOall -v7

% (the thresholds n=170 and S=0.11 provided the best performance)


n=170;
S=0.11;

high_det_thres=HPOprop_GO_mapping_sematic_sim_high_det(HPOprop_GO_mapping_sematic_sim_high_det(:,2)>=S,:);
high_det_thres=high_det_thres(high_det_thres(:,3)>=n,:);
HPOprop2GOall_Finalized_Mapping=HPOprop_GO_mapping_terms_unique(high_det_thres(:,1),:);
HPOprop2GOall_Finalized_Mapping(:,3:4)=num2cell(HPOprop_GO_mapping_sematic_sim_high_det(high_det_thres(:,1),2:3));

save HPO2GO_Files/HPOprop2GOall_Finalized_Mapping.mat HPOprop2GOall_Finalized_Mapping -v7
length(unique(HPOprop2GOall_Finalized_Mapping(:,1)))
length(unique(HPOprop2GOall_Finalized_Mapping(:,2)))

mapGO_unique=unique(HPOprop2GOall_Finalized_Mapping(:,2));
[~,Locb0]=ismember(HPOprop2GOall_Finalized_Mapping(:,2),mapGO_unique);
[Lia,Locb]=ismember(GO_annot_manual_human_all_CAFA2hpobench(:,2),mapGO_unique);
to=0;
CAFA2_HPOprop_predictions=cell(15000000,3);
for i=1:length(mapGO_unique);
    disp(['Predicting HPOprop terms for the mapped GO term #: ', num2str(i), ' / ', num2str(length(mapGO_unique))])
    genesymbol_temp=GO_annot_manual_human_all_CAFA2hpobench(Locb==i,1);
    HPOprop_temp=HPOprop2GOall_Finalized_Mapping(Locb0==i,1);
    score_temp=HPOprop2GOall_Finalized_Mapping(Locb0==i,3);
    if sum(Lia)>0
        siz=length(genesymbol_temp)*length(HPOprop_temp);
        [A,B]=meshgrid(1:length(genesymbol_temp),1:length(HPOprop_temp));
        c=cat(2,A',B');
        d=reshape(c,[],2);
        CAFA2_HPOprop_predictions(to+1:to+siz,1)=genesymbol_temp(d(:,1),1);
        CAFA2_HPOprop_predictions(to+1:to+siz,2)=HPOprop_temp(d(:,2),1);
        CAFA2_HPOprop_predictions(to+1:to+siz,3)=score_temp(d(:,2),1);
        to=to+siz;
    end
end
CAFA2_HPOprop_predictions(cellfun(@isempty,CAFA2_HPOprop_predictions(:,1))==1,:)=[];
[CAFA2_HPOprop_predictions_un,I,C]=uniqueRowsCA(CAFA2_HPOprop_predictions(:,1:2));
score=cell2mat(CAFA2_HPOprop_predictions(:,3));
score_unique=splitapply(@max,score,C);  % Alternative 3) Very fast but requires Matlab 2017b
CAFA2_HPOprop_predictions_un(:,3)=num2cell(score_unique);
CAFA2_HPOprop_predictions=CAFA2_HPOprop_predictions_un;
CAFA2_HPOprop_predictions(ismember(CAFA2_HPOprop_predictions(:,1),'2xchrna4-3xchrnb2_human')==1,1)=cellstr('CHRNA4');
CAFA2_HPOprop_predictions(ismember(CAFA2_HPOprop_predictions(:,1),'3xchrna4-2xchrnb2_human')==1,1)=cellstr('CHRNA4');


mapsymbol_unique=unique(CAFA2_HPOprop_predictions(:,1));
[Lia0,Locb0]=ismember(CAFA2_HPOprop_predictions(:,1),mapsymbol_unique);
[Lia,Locb]=ismember(mapsymbol_unique,CAFA2_Targets_human_all_mappings(:,4));
Locb(Locb==0,1)=length(CAFA2_Targets_human_all_mappings)+1;
CAFA2_Targets_human_all_mappings(end+1,:)=cellstr(' ');
mapCAFAid_unique=CAFA2_Targets_human_all_mappings(Locb,1);
CAFA2_HPOprop_predictions_CAFAid_genesymbol_HPOpropid_sco=mapCAFAid_unique(Locb0,1);
CAFA2_HPOprop_predictions_CAFAid_genesymbol_HPOpropid_sco(:,2:4)=CAFA2_HPOprop_predictions;
%save HPO2GO_Files/HPOprop2GOall_CAFA2_HPOprop_predictions_CAFAid_genesymbol_HPOpropid_sco.mat CAFA2_HPOprop_predictions_CAFAid_genesymbol_HPOpropid_sco -v7
CAFA2_pred_eval_save=CAFA2_HPOprop_predictions_CAFAid_genesymbol_HPOpropid_sco;
CAFA2_pred_eval_save(cellfun(@isempty,CAFA2_pred_eval_save(:,1))==1,:)=[];
CAFA2_HPOprop_target_predictions=CAFA2_pred_eval_save(:,[1 3 4]);
ind_greaterone=find(cell2mat(CAFA2_HPOprop_target_predictions(:,3))>1);
CAFA2_HPOprop_target_predictions(ind_greaterone,3)=num2cell(1.0000);
%save HPO2GO_Files/HPOprop2GOall_CAFA2_HPOprop_target_predictions.mat CAFA2_HPOprop_target_predictions

CAFA2_HPOprop_target_predictions_txt=CAFA2_HPOprop_target_predictions';
fid=fopen('HPO2GO_Files/CAFA2_HPO_benchmark_predictions_HPOprop2GOall.txt', 'w');
fprintf(fid, '%s\t%s\t%.4f\n', CAFA2_HPOprop_target_predictions_txt{:});
fclose(fid);


% Calculating performance measures using CAFA scripts (CAFA scripts and datasets were downloaded from CAFA2 repository):

load CAFA2/CAFA2-master/ontology/HPO.mat
load CAFA2/CAFA2-master/benchmark/groundtruth/hpoa.mat
benchmark=pfp_loaditem('CAFA2/CAFA2-master/benchmark/lists/hpo_HUMAN_type1.txt', 'char');
pred=cafa_import('HPO2GO_Files/CAFA2_HPO_benchmark_predictions_HPOprop2GOall.txt', HPO, false);
pred2=cafa_import('HPO2GO_Files/CAFA2_HPO_target_predictions_n1.txt', HPO, false);
fmax=pfp_seqmetric(benchmark, pred, oa, 'fmax');
wfmax=pfp_seqmetric(benchmark, pred, oa, 'wfmax', 'w', eia);
smin=pfp_seqmetric(benchmark, pred, oa, 'smin', 'w', eia);
pr=pfp_seqmetric(benchmark, pred, oa, 'pr');
fmax_all=(2*pr(:,1).*pr(:,2))./(pr(:,1)+pr(:,2));
wpr=pfp_seqmetric(benchmark, pred, oa, 'wpr', 'w', eia);
wfmax_all=(2*wpr(:,1).*wpr(:,2))./(wpr(:,1)+wpr(:,2));
rm=pfp_seqmetric(benchmark, pred, oa, 'rm', 'w', eia);
cm=pfp_termcm(benchmark, pred, oa, 'full');
ev=cafa_eval_term_auc('HPO2GO', benchmark, pred2, oa,'full');
auc_mean=mean(ev.auc(isnan(ev.auc)==0));
save HPO2GO_Files/HPO2GO_CAFA2_perf_results.mat pred benchmark oa eia pred2 fmax wfmax smin pr rm cm ev auc_mean -v7

% (fmax=0.35, wfmax=0.29, smin=57.2 auc_mean=0.59)


% Plotting the weighted precision - recall curve for HPO2GO CAFA2 prediction test results:

thres=0:0.01:0.58;
figure;
hold on;
F=plot(wpr(1:59,2),wpr(1:59,1));
set(F,'Color',[0 0 0],'LineWidth',4)
x_lab=xlabel('Recall');
y_lab=ylabel('Precision');
set(x_lab,'FontSize',14);
set(y_lab,'FontSize',14);
title('Precision-recall Curve for CAFA2 HPO Prediction','FontSize',16);
%grid on
axis([0 1 0 1])
hold off;
S_mat=0:0.01:1;
indthres=find(wfmax_all==max(wfmax_all));
max_perf_thres=S_mat(max(indthres));





