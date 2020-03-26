% HPO2GO v1.1 - CAFA4 Run
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
% Copyright (C) 2020 CanSyL
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

[HPO_ID,~,HPO_gene_id,HPO_gene_symbol]=textread('ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes_11_2019.txt', '%s %s %s %s', 'delimiter', '\t','headerlines',1);
HPO_gene_annotation=[HPO_gene_id HPO_ID];
HPO_gene_annotation=uniqueRowsCA(HPO_gene_annotation);
save HPO_test_HPO_gene_annotation.mat HPO_gene_annotation -v7
% HPO_gene_symbol_unique=unique(HPO_gene_annotation(:,1));
% dlmcell('HPO_gene_symbol_unique.txt',HPO_gene_symbol_unique);
HPO_gene_id_unique=unique(HPO_gene_annotation(:,1));
dlmcell('HPO_gene_id_unique.txt',HPO_gene_id_unique);

% Loading all GO annoatation data for human proteins with manual evidence code:
% (download manual_GOA_20191024_propagated.tsv from: https://drive.google.com/open?id=1mz6VJjUmQN_SuAuLdLNdRFzRWMKNURcp)

[GO_term_id_all,GO_UniProt_acc_all,~,~,~,GO_tax_id_all,~]=textread('manual_GOA_20191024_propagated.tsv', '%s %s %s %s %s %s %s', 'delimiter', '\t','headerlines',1);
ind=find(ismember(GO_tax_id_all,'9606')==1);
GO_annot_manual_human_all=[GO_UniProt_acc_all(ind) GO_term_id_all(ind)];
GO_annot_manual_human_all=uniqueRowsCA(GO_annot_manual_human_all);
clear GO_term_id_all GO_UniProt_acc_all GO_tax_id_all

[Mapping_UniProtAcc,Mapping_GeneID]=textread('Mapping_UniProtAcc_GeneID.tsv', '%s %s', 'delimiter', '\t','headerlines',1);
[Lia,Locb]=ismember(GO_annot_manual_human_all(:,1),Mapping_UniProtAcc);
GO_annot_manual_human_all(:,3)=cell(length(GO_annot_manual_human_all),1);
GO_annot_manual_human_all(Lia==1,3)=Mapping_GeneID(Locb(Locb>0));
GO_annot_manual_human_all(cellfun(@ischar,GO_annot_manual_human_all(:,3))==0,3)=repmat({''},2115,1);
save HPO_test_GO_annotation_all.mat GO_annot_manual_human_all -v7

% GO_annot_manual_human_all_addcol=[GO_gene_symbol_all GO_term_id_all GO_term_name_all GO_term_aspect_all];
% GO_annot_manual_human_all_addcol=uniqueRowsCA(GO_annot_manual_human_all_addcol);
% save HPO_test_GO_annotation_all_addcol.mat GO_annot_manual_human_all_addcol -v7
% 
% % Saving the processed GO annotation file will additional columns as text:
% 
% GO_annot_manual_human_all_addcol_txt=GO_annot_manual_human_all_addcol';
% fid=fopen('GO_annot_human_proteins_UniProtGOA_01_2017.txt', 'w');
% fprintf(fid, '%s\t%s\t%s\t%s\n', GO_annot_manual_human_all_addcol_txt{:});
% fclose(fid);

% GO_annot_manual_human_all_uniprotacc=[GO_UniProt_acc_all GO_gene_symbol_all GO_term_id_all];
% GO_annot_manual_human_all_uniprotacc=uniqueRowsCA(GO_annot_manual_human_all_uniprotacc);
% save HPO_test_GO_annotation_all_uniprotacc.mat GO_annot_manual_human_all_uniprotacc -v7

% Selecting the GO annotations for the HPO annotated protein list:

[Lia,Locb]=ismember(GO_annot_manual_human_all(:,3),HPO_gene_id_unique);
GO_annot_manual_human=GO_annot_manual_human_all(Lia==1,:);
save HPO_test_GO_annotation.mat GO_annot_manual_human -v7

% Generating the HPO and GO term pairs:

to=0;
HPO_GO_mapping_ind=zeros(20000000,2);
for i=1:length(HPO_gene_id_unique)
    disp(['Mapping HPO terms with GO terms with gene #: ', num2str(i), ' / ', num2str(length(HPO_gene_id_unique))])
    [Lia,Locb]=ismember(HPO_gene_annotation(:,1),HPO_gene_id_unique(i,1));
    [Lia2,Locb2]=ismember(GO_annot_manual_human(:,3),HPO_gene_id_unique(i,1));
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
HPO_GO_mapping_terms=HPO_gene_annotation(HPO_GO_mapping_ind(:,1),2);
HPO_GO_mapping_terms(:,2)=GO_annot_manual_human(HPO_GO_mapping_ind(:,2),2);
save HPO_GO_mapping_terms.mat HPO_GO_mapping_terms -v7


% Calculating the co-occurrence similarity between ontology terms:

[HPO_GO_mapping_terms_unique,I,J]=uniqueRowsCA(HPO_GO_mapping_terms);
HPO_GO_mapping_terms_freq=sort(J);
HPO_GO_mapping_terms_pair_hist=hist(HPO_GO_mapping_terms_freq,0.5:1:max(HPO_GO_mapping_terms_freq)+0.5);
HPO_GO_mapping_terms_pair_hist=HPO_GO_mapping_terms_pair_hist(1,1:end-1)';
save HPO_GO_mapping_terms_unique.mat HPO_GO_mapping_terms_unique -v7
save HPO_GO_mapping_terms_pair_hist.mat HPO_GO_mapping_terms_pair_hist -v7

[HPO_terms_unique,~,N]=unique(HPO_gene_annotation(:,2));
HPO_terms_annot_hist=hist(N,0.5:1:max(N)+0.5);
HPO_terms_annot_hist=HPO_terms_annot_hist(1,1:end-1)';
save HPO_terms_unique.mat HPO_terms_unique -v7
save HPO_terms_annot_hist.mat HPO_terms_annot_hist -v7

[GO_terms_unique,~,N]=unique(GO_annot_manual_human(:,2));
GO_terms_annot_hist=hist(N,0.5:1:max(N)+0.5);
GO_terms_annot_hist=GO_terms_annot_hist(1,1:end-1)';
save GO_terms_unique.mat GO_terms_unique -v7
save GO_terms_annot_hist.mat GO_terms_annot_hist -v7

[~,HPO_GO_mapping_unique_array_ind]=ismember(HPO_GO_mapping_terms_unique(:,1),HPO_terms_unique);
[~,HPO_GO_mapping_unique_array_ind(:,2)]=ismember(HPO_GO_mapping_terms_unique(:,2),GO_terms_unique);

HPO_GO_mapping_sematic_sim_high_det=zeros(length(HPO_GO_mapping_terms_unique),5);
for i=1:length(HPO_GO_mapping_terms_unique)
    disp(['Calculating the semantic similarity between the HPO ang GO term pair #: ', num2str(i), ' / ', num2str(length(HPO_GO_mapping_terms_unique))])
    t1=HPO_GO_mapping_terms_pair_hist(i,1);
    t2=HPO_terms_annot_hist(HPO_GO_mapping_unique_array_ind(i,1),1);
    t3=GO_terms_annot_hist(HPO_GO_mapping_unique_array_ind(i,2),1);
    HPO_GO_mapping_sematic_sim_high_det(i,1:5)=[i (2*t1)/(t2+t3) t1 t2 t3]; % columns: indice of mapping, semantic similarity of the mapped terms, # of co-annotated genes, total # of annotation of the mapped HPO term on different genes, total # of annotation of the mapped GO term on different genes
end
save HPO_GO_mapping_sematic_sim_high_det.mat HPO_GO_mapping_sematic_sim_high_det -v7

HPO_GO_Raw_Mapping_merged=[HPO_GO_mapping_terms_unique num2cell(HPO_GO_mapping_sematic_sim_high_det(:,2:end))];
save HPO_GO_Raw_Mapping_merged.mat HPO_GO_Raw_Mapping_merged -v7

% Saving the original raw mapping file as text:

HPO_GO_Raw_Mapping_merged_txt=HPO_GO_Raw_Mapping_merged';
fid=fopen('HPO_GO_Raw_Original_Mapping.txt', 'w');
fprintf(fid, '%s\t%s\t%.4f\t%d\t%d\t%d\n', HPO_GO_Raw_Mapping_merged_txt{:});
fclose(fid);

% Saving the finalized HPO2GO mapping file:

n=5;
S=0.1;

high_det_thres=HPO_GO_mapping_sematic_sim_high_det(HPO_GO_mapping_sematic_sim_high_det(:,2)>=S,:);
high_det_thres=high_det_thres(high_det_thres(:,3)>=n,:);
HPO2GO_Finalized_Mapping=HPO_GO_mapping_terms_unique(high_det_thres(:,1),:);
HPO2GO_Finalized_Mapping(:,3:4)=num2cell(HPO_GO_mapping_sematic_sim_high_det(high_det_thres(:,1),2:3));

% (deleting the mappings to HPO terms that are not belong to phenotypic abnormality sub-ontology)

[HPO_ID_not_abnor,HPO_Name_not_abnor]=textread('HPO_terms_not_abnormality.txt', '%s %s', 'delimiter', '\t');
[Lia,Locb]=ismember(HPO2GO_Finalized_Mapping(:,1),HPO_ID_not_abnor);
HPO2GO_Finalized_Mapping(Lia==1,:)=[];

save HPO2GO_Finalized_Mapping.mat HPO2GO_Finalized_Mapping -v7
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






% Predicting HPO terms for CAFA4 human target proteins:


% Generating CAFA id vs. gene symbol vs. entry name vs. UniProt acc. mapping file:

[Humprot_UniProt_acc,Humprot_gene_symbol,Humprot_entry_name,Humprot_gene_id]=textread('Protein_acc_name_mapping.tab', '%s %s %s %s', 'delimiter', '\t','headerlines',1);
Humprot_gene_symbol=strtok(Humprot_gene_symbol);
save Humprot_mapping_variables.mat Humprot_UniProt_acc Humprot_gene_symbol Humprot_entry_name Humprot_gene_id -v7

% (download CAFA4-export files from: https://www.biofunctionprediction.org/cafa-targets/CAFA4-export.tgz)

CAFA4_Targets_human_mappings=cell(20430,2);
[CAFA4_Targets_human_mappings(:,1),CAFA4_Targets_human_mappings(:,2)]=textread('CAFA4-export/MappingFiles/mapping.9606.map', '%s %s', 'delimiter', '\t');
save CAFA4_HPO_target_variables.mat CAFA4_Targets_human_mappings -v7

[~,Locb]=ismember(CAFA4_Targets_human_mappings(:,2),Humprot_entry_name);
Locb(Locb==0,1)=length(Humprot_entry_name)+1;
Humprot_entry_name(end+1,1)=cellstr(' ');
Humprot_UniProt_acc(end+1,1)=cellstr(' ');
Humprot_gene_symbol(end+1,1)=cellstr(' ');
Humprot_gene_id(end+1,1)=cellstr(' ');
CAFA4_Targets_human_all_mappings=CAFA4_Targets_human_mappings;
CAFA4_Targets_human_all_mappings(:,3)=Humprot_UniProt_acc(Locb,1);
CAFA4_Targets_human_all_mappings(:,4)=Humprot_gene_symbol(Locb,1);
CAFA4_Targets_human_all_mappings(:,5)=Humprot_gene_id(Locb,1);
save CAFA4_Targets_human_all_mappings.mat CAFA4_Targets_human_all_mappings -v7


load HPO_GO_mapping_sematic_sim_high_det.mat
load HPO_GO_mapping_terms_unique.mat
load CAFA4_Targets_human_all_mappings.mat
load CAFA4_HPO_target_variables.mat
load Humprot_mapping_variables.mat
load HPO_GO_mapping_terms.mat
load HPO_test_GO_annotation.mat
load HPO_test_HPO_gene_annotation.mat
load HPO_test_GO_annotation_all.mat
load HPO2GO_Finalized_Mapping.mat


% Prediction set 1:

mapGO_unique=unique(HPO2GO_Finalized_Mapping(:,2));
[~,Locb0]=ismember(HPO2GO_Finalized_Mapping(:,2),mapGO_unique);
[Lia,Locb]=ismember(GO_annot_manual_human_all(:,2),mapGO_unique);

to=0;
CAFA4_HPO_predictions=cell(20000000,3);
for i=1:1600
    disp(['Predicting HPO terms for the mapped GO term #: ', num2str(i), ' / ', num2str(length(mapGO_unique))])
    genesymbol_temp=GO_annot_manual_human_all(Locb==i,1);
    HPO_temp=HPO2GO_Finalized_Mapping(Locb0==i,1);
    score_temp=HPO2GO_Finalized_Mapping(Locb0==i,3);
    if sum(Lia)>0
        siz=length(genesymbol_temp)*length(HPO_temp);
        [A,B]=meshgrid(1:length(genesymbol_temp),1:length(HPO_temp));
        c=cat(2,A',B');
        d=reshape(c,[],2);
        CAFA4_HPO_predictions(to+1:to+siz,1)=genesymbol_temp(d(:,1),1);
        CAFA4_HPO_predictions(to+1:to+siz,2)=HPO_temp(d(:,2),1);
        CAFA4_HPO_predictions(to+1:to+siz,3)=score_temp(d(:,2),1);
        to=to+siz;
    end
end
CAFA4_HPO_predictions(cellfun(@isempty,CAFA4_HPO_predictions(:,1))==1,:)=[];


% Prediction set 2:

mapGO_unique=unique(HPO2GO_Finalized_Mapping(:,2));
[~,Locb0]=ismember(HPO2GO_Finalized_Mapping(:,2),mapGO_unique);
[Lia,Locb]=ismember(GO_annot_manual_human_all(:,2),mapGO_unique);

to=0;
CAFA4_HPO_predictions_2=cell(20000000,3);
for i=1601:length(mapGO_unique)
    disp(['Predicting HPO terms for the mapped GO term #: ', num2str(i), ' / ', num2str(length(mapGO_unique))])
    genesymbol_temp=GO_annot_manual_human_all(Locb==i,1);
    HPO_temp=HPO2GO_Finalized_Mapping(Locb0==i,1);
    score_temp=HPO2GO_Finalized_Mapping(Locb0==i,3);
    if sum(Lia)>0
        siz=length(genesymbol_temp)*length(HPO_temp);
        [A,B]=meshgrid(1:length(genesymbol_temp),1:length(HPO_temp));
        c=cat(2,A',B');
        d=reshape(c,[],2);
        CAFA4_HPO_predictions_2(to+1:to+siz,1)=genesymbol_temp(d(:,1),1);
        CAFA4_HPO_predictions_2(to+1:to+siz,2)=HPO_temp(d(:,2),1);
        CAFA4_HPO_predictions_2(to+1:to+siz,3)=score_temp(d(:,2),1);
        to=to+siz;
    end
end
CAFA4_HPO_predictions_2(cellfun(@isempty,CAFA4_HPO_predictions_2(:,1))==1,:)=[];



CAFA4_HPO_predictions=[CAFA4_HPO_predictions;CAFA4_HPO_predictions_2];

[CAFA4_HPO_predictions_un,I,C]=uniqueRowsCA(CAFA4_HPO_predictions(:,1:2));

score=cell2mat(CAFA4_HPO_predictions(:,3));
score_unique=splitapply(@max,score,C);  % Alternative 3) Very fast but requires Matlab 2017b
CAFA4_HPO_predictions_un(:,3)=num2cell(score_unique);
CAFA4_HPO_predictions=CAFA4_HPO_predictions_un;
CAFA4_HPO_predictions(ismember(CAFA4_HPO_predictions(:,1),'2xchrna4-3xchrnb2_human')==1,1)=cellstr('CHRNA4');
CAFA4_HPO_predictions(ismember(CAFA4_HPO_predictions(:,1),'3xchrna4-2xchrnb2_human')==1,1)=cellstr('CHRNA4');

mapsymbol_unique=unique(CAFA4_HPO_predictions(:,1));
[Lia0,Locb0]=ismember(CAFA4_HPO_predictions(:,1),mapsymbol_unique);
[Lia,Locb]=ismember(mapsymbol_unique,CAFA4_Targets_human_all_mappings(:,3));
Locb(Locb==0,1)=length(CAFA4_Targets_human_all_mappings)+1;
CAFA4_Targets_human_all_mappings(end+1,:)=cellstr(' ');
mapCAFAid_unique=CAFA4_Targets_human_all_mappings(Locb,1);
CAFA4_HPO_predictions_CAFAid_genesymbol_HPOid_sco=mapCAFAid_unique(Locb0,1);
CAFA4_HPO_predictions_CAFAid_genesymbol_HPOid_sco(:,2:4)=CAFA4_HPO_predictions;
save CAFA4_HPO_predictions_CAFAid_genesymbol_HPOid_sco.mat CAFA4_HPO_predictions_CAFAid_genesymbol_HPOid_sco -v7
CAFA4_pred_eval_save=CAFA4_HPO_predictions_CAFAid_genesymbol_HPOid_sco;
CAFA4_pred_eval_save(cellfun(@isempty,CAFA4_pred_eval_save(:,1))==1,:)=[];
CAFA4_pred_eval_save=uniqueRowsCA(CAFA4_pred_eval_save);
save CAFA4_HPO_predictions_semantic.mat CAFA4_pred_eval_save -v7
CAFA4_HPO_target_predictions=CAFA4_pred_eval_save(:,[1 3 4]);
save CAFA4_HPO_target_predictions.mat CAFA4_HPO_target_predictions -v7
length(CAFA4_HPO_target_predictions)
length(unique(CAFA4_HPO_target_predictions(:,1)))
length(unique(CAFA4_HPO_target_predictions(:,2)))

% Saving CAFA4 HPO target predictions in a text file:

CAFA4_pred_eval_save_txt=CAFA4_pred_eval_save(:,[1 2 3 4])';
fid=fopen('CAFA4_HPO_target_predictions.txt', 'w');
fprintf(fid, '%s\t%s\t%s\t%.2f\n', CAFA4_pred_eval_save_txt{:});
fclose(fid);




