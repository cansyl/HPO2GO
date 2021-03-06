## HPO2GO v1.1

Mapping between Human Phenotype Ontology (HPO) and Gene Ontology (GO) terms for the prediction of gene/protein - function - phenotype - disease associations.

## Citation

If you find HPO2GO useful, please consider citing this publication:

Doğan, T. (2018). HPO2GO: Prediction of Human Phenotype Ontology Term Associations Using Cross Ontology Annotation Co-occurrences. PeerJ 6:e5298.

## Title, Author Info and Abstract

HPO2GO: Prediction of Human Phenotype Ontology Term Associations Using Cross Ontology Annotation Co-occurrences

Tunca Doğan1,2,3,*

1 Cancer Systems Biology Laboratory (CanSyL), Graduate School of Informatics, METU, Ankara, 06800, Turkey  
2 Department of Health Informatics, Graduate School of Informatics, METU, Ankara, 06800, Turkey  
3 European Molecular Biology Laboratory, European Bioinformatics Institute (EMBL-EBI), Hinxton, Cambridge, CB10 1SD, UK

\* corresponding author email addresses: tdogan@metu.edu.tr, tdogan@ebi.ac.uk

**ABSTRACT**

Analysing the relationships between biomolecules and the genetic diseases is a highly active area of research, where the aim is to identify the genes and their products that cause a particular disease due to functional changes originated from mutations. Biological ontologies are frequently employed in these studies, which provides researchers with extensive opportunities for knowledge discovery through computational data analysis.

In this study, a novel approach is proposed for the identification of relationships between biomedical entities by automatically mapping phenotypic abnormality defining HPO terms with biomolecular function defining GO terms, where each association indicates the occurrence of the abnormality due to the loss of the biomolecular function expressed by the corresponding GO term. The proposed HPO2GO mappings were extracted by calculating the frequency of the co-annotations of the terms on the same genes/proteins, using already existing curated HPO and GO annotation sets. This was followed by the filtering of the unreliable mappings that could be observed due to chance, by statistical resampling of the co-occurrence similarity distributions. Furthermore, the biological relevance of the finalized mappings were discussed over selected cases, using the literature.

The resulting HPO2GO mappings can be employed in different settings to predict and to analyse novel gene/protein—ontology term—disease relations. As an application of the proposed approach, HPO term—protein associations (i.e., HPO2protein) were predicted. In order to test the predictive performance of the method on a quantitative basis, and to compare it with the state-of-the-art, CAFA2 challenge HPO prediction target protein set was employed. The results of the benchmark indicated the potential of the proposed approach, as HPO2GO performance was among the best (Fmax = 0.35). The automated cross ontology mapping approach developed in this work may be extended to other ontologies as well, to identify unexplored relation patterns at the systemic level. The datasets, results and the source code of HPO2GO are available for download at: https://github.com/cansyl/HPO2GO.

## Programming Environment

All analyses and computation in this study (available through: 'HPO2GO_Source_Code.m') were carried out in MATLAB Release 2017b by The MathWorks.

## License
HPO2GO: Prediction of Human Phenotype Ontology Term Associations with Cross Ontology Annotation Co-occurrences

Copyright (C) 2018 CanSyL

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
