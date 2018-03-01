## HPO2GO v1.0

Mapping between Human Phenotype Ontology (HPO) and Gene Ontology (GO) terms for the prediction of gene/protein - function - phenotype - disease associations.

## Title, Author Info and Abstract

HPO2GO: Prediction of Human Phenotype Ontology Term Associations with Cross Ontology Annotation Co-occurrences

Tunca Doğan1,2,3,*

1	Cancer Systems Biology Laboratory (CanSyL), Graduate School of Informatics, METU, Ankara, 06800, Turkey

2	Department of Health Informatics, Graduate School of Informatics, METU, Ankara, 06800, Turkey

3	European Molecular Biology Laboratory, European Bioinformatics Institute (EMBL-EBI), Hinxton, Cambridge, CB10 1SD, UK

	* corresponding author email addresses: tdogan@metu.edu.tr, tdogan@ebi.ac.uk

ABSTRACT

Analysing the relations between diseases with genetic origin and biomolecules is a highly active area of research, where the aim is to identify the genes and gene products that cause a particular disease due to functional changes originated from mutations. Biological ontologies are frequently employed for this purpose. The Gene Ontology (GO) systematically defines the functions of genes/proteins using uniquely defined terms. Another system, the Human Phenotype Ontology (HPO) defines the phenotypic abnormalities with a controlled vocabulary that maps to disease entries in data resources such as OMIM and Orphanet. The employment of ontologies provided researchers with opportunities for knowledge discovery through computational data analysis, which significantly advanced the field.

In this study, a novel approach is proposed for the identification of relations between biomedical entities by semantically mapping phenotypic abnormality defining HPO terms with biomolecular function defining GO terms, where each association indicates the occurrence of the abnormality due to the loss of the biomolecular function expressed by the corresponding GO term. The proposed HPO2GO mappings were extracted by calculating the frequency of the co-annotations of the terms on the same genes/proteins, using already existing curated HPO and GO annotation sets. This was followed by the filtering of the unreliable mappings that could be observed due to chance by statistical resampling of the co-occurrence similarity distributions using the Kolmogorov-Smirnov (KS) statistic. Furthermore, biological relevance of the finalized mappings were discussed over selected cases, using the literature based information.

The resulting HPO2GO mappings can be employed in different settings to predict and to analyse novel gene/protein - ontological term - disease relations. As an application of the proposed approach, HPO term - protein associations (i.e., HPO2protein) are predicted using the HPO2GO mappings. In order to test the predictive performance of the method on a quantitative basis, and to compare it with the state-of-the-art, Critical Assessment of Functional Annotation 2 (CAFA2) challenge HPO prediction track target protein set was employed. The results of the benchmark experiments showed that the proposed approach is effective, as HPO2GO was among the top three performers out of 38 CAFA2 participating groups (Fmax=0.35). It is important to note that, HPO2GO was not proposed to replace, but to complement the conventional approaches used in the field of biomedical relation discovery. The cross ontology mapping approach developed in this work can easily be extended to other ontologies as well, to identify unexplored relation patterns at the systemic level. The proposed approach will be especially effective when combined with powerful techniques such as text/literature mining.

## License
HPO2GO: Prediction of Human Phenotype Ontology Term Associations with Cross Ontology Annotation Co-occurrences

Copyright (C) 2018 CanSyL

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
