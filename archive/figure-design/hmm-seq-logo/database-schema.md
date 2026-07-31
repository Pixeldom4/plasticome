# Database Schema

## Overview

PETadex SQL database is an ordered and accessible way to access the petadex project data via pgadmin4 or your favourite PostgreSQL program.

This documentation is automatically generated from the live database schema.

---

Auto-generated on Sun Jul 26 05:38:09 UTC 2026


## Tables

### 1.ENZYMATIC_ACTIVITIES Table

**Description:** 

**Estimated rows:** _not analyzed_

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|


### 2.ASSAY_CONDITIONS Table

**Description:** 

**Estimated rows:** _not analyzed_

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|


### 3.EXPERIMENT_METADATA Table

**Description:** 

**Estimated rows:** _not analyzed_

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|


### 30PID_SUPERFAMILY_CLUSTERS Table

**Description:** DIAMOND clustering of PETadex 60pid family centroids by 30% shared amino acid identity. Clustering is forward compatible: new sequences are put into existing clusters, and the ones that do not align are then clustered. The order of database clustering was hierarchical: PAZy then NR, then Logan. 

**Estimated rows:** ~22,235

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | 30pid_superfamily_id | integer | No |  | primary.key  The integer accession of the 30 percent identity clusters within the PETadex. |
 | centroid_orf_id | bigint | No |  | foreign.key refers to "orf_origins"  The PETadex ORF integer accession of 30 percent identity cluster centroid. |
 | date_clustered | date | No |  | The date (ISO-8601 format) the clustering algorithm was performed. |


### 4.BIBLIOGRAPHY Table

**Description:** 

**Estimated rows:** _not analyzed_

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|


### 60PID_FAMILY_CLUSTERS Table

**Description:** DIAMOND clustering of PETadex 90pid centroids by 60% shared amino acid identity. Clustering is forward compatible: new sequences are put into existing clusters, and the ones that do not align are then clustered. The order of database clustering was hierarchical: PAZy then NR, then Logan. 

**Estimated rows:** ~1,814,296

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | 60pid_family_id | bigint | No |  | primary.key  The integer accession of the 60 percent identity clusters within the PETadex. |
 | centroid_orf_id | bigint | No |  | foreign.key refers to "orf_origins"  The PETadex ORF integer accession of 60 percent identity cluster centroid. |
 | date_clustered | date | No |  | The date (ISO-8601 format) the clustering algorithm was performed. |


### 90PID_ENZYME_CLUSTERS Table

**Description:** DIAMOND clustering of PETadex orfs by 90% shared amino acid identity. Clustering is forward compatible: new sequences are put into existing clusters, and the ones that do not align are then clustered. The order of database clustering was hierarchical: PAZy then NR, then Logan. 

**Estimated rows:** ~18,173,654

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | 90pid_enzyme_id | bigint | No |  | primary.key  The integer accession of the 90 percent identity clusters within the PETadex. |
 | centroid_orf_id | integer | No |  | foreign.key refers to "orf_origins"  The PETadex ORF integer accession of 90 percent identity cluster centroid. |
 | date_clustered | date | No |  | The date (ISO-8601 format) the clustering algorithm was performed. |


### BLAST_NR_METADATA Table

**Description:** No description available

**Estimated rows:** ~2,736,843

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | genbank_accession_id | character varying | No |  |  |
 | organism | text | Yes |  |  |
 | protein_id | integer | Yes |  |  |
 | definition | text | Yes |  |  |
 | taxonomy | text | Yes |  |  |
 | journal | text | Yes |  |  |
 | collection_date | text | Yes |  |  |
 | country | text | Yes |  |  |


### COMPONENT_CATH_DICTIONARY Table

**Description:** No description available

**Estimated rows:** _not analyzed_

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | component | integer | Yes |  |  |
 | gene | text | Yes |  |  |
 | domain | text | Yes |  |  |
 | cath | text | Yes |  |  |


### DOCKING_ANALYSIS Table

**Description:** No description available

**Estimated rows:** ~4,385,957

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | docking_job_id | character varying(12) | No |  | primary.key + foreign.key  Refers to 'docking_data'. |
 | conformation | real | No |  | primary.key  The number of the conformation within the outputted docking pdbqt. |
 | spatial_properties | ARRAY | Yes |  |  |
 | binding_energy | real | Yes |  |  |
 | properties_script | text | Yes |  |  |
 | date_entered | timestamp with time zone | Yes | now() |  |


### DOCKING_DATA Table

**Description:** The table will contain the predicted conformation and gibbs free energy of binding for each run within a job.

**Estimated rows:** ~5,069

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | docking_job_id | character varying | No |  | Accession for the docking conformation | letters refer to the software used, first set of numbers is the batch of conformations, second set of numbers refers to the conformation within the job |
 | docking_batch_id | character varying | Yes |  | Accession of the docking job | refers to docking_metadata |
 | date_entered | timestamp with time zone | Yes | now() | now() |
 | pdb_id | character varying | Yes |  | The structure used as the receptor | refers to pdb_accessions |
 | ligand_id | character varying | Yes |  | The ligand used for docking | refers to ligands |


### DOCKING_METADATA Table

**Description:** No description available

**Estimated rows:** _not analyzed_

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | docking_batch_id | character varying | No |  | Unique docking job identifier with 2 letters (refers to the assay name) and 5 numbers (refers to order of plate creation) ex. LS00001 |
 | prepare_receptor | character varying | Yes |  | The structure used as the receptor | refers to pdb_accessions |
 | prepare_ligand | character varying | Yes |  | The ligand used for docking | refers to ligands |
 | gridbox_center | ARRAY | Yes |  | [x, y, z] of the gridbox used as the area for docking |
 | gridbox_size | ARRAY | Yes |  | [x, y, z] of the size of the gridbox |
 | rigid_receptor | boolean | Yes |  | Boolean if the receptor is rigid or includes flexible residues |
 | software | character varying | Yes |  | The software used to calculate the docking |
 | algorithm | character varying | Yes |  | Algorithm used to perform the docking |
 | exhaustiveness | integer | Yes |  |  |
 | procedure_script | character varying | Yes |  |  |
 | date_entered | timestamp with time zone | No | now() |  |
 | date_ran | timestamp with time zone | Yes |  |  |
 | experiment_id | character varying | Yes |  |  |
 | exp_description | character varying | Yes |  |  |


### ENZYME_FASTAA Table

**Description:** No description available

**Estimated rows:** ~1,048,585

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | enzyme_id | integer | No |  | primary.key  Enzymes are the PETadex sequences amino acid sequences clustered to 90%. The sequences that compose an 'enzyme' cluster are called 'variants'. The 'variant' that is the 90% clustering centroid is the 'enzyme' representative and is given a unique numeric identifier. |
 | contig_id | integer | Yes |  | foreign.key  The assembled contig from Logan that contains the open reading frame. |
 | orf_start | integer | Yes |  | The start nucleotide of the open reading frame. If the ORF is partial, where the start codon is missing, then the orf_start is 0. |
 | orf_end | integer | Yes |  | The end nucleotide of the open reading frame. If the ORF is partial, where the stop codon is missing, then the orf_end is -1. |
 | translated_sequence | text | Yes |  | The nucleotide sequence of the open reading frame translated into amino acids. |
 | genbank_accession_id | character varying(16) | Yes |  | foreign.key  The genbank protein accession from PETadex NR. |
 | orf_type | smallint | Yes |  | Describes the completeness of the ORF. 0 = complete + positive sense 1 = 5' partial + positive sense 2 = 3' partial + positive sense 3 = 5' and 3' partial + positive sense 4 = complete + negative sense 5 = 5' partial + negative sense 6 = 3' partial + negative sense 7 = 5' and 3' partial + negative sense |
 | library_id | character varying | Yes |  | foreign.key  The SRA library that contains the assembled contig from Logan. |


### ENZYME_TAXONOMY Table

**Description:** No description available

**Estimated rows:** ~1,046,402

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | enzyme_id | integer | No |  |  |
 | family | integer | Yes |  |  |
 | family_pid | numeric | Yes |  |  |
 | component | integer | Yes |  |  |
 | cath_domain | text | Yes |  |  |
 | date_entered | timestamp with time zone | Yes | now() |  |
 | gene | text | Yes |  |  |
 | domain_name | text | Yes |  |  |


### FAMILY_UMAP_COORDINATES Table

**Description:** No description available

**Estimated rows:** ~64,730

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | family_id | integer | No |  |  |
 | umap_x | double precision | No |  |  |
 | umap_y | double precision | No |  |  |
 | family_size | integer | No |  |  |


### FASTAA Table

**Description:** Primary table for storing biological amino acid sequence data of wildtype and mutant enzymes.

**Estimated rows:** ~602

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | accession | character varying | No |  | Unique sequence identifier |
 | aa_sequence | character varying | Yes |  | Full ORF amino acid sequence |
 | source | character varying | Yes |  | Database containing the accession or method of synthesis |
 | date_entered | timestamp with time zone | Yes | now() | Auto |
 | genotype | character varying | Yes |  | HGVS nomenclature describing changes applied to the sequence |
 | genotype_description | character varying | Yes |  | Brief description of the genotype / modifications applied to the sequence |
 | synthetic | boolean | Yes |  | Is the sequence found within nature or synthetic |
 | parent_accessions | ARRAY | Yes |  | The accessions used to synthesize the sequence |
 | parent_genes | ARRAY | Yes |  | The genes used to synthesize the accession  |
 | synonyms | ARRAY | Yes |  | Common names of the enzyme |
 | in_gene_metadata | boolean | Yes |  |  |


### GENE_METADATA Table

**Description:** This table contains all of the metadata that is tied to a synthesized gene.

**Estimated rows:** ~259

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | gene | character varying | No |  | Unique gene identifier that is the common name of the gene or accession consisting of 2 letters (refers to the accession source) and 5 numbers (refers to order of synthesis) ex. LG00001 |
 | nickname | character varying | Yes |  | A common name for the gene that is randomly generated in the format of: adjective + noun. Automatically filled from nicknames_list. ex. DaffyDuck |
 | accession | character varying | Yes |  | The accession for the amino acid sequence of the gene |
 | orf_nt_sequence | character varying | Yes |  | The ORF nucleotide sequence |
 | left_homology_arm | character varying | Yes |  | The left homology arm used to insert the gene into the plasmid. |
 | right_homology_arm | character varying | Yes |  | The right homology arm used to insert the gene into the plasmid. |
 | batch | character varying | Yes |  | An identifier to the batch that the gene was ordered in. |
 | date_entered | timestamp with time zone | Yes | now() | Auto |
 | genetic_code | character varying | Yes |  | The organism that the gene was optimized for |


### INTERPRO_SCAN_RESULTS Table

**Description:** Stores relevant domains for the corresponding protein.

**Estimated rows:** ~3,771

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | id | integer | No | nextval('interpro_scan_results_id_seq'::regclass) |  |
 | protein_id | character varying(255) | No |  |  |
 | sequence_md5 | character varying(32) | Yes |  |  |
 | sequence_length | integer | Yes |  |  |
 | library_name | character varying(100) | No |  |  |
 | signature_accession | character varying(100) | Yes |  |  |
 | signature_name | character varying(255) | Yes |  |  |
 | signature_description | text | Yes |  |  |
 | match_start | integer | Yes |  |  |
 | match_end | integer | Yes |  |  |
 | evalue | double precision | Yes |  |  |
 | score | double precision | Yes |  |  |
 | is_representative | boolean | Yes | false |  |
 | interpro_accession | character varying(50) | Yes |  |  |
 | interpro_name | character varying(255) | Yes |  |  |
 | interpro_description | text | Yes |  |  |
 | interpro_type | character varying(100) | Yes |  |  |
 | signature_type | character varying(100) | Yes |  |  |
 | model_accession | character varying(100) | Yes |  |  |
 | hmm_start | integer | Yes |  |  |
 | hmm_end | integer | Yes |  |  |
 | hmm_length | integer | Yes |  |  |
 | hmm_bounds | character varying(50) | Yes |  |  |
 | envelope_start | integer | Yes |  |  |
 | envelope_end | integer | Yes |  |  |
 | fragment_start | integer | Yes |  |  |
 | fragment_end | integer | Yes |  |  |
 | dc_status | character varying(50) | Yes | 'CONTINUOUS'::character varying |  |
 | binding_site_description | text | Yes |  |  |
 | binding_site_type | character varying(100) | Yes |  |  |
 | binding_site_num_locations | integer | Yes |  |  |
 | binding_position_start | integer | Yes |  |  |
 | binding_position_end | integer | Yes |  |  |
 | binding_residue | character(1) | Yes |  |  |
 | binding_confidence_score | double precision | Yes |  |  |
 | go_accessions | jsonb | Yes |  |  |
 | go_terms | jsonb | Yes |  |  |
 | pathway_database | character varying(100) | Yes |  |  |
 | pathway_accession | character varying(100) | Yes |  |  |
 | pathway_name | character varying(255) | Yes |  |  |
 | xref_name | character varying(100) | Yes |  |  |
 | xref_id_value | character varying(255) | Yes |  |  |
 | library_version | character varying(50) | Yes |  |  |
 | created_at | timestamp without time zone | Yes | CURRENT_TIMESTAMP |  |
 | updated_at | timestamp without time zone | Yes | CURRENT_TIMESTAMP |  |


### LOGAN_CATALYTIC_ORFS Table

**Description:** This table contains all the PAZy/BlastNR homologs retrieved from Logan using DIAMOND, which contain a complete catalytic domain, described as matching to all of the catalytic sites of a PAZy HMM.

**Estimated rows:** ~302,432,416

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | orf_id | bigint | No |  | foreign.key refers to "orf_origins"  The PETadex ORF integer accession of the amino acid sequence. |
 | library_id | text | No |  | foreign.key  The SRA library that contains the assembled contig from Logan. |
 | contig | bigint | No |  | The assembled contig from Logan that contains the ORF. |
 | orf_start | integer | No |  | The start nucleotide of the open reading frame. |
 | orf_end | integer | No |  | The end nucleotide of the open reading frame. |
 | orf_type | smallint | No |  | Describes the completeness and sense of the ORF. 0 = complete & sense strand 1 = 3' partial & sense strand 2 = complete & antisense strand 3 = 3' partial & antisense strand |


### LOGAN_PAZY_ACCESSIONS Table

**Description:** No description available

**Estimated rows:** _not analyzed_

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | protein_id | text | No |  |  |
 | component | real | Yes |  |  |
 | family | real | Yes |  |  |
 | enzyme | real | Yes |  |  |
 | mutant | real | Yes |  |  |
 | nucleotide_accession | text | Yes |  |  |


### NCBI_METADATA Table

**Description:** No description available

**Estimated rows:** ~64,730

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | protein_id | integer | No |  |  |
 | genbank_accession_id | text | No |  |  |
 | definition | text | Yes |  |  |
 | organism | text | Yes |  |  |
 | taxonomy | text | Yes |  |  |
 | journal | text | Yes |  |  |
 | collection_date | text | Yes |  |  |
 | country | text | Yes |  |  |


### NICKNAMES_LIST Table

**Description:** Unique Adjustive-Noun Nickname for each protein.primary.key | Adjective + NounFrom: https://gist.github.com/hugsy/8910dc78d208e40de42deb29e62df913#file-english-adjectives-txt

**Estimated rows:** ~2,040,636

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | nickname | character varying | No |  |  |


### NR_CATALYTIC_ORFS Table

**Description:** PAZy homologs retrieved from the NCBI non-redundant protein database (NR) using DIAMOND2 (--very-sensitive), which contain a complete catalytic domain, described as matching to all of the catalytic sites of a PAZy HMM.

**Estimated rows:** ~4,723,128

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | orf_id | bigint | No |  | foreign.key refers to "orf_origins"  The PETadex ORF integer accession of the amino acid sequence. |
 | genbank_accession_id | text | No |  | foreign.key  The GenBank protein accession of the amino acid sequence. |


### ORF_OFFSET Table

**Description:** Index table used by frontend for finding sequence from the fasta file quickly. Byte-offset and length references the full 300M fasta file in s3.

**Estimated rows:** ~307,155,744

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | orf_id | bigint | No |  |  |
 | byte_offset | bigint | No |  |  |
 | byte_length | integer | No |  |  |


### ORF_ORIGINS Table

**Description:** The database that the translated ORF was retrieved from, or the script the translated ORF was generated from. The origin in this table will direct which table to query for retrieval/generation information.

**Estimated rows:** ~307,155,744

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | orf_id | bigint | No |  | primary.key  The PETadex ORF integer accession of the amino acid sequence. |
 | orf_origin | smallint | No |  | The origin database/script the translated ORF was retrieved/generated from: 0 = PAZy 1 = NR 2 = Logan |
 | date_retrieval | date | No |  | The date (ISO-8601 format) the origin database was accessed or origin script was performed. |


### PAZY_CATALYTIC_ORFS Table

**Description:** Sequences retrieved from The Plastics-Active Enzymes Database (PAZy), which contain a complete catalytic domain, described as matching to all of the catalytic sites of a PAZy HMM.

**Estimated rows:** ~211

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | orf_id | bigint | No |  | foreign.key refers to "orf_origins"  The PETadex ORF integer accession of the amino acid sequence. |
 | genbank_accession_id | text | No |  | foreign.key  The GenBank protein accession of the amino acid sequence. |


### PAZY_HMMS Table

**Description:** PHMMs created through iterative cycles of adding catalytic domains to the profile (similar to HMMER). The seed profile is the PAZy sequences for the given component, while the full profile is the BLASTnr sequences.

**Estimated rows:** _not analyzed_

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | pazy_hmm_id | text | No |  | primary.key  PHMMs created through iterative cycles of adding catalytic domains to the profile (similar to HMMER). The seed profile is the PAZy sequences for the given component, while the full profile is the BLASTnr sequences. |
 | domain | text | No |  | The common name for the catalytic domain represented by the PAZy HMM. |
 | catalytic_residues | ARRAY | No |  | A functional descriptor for the catalytic residues that were used to determine if an ORF was catalytically complete. The order of this list is the same to "catalytic_match_states" |
 | catalytic_match_states | ARRAY | No |  | The match states of the catalytic residues that were used to determine if an ORF was catalytically complete. The order of this list is the same to "catalytic_residues". |


### PDB_ACCESSIONS Table

**Description:** Describes the pdb_id for each accession, which can be linked in the s3 bucket to the pdb file. It also contains relevant metadata. 

**Estimated rows:** ~358

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | pdb_id | character varying | No |  | Accession for the pdb | follows the format of two letters for technique, 8 digits for the batch, then _ and then the rank of the model ex. CF00000001_1 |
 | accession | character varying | Yes |  | accession of the aa sequence | refers to fastaa |
 | technique | character varying | Yes |  | The technique / model used to generate the structure ex. X-ray crystallography or AlphaFold3 |
 | relaxed | boolean | Yes |  | Relaxation of side chains to reduce steric clashes |
 | date_created | timestamp with time zone | Yes |  |  |
 | date_entered | timestamp with time zone | No | now() |  |
 | alignment | character varying | Yes |  | Script used for alignment to a reference structure |


### PETADEX_ANCESTRAL_RECONSTRUCTIONS Table

**Description:** No description available

**Estimated rows:** _not analyzed_

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | protein_id | integer | No |  |  |
 | script_used | text | No |  |  |
 | daughter_nodes | ARRAY | No |  |  |
 | date_entered | timestamp with time zone | Yes |  |  |


### PETADEX_CATALYTIC_DOMAINS Table

**Description:** HMMsearch (e-value < 1e-5) top hits using PAZy HMMs, which contain a complete catalytic domain, described as matching to all of the catalytic sites of the PAZy HMM.

**Estimated rows:** ~307,175,776

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | orf_id | bigint | No |  | foreign.key refers to "orf_origins"  The PETadex ORF integer accession of the amino acid sequence. |
 | pazy_hmm_id | text | No |  | foreign.key refers to "pazy_hmms"  PAZy HMM used for the HMMsearch. |
 | domain_start | integer | No |  | First amino acid in the ORF hit by the PAZy HMM. |
 | domain_end | integer | No |  | Last amino acid in the ORF hit by the PAZy HMM. |
 | catalytic_residues | ARRAY | No |  | Amino acid identities corresponding the catalytic match states of the PAZy HMM. Order of the catalytic match states determined by "pazy_hmms". |
 | date_performed | date | No | '2026-05-09'::date | The date (ISO-8601) HMMsearch was ran. |


### PETADEX_CLUSTERING Table

**Description:** No description available

**Estimated rows:** ~307,161,920

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | orf_id | bigint | No |  |  |
 | 90pid_enzyme_id | bigint | No |  |  |
 | 60pid_family_id | integer | No |  |  |
 | 30pid_superfamily_id | integer | No |  |  |


### PETADEX_LOGAN_MATCHES Table

**Description:** No description available

**Estimated rows:** ~448,790

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | parent | text | Yes |  |  |
 | child | text | Yes |  |  |
 | identity | numeric | Yes |  |  |
 | pval | numeric | Yes |  |  |


### PETADEX_SYNTHETIC Table

**Description:** Contains the amino acid sequences of proteins (synthetic proteins) derived from natural proeins.

**Estimated rows:** ~64

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | protein_id | integer | No |  | primary.key  A protein is an amino acid sequence that is either found from natural diversity: Logan, NCBI NR, or is engineered. |
 | wildtype_protein_id | ARRAY | Yes |  | The natural protein sequences that the synthetic protein was derived from. The first listed protein is the amino acid sequence that the genotype is applied to. |
 | genotype | text | No |  | HGVS nomenclature describing changes applied to a natural amino acid sequence to produce the synthetic protein amino acid sequence. |
 | genotype_description | text | No |  | Brief description of the genotype / mutations made to the natural amino acid sequence. |
 | engineering_method | text | No |  | The method/script used to create/engineer the amino acid sequence. |
 | date_entered | timestamp with time zone | No | now() |  |
 | parent_accessions | ARRAY | Yes |  |  |


### PLASTIC_KINETICS_PUBLISHED Table

**Description:** No description available

**Estimated rows:** ~464

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | Class | text | Yes |  |  |
 | Experiment_# | text | Yes |  |  |
 | Enzyme | text | Yes |  |  |
 | Species | text | Yes |  |  |
 | Uniprot | text | Yes |  |  |
 | GenBank | text | Yes |  |  |
 | AA_SEQ | text | Yes |  |  |
 | Substrate | text | Yes |  |  |
 | Substrate_SMILES | text | Yes |  |  |
 | Temperature | text | Yes |  |  |
 | pH | text | Yes |  |  |
 | Kcat_(/s) | text | Yes |  |  |
 | Km_(M) | text | Yes |  |  |
 | Kcat/Km_(/s/M) | text | Yes |  |  |
 | Paper_# | text | Yes |  |  |
 | Pubmed | text | Yes |  |  |
 | DOI | text | Yes |  |  |
 | Supplemental | text | Yes |  |  |


### PLATE_DATA Table

**Description:** This table contains all of the output readings from plate assays.

**Estimated rows:** ~27,648

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | id | integer | No | nextval('plate_data_id_seq'::regclass) | primary.key |
 | gene | character varying | No |  | Refers to the gene that was measured |
 | plate | character varying | No |  | Refers to the plate identifier |
 | plasmid | character varying | Yes |  | Plasmid expressing the gene |
 | column | integer | Yes |  | Column of the well |
 | row | integer | Yes |  | Row of the well |
 | normalization_method | character varying | No | 'raw'::character varying | The script name that normalized the data: default = "raw" |
 | readout_value | real | Yes |  | Measured value of the well |
 | colony_size | real | Yes |  | Size of the colony |
 | date_entered | timestamp with time zone | No | now() | Auto |
 | measurement_type | character varying | Yes |  | The type of measurement the readout value corresponds to |


### PLATE_METADATA Table

**Description:** Relevant metadata for each plate assay, corresponds to plate_data

**Estimated rows:** ~72

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | plate | character varying | No |  | Unique plate identifier with 2 letters (refers to the assay name) and 5 numbers (refers to order of plate creation) ex. HA00001 |
 | exp_id | character varying | Yes |  | Name of the experiment the plate was included in |
 | exp_description | character varying | Yes |  | A brief description of the experiment |
 | media | character varying | Yes |  | Unique media identifier |
 | timepoint_hours | real | Yes |  | Time read after pinning |
 | temp_celsius | real | Yes |  | Temperature of plate, if constant |
 | ph | real | Yes |  | pH of the media, if constant |
 | organism | character varying | Yes |  | Organism expressing the plasmid |
 | control_genes | ARRAY | Yes |  | List of genes acting as the controls |
 | operator | character varying | Yes |  | Who was in charge of the plate |
 | date_created | timestamp with time zone | Yes |  | Plate pinning date |
 | date_read | timestamp with time zone | Yes |  | Plate readout date |
 | date_entered | timestamp with time zone | No | now() | Auto |


### PREDICTED_SIGNAL_SEQUENCE Table

**Description:** Signal peptide cleavage site of sequences from fastaa table predicted using SignalP 6.0

**Estimated rows:** ~516

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | id | integer | No | nextval('predicted_signal_sequence_id_seq'::regclass) |  |
 | accession | text | Yes |  | Unique accession of the sequence |
 | organism | text | Yes |  | Organism type selected in SignalP 6.0 to guide prediction (Eukarya/Other) |
 | signal | text | Yes |  | Type of signal sequence predited by SignalP 6.0 Possible signal for "Eukarya": SP(Sec/SPI), OTHER Possible signal for "Other": SP(Sec/SPI), LIPO(Sec/SPII), TAT(Tat/SPI), TATLIPO(Tat/SPII), PILIN(Sec/SPIII) |
 | cleave_after | text | Yes |  | The site after which cleavage occurs |
 | p(signal) | double precision | Yes |  | Probability of the predicted signal type |
 | p(cleavage) | double precision | Yes |  | Probability of the predicted cleavage site |
 | script | text | Yes |  | Bash script used to automatically run the signal sequence prediction |
 | date_entered | timestamp with time zone | Yes | now() | Automatically generated |


### SARA_DOMAINS Table

**Description:** No description available

**Estimated rows:** _not analyzed_

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | orf_id | integer | No |  |  |
 | domain | text | No |  |  |
 | start_aa | integer | No |  |  |
 | end_aa | integer | No |  |  |


### SARA_IMPORTANT_MOTFIS Table

**Description:** No description available

**Estimated rows:** _not analyzed_

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | orf_id | integer | No |  |  |
 | motif | text | No |  |  |
 | start_aa | integer | No |  |  |
 | end_aa | integer | No |  |  |


### SARA_SIGNAL_SEQUENCES Table

**Description:** No description available

**Estimated rows:** _not analyzed_

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | orf_id | integer | No |  |  |
 | signal_sequence | text | No |  |  |
 | cleavage_site | integer | No |  |  |


### SEQ_COLLECTIONS Table

**Description:** Contains the dataset that each sequence is associated with [pulled from]

**Estimated rows:** ~4

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | collection | text | No |  | primary.key | Experiment / collection name |
 | accession_list | ARRAY | No |  | All accessions within the collection ['{name}', '{name}', {...}] |
 | created_at | timestamp with time zone | No | now() | Auto |


### SIGNALP6_ORF_PREDICTIONS Table

**Description:** SignalP 6.0 positive signal peptide predictions for PETadex ORFs. organism_mode=other, model_mode=fast. Contains only ORFs with detected signal types; OTHER rows are intentionally omitted.

**Estimated rows:** ~42,530,188

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | orf_id | bigint | No |  | The PETadex ORF integer accession. Refers to orf_origins.orf_id. |
 | top_signal | smallint | No |  | SignalP 6.0 predicted signal type, encoded as SMALLINT:     1 = SP      (Sec/SPI classical signal peptide)     2 = LIPO    (Sec/SPII lipoprotein signal peptide)     3 = TAT     (Tat/SPI twin-arginine signal peptide)     4 = TATLIPO (Tat/SPII twin-arginine lipoprotein)     5 = PILIN   (Sec/SPIII pilin-type signal peptide).     OTHER=0 rows are intentionally not stored in this positive-hit table. |
 | signal_prob | real | No |  | Probability of the winning SignalP 6.0 signal class, 0–1. |
 | cleavage_pos | smallint | Yes |  | Residue position after which cleavage is predicted, e.g. 21 = cut after residue 21. Also the signal peptide length in residues. |
 | cleavage_prob | real | Yes |  | SignalP 6.0 cleavage-site probability, 0–1. May be NULL if unavailable for rare positive-label edge cases. |


### SPATIAL_REF_SYS Table

**Description:** No description available

**Estimated rows:** ~8,500

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | srid | integer | No |  |  |
 | auth_name | character varying(256) | Yes |  |  |
 | auth_srid | integer | Yes |  |  |
 | srtext | character varying(2048) | Yes |  |  |
 | proj4text | character varying(2048) | Yes |  |  |


### SRA_METADATA Table

**Description:** No description available

**Estimated rows:** ~8,273,489

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | acc | text | No |  |  |
 | assay_type | text | Yes |  |  |
 | center_name | text | Yes |  |  |
 | consent | text | Yes |  |  |
 | experiment | text | Yes |  |  |
 | sample_name | text | Yes |  |  |
 | instrument | text | Yes |  |  |
 | librarylayout | text | Yes |  |  |
 | libraryselection | text | Yes |  |  |
 | librarysource | text | Yes |  |  |
 | platform | text | Yes |  |  |
 | sample_acc | text | Yes |  |  |
 | biosample | text | Yes |  |  |
 | organism | text | Yes |  |  |
 | sra_study | text | Yes |  |  |
 | releasedate | timestamp without time zone | Yes |  |  |
 | bioproject | text | Yes |  |  |
 | mbytes | integer | Yes |  |  |
 | avgspotlen | integer | Yes |  |  |
 | mbases | integer | Yes |  |  |
 | library_name | text | Yes |  |  |
 | biosamplemodel_sam | text | Yes |  |  |
 | collection_date_sam | text | Yes |  |  |
 | geo_loc_name_country_calc | text | Yes |  |  |
 | geo_loc_name_country_continent_calc | text | Yes |  |  |
 | geo_loc_name_sam | text | Yes |  |  |
 | latitude | double precision | Yes |  |  |
 | longitude | double precision | Yes |  |  |
 | elevation | numeric | Yes |  |  |
 | country | text | Yes |  |  |
 | biome | text | Yes |  |  |
 | confidence | numeric | Yes |  |  |


### SYNTHESIZED_CONSTRUCTS Table

**Description:** Contains the vectors used to produce a synthetic protein, whose activity has been validated in a wet lab experiment.

**Estimated rows:** ~516

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | construct_id | integer | No | nextval('construct_id_seq'::regclass) | Primary key. Unique identifier for plasmids cloned with the genes of synthetic proteins. The identifier is an integer. |
 | plasmid | character varying | Yes |  | The name of the plasmid backbone that the gene for the synthetic protein was cloned into.  |
 | protein_id | integer | Yes |  | Foreign key to table "protein_sources". The protein that corresponds to the translated open reading frame of the plasmid. |
 | codon_optimized_organism | character varying | Yes |  | The organism that the gene of the synthetic protein was optimized for. |
 | nt_gene_orf | character varying | Yes |  | The main nucleotide sequence of the open reading frame of the plasmid. This is the reverse translated sequence of the synthetic protein. The nucleotide sequence may include a portion of the backbone if that sequence is present within the open reading frame. |
 | date_entered | timestamp with time zone | Yes | now() | Automatically enters the current date. |
 | temp_accession | text | Yes |  |  |
 | temp_gene | text | Yes |  |  |


### TEMP_LOGAN_ANCESTRAL Table

**Description:** No description available

**Estimated rows:** _not analyzed_

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | accession | character varying | Yes |  |  |
 | protein_id | bigint | Yes |  |  |


### TEMP_LOGAN_GENBANK Table

**Description:** No description available

**Estimated rows:** ~90

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | accession | character varying | Yes |  |  |
 | protein_id | integer | Yes |  |  |


### TEMP_LOGAN_LOGAN Table

**Description:** No description available

**Estimated rows:** ~180

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | accession | character varying | Yes |  |  |
 | aa_sequence | character varying | Yes |  |  |


### TEMP_LOGAN_SYNTHETIC Table

**Description:** No description available

**Estimated rows:** ~64

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | protein_id | bigint | Yes |  |  |
 | accession | character varying | Yes |  |  |
 | genotype | character varying | Yes |  |  |
 | genotype_description | character varying | Yes |  |  |
 | parent_accessions | ARRAY | Yes |  |  |


### TEMP_LOGAN_UNIPROT Table

**Description:** No description available

**Estimated rows:** _not analyzed_

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | accession | character varying | Yes |  |  |
 | aa_sequence | character varying | Yes |  |  |
 | source | character varying | Yes |  |  |
 | date_entered | timestamp with time zone | Yes |  |  |
 | genotype | character varying | Yes |  |  |
 | genotype_description | character varying | Yes |  |  |
 | synthetic | boolean | Yes |  |  |
 | parent_accessions | ARRAY | Yes |  |  |
 | parent_genes | ARRAY | Yes |  |  |
 | synonyms | ARRAY | Yes |  |  |
 | in_gene_metadata | boolean | Yes |  |  |


### UPDATE_PLATE_DATA Table

**Description:** This table contains all of the output readings from plate assays.

**Estimated rows:** _not analyzed_

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | well_id | integer | No | nextval('well_id_seq'::regclass) | Primary key. A unique identifier for the measurement taken/calculated from a plate assay. |
 | construct_id | integer | No |  | Foreign key to table "synthesized_constructs". Unique identifier for plasmids cloned with the genes of synthetic proteins. |
 | plate_id | character varying(7) | No |  | Foreign key to table "plate_metadata". The unique identifier for a plate measurement at a give time. |
 | column | integer | Yes |  | The physical column of the well within the plate. |
 | row | integer | Yes |  | The physical row of the well within the plate. |
 | readout_value | real | No |  | The measured value of the given well. The datatype being measured is held in the column "readout_type". |
 | readout_type | character varying | No |  | The type of data that is being measured, specifying any modifications made to the data. ex. melting temperature |
 | modification_script | character varying | Yes |  | The script used to modify the measurement, held within the GitHub repository "petadex/root/scripts". |
 | date_entered | timestamp with time zone | No | now() | Automatically enters the current date. |


### VARIANT_DICTIONARY Table

**Description:** No description available

**Estimated rows:** ~2,735,959

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | enzyme_id | integer | Yes |  |  |
 | library_id | integer | Yes |  |  |
 | contig_id | integer | Yes |  |  |
 | genbank_accession_id | character varying(16) | Yes |  |  |
 | variant_id | integer | No | nextval('variant_id_seq'::regclass) |  |
 | enzyme_pid | numeric | Yes |  |  |


## Materialized Views

### ACCESSION_ACTIVITY_VIEW Materialized View

**Description:** No description available

**Estimated rows:** ~27,648

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | plasmid | character varying | Yes |  |  |
 | readout_value | real | Yes |  |  |
 | timepoint_hours | real | Yes |  |  |
 | gene | character varying | Yes |  |  |
 | plate | character varying | Yes |  |  |
 | column | integer | Yes |  |  |
 | row | integer | Yes |  |  |
 | media | character varying | Yes |  |  |
 | exp_name | character varying | Yes |  |  |
 | accession | character varying | Yes |  |  |
 | source | character varying | Yes |  |  |
 | synonyms | character varying[] | Yes |  |  |


### BLOCK_30PID Materialized View

**Description:** No description available

**Estimated rows:** ~22,235

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | level | text | Yes |  |  |
 | cluster_id | integer | Yes |  |  |
 | centroid_orf_id | bigint | Yes |  |  |
 | centroid_accession | character varying | Yes |  |  |
 | centroid_source | text | Yes |  |  |
 | member_count | bigint | Yes |  |  |
 | child_count | bigint | Yes |  |  |
 | n_pazy | bigint | Yes |  |  |
 | n_nr | bigint | Yes |  |  |
 | n_sra | bigint | Yes |  |  |
 | dominant_organism | text | Yes |  |  |
 | distinct_organism_count | bigint | Yes |  |  |
 | centroid_cath_domain | text | Yes |  |  |
 | centroid_domain_name | text | Yes |  |  |
 | centroid_component | integer | Yes |  |  |


### BLOCK_60PID Materialized View

**Description:** No description available

**Estimated rows:** ~1,814,166

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | level | text | Yes |  |  |
 | cluster_id | bigint | Yes |  |  |
 | centroid_orf_id | bigint | Yes |  |  |
 | centroid_accession | character varying | Yes |  |  |
 | centroid_source | text | Yes |  |  |
 | member_count | bigint | Yes |  |  |
 | child_count | bigint | Yes |  |  |
 | n_pazy | bigint | Yes |  |  |
 | n_nr | bigint | Yes |  |  |
 | n_sra | bigint | Yes |  |  |
 | dominant_organism | text | Yes |  |  |
 | distinct_organism_count | bigint | Yes |  |  |
 | centroid_cath_domain | text | Yes |  |  |
 | centroid_domain_name | text | Yes |  |  |
 | centroid_component | integer | Yes |  |  |


### BLOCK_90PID Materialized View

**Description:** No description available

**Estimated rows:** ~18,172,364

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | level | text | Yes |  |  |
 | cluster_id | bigint | Yes |  |  |
 | centroid_orf_id | integer | Yes |  |  |
 | centroid_accession | character varying | Yes |  |  |
 | centroid_source | text | Yes |  |  |
 | member_count | bigint | Yes |  |  |
 | child_count | bigint | Yes |  |  |
 | n_pazy | bigint | Yes |  |  |
 | n_nr | bigint | Yes |  |  |
 | n_sra | bigint | Yes |  |  |
 | dominant_organism | text | Yes |  |  |
 | distinct_organism_count | bigint | Yes |  |  |
 | centroid_cath_domain | text | Yes |  |  |
 | centroid_domain_name | text | Yes |  |  |
 | centroid_component | integer | Yes |  |  |


### CORPUS_SUMMARY Materialized View

**Description:** No description available

**Estimated rows:** ~1

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | id | integer | Yes |  |  |
 | total_enzymes | bigint | Yes |  |  |
 | total_families | bigint | Yes |  |  |
 | total_components | bigint | Yes |  |  |
 | total_variants | bigint | Yes |  |  |
 | pazy_total | bigint | Yes |  |  |
 | nr_total | bigint | Yes |  |  |
 | sra_total | bigint | Yes |  |  |
 | catalytic_core_total | bigint | Yes |  |  |
 | clusters_90pid | bigint | Yes |  |  |
 | clusters_60pid | bigint | Yes |  |  |
 | clusters_30pid | bigint | Yes |  |  |


### ENZYME_FAMILY_SUMMARY Materialized View

**Description:** No description available

**Estimated rows:** ~64,730

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | family_id | integer | Yes |  |  |
 | centroid_accession | character varying(16) | Yes |  |  |
 | variant_count | bigint | Yes |  |  |
 | component_count | bigint | Yes |  |  |
 | avg_identity | numeric | Yes |  |  |


### ENZYME_METADATA_EXPLORER Materialized View

**Description:** No description available

**Estimated rows:** ~1,047,216

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | enzyme_id | integer | Yes |  |  |
 | genbank_accession_id | character varying(16) | Yes |  |  |
 | translated_sequence | text | Yes |  |  |
 | library_id | character varying | Yes |  |  |
 | contig_id | integer | Yes |  |  |
 | orf_start | integer | Yes |  |  |
 | orf_end | integer | Yes |  |  |
 | orf_type | smallint | Yes |  |  |
 | orf_type_description | text | Yes |  |  |
 | family | integer | Yes |  |  |
 | family_pid | numeric | Yes |  |  |
 | component | integer | Yes |  |  |
 | logan_protein_id | text | Yes |  |  |
 | logan_nucleotide_accession | text | Yes |  |  |
 | biosample | text | Yes |  |  |
 | organism | text | Yes |  |  |
 | bioproject | text | Yes |  |  |
 | country | text | Yes |  |  |
 | continent | text | Yes |  |  |
 | geographic_location | text | Yes |  |  |
 | coordinates | geometry | Yes |  |  |
 | elevation | numeric | Yes |  |  |
 | biome | text | Yes |  |  |
 | location_confidence | numeric | Yes |  |  |


### ENZYME_STATS_OVERVIEW Materialized View

**Description:** No description available

**Estimated rows:** _not analyzed_

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | total_enzymes | bigint | Yes |  |  |
 | total_families | bigint | Yes |  |  |
 | total_components | bigint | Yes |  |  |
 | total_variants | bigint | Yes |  |  |


### FAMILY_ATLAS Materialized View

**Description:** No description available

**Estimated rows:** ~64,730

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | family_id | integer | Yes |  |  |
 | umap_x | double precision | Yes |  |  |
 | umap_y | double precision | Yes |  |  |
 | family_size | integer | Yes |  |  |
 | genbank_accession_id | text | Yes |  |  |
 | definition | text | Yes |  |  |
 | organism | text | Yes |  |  |
 | taxonomy | text | Yes |  |  |
 | journal | text | Yes |  |  |
 | collection_date | text | Yes |  |  |
 | country | text | Yes |  |  |
 | component | integer | Yes |  |  |
 | cath_domain | text | Yes |  |  |
 | domain_name | text | Yes |  |  |


### FAMILY_REPRESENTATIVES Materialized View

**Description:** No description available

**Estimated rows:** ~64,730

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | family_id | integer | Yes |  |  |
 | accession | character varying(16) | Yes |  |  |
 | sequence | text | Yes |  |  |


### LOGAN_MATCH_NODES_UNIQUE Materialized View

**Description:** No description available

**Estimated rows:** ~367,104

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | parent | text | Yes |  |  |
 | child | text | Yes |  |  |
 | identity | numeric | Yes |  |  |
 | pval | numeric | Yes |  |  |


### PLASTIC_KINETICS_PUBLISHED_SUMMARY Materialized View

**Description:** No description available

**Estimated rows:** ~464

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | Class | text | Yes |  |  |
 | Enzyme | text | Yes |  |  |
 | Species | text | Yes |  |  |
 | Uniprot | text | Yes |  |  |
 | GenBank | text | Yes |  |  |
 | AA_SEQ | text | Yes |  |  |
 | Substrate | text | Yes |  |  |
 | Substrate_SMILES | text | Yes |  |  |
 | Temperature | text | Yes |  |  |
 | pH | text | Yes |  |  |
 | Kcat_(/s) | text | Yes |  |  |
 | Km_(M) | text | Yes |  |  |
 | Pubmed | text | Yes |  |  |
 | DOI | text | Yes |  |  |
 | Supplemental | text | Yes |  |  |


### PLATE_ACTIVITY_VIEW Materialized View

**Description:** No description available

**Estimated rows:** _not analyzed_

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | id | integer | Yes |  |  |
 | gene | character varying | Yes |  |  |
 | plate | character varying | Yes |  |  |
 | readout_value | real | Yes |  |  |
 | measurement_type | character varying | Yes |  |  |
 | normalization_method | character varying | Yes |  |  |
 | colony_size | real | Yes |  |  |
 | row | integer | Yes |  |  |
 | column | integer | Yes |  |  |
 | plasmid | character varying | Yes |  |  |
 | exp_id | character varying | Yes |  |  |
 | exp_description | character varying | Yes |  |  |
 | media | character varying | Yes |  |  |
 | timepoint_hours | real | Yes |  |  |
 | temp_celsius | real | Yes |  |  |
 | ph | real | Yes |  |  |
 | organism | character varying | Yes |  |  |
 | control_genes | character varying[] | Yes |  |  |
 | operator | character varying | Yes |  |  |
 | date_created | timestamp with time zone | Yes |  |  |
 | date_read | timestamp with time zone | Yes |  |  |
 | data_date_entered | timestamp with time zone | Yes |  |  |
 | metadata_date_entered | timestamp with time zone | Yes |  |  |


### SEARCH_INDEX Materialized View

**Description:** No description available

**Estimated rows:** ~319,513,664

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | match_value | text | Yes |  |  |
 | match_type | text | Yes |  |  |
 | result_kind | text | Yes |  |  |
 | orf_id | bigint | Yes |  |  |
 | c90_id | bigint | Yes |  |  |
 | c60_id | integer | Yes |  |  |
 | c30_id | integer | Yes |  |  |


### WITH_SRA_AND_BIOSAMPLE_LOC_METADATA Materialized View

**Description:** No description available

**Estimated rows:** ~394

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | accession | character varying | Yes |  |  |
 | biosample | text | Yes |  |  |
 | organism | text | Yes |  |  |
 | bioproject | text | Yes |  |  |
 | geo_loc_name_country_calc | text | Yes |  |  |
 | geo_loc_name_country_continent_calc | text | Yes |  |  |
 | geo_loc_name_sam | text | Yes |  |  |
 | lat_lon | geometry | Yes |  |  |
 | elevation | numeric | Yes |  |  |
 | country | text | Yes |  |  |
 | biome | text | Yes |  |  |
 | confidence | numeric | Yes |  |  |


### WITH_SRA_METADATA Materialized View

**Description:** No description available

**Estimated rows:** ~192

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | accession | character varying | Yes |  |  |
 | aa_sequence | character varying | Yes |  |  |
 | source | character varying | Yes |  |  |
 | date_entered | timestamp with time zone | Yes |  |  |
 | genotype | character varying | Yes |  |  |
 | genotype_description | character varying | Yes |  |  |
 | synthetic | boolean | Yes |  |  |
 | parent_accessions | character varying[] | Yes |  |  |
 | parent_genes | character varying[] | Yes |  |  |
 | synonyms | character varying[] | Yes |  |  |
 | acc | text | Yes |  |  |
 | sample_acc | text | Yes |  |  |
 | biosample | text | Yes |  |  |
 | organism | text | Yes |  |  |
 | sra_study | text | Yes |  |  |
 | release_date | timestamp without time zone | Yes |  |  |
 | bioproject | text | Yes |  |  |
 | biosamplemodel_sam | text | Yes |  |  |
 | collection_date_sam | date | Yes |  |  |
 | geo_loc_name_country_calc | text | Yes |  |  |
 | geo_loc_name_country_continent_calc | text | Yes |  |  |
 | geo_loc_name_sam | text | Yes |  |  |


### WITH_STAT_DATA Materialized View

**Description:** petadex sequences labeled with sra stat data

**Estimated rows:** ~162

| Column | Type | Nullable | Default | Comment |
|--------|------|----------|---------|---------|
 | accession | character varying | Yes |  |  |
 | names | text[] | Yes |  |  |
 | kmer_percs | numeric[] | Yes |  |  |


