#!/usr/bin/env python3
"""Curated manifest content fields for the included records (spec v1.0 section 3).

Every value below was read off text returned by the Europe PMC API (title,
abstractText, or fullTextXML) for that record. `src` names which of those the
row's snippet came from: "abstract", "fullText", or "abstract_only" for records
whose full text could not be retrieved. NOT_FOUND means the retrieved text does
not state it -- never that the value is absent in reality.

Field keys map 1:1 onto the spec schema. Tissue flags are Y/N; where a Y is
recorded the snippet evidences that tissue type.
"""

# key -> dict of curated fields
M = {
 # ============================== 2020 =========================================
 "32061950": dict(
   src="abstract_only", modality="scRNA;scTCR",
   primary="Y", ln="N", dist="N", normal="N", blood="N",
   n_pat="3", n_samp="3", n_cells="NOT_FOUND",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="NOT_FOUND", nk_claim="NOT_FOUND", ebv="none",
   acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="we performed single cell RNA-seq and analyzed tumor cells together with the infiltrating immune cells from three NPC tumor tissues",
   note="full text not retrievable from Europe PMC; NK not assessable from abstract"),

 "32686767": dict(
   src="fullText", modality="scRNA",
   primary="Y", ln="N", dist="N", normal="Y", blood="N",
   n_pat="15", n_samp="16", n_cells="47866",
   treat="naive", longit="N", mstage="M0_only",
   nk_depth="subclustered",
   nk_claim="NK cell signature significantly associated with improved survival outcomes in NPC",
   ebv="continuous", acc="CNP0000428;GSE150430", acc_blood="N",
   ebv_kind="plasma_EBV_DNA", ebv_snip='We further explored their correlations with plasma EBV DNA concentrations.',
   snip="Workflow diagram showing the collection and processing of fresh biopsy samples from 15 primary NPC tumors and one normal tissue for scRNA-seq",
   note="NK trap check: 'non-keratinizing' here is a histology descriptor of the enrolled patients, not NK cells; NK cells are separately analysed as natural killer cells in T/NK clusters"),

 "32901110": dict(
   src="abstract_only", modality="scRNA",
   primary="Y", ln="N", dist="N", normal="Y", blood="N",
   n_pat="NOT_FOUND", n_samp="26", n_cells="~104000",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="NOT_FOUND", nk_claim="NOT_FOUND", ebv="categorical",
   acc="NOT_FOUND", acc_blood="UNKNOWN",
   ebv_kind="tumour_EBV_status", ebv_snip='we performed single-cell RNA sequencing on ~104,000 cells from 19 EBV + NPCs and 7 nonmalignant nasopharyngeal biopsies',
   snip="Here we performed single-cell RNA sequencing on ~104,000 cells from 19 EBV + NPCs and 7 nonmalignant nasopharyngeal biopsies",
   note="19 EBV+ NPC + 7 nonmalignant biopsies = 26 samples; patient count not stated in abstract"),

 # ============================== 2021 =========================================
 "33531485": dict(
   src="fullText", modality="scRNA;scTCR",
   primary="Y", ln="N", dist="N", normal="N", blood="Y",
   n_pat="10", n_samp="20", n_cells="176447",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="subclustered",
   nk_claim="NK cells were decreased in the tumours compared to the PBMC",
   ebv="categorical", acc="NOT_FOUND", acc_blood="UNKNOWN",
   ebv_kind="tumour_EBV_status", ebv_snip='malignant cells with different Epstein-Barr virus infection status',
   snip="a total of 176,447 cells were identified from the 10 patients (including 82,622 and 93,825 for tumours and PBMC, respectively)",
   note="10 tumour-blood pairs = 20 samples; NK subclusters named (e.g. NK_C2_FCER1G). GSE102349/GSE121600 cited are external datasets, not this study's deposit"),

 "33545035": dict(
   src="abstract_only", modality="scRNA",
   primary="NOT_FOUND", ln="NOT_FOUND", dist="NOT_FOUND", normal="NOT_FOUND", blood="NOT_FOUND",
   n_pat="210", n_samp="NOT_FOUND", n_cells="NOT_FOUND",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="NOT_FOUND", nk_claim="NOT_FOUND", ebv="none",
   acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="by performing a pan-cancer analysis of single myeloid cells from 210 patients across 15 human cancer types ... Mast cells in nasopharyngeal cancer were found to be associated with better prognosis",
   note="pan-cancer myeloid atlas; 210 patients is the whole pan-cancer cohort, NPC-specific counts not stated in retrieved text"),

 "33750785": dict(
   src="fullText", modality="scRNA;scTCR;scBCR",
   primary="Y", ln="N", dist="N", normal="Y", blood="N",
   n_pat="14", n_samp="NOT_FOUND", n_cells="66627",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="subclustered",
   nk_claim="NK cells in the NPC microenvironment remained immune activated",
   ebv="none", acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="we performed 5' single-cell RNA sequencing integrated with V(D)J profiling on 14 patients with either nasopharyngeal carcinoma (NPC) or nasopharyngeal lymphatic hyperplasia (NLH)",
   note="normal comparator is nasopharyngeal lymphatic hyperplasia (NLH); dedicated NK results section with MAST profiling and exhaustion-marker assessment"),

 "34253636": dict(
   src="fullText", modality="scRNA",
   primary="Y", ln="N", dist="N", normal="Y", blood="Y",
   n_pat="NOT_FOUND", n_samp="NOT_FOUND", n_cells="NOT_FOUND",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="absent", nk_claim="NOT_FOUND", ebv="none",
   acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="Mononuclear cells were isolated from the blood and tissues as described previously. Peripheral blood mononuclear cells (PBMCs) were isolated using Ficoll density gradient centrifugation.",
   note="blood used for FACS and in vitro differentiation; retrieved text does not state whether scRNA-seq itself covered PBMC"),

 # ============================== 2022 =========================================
 "35096485": dict(
   src="fullText", modality="scRNA",
   primary="Y", ln="N", dist="N", normal="N", blood="N",
   n_pat="6", n_samp="6", n_cells="62164",
   treat="mixed", longit="N", mstage="unstated",
   nk_depth="subclustered",
   nk_claim="NK cells grouped into CD56brightCD16- and CD56dimCD16+ cytotoxic subtypes",
   ebv="none", acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="we generated scRNA-seq profiles for nasopharyngeal tumors from the six patients with NPC, four treatment-naive, and two recurrent samples",
   note="abstract states ~60,000 cells; full text reports 62,164 cells after QC. NK further grouped into five clusters"),

 # ============================== 2023 =========================================
 "36739462": dict(
   src="fullText", modality="scRNA",
   primary="Y", ln="Y", dist="Y (site: liver/lung/bone/other - de novo metastatic)", normal="N", blood="Y",
   n_pat="2", n_samp="11", n_cells="NOT_FOUND",
   treat="mixed", longit="N", mstage="M1_included",
   nk_depth="absent", nk_claim="NOT_FOUND", ebv="none",
   acc="HRA000034;HRA000035;HRA000036", acc_blood="N",
   snip="we perform an integrative genomic analysis of 163 matched blood and primary, regional lymph node metastasis and distant metastasis tumour samples, combined with single-cell RNA-seq on 11 samples from two patients",
   note="163 matched blood+tissue samples are genomic (WES/WGS); the single-cell component is 11 samples from 2 patients. HRA000036 is the single-cell accession"),

 "37280275": dict(
   src="abstract_only", modality="scRNA;scTCR;scBCR",
   primary="Y", ln="N", dist="N", normal="N", blood="N",
   n_pat="15", n_samp="30", n_cells="NOT_FOUND",
   treat="mixed", longit="Y (treatment-naive and post-GP chemotherapy, matched pairs)", mstage="unstated",
   nk_depth="NOT_FOUND", nk_claim="NOT_FOUND", ebv="none",
   acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="using single-cell RNA sequencing and T cell and B cell receptor sequencing of matched, treatment-naive and post-GP chemotherapy NPC samples (n = 15 pairs)",
   note="15 matched pairs = 30 samples; full text not retrievable from Europe PMC"),

 "37349991": dict(
   src="fullText", modality="scRNA;snRNA;CyTOF",
   primary="Y", ln="N", dist="N", normal="Y", blood="N",
   n_pat="9", n_samp="10", n_cells="112442",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="NK_is_theme",
   nk_claim="a CD8+ natural killer (NK) cell cluster that is specific to the NPC TME",
   ebv="none", acc="HRA003340", acc_blood="N",
   snip="we collected three non-malignant tissues and seven NPC samples from nine donors via clinical biopsy",
   note="50,439 single cells + 62,003 single nuclei = 112,442; plus 670,044 cells from imaging mass cytometry. CyTOF here is imaging mass cytometry (IMC) on tissue"),

 # ============================== 2024 =========================================
 "38547605": dict(
   src="abstract_only", modality="scRNA;scTCR;scBCR",
   primary="Y", ln="N", dist="N", normal="Y", blood="N",
   n_pat="10", n_samp="10", n_cells="NOT_FOUND",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="NOT_FOUND", nk_claim="NOT_FOUND", ebv="none",
   acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="we analyzed data from 7 patients with NPC and 3 patients with nasopharyngeal lymphatic hyperplasia (NLH)",
   note="7 NPC + 3 NLH; provenance of the data not stated in the retrieved abstract"),

 "39197193": dict(
   src="abstract_only", modality="scRNA",
   primary="Y", ln="N", dist="N", normal="N", blood="N",
   n_pat="NOT_FOUND", n_samp="NOT_FOUND", n_cells="NOT_FOUND",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="NOT_FOUND", nk_claim="NOT_FOUND", ebv="none",
   acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="Single-cell transcriptome sequencing analysis of human nasopharyngeal carcinoma (NPC) revealed a significant enrichment of B cell subset characterized by high expression of EGR1 and EGR3",
   note="cohort size not stated in the retrieved abstract; full text not retrievable"),

 "39205595": dict(
   src="abstract_only", modality="scRNA",
   primary="Y", ln="N", dist="N", normal="N", blood="N",
   n_pat="NOT_FOUND", n_samp="NOT_FOUND", n_cells="NOT_FOUND",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="NOT_FOUND", nk_claim="NOT_FOUND", ebv="none",
   acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="using single-cell RNA sequencing techniques, we investigated the cellular landscape in NPC and oral cancers",
   note="NPC and oral cancer compared; NPC-specific cohort size not stated in the retrieved abstract"),

 "39231979": dict(
   src="fullText", modality="scRNA;scTCR;scBCR;spatial",
   primary="Y", ln="N", dist="N", normal="Y", blood="Y",
   n_pat="56", n_samp="77", n_cells="343829",
   treat="mixed", longit="N", mstage="unstated",
   nk_depth="cluster_only",
   nk_claim="20,051 NK cells identified among the annotated clusters",
   ebv="categorical", acc="HRA006885;zenodo.12728503", acc_blood="Y",
   ebv_kind="tumour_EBV_status", ebv_snip='All NPC tumours were EBV positive, as confirmed using in situ hybridisation of EBV encoded small RNAs (EBERs) in tumour biopsy.',
   snip="we performed single-cell RNA sequencing (scRNA-seq) analysis of transcriptome and immune cell receptor profiles for 77 samples, including 56 tumours and 10 peripheral blood mononuclear cells (PBMC) from 56 patients with NPC and 11 nasopharyngeal non-cancerous (NPN) tissues",
   note="73,477 of the 343,829 cells are PBMC; spatial transcriptomics covers 31,316 spots from 15 tumours"),

 "39415331": dict(
   src="fullText", modality="scRNA",
   primary="N", ln="Y", dist="N", normal="N", blood="N",
   n_pat="6", n_samp="11", n_cells="87191",
   treat="mixed", longit="Y (paired before-treatment and on-treatment)", mstage="unstated",
   nk_depth="cluster_only",
   nk_claim="rising trend of T/NK cell ratio in the majority of on-treatment samples",
   ebv="categorical", acc="zenodo.13736891", acc_blood="N",
   ebv_kind="plasma_EBV_DNA", ebv_snip='one of the patients with negative plasma EBV DNA at the baseline (P5) showed decreased CD8T_IFN cells on the treatment with induction therapy, which was different from other plasma EBV DNA positive samples',
   snip="11 lymph node samples from six patients were fresh-processed for scRNA-seq and another 18 lymph node samples from nine patients for bulk RNA-seq",
   note="material is metastatic lymph node biopsies from high-risk metastatic locally advanced NPC; abstract's 11 scRNA-seq refers to samples not patients"),

 "39438442": dict(
   src="fullText", modality="CyTOF",
   primary="N", ln="N", dist="N", normal="N", blood="Y",
   n_pat="24", n_samp="NOT_FOUND", n_cells="NOT_FOUND",
   treat="post-ICI", longit="Y (baseline and during treatment)", mstage="M0_only",
   nk_depth="subclustered",
   nk_claim="CD56+ NK cells resolved into phenotypic subpopulations by CyTOF",
   ebv="none", acc="NOT_FOUND", acc_blood="N",
   snip="A dynamic single-cell atlas was profiled using mass cytometry on peripheral blood mononuclear cell samples from 12 pairs of matched relapsing and non-relapsing patients in the aPD1-CRT arm",
   note="12 matched pairs = 24 patients profiled by CyTOF; locoregionally advanced NPC (CONTINUUM phase 3). GEO/EGA accessions in the paper are external validation datasets, not this study's CyTOF deposit"),

 # ============================== 2025 =========================================
 "40059116": dict(
   src="fullText", modality="scRNA",
   primary="Y", ln="N", dist="N", normal="N", blood="N",
   n_pat="10", n_samp="10", n_cells="NOT_FOUND",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="absent", nk_claim="NOT_FOUND", ebv="categorical",
   acc="HRA010229", acc_blood="N",
   ebv_kind="tumour_EBV_status", ebv_snip='NPC scRNA-seq data was retrieved from a published study, which classified NPC cells into EBV-high and EBV-low',
   snip="Clinical sample collection Ten NPC samples for single-cell RNA sequencing (scRNA-seq) wer[e collected]",
   note="the retrieved text also says 'NPC scRNA-seq data was retrieved from a published study' for the EBV-high/low classification; both statements are recorded"),

 "40265092": dict(
   src="fullText", modality="scRNA",
   primary="N", ln="Y", dist="Y (site: cervical lymph node metastasis)", normal="N", blood="N",
   n_pat="1", n_samp="1", n_cells="NOT_FOUND",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="cluster_only", nk_claim="NOT_FOUND", ebv="none",
   acc="HRA005638", acc_blood="N",
   snip="a tumor metastatic lymph node (TM_LN) of nasopharyngeal carcinoma",
   note="primarily an HIV-associated lymphoma study; exactly one NPC metastatic lymph node sample among 17 lymph node tissues profiled. NPC-specific cell count not stated"),

 "40315843": dict(
   src="fullText", modality="scRNA;scTCR",
   primary="Y", ln="N", dist="N", normal="N", blood="Y",
   n_pat="22", n_samp="NOT_FOUND", n_cells="NOT_FOUND",
   treat="mixed", longit="N", mstage="M0_only",
   nk_depth="cluster_only",
   nk_claim="NK cell marker genes profiled across CD56- and CD56+ double-negative TIL clusters",
   ebv="none", acc="HRA008590", acc_blood="UNKNOWN",
   snip="The generated bulk RNA-seq, scRNA-seq, and scTCR-seq data in this study have been deposited to Genome Sequence Archive (GSA) ... HRA008590",
   note="47 TIL infusion products and 62 pretreatment TMEs analysed overall; 22 patients received TIL plus CCRT. PBMCs isolated from peripheral blood for analyses"),

 "40389670": dict(
   src="fullText", modality="scRNA",
   primary="Y", ln="N", dist="N", normal="Y", blood="N",
   n_pat="NOT_FOUND", n_samp="NOT_FOUND", n_cells="359566",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="cluster_only",
   nk_claim="CD8+ effector T cells and NK cells were markedly reduced in tumors",
   ebv="none", acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="Single-cell RNA library preparation used the 10 x Genomics Chromium platform and Single Cell 3' v3 kit, strictly following the manufacturer's protocol.",
   note="GSE102349/GSE53819/GSE68799 named in the paper are bulk validation cohorts, not the single-cell source"),

 "40392335": dict(
   src="fullText", modality="scRNA",
   primary="Y", ln="N", dist="N", normal="Y", blood="N",
   n_pat="NOT_FOUND", n_samp="NOT_FOUND", n_cells="NOT_FOUND",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="absent", nk_claim="NOT_FOUND", ebv="none",
   acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="The 10 Chromium Single Cell platform was employed to barcode-label cell suspensions, enabling a loading range of 300-500,000 cells.",
   note="data availability states the generated datasets are not publicly available; TCGA-NPC used for bulk analyses"),

 "40530955": dict(
   src="fullText", modality="scRNA",
   primary="Y", ln="N", dist="N", normal="N", blood="N",
   n_pat="4", n_samp="4", n_cells="28957",
   treat="post-CRT", longit="N", mstage="unstated",
   nk_depth="cluster_only",
   nk_claim="activated NK cells elevated in radiosensitive versus radioresistant samples",
   ebv="none", acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="Following single-cell suspension preparation, library construction, sequencing analysis, and quality control of the obtained data, we acquired transcriptome data for a total of 28,957 cells.",
   note="supplementary file titles indicate 34 NPC bulk samples and 4 NPC single-cell samples; GSE32389 is an external bulk validation set"),

 "40543508": dict(
   src="fullText", modality="spatial;scRNA",
   primary="Y", ln="N", dist="N", normal="Y", blood="N",
   n_pat="135", n_samp="NOT_FOUND", n_cells="105784",
   treat="post-ICI", longit="N", mstage="unstated",
   nk_depth="absent",
   nk_claim="NOT_FOUND",
   ebv="none", acc="NOT_FOUND", acc_blood="N",
   snip="Raw sequencing data, including GeoMx DSP data, HEV RNA-seq data; and RNA-seq data of HUVECs after cytokine induction, have been deposited in the Genome Sequence Archive",
   note="new data are GeoMx DSP spatial (region-level, not single-cell resolution) plus bulk; the 105,784 vascular endothelial cells come from an integrated pan-cancer single-cell atlas. Cohort 1 = 135 patients"),

 "40691404": dict(
   src="abstract_only", modality="scRNA;spatial",
   primary="Y", ln="N", dist="N", normal="N", blood="N",
   n_pat="24", n_samp="39", n_cells="NOT_FOUND",
   treat="post-CRT", longit="N", mstage="unstated",
   nk_depth="NOT_FOUND", nk_claim="NOT_FOUND", ebv="none",
   acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="we conducted single-cell and spatial transcriptomics analysis of 39 tumors from 24 patients to reveal the microenvironmental differences between primary and rNPC",
   note="primary versus recurrent NPC; full text not retrievable from Europe PMC"),

 "40746726": dict(
   src="fullText", modality="scRNA",
   primary="Y", ln="N", dist="N", normal="N", blood="N",
   n_pat="15", n_samp="15", n_cells="NOT_FOUND",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="cluster_only",
   nk_claim="key regulatory genes predominantly localized within the T & NK cell population",
   ebv="none", acc="CNP0001341;CNP0001503", acc_blood="N",
   snip="A total of 15 primary NPC tumor samples were collected for single-cell RNA sequencing (scRNA-seq).",
   note="each well designed to capture 8,000-14,000 cells; total cell count not stated"),

 "40818459": dict(
   src="fullText", modality="spatial;scRNA",
   primary="Y", ln="N", dist="N", normal="Y", blood="N",
   n_pat="NOT_FOUND", n_samp="NOT_FOUND", n_cells="355955",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="absent", nk_claim="NOT_FOUND", ebv="none",
   acc="GSE150825;GSE150430;GSE162025", acc_blood="Y",
   snip="we generated an extensive scRNA-seq profile by integrating datasets from three independent studies. We processed the dataset and clustered a total of 355,955 single cells into major cell lineages",
   note="the single-cell component is an integration of three published datasets; the study's own new data are in-house Visium spatial transcriptomics (spot-level). GSE162025 contains paired blood"),

 "41067232": dict(
   src="abstract_only", modality="spatial",
   primary="Y", ln="N", dist="N", normal="N", blood="N",
   n_pat="25", n_samp="NOT_FOUND", n_cells="NOT_FOUND",
   treat="post-ICI", longit="N", mstage="unstated",
   nk_depth="NOT_FOUND", nk_claim="NOT_FOUND", ebv="none",
   acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="Integrated genomic and spatial transcriptomic analyses were performed to characterize the patient population benefitting from this combination therapy.",
   note="25 patients with unresectable recurrent NPC; spatial platform not named in the retrieved abstract, so single-cell resolution is unconfirmed"),

 "41395113": dict(
   src="fullText", modality="scRNA",
   primary="Y", ln="N", dist="N", normal="Y", blood="N",
   n_pat="NOT_FOUND", n_samp="NOT_FOUND", n_cells="NOT_FOUND",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="cluster_only", nk_claim="NOT_FOUND", ebv="none",
   acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="Researchers performed single-cell RNA sequencing on NPC and adjacent normal tissue samples.",
   note="cohort size not stated in retrieved text; data availability says data obtainable from the corresponding author"),

 "41619722": dict(
   src="fullText", modality="scRNA;CyTOF",
   primary="N", ln="N", dist="Y (site: bone/spine)", normal="N", blood="N",
   n_pat="3", n_samp="NOT_FOUND", n_cells="NOT_FOUND",
   treat="unstated", longit="N", mstage="M1_included",
   nk_depth="cluster_only",
   nk_claim="T/NK cells were more abundant in NPC, COAD, HCC, RCC bone metastases",
   ebv="none", acc="OEP005136", acc_blood="N",
   snip="The samples were collected from 52 patients diagnosed with one of the 13 cancer types, including ... nasopharyngeal carcinoma (NPC) (n = 3)",
   note="pan-cancer bone metastasis atlas; NPC contributes 3 patients. 1,304,048 cells is the whole atlas, NPC-specific count not stated"),

 "41986499": dict(
   src="abstract_only", modality="scRNA",
   primary="Y", ln="N", dist="N", normal="N", blood="N",
   n_pat="240", n_samp="NOT_FOUND", n_cells="NOT_FOUND",
   treat="mixed", longit="N", mstage="unstated",
   nk_depth="NOT_FOUND", nk_claim="NOT_FOUND", ebv="none",
   acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="we performed multi-omics profiling, including proteomics, phosphoproteomics, genomics and transcriptomics, on 240 patients with NPC who were receiving GP-IC or concurrent chemoradiotherapy (CCRT) alone",
   note="240 is the multi-omics cohort; how many of those had scRNA-seq is not stated in the abstract"),

 "41992060": dict(
   src="abstract_only", modality="scRNA",
   primary="Y", ln="N", dist="N", normal="N", blood="Y",
   n_pat="NOT_FOUND", n_samp="NOT_FOUND", n_cells="NOT_FOUND",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="NK_is_theme",
   nk_claim="Higher levels of CD16+CD57+ NK cells in blood correlated with better patient outcomes",
   ebv="none", acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="This study investigates changes in NK cell subsets in the blood of NPC patients ... Dimensionality reduction and clustering analyses were conducted on paired single-cell RNA sequencing data to explore differences in NK cell subsets",
   note="NK trap check: 'NK' here is unambiguously natural killer cells (CD16+CD57+ subsets, cytotoxicity assays), not non-keratinizing"),

 # ============================== 2026 =========================================
 "42230534": dict(
   src="abstract_only", modality="scRNA",
   primary="Y", ln="N", dist="Y (site: not specified - metastatic vs non-metastatic NPC compared)", normal="N", blood="N",
   n_pat="NOT_FOUND", n_samp="NOT_FOUND", n_cells="NOT_FOUND",
   treat="unstated", longit="N", mstage="M1_included",
   nk_depth="NOT_FOUND", nk_claim="NOT_FOUND", ebv="none",
   acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="single-cell RNA sequencing combined with multiplex immunofluorescence revealed that B cells were highly enriched in non-metastatic NPC and formed complete tertiary lymphoid structures (TLS), whereas TLS was absent in metastatic NPC",
   note="metastatic site not named in the abstract; full text not retrievable"),

 "42245668": dict(
   src="fullText", modality="scRNA;spatial",
   primary="Y", ln="N", dist="N", normal="Y", blood="N",
   n_pat="20", n_samp="20", n_cells="89875",
   treat="naive", longit="N", mstage="unstated",
   nk_depth="subclustered",
   nk_claim="T and NK cell subsets resolved with distinct metabolic characteristics",
   ebv="categorical", acc="HRA003609;OMIX013224", acc_blood="N",
   ebv_kind="tumour_EBV_status", ebv_snip='EBER in situ hybridization (ISH) was performed on formalin-fixed, paraffin-embedded (FFPE) tissue specimens',
   snip="We collected 20 untreated samples (15 nasopharyngeal carcinoma (NPC) and 5 [chronic inflammation]) ... Of these, 9 samples underwent single-cell transcriptome sequencing, 8 samples were processed for 10x spatial transcriptomics",
   note="9 of the 20 samples went to scRNA-seq; 10x Visium spatial is spot-level, not single-cell resolution"),

 "42456326": dict(
   src="abstract_only", modality="scRNA",
   primary="Y", ln="N", dist="Y (site: liver)", normal="N", blood="N",
   n_pat="NOT_FOUND", n_samp="NOT_FOUND", n_cells="NOT_FOUND",
   treat="unstated", longit="N", mstage="M1_included",
   nk_depth="NOT_FOUND", nk_claim="NOT_FOUND", ebv="none",
   acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="the construction of a comprehensive single-cell transcriptomic atlas of primary NPC and NPCLM [NPC liver metastasis]",
   note="cohort size and data provenance not stated in the retrieved abstract; full text not retrievable"),

 "42625209": dict(
   src="abstract_only", modality="scRNA",
   primary="Y", ln="Y", dist="N", normal="Y", blood="N",
   n_pat="NOT_FOUND", n_samp="NOT_FOUND", n_cells="27330",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="NOT_FOUND", nk_claim="NOT_FOUND", ebv="none",
   acc="GSE150825", acc_blood="N",
   snip="Fresh normal nasopharyngeal tissue, primary NPC, and nodal metastases were profiled by scRNA-seq, yielding 27,330 cells.",
   note="GSE150825 is the public validation cohort re-analysed, not this study's own deposit; own accession not stated in the abstract"),

 "PPR1118410": dict(
   src="abstract_only", modality="scRNA",
   primary="N", ln="N", dist="N", normal="N", blood="Y",
   n_pat="NOT_FOUND", n_samp="NOT_FOUND", n_cells="NOT_FOUND",
   treat="post-CRT", longit="Y (before and after 18 Gy RT; leukocyte counts at week 1-2, 3-4, 5-6)", mstage="M0_only",
   nk_depth="NOT_FOUND", nk_claim="NOT_FOUND", ebv="none",
   acc="NOT_FOUND", acc_blood="UNKNOWN",
   snip="Single-cell RNA sequencing (scRNA-seq) was performed on peripheral blood mononuclear cell (PBMC) samples from stage I NPC patients (T1N0M0) treated with RT alone. Blood samples were collected before and after 18 Gy RT",
   note="preprint (bioRxiv/medRxiv), no PMID; T1N0M0 = M0 only. Number of patients with scRNA-seq not stated"),

 "37607536": dict(
   src="abstract+external_repository", modality="scRNA",
   primary="NOT_FOUND", ln="N", dist="N", normal="N", blood="Y",
   n_pat="10", n_samp="10", n_cells="9512",
   treat="unstated", longit="N", mstage="unstated",
   nk_depth="NK_is_theme",
   nk_claim="tumor-associated NK cells enriched in tumors show impaired anti-tumor functions",
   ebv="none", acc="zenodo.8275845 (via GSE162025)", acc_blood="Y",
   snip="we perform integrative single-cell RNA sequencing analyses on NK cells from 716 patients with cancer, covering 24 cancer types",
   note="NPC counts are from the deposited object comb_CD56_CD16_NK_blood.h5ad (Zenodo record 8275845): meta_histology 'Nasopharyngeal Carcinoma(NPC)' = 9,512 blood NK cells from 10 patients / 10 samples, all sourced from dataset GSE162025. Cell counts here are circulating NK cells only, not the whole NPC sample"),
}
