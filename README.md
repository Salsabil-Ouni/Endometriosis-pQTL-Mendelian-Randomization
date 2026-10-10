# Endometriosis pQTL Mendelian Randomization (subtype-resolved, two-platform)

## Question
Which circulating proteins causally influence endometriosis risk, do the causal proteins differ by
subtype (ovarian, deep, peritoneal), and in which endometrial/lesion cell types do they act?

## Why this design (novelty check done 2026-09-30)
Plain plasma-proteome MR on endometriosis is already published many times (UKB-PPP vs FinnGen;
recurring hits RSPO3, FSHB). Our added value:
1. Two independent pQTL platforms (UKB-PPP Olink + deCODE SomaScan) for replication.
2. Subtype-specific MR (FinnGen ovary / deep / pelvic peritoneum endpoints).
3. Cell-type mapping of causal proteins using a public endometrial single-cell atlas.
4. Druggability annotation (Open Targets, DGIdb).
Not yet verified: that no paper already combines 1-3. Re-check before writing up.

## Data sources (all public; verify terms at download)
- Endometriosis GWAS (overall): GWAS Catalog GCST90205183 (Rahmioglu 2023, 23andMe excluded, GRCh37)
- Subtype GWAS: FinnGen public results (endometriosis of ovary, deep, pelvic peritoneum)
- pQTL: UKB-PPP (Sun 2023; AWS open data / Synapse syn51364943), deCODE (Ferkingstad 2021; decode.com/summarydata)
- Tissue eQTL (optional): GTEx uterus
- Single cell: Human Endometrial Cell Atlas

## Pipeline
1. Instruments: cis-pQTLs (within 1 Mb of gene), p < 5e-8, LD clumped r2 < 0.001 (1000G EUR), F > 10.
2. cis-MR per protein (Wald ratio / IVW), FDR correction across proteins.
3. Sensitivity: weighted median, MR-Egger, Steiger directionality, HEIDI/colocalization (PP.H4 > 0.8).
4. Replication across platforms and across GWAS (Rahmioglu vs FinnGen).
5. Subtype MR (exploratory for deep, lower power).
6. Single-cell expression mapping of surviving proteins.
7. Druggability annotation.

## Known risks
- Genome build: Rahmioglu is GRCh37; check each pQTL set and FinnGen build before harmonizing.
- Sample overlap: UKB-PPP participants vs FinnGen is minimal; UKB vs Rahmioglu may overlap (check cohorts).
- Weak instruments, pleiotropy, Olink vs SomaScan epitope effects.
- Deep endometriosis has ~3,800 cases: report as exploratory.

## Environment
R 4.6.1 (TwoSampleMR, coloc, MendelianRandomization) + Python 3.11 (scanpy for single cell).
Packages not yet installed.

## Status (updated 2026-09-30)
- [x] Folder set up
- [x] Novelty re-check: cross-platform replication (deCODE + INTERVAL) and MR + scRNA-seq endometriosis papers already exist. Subtype-specific pQTL MR with cell-type mapping of causal proteins not found, but the space is crowded. Treat subtype resolution as the main differentiator and re-check before submission.
- [x] R packages installed (TwoSampleMR, coloc, MendelianRandomization, data.table etc.) in the user R library
- [~] GWAS download running (script 01): Rahmioglu overall + FinnGen R12 overall, ovary, deep, pelvic peritoneum. FinnGen R12 manifest saved in data/gwas.
- [ ] pQTL data: BLOCKED on access.
  - UKB-PPP: Synapse project syn51364943 needs a free Synapse login and terms acceptance.
  - deCODE: download site is a JavaScript form app, not scriptable with plain HTTP; files must be fetched through a browser (deCODE proteomics folder) after accepting the form.
- [ ] Python single-cell packages: installed into the project venv (.venv), NOT the shared Python 3.11. A global scanpy install upgraded numpy and broke pandas there on 2026-10-01; versions restored (numpy 1.26.4, pandas 2.0.3, scipy 1.10.1). Always use the venv.

## Results so far (2026-10-01)
- Design change: full UKB-PPP archives are about 1.6 TB, so the screen uses the lead cis-pQTL list from Sun 2023 ST15 (1,928 non-MHC cis leads, 1,927 proteins, min F 45.5) with Wald-ratio MR. Full archives are downloaded only for hits, for colocalization.
- Screen across 5 outcomes (Rahmioglu overall, FinnGen overall, ovary, deep, pelvic peritoneum). FDR<0.05 hits: 5, 11, 12, 1, 4 respectively; 19 distinct proteins in the first three outcomes (hit list in results).
- Sanity check passed: published hits FSHB, RSPO3, WASHC3, EPHB4, HSPG2, ITGA11, TFPI reappear.
- Cross-GWAS replication (FinnGen and Rahmioglu both p<0.05, same direction): FSHB, WASHC3, RSPO3, TFPI, ITGA11, ABO, PCDH7, FN1.
- Subtype heterogeneity (Cochran Q across ovary/deep/peritoneum): weak. Only FSHB (p=0.049) and SHBG (p=0.012) nominal, none survive correction for 19 tests. Control sets overlap, so exploratory. This weakens the subtype differentiator; do not oversell it.
- Coloc (Rahmioglu overall, +/-500 kb): pipeline works. pGWAS GENPOS is GRCh38, the GRCh37 position is inside the variant ID. Early values: ABO PP.H4 0.61, ARSB 0.005. Remaining proteins pending download.
- Steiger directionality (overall Rahmioglu and FinnGen, prevalence 5% assumed): correct direction for FSHB, RSPO3, WASHC3, EPHB4, ERBB2, HSPG2, ITGA11, OGA, PCDH7, PPP1R14D, SPON2, CRIM1. Not supported (direction false, non-significant): ABO, ARSB, CXCL5, FN1, SHBG, TFPI, TIMP3. Treat the latter as weak or possibly pleiotropic until coloc says otherwise. Exposure sample size uses 33,985 (UKB-PPP European discovery N).
- Druggability (Open Targets, coarse): clinical candidates exist for ERBB2 (48), EPHB4 (4), FN1 (5), TFPI (4), RSPO3 (1); most other hits have tractability flags only. FSHB shows 0 because drugs act on the receptor, not the beta subunit; interpret with care.
- Coloc vs Rahmioglu overall (all 19 hit proteins, +/-500 kb, default priors, single-causal-variant assumption): PP.H4 > 0.8 for FN1 (0.98) and TFPI (0.89). Suggestive: RSPO3 (0.76), ABO (0.61), ITGA11 (0.60). Strong evidence of DISTINCT causal variants (PP.H3 > 0.9) for FSHB (0.94), WASHC3 (0.98) and HSPG2 (1.00), despite strong MR p-values; likely multi-signal loci where the single-variant assumption fails, so these need conditional/SuSiE-style coloc before any conclusion. ERBB2, EPHB4, SPON2, CRIM1, OGA, PPP1R14D, SHBG, ARSB, CXCL5 show no colocalization.
- Tension to resolve: FN1, TFPI and ABO fail Steiger but FN1 and TFPI colocalize. Do not report a protein as causal on one criterion alone. Current tiers: (A) MR FDR + replication + coloc: FN1, TFPI; (B) MR FDR + replication + Steiger OK, coloc suggestive or multi-signal: RSPO3, ITGA11, FSHB, WASHC3; (C) MR only: the rest.
- Single cell (HECA subset, 90,001 cells, 36 cell types; raw counts normalised to CP10k + log1p; an earlier un-normalised run was discarded as a depth artefact). 16 of 19 hit genes present (FSHB, OGA, PPP1R14D absent from the atlas; FSHB is pituitary-expressed). Most specific cell type by gene: stromal/mesenchymal for RSPO3, SPON2, TIMP3, WASHC3, FN1 (early decidualising stroma), ITGA11 (basal fibroblasts), CRIM1 (perivascular); endothelial/lymphatic for TFPI, EPHB4, HSPG2; epithelial for CXCL5, ARSB, ABO, PCDH7, ERBB2. Caution: the "dHormones" group is hormone-treated stroma, so its top rankings partly reflect treatment, not disease biology.
- Endometriosis vs control in eutopic endometrium (sample-level pseudobulk, 288 tests): essentially null. One FDR<0.05 result (ITGA11 in lymphoid cells) at near-zero expression, treat as noise. Atlas has no lesion tissue, so this does not address lesion biology.
- Multi-signal coloc (SuSiE + coloc.susie, 1000G EUR LD, Rahmioglu overall, +/-500 kb): FSHB has 3 independent pQTL signals and 1 GWAS signal, max PP.H4 = 0.034, so no pair of signals colocalizes; the earlier distinct-variant verdict stands even after allowing for multiple signals. WASHC3 has 1 pQTL signal but SuSiE finds no GWAS credible set in the region, so coloc cannot be run (GWAS signal at this locus too weak to resolve); WASHC3 stays MR-only. HSPG2 pending (needs chromosome 1 LD reference, still downloading). Updated tiers: FSHB and WASHC3 move out of tier B (MR signal not supported by colocalization; this does not refute FSH biology, since plasma FSHB protein here is a pituitary product and the instrument may tag a different variant).
- Multi-signal coloc, all 8 loci run (2026-10-06; results in part1 and part2 tables, signal pairs in matching files). FN1: PP.H4 = 0.96 (the single GWAS signal colocalizes with one of 9 pQTL signals), confirming tier A. HSPG2: GWAS signal present but PP.H4 = 6e-5 vs every pQTL signal, so not supported. FSHB: not supported (above). For TFPI, ITGA11, ABO, RSPO3 and WASHC3 SuSiE finds NO GWAS credible set in the region, so signal-level coloc cannot be run; this is not evidence against them (the earlier single-variant coloc gave TFPI 0.89, RSPO3 0.76, ABO 0.61, ITGA11 0.60) but they stay unconfirmed. Caution: SuSiE returned 7 to 10 pQTL credible sets for most proteins, which is probably over-splitting from the mismatch between 1000G LD and the UKB-PPP LD; treat the pQTL-side signal count as unreliable and do not report it as biology. A fairer check would use an in-sample LD reference.
- Fixed along the way: two of the 1000G files (chromosomes 1 and 6) were corrupted by stop-and-resume downloads (file size matched but gzip CRC failed), so always gzip -t after downloading; the genotype extractor now reads only the position field before the region.
- FINAL INTEGRATED RESULTS (2026-10-10), supersede earlier tier statements. FinnGen coloc added (script 14; both pGWAS and FinnGen R12 are GRCh38). CORRECTION: the earlier "FSHB not supported" verdict came from Rahmioglu data only. Diagnostics (script 15): the Rahmioglu public file has about 2,600 variants in the FSHB 1 Mb window vs about 6,700 in FinnGen, so a missing causal variant biases coloc toward H3. In FinnGen the top disease variant is 0.2 kb from the top pQTL variant and PP.H4 = 0.999 in all four outcomes. Rahmioglu-based coloc must therefore be treated as lower-resolution, not as contradicting evidence.
- Tiering rule (script 16; thresholds were fixed after seeing results, so describe as conventional-but-post-hoc): Tier 1 = MR replicated across the two overall GWAS + Steiger OK in both + coloc PP.H4 >= 0.7 in both GWAS; Tier 2 = replicated + PP.H4 >= 0.7 in at least one GWAS; Tier 3 = rest. Result: Tier 1 = RSPO3 only (Rahm 0.755, FinnGen 0.778). Tier 2 = FSHB, WASHC3, TFPI, ITGA11, FN1 (TFPI and FN1 fail Steiger; FN1 is 0.98 in Rahmioglu but 0.31 in FinnGen, so GWAS-dependent). Tier 3 = the other 13. WARNING: at the more conventional PP.H4 >= 0.8 threshold nothing is Tier 1 (RSPO3 falls below). Do not call any protein causal; call them prioritised candidates.
- Subtype colocalization (FinnGen): FSHB colocalizes in ovary, deep and overall alike (about 0.998), ITGA11 in ovary and deep (0.98, 0.93), RSPO3 and WASHC3 only in ovary (0.85); this is the only subtype signal that survives, and deep has few cases, so exploratory.
- Full table: results/integrated_evidence_table.tsv.
- deCODE (SomaScan) replication, first file only (FSHB, aptamer 8036_75; file verified against the published md5 of the uncompressed data; script 17; an earlier run was discarded because a variable-shadowing bug analysed the wrong region, caught by checking the lead position). RESULT: the UKB-PPP FSHB cis-pQTL does NOT replicate on SomaScan. Best deCODE variant in the 1 Mb region has p = 9e-4 (no genome-wide significant cis-pQTL), the UKB-PPP lead variant has p = 0.73 in deCODE, and coloc with all four FinnGen outcomes is H4 < 0.01 (H0 dominant, i.e. no protein signal to colocalize). So FSHB's MR evidence rests on the Olink assay alone; add a flag "not replicated on second proteomic platform" to FSHB in the paper, keep it as a hypothesis, and do not describe it as cross-platform supported. Possible reasons (cannot be resolved with these data): different assay targets (FSHB is one subunit of a heterodimer, low plasma abundance), epitope effects, or Olink-specific artefact. The same variant is strongly associated with endometriosis (FinnGen p 4e-40), so the disease signal itself is not in question.
- 2026-10-10 progress: figures made (scripts 19; figure 1 forest plot, figure 2 evidence matrix, both with FSHB flagged as not replicated on deCODE); manuscript draft v0.1 written (docs folder; methods and results drafted, discussion outlined, numbers marked to be verified); sensitivity analyses (script 18: coloc vs window 250 kb/500 kb/1 Mb and prior p12 1e-5/1e-6; Steiger vs prevalence 2/5/10%) were still running when this was written. All deep and peritoneum FDR hits were among the 19 carried-forward proteins, so none were omitted.
- SENSITIVITY RESULTS (script 18; results/sensitivity_*): Steiger verdicts identical at prevalence 2%, 5%, 10% in both GWAS (robust). Coloc insensitive to window (250 kb, 500 kb, 1 Mb) but STRONGLY sensitive to the p12 prior: with p12 = 1e-6 (instead of 1e-5) RSPO3 falls from 0.76/0.78 to 0.24/0.26, TFPI 0.89/0.87 to 0.44/0.39, WASHC3 0.71 to 0.20 (FinnGen), ITGA11 0.60/0.94 to 0.13/0.63. Only FSHB (FinnGen, 0.99) and FN1 (Rahmioglu, 0.82) stay above 0.8 and both fail another check. So RSPO3 as Tier 1 is TENTATIVE (prior-dependent); no protein has robust joint support across both GWAS and both priors. Write the paper as prioritisation of candidates with this stated up front, not as a list of causal proteins.
- Still needed before any preprint: remaining deCODE files, apply deCODE excluded-variants list, final verification of every number in the draft against the result tables, independent review and from-scratch rerun by someone experienced, co-author/supervisor agreement, data-terms check, STROBE-MR checklist.
- deCODE files still to obtain (user approved all 7): RSPO3 (13094_75, 8427_118), TFPI (3336_50), FN1 (3434_34, 3435_53, 4131_72). ITGA11 and WASHC3 are not on the deCODE aptamer list, so they cannot be replicated on that platform. Files downloaded so far are kept in D:\data D\Downloads D\ and copied into data/pqtl_decode. deCODE terms: scientific use only, no redistribution (never include these files in a public deposit).
- LD approximation to state in the methods: 1000 Genomes phase 3 European LD (n=503) applied to UKB-PPP and Rahmioglu summary data.
- Not yet done: FinnGen-based coloc (needs rsID mapping), deCODE replication (needs terms form accepted in browser), Steiger, single-cell mapping, druggability.
