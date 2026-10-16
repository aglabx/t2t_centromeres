# t2t_centromeres

Centromeric satellite arrays, one locus per chromosome, for **121 vertebrate assemblies of 68 species** (and
three interspecific fish hybrids): the 32 T2T assemblies of 14 mammals of v1.2 (human CHM13, chimpanzee, bonobo,
gorilla, Sumatran and Bornean orangutan, siamang, crab-eating and rhesus macaque, pig, sheep, goat, horse,
giraffe, with their alternate and parental haplotype assemblies) and 89 more chosen from NCBI by their sequence
statistics, not their names: 63 ray-finned fishes, 21 mammals (among them mouse, cat, cattle, donkey, kiang,
water deer and Yarkand hare), duck, red-legged partridge, sand martin, Chinese softshell turtle and Gila monster.
Every assembly carries a completeness tier measured from its sequence ([below](#which-assemblies-are-t2t)).

The calls are made from [satellome](https://github.com/aglabx/satellome) (FasTAN engine) tandem-repeat
catalogues. Where the assembly's authors published a centromere annotation, it sits in the same GFF and
every chromosome says **where and how the two intersect**. Everything here is rebuilt from public
inputs by one script (`workflow/reproduce.sh`) and checked by md5 (`workflow/verify.sh`).

## What is in `annotations/<Organism>__<assembly>__<accession>/`

One directory per assembly, named so that it sorts by species, e.g.
`Homo_sapiens__T2T-CHM13v2.0__GCA_009914755.4`, `Pan_paniscus__NHGRI_mPanPan1-v2.0_mat__GCA_028858845.2`
(`workflow/asm_dir.sh ACC` prints the name).

| file | content |
|---|---|
| `karyotype.svg` | open in a browser. **a**: every chromosome to scale with its satellite arrays coloured by superfamily, our locus boxed, each published interval bracketed beside it, one dot per published set (overlap / no overlap / no call). **b**: the centromeric region of each chromosome close up, our locus and the published intervals side by side, one bracket per interval. Hover for coordinates |
| `centromeres.gff3` | our locus (`centromere`, source `satellome_cen`); its arrays (`satellite_DNA`, `largest=true` on the longest one, a size fact only); every published interval (source = publisher); every overlap between the two (`region`, source `intersection`). seqid = GenBank accession |
| `centromeres.tsv` | one row per chromosome: locus, family, bp, largest array, monomer / HOR period, monomer consensus, co-located families |
| `centromeric_arrays.tsv` / `.bed` | every array of every locus with FasTAN period, monomer and HOR period, monomer and HOR consensus, and how it joined the locus |
| `official_overlap.tsv` | chromosome × published set: relation, overlap bp, fractions, gap |
| `satellite_arrays.tsv` | every satellite array ≥ 10 kb on a chromosome with its superfamily and monomer / HOR period (what the karyotype draws) |
| `family_breadth.tsv` | every satellite superfamily: chromosomes carrying it (≥ 10 kb; and ≥ 100 kb), bp, base monomer period, chromosomes with it at both ends |
| `outliers.tsv` | chromosomes with no array of the centromeric family; arrays without a usable consensus |
| `exceptions.tsv` | near ties, broad short-consensus (not eligible) superfamilies, co-families, cap mode, skipped subterminal caps, terminal (acro/telocentric) loci, a second comparable locus |
| `MANIFEST.md5` | md5 of the above |

Panel tables: `annotations/SUMMARY.tsv`, `TIERS.tsv`, `ALL_centromeres.tsv`, `ALL_official_overlap.tsv`, `STATUS.tsv`.
Every GFF attribute and relation is defined in [`docs/OUTPUT_FORMAT.md`](docs/OUTPUT_FORMAT.md).

**The locus is the centromeric satellite zone, not the centromere itself.** The functional centromere (the
active HOR in primates, the CENP-A domain) is usually a small part of that zone, and the size of an array
does not identify it: on the 15 ape chromosomes where CenSat placed the active array from CENP-A and CDR,
our locus overlaps it on 14, but the largest array of our locus covers at least half of it on only 6. No
array is put forward as the active one; the comparison with every published set is made at the level of
the locus.

![CHM13 karyotype: (a) alpha satellite arrays in blue, our locus boxed, the T2T Consortium live HOR bracketed in magenta; (b) each centromeric region close up](annotations/Homo_sapiens__T2T-CHM13v2.0__GCA_009914755.4/karyotype.svg)

## Which assemblies are T2T

An assembly is called T2T here by what its sequence shows, never by its name or by the NCBI assembly level (which
the submitter declares; "Complete Genome" also labels mitochondria, endosymbiont genomes and test uploads).

| tier | rule |
|---|---|
| real T2T | no gap (run of N), no unplaced or unlocalized nuclear sequence, the telomeric repeat (TTAGGG) at both ends of every chromosome |
| gap-free, telomeres incomplete | no gap, no unplaced sequence, the telomeric repeat missing or below density at one or more chromosome ends |
| gaps | at least one gap or one unplaced / unlocalized sequence |

Gaps come from satellome's `gaps.bed`, unplaced sequences from the NCBI assembly report, telomeres from satellome's
telomere check (TTAGGG / CCCTAA in the terminal 10 kb, density ≥ 0.15); `src/t2t_centromeres/tier.py` writes
`annotations/TIERS.tsv`, and `SUMMARY.tsv` carries the tier. Chromosomes are the assembled molecules of type
Chromosome or Linkage Group (some fish reports use the latter).

The 89 assemblies added in v1.3 were selected from all 4,957 GenBank vertebrate assemblies at chromosome or
complete-genome level (NCBI Datasets, 2026-10-16) whose statistics show no gap (contigs = scaffolds) and no
unplaced sequence, with more than one chromosome and no NCBI "atypical" flag: 99 assemblies, 10 of them already
in v1.2. Names do not predict this: 21 of the 54 vertebrate assemblies with "T2T" in their name have gaps, and most
of the 99 have no "T2T" in theirs. The sequence confirmed every one of the 89 as gap-free, but only 13 also carry the
telomeric repeat at both ends of every chromosome, and in 11 not one chromosome has it at both ends (among them both
Najdi sheep assemblies named T2T). Of the 32 v1.2 assemblies, 7 are real T2T (CHM13, both macaques,
goat 1.0 and 2.0, giraffe primary and paternal), 3 are gap-free with incomplete telomeres (pig, sheep 1.0, horse),
and 22 have gaps (all NHGRI ape assemblies, the goat and sheep haplotype assemblies, giraffe alternate and
maternal); they stay in the release with their tier.

## Results

<!-- BEGIN TABLES (generated by src/t2t_centromeres/readme_tables.py; do not edit) -->

| completeness tier | assemblies |
|---|---:|
| real T2T | 20 |
| gap-free, telomeres incomplete | 79 |
| gaps | 22 |

| organism | assembly | accession | completeness | gaps | both telomeres | chromosomes with a call | centromeric family (monomer bp) | karyotype |
|---|---|---|---|---:|---:|---:|---|---|
| *Acanthopagrus latus* | [ASM5614952v1](annotations/Acanthopagrus_latus__ASM5614952v1__GCA_056149525.1) | GCA_056149525.1 | gap-free, telomeres incomplete | 0 | 0 / 24 | 24 / 24 | 186 | [svg](annotations/Acanthopagrus_latus__ASM5614952v1__GCA_056149525.1/karyotype.svg) |
| *Acanthopagrus latus* | [ASM5813089v1](annotations/Acanthopagrus_latus__ASM5813089v1__GCA_058130895.1) | GCA_058130895.1 | real T2T | 0 | 24 / 24 | 24 / 24 | 186 | [svg](annotations/Acanthopagrus_latus__ASM5813089v1__GCA_058130895.1/karyotype.svg) |
| *Acanthopagrus schlegelii* | [ASM4175387v1](annotations/Acanthopagrus_schlegelii__ASM4175387v1__GCA_041753875.1) | GCA_041753875.1 | gap-free, telomeres incomplete | 0 | 19 / 24 | 24 / 24 | 186 | [svg](annotations/Acanthopagrus_schlegelii__ASM4175387v1__GCA_041753875.1/karyotype.svg) |
| *Alectoris rufa* | [Arufa_2401F](annotations/Alectoris_rufa__Arufa_2401F__GCA_051903495.1) | GCA_051903495.1 | gap-free, telomeres incomplete | 0 | 0 / 40 | 20 / 40 | 996 | [svg](annotations/Alectoris_rufa__Arufa_2401F__GCA_051903495.1/karyotype.svg) |
| *Anas platyrhynchos* | [IASCAAS_PekinDuck_T2T](annotations/Anas_platyrhynchos__IASCAAS_PekinDuck_T2T__GCA_047663525.1) | GCA_047663525.1 | gap-free, telomeres incomplete | 0 | 40 / 41 | 41 / 41 | 190 | [svg](annotations/Anas_platyrhynchos__IASCAAS_PekinDuck_T2T__GCA_047663525.1/karyotype.svg) |
| *Arapaima gigas* | [ASM5488337v1](annotations/Arapaima_gigas__ASM5488337v1__GCA_054883375.1) | GCA_054883375.1 | gap-free, telomeres incomplete | 0 | 13 / 28 | 8 / 28 | 2825 | [svg](annotations/Arapaima_gigas__ASM5488337v1__GCA_054883375.1/karyotype.svg) |
| *Astatotilapia calliptera* | [fAstCal671.pri](annotations/Astatotilapia_calliptera__fAstCal671.pri__GCA_988222555.1) | GCA_988222555.1 | gap-free, telomeres incomplete | 0 | 9 / 22 | 22 / 22 | 237 | [svg](annotations/Astatotilapia_calliptera__fAstCal671.pri__GCA_988222555.1/karyotype.svg) |
| *Bidyanus bidyanus* | [Bidyanus_bidyanus_T2T_v1](annotations/Bidyanus_bidyanus__Bidyanus_bidyanus_T2T_v1__GCA_986341675.1) | GCA_986341675.1 | gap-free, telomeres incomplete | 0 | 21 / 24 | 22 / 24 | 199 | [svg](annotations/Bidyanus_bidyanus__Bidyanus_bidyanus_T2T_v1__GCA_986341675.1/karyotype.svg) |
| *Bos taurus* | [ASM5933614v1](annotations/Bos_taurus__ASM5933614v1__GCA_059336145.1) | GCA_059336145.1 | gap-free, telomeres incomplete | 0 | 0 / 30 | 13 / 30 | 3539 | [svg](annotations/Bos_taurus__ASM5933614v1__GCA_059336145.1/karyotype.svg) |
| *Capra hircus* | [T2T-goat1.0](annotations/Capra_hircus__T2T-goat1.0__GCA_040806595.1) | GCA_040806595.1 | real T2T | 0 | 31 / 31 | 30 / 31 | 703 + co-family | [svg](annotations/Capra_hircus__T2T-goat1.0__GCA_040806595.1/karyotype.svg) |
| *Capra hircus* | [T2T-goat1.0M](annotations/Capra_hircus__T2T-goat1.0M__GCA_041736475.1) | GCA_041736475.1 | gaps | 6 | 21 / 30 | 28 / 30 | 816 + co-family | [svg](annotations/Capra_hircus__T2T-goat1.0M__GCA_041736475.1/karyotype.svg) |
| *Capra hircus* | [T2T-goat1.0P](annotations/Capra_hircus__T2T-goat1.0P__GCA_041920685.1) | GCA_041920685.1 | gaps | 11 | 20 / 30 | 28 / 30 | 816 + co-family | [svg](annotations/Capra_hircus__T2T-goat1.0P__GCA_041920685.1/karyotype.svg) |
| *Capra hircus* | [T2T-goat2.0](annotations/Capra_hircus__T2T-goat2.0__GCA_041735815.1) | GCA_041735815.1 | real T2T | 0 | 31 / 31 | 30 / 31 | 703 + co-family | [svg](annotations/Capra_hircus__T2T-goat2.0__GCA_041735815.1/karyotype.svg) |
| *Capra hircus* | [T2T-goat2.0M](annotations/Capra_hircus__T2T-goat2.0M__GCA_040970035.1) | GCA_040970035.1 | gaps | 6 | 21 / 30 | 28 / 30 | 816 + co-family | [svg](annotations/Capra_hircus__T2T-goat2.0M__GCA_040970035.1/karyotype.svg) |
| *Capra hircus* | [T2T-goat2.0P](annotations/Capra_hircus__T2T-goat2.0P__GCA_041053135.1) | GCA_041053135.1 | gaps | 11 | 20 / 30 | 28 / 30 | 816 + co-family | [svg](annotations/Capra_hircus__T2T-goat2.0P__GCA_041053135.1/karyotype.svg) |
| *Cephalopholis sonnerati* | [ASM4338838v1](annotations/Cephalopholis_sonnerati__ASM4338838v1__GCA_043388385.1) | GCA_043388385.1 | real T2T | 0 | 24 / 24 | 24 / 24 | 38 | [svg](annotations/Cephalopholis_sonnerati__ASM4338838v1__GCA_043388385.1/karyotype.svg) |
| *Cephalopholis sonnerati* | [ASM4338842v1](annotations/Cephalopholis_sonnerati__ASM4338842v1__GCA_043388425.1) | GCA_043388425.1 | real T2T | 0 | 24 / 24 | 24 / 24 | 38 | [svg](annotations/Cephalopholis_sonnerati__ASM4338842v1__GCA_043388425.1/karyotype.svg) |
| *Channa maculata* | [ASM5991028v1](annotations/Channa_maculata__ASM5991028v1__GCA_059910285.1) | GCA_059910285.1 | gap-free, telomeres incomplete | 0 | 19 / 21 | 21 / 21 | 178 | [svg](annotations/Channa_maculata__ASM5991028v1__GCA_059910285.1/karyotype.svg) |
| *Cheilinus undulatus* | [ASM4981552v1](annotations/Cheilinus_undulatus__ASM4981552v1__GCA_049815525.1) | GCA_049815525.1 | real T2T | 0 | 24 / 24 | 23 / 24 | 204 | [svg](annotations/Cheilinus_undulatus__ASM4981552v1__GCA_049815525.1/karyotype.svg) |
| *Chromileptes altivelis (female) x Epinephelus fuscoguttatus (male)* | [ASM6075322v1](annotations/Chromileptes_altivelis_(female)_x_Epinephelus_fuscoguttatus_(male)__ASM6075322v1__GCA_060753225.1) | GCA_060753225.1 | gap-free, telomeres incomplete | 0 | 11 / 24 | 21 / 24 | 186 | [svg](annotations/Chromileptes_altivelis_(female)_x_Epinephelus_fuscoguttatus_(male)__ASM6075322v1__GCA_060753225.1/karyotype.svg) |
| *Chromileptes altivelis (female) x Epinephelus fuscoguttatus (male)* | [ASM6075323v1](annotations/Chromileptes_altivelis_(female)_x_Epinephelus_fuscoguttatus_(male)__ASM6075323v1__GCA_060753235.1) | GCA_060753235.1 | gap-free, telomeres incomplete | 0 | 13 / 24 | 24 / 24 | 38 | [svg](annotations/Chromileptes_altivelis_(female)_x_Epinephelus_fuscoguttatus_(male)__ASM6075323v1__GCA_060753235.1/karyotype.svg) |
| *Cromileptes altivelis x Epinephelus fuscoguttatus* | [YSFRI_Yushuban_mat1.0](annotations/Cromileptes_altivelis_x_Epinephelus_fuscoguttatus__YSFRI_Yushuban_mat1.0__GCA_061495365.1) | GCA_061495365.1 | gap-free, telomeres incomplete | 0 | 15 / 24 | 24 / 24 | 186 | [svg](annotations/Cromileptes_altivelis_x_Epinephelus_fuscoguttatus__YSFRI_Yushuban_mat1.0__GCA_061495365.1/karyotype.svg) |
| *Cromileptes altivelis x Epinephelus fuscoguttatus* | [YSFRI_Yushuban_pat1.0](annotations/Cromileptes_altivelis_x_Epinephelus_fuscoguttatus__YSFRI_Yushuban_pat1.0__GCA_061495325.1) | GCA_061495325.1 | gap-free, telomeres incomplete | 0 | 12 / 24 | 22 / 24 | 186 | [svg](annotations/Cromileptes_altivelis_x_Epinephelus_fuscoguttatus__YSFRI_Yushuban_pat1.0__GCA_061495325.1/karyotype.svg) |
| *Ctenopharyngodon idella* | [ASM4991979v1](annotations/Ctenopharyngodon_idella__ASM4991979v1__GCA_049919795.1) | GCA_049919795.1 | real T2T | 0 | 24 / 24 | 24 / 24 | 1511 | [svg](annotations/Ctenopharyngodon_idella__ASM4991979v1__GCA_049919795.1/karyotype.svg) |
| *Ctenopharyngodon idella* | [ASM5190354v1](annotations/Ctenopharyngodon_idella__ASM5190354v1__GCA_051903545.1) | GCA_051903545.1 | gap-free, telomeres incomplete | 0 | 23 / 24 | 24 / 24 | 1281 | [svg](annotations/Ctenopharyngodon_idella__ASM5190354v1__GCA_051903545.1/karyotype.svg) |
| *Cyprinus carpio* | [ASM5088196v1](annotations/Cyprinus_carpio__ASM5088196v1__GCA_050881965.1) | GCA_050881965.1 | gap-free, telomeres incomplete | 0 | 26 / 50 | 49 / 50 | 260 | [svg](annotations/Cyprinus_carpio__ASM5088196v1__GCA_050881965.1/karyotype.svg) |
| *Cyprinus carpio var. Quanzhounensis* | [ASM6041419v1](annotations/Cyprinus_carpio_var._Quanzhounensis__ASM6041419v1__GCA_060414195.1) | GCA_060414195.1 | gap-free, telomeres incomplete | 0 | 24 / 50 | 50 / 50 | 260 | [svg](annotations/Cyprinus_carpio_var._Quanzhounensis__ASM6041419v1__GCA_060414195.1/karyotype.svg) |
| *Danio rerio* | [ASM3317019v3](annotations/Danio_rerio__ASM3317019v3__GCA_033170195.3) | GCA_033170195.3 | gap-free, telomeres incomplete | 0 | 19 / 25 | 25 / 25 | 92 | [svg](annotations/Danio_rerio__ASM3317019v3__GCA_033170195.3/karyotype.svg) |
| *Danio rerio* | [GRCz12tu](annotations/Danio_rerio__GRCz12tu__GCA_049306965.1) | GCA_049306965.1 | gap-free, telomeres incomplete | 0 | 21 / 25 | 25 / 25 | 92 | [svg](annotations/Danio_rerio__GRCz12tu__GCA_049306965.1/karyotype.svg) |
| *Danio rerio* | [GRCz13ab](annotations/Danio_rerio__GRCz13ab__GCA_052040795.2) | GCA_052040795.2 | real T2T | 0 | 25 / 25 | 25 / 25 | 92 | [svg](annotations/Danio_rerio__GRCz13ab__GCA_052040795.2/karyotype.svg) |
| *Epinephelus cyanopodus* | [ASM5856915v1](annotations/Epinephelus_cyanopodus__ASM5856915v1__GCA_058569155.1) | GCA_058569155.1 | gap-free, telomeres incomplete | 0 | 23 / 24 | 23 / 24 | 205 | [svg](annotations/Epinephelus_cyanopodus__ASM5856915v1__GCA_058569155.1/karyotype.svg) |
| *Epinephelus fasciatus* | [ASM5964275v1](annotations/Epinephelus_fasciatus__ASM5964275v1__GCA_059642755.1) | GCA_059642755.1 | gap-free, telomeres incomplete | 0 | 15 / 24 | 24 / 24 | 186 | [svg](annotations/Epinephelus_fasciatus__ASM5964275v1__GCA_059642755.1/karyotype.svg) |
| *Epinephelus fuscoguttatus (female) x Cephalopholis sonnerati (male)* | [ASM62009623v1](annotations/Epinephelus_fuscoguttatus_(female)_x_Cephalopholis_sonnerati_(male)__ASM62009623v1__GCA_062009623.1) | GCA_062009623.1 | gap-free, telomeres incomplete | 0 | 2 / 24 | 21 / 24 | 186 | [svg](annotations/Epinephelus_fuscoguttatus_(female)_x_Cephalopholis_sonnerati_(male)__ASM62009623v1__GCA_062009623.1/karyotype.svg) |
| *Epinephelus fuscoguttatus (female) x Cephalopholis sonnerati (male)* | [ASM62009624v1](annotations/Epinephelus_fuscoguttatus_(female)_x_Cephalopholis_sonnerati_(male)__ASM62009624v1__GCA_062009624.1) | GCA_062009624.1 | gap-free, telomeres incomplete | 0 | 17 / 24 | 24 / 24 | 38 | [svg](annotations/Epinephelus_fuscoguttatus_(female)_x_Cephalopholis_sonnerati_(male)__ASM62009624v1__GCA_062009624.1/karyotype.svg) |
| *Equus asinus* | [CAU_T2T_donkey1.0](annotations/Equus_asinus__CAU_T2T_donkey1.0__GCA_054375355.1) | GCA_054375355.1 | gap-free, telomeres incomplete | 0 | 21 / 32 | 22 / 32 | 4678 + co-family | [svg](annotations/Equus_asinus__CAU_T2T_donkey1.0__GCA_054375355.1/karyotype.svg) |
| *Equus caballus* | [ASM5232381v1](annotations/Equus_caballus__ASM5232381v1__GCA_052323815.1) | GCA_052323815.1 | real T2T | 0 | 33 / 33 | 32 / 33 | 221 | [svg](annotations/Equus_caballus__ASM5232381v1__GCA_052323815.1/karyotype.svg) |
| *Equus caballus* | [T2T-horse1.0](annotations/Equus_caballus__T2T-horse1.0__GCA_052818215.1) | GCA_052818215.1 | gap-free, telomeres incomplete | 0 | 23 / 33 | 32 / 33 | 221 | [svg](annotations/Equus_caballus__T2T-horse1.0__GCA_052818215.1/karyotype.svg) |
| *Equus kiang* | [ASM5215063v1](annotations/Equus_kiang__ASM5215063v1__GCA_052150635.1) | GCA_052150635.1 | gap-free, telomeres incomplete | 0 | 26 / 27 | 18 / 27 | 117 | [svg](annotations/Equus_kiang__ASM5215063v1__GCA_052150635.1/karyotype.svg) |
| *Felis catus* | [F.catus_OcraZ1_1.0](annotations/Felis_catus__F.catus_OcraZ1_1.0__GCA_051911835.1) | GCA_051911835.1 | gap-free, telomeres incomplete | 0 | 19 / 20 | 16 / 20 | 483 | [svg](annotations/Felis_catus__F.catus_OcraZ1_1.0__GCA_051911835.1/karyotype.svg) |
| *Felis catus* | [F.catus_OcraZ1_mat1.0](annotations/Felis_catus__F.catus_OcraZ1_mat1.0__GCA_051911995.1) | GCA_051911995.1 | gap-free, telomeres incomplete | 0 | 18 / 19 | 15 / 19 | 483 | [svg](annotations/Felis_catus__F.catus_OcraZ1_mat1.0__GCA_051911995.1/karyotype.svg) |
| *Felis catus* | [F.catus_OcraZ1_pat1.0](annotations/Felis_catus__F.catus_OcraZ1_pat1.0__GCA_051911875.1) | GCA_051911875.1 | real T2T | 0 | 19 / 19 | 15 / 19 | 483 | [svg](annotations/Felis_catus__F.catus_OcraZ1_pat1.0__GCA_051911875.1/karyotype.svg) |
| *Giraffa tippelskirchi* | [T2T-mGirTip1v1.2.alt](annotations/Giraffa_tippelskirchi__T2T-mGirTip1v1.2.alt__GCA_054371575.1) | GCA_054371575.1 | gaps | 1 | 14 / 14 | 14 / 14 | 328 | [svg](annotations/Giraffa_tippelskirchi__T2T-mGirTip1v1.2.alt__GCA_054371575.1/karyotype.svg) |
| *Giraffa tippelskirchi* | [T2T-mGirTip1v1.2.mat](annotations/Giraffa_tippelskirchi__T2T-mGirTip1v1.2.mat__GCA_054371825.1) | GCA_054371825.1 | gaps | 1 | 15 / 15 | 15 / 15 | 328 | [svg](annotations/Giraffa_tippelskirchi__T2T-mGirTip1v1.2.mat__GCA_054371825.1/karyotype.svg) |
| *Giraffa tippelskirchi* | [T2T-mGirTip1v1.2.pat](annotations/Giraffa_tippelskirchi__T2T-mGirTip1v1.2.pat__GCA_054371615.1) | GCA_054371615.1 | real T2T | 0 | 15 / 15 | 14 / 15 | 329 | [svg](annotations/Giraffa_tippelskirchi__T2T-mGirTip1v1.2.pat__GCA_054371615.1/karyotype.svg) |
| *Giraffa tippelskirchi* | [T2T-mGirTip1v1.2.pri](annotations/Giraffa_tippelskirchi__T2T-mGirTip1v1.2.pri__GCA_054371585.1) | GCA_054371585.1 | real T2T | 0 | 16 / 16 | 15 / 16 | 328 | [svg](annotations/Giraffa_tippelskirchi__T2T-mGirTip1v1.2.pri__GCA_054371585.1/karyotype.svg) |
| *Gorilla gorilla gorilla* | [NHGRI_mGorGor1-v2.0_mat](annotations/Gorilla_gorilla_gorilla__NHGRI_mGorGor1-v2.0_mat__GCA_028885495.2) | GCA_028885495.2 | gaps | 5 | 23 / 24 | 24 / 24 | 171 | [svg](annotations/Gorilla_gorilla_gorilla__NHGRI_mGorGor1-v2.0_mat__GCA_028885495.2/karyotype.svg) |
| *Gorilla gorilla gorilla* | [NHGRI_mGorGor1-v2.0_pat](annotations/Gorilla_gorilla_gorilla__NHGRI_mGorGor1-v2.0_pat__GCA_028885475.2) | GCA_028885475.2 | gaps | 2 | 24 / 24 | 24 / 24 | 171 | [svg](annotations/Gorilla_gorilla_gorilla__NHGRI_mGorGor1-v2.0_pat__GCA_028885475.2/karyotype.svg) |
| *Gorilla gorilla gorilla* | [NHGRI_mGorGor1-v2.1_pri](annotations/Gorilla_gorilla_gorilla__NHGRI_mGorGor1-v2.1_pri__GCA_029281585.3) | GCA_029281585.3 | gaps | 2 | 25 / 25 | 25 / 25 | 171 | [svg](annotations/Gorilla_gorilla_gorilla__NHGRI_mGorGor1-v2.1_pri__GCA_029281585.3/karyotype.svg) |
| *Heloderma suspectum* | [16-S_hap1](annotations/Heloderma_suspectum__16-S_hap1__GCA_057233065.1) | GCA_057233065.1 | real T2T | 0 | 18 / 18 | 10 / 18 | 3803 | [svg](annotations/Heloderma_suspectum__16-S_hap1__GCA_057233065.1/karyotype.svg) |
| *Hemibarbus labeo* | [ASM5663265v1](annotations/Hemibarbus_labeo__ASM5663265v1__GCA_056632655.1) | GCA_056632655.1 | gap-free, telomeres incomplete | 0 | 1 / 25 | 20 / 25 | 205 | [svg](annotations/Hemibarbus_labeo__ASM5663265v1__GCA_056632655.1/karyotype.svg) |
| *Homo sapiens* | [HuRef2.0](annotations/Homo_sapiens__HuRef2.0__GCA_059792245.1) | GCA_059792245.1 | real T2T | 0 | 24 / 24 | 24 / 24 | 171 | [svg](annotations/Homo_sapiens__HuRef2.0__GCA_059792245.1/karyotype.svg) |
| *Homo sapiens* | [T2T-CHM13v2.0](annotations/Homo_sapiens__T2T-CHM13v2.0__GCA_009914755.4) | GCA_009914755.4 | real T2T | 0 | 24 / 24 | 24 / 24 | 171 | [svg](annotations/Homo_sapiens__T2T-CHM13v2.0__GCA_009914755.4/karyotype.svg) |
| *Hucho bleekeri* | [ASM5731887v1](annotations/Hucho_bleekeri__ASM5731887v1__GCA_057318875.1) | GCA_057318875.1 | gap-free, telomeres incomplete | 0 | 11 / 45 | 45 / 45 | 35 | [svg](annotations/Hucho_bleekeri__ASM5731887v1__GCA_057318875.1/karyotype.svg) |
| *Hydropotes inermis* | [ASM5624782v2](annotations/Hydropotes_inermis__ASM5624782v2__GCA_056247825.2) | GCA_056247825.2 | gap-free, telomeres incomplete | 0 | 15 / 35 | 34 / 35 | 980 | [svg](annotations/Hydropotes_inermis__ASM5624782v2__GCA_056247825.2/karyotype.svg) |
| *Hypomesus nipponensis* | [ASM5449105v1](annotations/Hypomesus_nipponensis__ASM5449105v1__GCA_054491055.1) | GCA_054491055.1 | real T2T | 0 | 28 / 28 | 28 / 28 | 196 | [svg](annotations/Hypomesus_nipponensis__ASM5449105v1__GCA_054491055.1/karyotype.svg) |
| *Ilisha elongata* | [ASM4812638v1](annotations/Ilisha_elongata__ASM4812638v1__GCA_048126385.1) | GCA_048126385.1 | gap-free, telomeres incomplete | 0 | 0 / 24 | 8 / 24 | 7065 | [svg](annotations/Ilisha_elongata__ASM4812638v1__GCA_048126385.1/karyotype.svg) |
| *Larimichthys crocea* | [ASM5180228v1](annotations/Larimichthys_crocea__ASM5180228v1__GCA_051802285.1) | GCA_051802285.1 | gap-free, telomeres incomplete | 0 | 23 / 24 | 22 / 24 | 42 | [svg](annotations/Larimichthys_crocea__ASM5180228v1__GCA_051802285.1/karyotype.svg) |
| *Larimichthys crocea* | [ASM5180229v1](annotations/Larimichthys_crocea__ASM5180229v1__GCA_051802295.1) | GCA_051802295.1 | gap-free, telomeres incomplete | 0 | 19 / 24 | 21 / 24 | 42 | [svg](annotations/Larimichthys_crocea__ASM5180229v1__GCA_051802295.1/karyotype.svg) |
| *Larimichthys crocea* | [Lc-T2T](annotations/Larimichthys_crocea__Lc-T2T__GCA_059630365.1) | GCA_059630365.1 | gap-free, telomeres incomplete | 0 | 12 / 25 | 22 / 25 | 42 | [svg](annotations/Larimichthys_crocea__Lc-T2T__GCA_059630365.1/karyotype.svg) |
| *Lates calcarifer* | [ASM5102725v1](annotations/Lates_calcarifer__ASM5102725v1__GCA_051027255.1) | GCA_051027255.1 | gap-free, telomeres incomplete | 0 | 22 / 24 | 22 / 24 | 38 | [svg](annotations/Lates_calcarifer__ASM5102725v1__GCA_051027255.1/karyotype.svg) |
| *Lepus yarkandensis* | [LepYark_1.0](annotations/Lepus_yarkandensis__LepYark_1.0__GCA_047496845.1) | GCA_047496845.1 | gap-free, telomeres incomplete | 0 | 7 / 24 | 22 / 24 | 372 + co-family | [svg](annotations/Lepus_yarkandensis__LepYark_1.0__GCA_047496845.1/karyotype.svg) |
| *Lutjanus argentimaculatus* | [GDOU_Ler_1.0](annotations/Lutjanus_argentimaculatus__GDOU_Ler_1.0__GCA_059048685.1) | GCA_059048685.1 | gap-free, telomeres incomplete | 0 | 0 / 24 | 24 / 24 | 88 | [svg](annotations/Lutjanus_argentimaculatus__GDOU_Ler_1.0__GCA_059048685.1/karyotype.svg) |
| *Macaca fascicularis* | [T2T-MFA8v1.1](annotations/Macaca_fascicularis__T2T-MFA8v1.1__GCA_037993035.2) | GCA_037993035.2 | real T2T | 0 | 22 / 22 | 22 / 22 | 171 | [svg](annotations/Macaca_fascicularis__T2T-MFA8v1.1__GCA_037993035.2/karyotype.svg) |
| *Macaca mulatta* | [T2T-MMU8v2.0](annotations/Macaca_mulatta__T2T-MMU8v2.0__GCA_049350105.2) | GCA_049350105.2 | real T2T | 0 | 22 / 22 | 22 / 22 | 171 | [svg](annotations/Macaca_mulatta__T2T-MMU8v2.0__GCA_049350105.2/karyotype.svg) |
| *Micropterus salmoides* | [T2T_MSal_F1.0](annotations/Micropterus_salmoides__T2T_MSal_F1.0__GCA_056957265.1) | GCA_056957265.1 | gap-free, telomeres incomplete | 0 | 12 / 23 | 23 / 23 | 608 | [svg](annotations/Micropterus_salmoides__T2T_MSal_F1.0__GCA_056957265.1/karyotype.svg) |
| *Micropterus salmoides* | [T2T_MSal_M1.0](annotations/Micropterus_salmoides__T2T_MSal_M1.0__GCA_056957295.1) | GCA_056957295.1 | gap-free, telomeres incomplete | 0 | 17 / 23 | 23 / 23 | 608 | [svg](annotations/Micropterus_salmoides__T2T_MSal_M1.0__GCA_056957295.1/karyotype.svg) |
| *Monopterus albus* | [Monopterus_albus_T2T_assembly](annotations/Monopterus_albus__Monopterus_albus_T2T_assembly__GCA_986269415.1) | GCA_986269415.1 | gap-free, telomeres incomplete | 0 | 2 / 12 | 12 / 12 | 237 | [svg](annotations/Monopterus_albus__Monopterus_albus_T2T_assembly__GCA_986269415.1/karyotype.svg) |
| *Mus musculus* | [T2T_mhaESC_v1.0](annotations/Mus_musculus__T2T_mhaESC_v1.0__GCA_050437135.1) | GCA_050437135.1 | gap-free, telomeres incomplete | 0 | 18 / 20 | 20 / 20 | 58 | [svg](annotations/Mus_musculus__T2T_mhaESC_v1.0__GCA_050437135.1/karyotype.svg) |
| *Mus musculus domesticus* | [T2T_mhaESC_v1.5](annotations/Mus_musculus_domesticus__T2T_mhaESC_v1.5__GCA_056825265.1) | GCA_056825265.1 | gap-free, telomeres incomplete | 0 | 19 / 21 | 20 / 21 | 58 | [svg](annotations/Mus_musculus_domesticus__T2T_mhaESC_v1.5__GCA_056825265.1/karyotype.svg) |
| *Opsariichthys bidens* | [HZAU_OBidXX_1.0](annotations/Opsariichthys_bidens__HZAU_OBidXX_1.0__GCA_060511635.1) | GCA_060511635.1 | gap-free, telomeres incomplete | 0 | 26 / 39 | 34 / 39 | 188 + co-family | [svg](annotations/Opsariichthys_bidens__HZAU_OBidXX_1.0__GCA_060511635.1/karyotype.svg) |
| *Opsariichthys bidens* | [HZAU_OBidYY_1.0](annotations/Opsariichthys_bidens__HZAU_OBidYY_1.0__GCA_060251745.1) | GCA_060251745.1 | gap-free, telomeres incomplete | 0 | 29 / 37 | 21 / 37 | 536 | [svg](annotations/Opsariichthys_bidens__HZAU_OBidYY_1.0__GCA_060251745.1/karyotype.svg) |
| *Ovis aries* | [NLFDP_OvArieNajdiFemle_T2T_1.0](annotations/Ovis_aries__NLFDP_OvArieNajdiFemle_T2T_1.0__GCA_057332215.1) | GCA_057332215.1 | gap-free, telomeres incomplete | 0 | 0 / 27 | 23 / 27 | 816 | [svg](annotations/Ovis_aries__NLFDP_OvArieNajdiFemle_T2T_1.0__GCA_057332215.1/karyotype.svg) |
| *Ovis aries* | [NLFDP_OvArieNajdiMale_T2T_1.0](annotations/Ovis_aries__NLFDP_OvArieNajdiMale_T2T_1.0__GCA_057332225.1) | GCA_057332225.1 | gap-free, telomeres incomplete | 0 | 0 / 28 | 24 / 28 | 816 | [svg](annotations/Ovis_aries__NLFDP_OvArieNajdiMale_T2T_1.0__GCA_057332225.1/karyotype.svg) |
| *Ovis aries* | [T2T-sheep1.0](annotations/Ovis_aries__T2T-sheep1.0__GCA_040805955.1) | GCA_040805955.1 | gap-free, telomeres incomplete | 0 | 26 / 28 | 28 / 28 | 816 | [svg](annotations/Ovis_aries__T2T-sheep1.0__GCA_040805955.1/karyotype.svg) |
| *Ovis aries* | [T2T-sheep1.0M](annotations/Ovis_aries__T2T-sheep1.0M__GCA_041121725.1) | GCA_041121725.1 | gaps | 10 | 4 / 27 | 25 / 27 | 816 | [svg](annotations/Ovis_aries__T2T-sheep1.0M__GCA_041121725.1/karyotype.svg) |
| *Ovis aries* | [T2T-sheep1.0P](annotations/Ovis_aries__T2T-sheep1.0P__GCA_040836155.1) | GCA_040836155.1 | gaps | 4 | 3 / 27 | 24 / 27 | 816 | [svg](annotations/Ovis_aries__T2T-sheep1.0P__GCA_040836155.1/karyotype.svg) |
| *Pan paniscus* | [NHGRI_mPanPan1-v2.0_mat](annotations/Pan_paniscus__NHGRI_mPanPan1-v2.0_mat__GCA_028858845.2) | GCA_028858845.2 | gaps | 5 | 24 / 24 | 24 / 24 | 171 | [svg](annotations/Pan_paniscus__NHGRI_mPanPan1-v2.0_mat__GCA_028858845.2/karyotype.svg) |
| *Pan paniscus* | [NHGRI_mPanPan1-v2.0_pat](annotations/Pan_paniscus__NHGRI_mPanPan1-v2.0_pat__GCA_028858825.2) | GCA_028858825.2 | gaps | 7 | 23 / 24 | 24 / 24 | 171 | [svg](annotations/Pan_paniscus__NHGRI_mPanPan1-v2.0_pat__GCA_028858825.2/karyotype.svg) |
| *Pan paniscus* | [NHGRI_mPanPan1-v2.1_pri](annotations/Pan_paniscus__NHGRI_mPanPan1-v2.1_pri__GCA_029289425.3) | GCA_029289425.3 | gaps | 6 | 25 / 25 | 25 / 25 | 171 | [svg](annotations/Pan_paniscus__NHGRI_mPanPan1-v2.1_pri__GCA_029289425.3/karyotype.svg) |
| *Pan troglodytes* | [NHGRI_mPanTro3-v2.0_alt](annotations/Pan_troglodytes__NHGRI_mPanTro3-v2.0_alt__GCA_028858805.2) | GCA_028858805.2 | gaps | 5 | 20 / 23 | 23 / 23 | 171 | [svg](annotations/Pan_troglodytes__NHGRI_mPanTro3-v2.0_alt__GCA_028858805.2/karyotype.svg) |
| *Pan troglodytes* | [NHGRI_mPanTro3-v2.1_pri](annotations/Pan_troglodytes__NHGRI_mPanTro3-v2.1_pri__GCA_028858775.3) | GCA_028858775.3 | gaps | 5 | 25 / 25 | 25 / 25 | 171 | [svg](annotations/Pan_troglodytes__NHGRI_mPanTro3-v2.1_pri__GCA_028858775.3/karyotype.svg) |
| *Pangasianodon hypophthalmus* | [ASM5136230v1](annotations/Pangasianodon_hypophthalmus__ASM5136230v1__GCA_051362305.1) | GCA_051362305.1 | gap-free, telomeres incomplete | 0 | 21 / 30 | 23 / 30 | 300 | [svg](annotations/Pangasianodon_hypophthalmus__ASM5136230v1__GCA_051362305.1/karyotype.svg) |
| *Parabramis pekinensis* | [BGI_Ppek_1.0](annotations/Parabramis_pekinensis__BGI_Ppek_1.0__GCA_057375095.1) | GCA_057375095.1 | gap-free, telomeres incomplete | 0 | 11 / 24 | 21 / 24 | 102 | [svg](annotations/Parabramis_pekinensis__BGI_Ppek_1.0__GCA_057375095.1/karyotype.svg) |
| *Pelodiscus sinensis* | [ASM5242630v1](annotations/Pelodiscus_sinensis__ASM5242630v1__GCA_052426305.1) | GCA_052426305.1 | gap-free, telomeres incomplete | 0 | 18 / 34 | 27 / 34 | 1131 | [svg](annotations/Pelodiscus_sinensis__ASM5242630v1__GCA_052426305.1/karyotype.svg) |
| *Platax teira* | [ASM5870766v1](annotations/Platax_teira__ASM5870766v1__GCA_058707665.1) | GCA_058707665.1 | gap-free, telomeres incomplete | 0 | 0 / 24 | 15 / 24 | 31 | [svg](annotations/Platax_teira__ASM5870766v1__GCA_058707665.1/karyotype.svg) |
| *Platax teira* | [ASM6109948v1](annotations/Platax_teira__ASM6109948v1__GCA_061099485.1) | GCA_061099485.1 | gap-free, telomeres incomplete | 0 | 0 / 24 | 15 / 24 | 31 | [svg](annotations/Platax_teira__ASM6109948v1__GCA_061099485.1/karyotype.svg) |
| *Plectropomus leopardus* | [Pleo_t2t_v1.0](annotations/Plectropomus_leopardus__Pleo_t2t_v1.0__GCA_060308525.1) | GCA_060308525.1 | gap-free, telomeres incomplete | 0 | 21 / 24 | 24 / 24 | 168 | [svg](annotations/Plectropomus_leopardus__Pleo_t2t_v1.0__GCA_060308525.1/karyotype.svg) |
| *Pongo abelii* | [NHGRI_mPonAbe1-v2.0_alt](annotations/Pongo_abelii__NHGRI_mPonAbe1-v2.0_alt__GCA_028885685.2) | GCA_028885685.2 | gaps | 10 | 19 / 23 | 23 / 23 | 171 | [svg](annotations/Pongo_abelii__NHGRI_mPonAbe1-v2.0_alt__GCA_028885685.2/karyotype.svg) |
| *Pongo abelii* | [NHGRI_mPonAbe1-v2.1_pri](annotations/Pongo_abelii__NHGRI_mPonAbe1-v2.1_pri__GCA_028885655.3) | GCA_028885655.3 | gaps | 11 | 25 / 25 | 25 / 25 | 171 | [svg](annotations/Pongo_abelii__NHGRI_mPonAbe1-v2.1_pri__GCA_028885655.3/karyotype.svg) |
| *Pongo pygmaeus* | [NHGRI_mPonPyg2-v2.0_alt](annotations/Pongo_pygmaeus__NHGRI_mPonPyg2-v2.0_alt__GCA_028885525.2) | GCA_028885525.2 | gaps | 10 | 19 / 23 | 23 / 23 | 171 | [svg](annotations/Pongo_pygmaeus__NHGRI_mPonPyg2-v2.0_alt__GCA_028885525.2/karyotype.svg) |
| *Pongo pygmaeus* | [NHGRI_mPonPyg2-v2.1_pri](annotations/Pongo_pygmaeus__NHGRI_mPonPyg2-v2.1_pri__GCA_028885625.3) | GCA_028885625.3 | gaps | 8 | 24 / 25 | 25 / 25 | 171 | [svg](annotations/Pongo_pygmaeus__NHGRI_mPonPyg2-v2.1_pri__GCA_028885625.3/karyotype.svg) |
| *Premnas biaculeatus* | [ASM5381358v1](annotations/Premnas_biaculeatus__ASM5381358v1__GCA_053813585.1) | GCA_053813585.1 | gap-free, telomeres incomplete | 0 | 13 / 24 | 16 / 24 | 188 | [svg](annotations/Premnas_biaculeatus__ASM5381358v1__GCA_053813585.1/karyotype.svg) |
| *Pseudopleuronectes yokohamae* | [ASM5131432v1](annotations/Pseudopleuronectes_yokohamae__ASM5131432v1__GCA_051314325.1) | GCA_051314325.1 | gap-free, telomeres incomplete | 0 | 8 / 24 | 20 / 24 | 625 + co-family | [svg](annotations/Pseudopleuronectes_yokohamae__ASM5131432v1__GCA_051314325.1/karyotype.svg) |
| *Pungitius pungitius* | [TVA-T2T](annotations/Pungitius_pungitius__TVA-T2T__GCA_056822235.1) | GCA_056822235.1 | gap-free, telomeres incomplete | 0 | 8 / 22 | 21 / 22 | 177 | [svg](annotations/Pungitius_pungitius__TVA-T2T__GCA_056822235.1/karyotype.svg) |
| *Riparia riparia* | [ICG_Riri_1](annotations/Riparia_riparia__ICG_Riri_1__GCA_055001455.1) | GCA_055001455.1 | gap-free, telomeres incomplete | 0 | 2 / 40 | 12 / 40 | 37 | [svg](annotations/Riparia_riparia__ICG_Riri_1__GCA_055001455.1/karyotype.svg) |
| *Scleropages jardinii* | [ASM5488336v1](annotations/Scleropages_jardinii__ASM5488336v1__GCA_054883365.1) | GCA_054883365.1 | gap-free, telomeres incomplete | 0 | 11 / 24 | 18 / 24 | 60 + co-family | [svg](annotations/Scleropages_jardinii__ASM5488336v1__GCA_054883365.1/karyotype.svg) |
| *Seriola dumerili* | [ASM5688231v1](annotations/Seriola_dumerili__ASM5688231v1__GCA_056882315.1) | GCA_056882315.1 | gap-free, telomeres incomplete | 0 | 5 / 24 | 23 / 24 | 226 | [svg](annotations/Seriola_dumerili__ASM5688231v1__GCA_056882315.1/karyotype.svg) |
| *Seriola dumerili* | [ASM5688233v1](annotations/Seriola_dumerili__ASM5688233v1__GCA_056882335.1) | GCA_056882335.1 | gap-free, telomeres incomplete | 0 | 8 / 24 | 24 / 24 | 226 | [svg](annotations/Seriola_dumerili__ASM5688233v1__GCA_056882335.1/karyotype.svg) |
| *Siniperca chuatsi* | [ASM5996270v1](annotations/Siniperca_chuatsi__ASM5996270v1__GCA_059962705.1) | GCA_059962705.1 | gap-free, telomeres incomplete | 0 | 15 / 24 | 16 / 24 | 3837 | [svg](annotations/Siniperca_chuatsi__ASM5996270v1__GCA_059962705.1/karyotype.svg) |
| *Sparus aurata* | [ASM6129181v1](annotations/Sparus_aurata__ASM6129181v1__GCA_061291815.1) | GCA_061291815.1 | gap-free, telomeres incomplete | 0 | 0 / 24 | 22 / 24 | 187 | [svg](annotations/Sparus_aurata__ASM6129181v1__GCA_061291815.1/karyotype.svg) |
| *Sus scrofa* | [ASM3070493v2](annotations/Sus_scrofa__ASM3070493v2__GCA_030704935.2) | GCA_030704935.2 | gap-free, telomeres incomplete | 0 | 16 / 20 | 16 / 20 | 336 | [svg](annotations/Sus_scrofa__ASM3070493v2__GCA_030704935.2/karyotype.svg) |
| *Sus scrofa* | [ASM5778240v1](annotations/Sus_scrofa__ASM5778240v1__GCA_057782405.1) | GCA_057782405.1 | gap-free, telomeres incomplete | 0 | 19 / 20 | 19 / 20 | 336 | [svg](annotations/Sus_scrofa__ASM5778240v1__GCA_057782405.1/karyotype.svg) |
| *Sus scrofa* | [ASM5999548v1](annotations/Sus_scrofa__ASM5999548v1__GCA_059995485.1) | GCA_059995485.1 | gap-free, telomeres incomplete | 0 | 16 / 20 | 16 / 20 | 336 | [svg](annotations/Sus_scrofa__ASM5999548v1__GCA_059995485.1/karyotype.svg) |
| *Sus scrofa* | [CAU_T2T_PIGHN1](annotations/Sus_scrofa__CAU_T2T_PIGHN1__GCA_059440985.1) | GCA_059440985.1 | gap-free, telomeres incomplete | 0 | 16 / 20 | 19 / 20 | 336 | [svg](annotations/Sus_scrofa__CAU_T2T_PIGHN1__GCA_059440985.1/karyotype.svg) |
| *Sus scrofa* | [CAU_T2T_PIGJH1](annotations/Sus_scrofa__CAU_T2T_PIGJH1__GCA_059440895.1) | GCA_059440895.1 | gap-free, telomeres incomplete | 0 | 17 / 20 | 19 / 20 | 336 | [svg](annotations/Sus_scrofa__CAU_T2T_PIGJH1__GCA_059440895.1/karyotype.svg) |
| *Sus scrofa* | [CAU_T2T_PIGTW1](annotations/Sus_scrofa__CAU_T2T_PIGTW1__GCA_059440865.1) | GCA_059440865.1 | gap-free, telomeres incomplete | 0 | 15 / 19 | 19 / 19 | 336 | [svg](annotations/Sus_scrofa__CAU_T2T_PIGTW1__GCA_059440865.1/karyotype.svg) |
| *Sus scrofa* | [T2T-Sscrofa](annotations/Sus_scrofa__T2T-Sscrofa__GCA_000003025.7) | GCA_000003025.7 | gap-free, telomeres incomplete | 0 | 19 / 20 | 19 / 20 | 336 | [svg](annotations/Sus_scrofa__T2T-Sscrofa__GCA_000003025.7/karyotype.svg) |
| *Sus scrofa* | [T2T_pig1.0](annotations/Sus_scrofa__T2T_pig1.0__GCA_048338725.1) | GCA_048338725.1 | real T2T | 0 | 20 / 20 | 20 / 20 | 336 | [svg](annotations/Sus_scrofa__T2T_pig1.0__GCA_048338725.1/karyotype.svg) |
| *Symphalangus syndactylus* | [NHGRI_mSymSyn1-v2.0_alt](annotations/Symphalangus_syndactylus__NHGRI_mSymSyn1-v2.0_alt__GCA_028878085.3) | GCA_028878085.3 | gaps | 3 | 20 / 24 | 24 / 24 | 171 | [svg](annotations/Symphalangus_syndactylus__NHGRI_mSymSyn1-v2.0_alt__GCA_028878085.3/karyotype.svg) |
| *Symphalangus syndactylus* | [NHGRI_mSymSyn1-v2.1_pri](annotations/Symphalangus_syndactylus__NHGRI_mSymSyn1-v2.1_pri__GCA_028878055.3) | GCA_028878055.3 | gaps | 1 | 26 / 26 | 26 / 26 | 171 | [svg](annotations/Symphalangus_syndactylus__NHGRI_mSymSyn1-v2.1_pri__GCA_028878055.3/karyotype.svg) |
| *Tachysurus ussuriensis* | [ASM5413144v1](annotations/Tachysurus_ussuriensis__ASM5413144v1__GCA_054131445.1) | GCA_054131445.1 | gap-free, telomeres incomplete | 0 | 23 / 26 | 24 / 26 | 2531 | [svg](annotations/Tachysurus_ussuriensis__ASM5413144v1__GCA_054131445.1/karyotype.svg) |
| *Tetraodon nigroviridis* | [TetNigT2Thap1](annotations/Tetraodon_nigroviridis__TetNigT2Thap1__GCA_051020865.1) | GCA_051020865.1 | gap-free, telomeres incomplete | 0 | 17 / 21 | 21 / 21 | 118 | [svg](annotations/Tetraodon_nigroviridis__TetNigT2Thap1__GCA_051020865.1/karyotype.svg) |
| *Tetraodon nigroviridis* | [TetNigT2Thap2](annotations/Tetraodon_nigroviridis__TetNigT2Thap2__GCA_051020875.1) | GCA_051020875.1 | gap-free, telomeres incomplete | 0 | 18 / 21 | 21 / 21 | 118 | [svg](annotations/Tetraodon_nigroviridis__TetNigT2Thap2__GCA_051020875.1/karyotype.svg) |
| *Thunnus orientalis* | [ASM61617105v1](annotations/Thunnus_orientalis__ASM61617105v1__GCA_061617105.1) | GCA_061617105.1 | gap-free, telomeres incomplete | 0 | 21 / 24 | 20 / 24 | 127 | [svg](annotations/Thunnus_orientalis__ASM61617105v1__GCA_061617105.1/karyotype.svg) |
| *Tinca tinca* | [ASM61535325v1](annotations/Tinca_tinca__ASM61535325v1__GCA_061535325.1) | GCA_061535325.1 | gap-free, telomeres incomplete | 0 | 18 / 24 | 22 / 24 | 106 | [svg](annotations/Tinca_tinca__ASM61535325v1__GCA_061535325.1/karyotype.svg) |
| *Trachinotus anak* | [Tana_hap1](annotations/Trachinotus_anak__Tana_hap1__GCA_977005155.1) | GCA_977005155.1 | gap-free, telomeres incomplete | 0 | 22 / 24 | 19 / 24 | 109 | [svg](annotations/Trachinotus_anak__Tana_hap1__GCA_977005155.1/karyotype.svg) |
| *Trachinotus anak* | [Tana_hap2](annotations/Trachinotus_anak__Tana_hap2__GCA_977005145.1) | GCA_977005145.1 | real T2T | 0 | 24 / 24 | 18 / 24 | 109 | [svg](annotations/Trachinotus_anak__Tana_hap2__GCA_977005145.1/karyotype.svg) |
| *Trachinotus ovatus* | [ASM6119389v1](annotations/Trachinotus_ovatus__ASM6119389v1__GCA_061193895.1) | GCA_061193895.1 | gap-free, telomeres incomplete | 0 | 0 / 24 | 7 / 24 | 109 | [svg](annotations/Trachinotus_ovatus__ASM6119389v1__GCA_061193895.1/karyotype.svg) |
| *Trachinotus ovatus* | [ASM6119397v1](annotations/Trachinotus_ovatus__ASM6119397v1__GCA_061193975.1) | GCA_061193975.1 | gap-free, telomeres incomplete | 0 | 5 / 24 | 8 / 24 | 118 | [svg](annotations/Trachinotus_ovatus__ASM6119397v1__GCA_061193975.1/karyotype.svg) |
| *Triplophysa stolickai* | [ASM5794468v1](annotations/Triplophysa_stolickai__ASM5794468v1__GCA_057944685.1) | GCA_057944685.1 | gap-free, telomeres incomplete | 0 | 11 / 25 | 25 / 25 | 128 | [svg](annotations/Triplophysa_stolickai__ASM5794468v1__GCA_057944685.1/karyotype.svg) |
| *Triplophysa yaopeizhii* | [HZAU_TYAO_T2T_v1](annotations/Triplophysa_yaopeizhii__HZAU_TYAO_T2T_v1__GCA_048296945.1) | GCA_048296945.1 | gap-free, telomeres incomplete | 0 | 24 / 25 | 23 / 25 | 128 | [svg](annotations/Triplophysa_yaopeizhii__HZAU_TYAO_T2T_v1__GCA_048296945.1/karyotype.svg) |

| assembly | published set | chromosomes compared | agrees | contested: no overlap | contested: no call |
|---|---|---:|---:|---:|---:|
| T2T-CHM13v2.0 | T2T_CHM13_censat_v2.1 (active_hor) | 24 | 24 | 0 | 0 |
| NHGRI_mPanTro3-v2.1_pri | T2T_primates_CenSat_v2.0 (active_hor) | 25 | 25 | 0 | 0 |
| NHGRI_mPanTro3-v2.0_alt | T2T_primates_CenSat_v2.0 (active_hor) | 23 | 22 | 1 | 0 |
| NHGRI_mPanPan1-v2.0_pat | T2T_primates_CenSat_v2.0 (active_hor) | 24 | 24 | 0 | 0 |
| NHGRI_mPanPan1-v2.0_mat | T2T_primates_CenSat_v2.0 (active_hor) | 24 | 23 | 1 | 0 |
| NHGRI_mSymSyn1-v2.1_pri | T2T_primates_CenSat_v2.0 (active_hor) | 26 | 23 | 3 | 0 |
| NHGRI_mSymSyn1-v2.0_alt | T2T_primates_CenSat_v2.0 (active_hor) | 24 | 21 | 3 | 0 |
| NHGRI_mGorGor1-v2.0_pat | T2T_primates_CenSat_v2.0 (active_hor) | 24 | 24 | 0 | 0 |
| NHGRI_mGorGor1-v2.0_mat | T2T_primates_CenSat_v2.0 (active_hor) | 24 | 24 | 0 | 0 |
| NHGRI_mPonPyg2-v2.0_alt | T2T_primates_CenSat_v2.0 (active_hor) | 23 | 23 | 0 | 0 |
| NHGRI_mPonPyg2-v2.1_pri | T2T_primates_CenSat_v2.0 (active_hor) | 25 | 25 | 0 | 0 |
| NHGRI_mPonAbe1-v2.1_pri | T2T_primates_CenSat_v2.0 (active_hor) | 24 | 24 | 0 | 0 |
| NHGRI_mPonAbe1-v2.0_alt | T2T_primates_CenSat_v2.0 (active_hor) | 22 | 22 | 0 | 0 |
| NHGRI_mGorGor1-v2.1_pri | T2T_primates_CenSat_v2.0 (active_hor) | 25 | 25 | 0 | 0 |
| NHGRI_mPanPan1-v2.1_pri | T2T_primates_CenSat_v2.0 (active_hor) | 25 | 24 | 1 | 0 |
| T2T-MFA8v1.1 | T2T-MFA8_Nature2025_S39 (presumed_active_array) | 21 | 21 | 0 | 0 |
| T2T-MFA8v1.1 | T2T-MFA8_release_v1.1 (centromere_region) | 22 | 22 | 0 | 0 |
| T2T-sheep1.0 | T2T-sheep_NatGenet2025_S8 (gene_in_centromeric_region) | 25 | 16 | 9 | 0 |
| T2T-MMU8v2.0 | T2T-MMU8_release_v2.0 (alpha_satellite_region) | 22 | 22 | 0 | 0 |

**Contested** — every chromosome where our call contradicts a published annotation (`agreement` in `official_overlap.tsv`, `published_agreement=contested` in the GFF):

- T2T-sheep1.0 vs T2T-sheep_NatGenet2025_S8, no overlap: 5, 7, 14, 16, 19, 21, 22, 24, X
- NHGRI_mPanPan1-v2.0_mat vs T2T_primates_CenSat_v2.0, no overlap: chr17_mat_hsa18
- NHGRI_mPanPan1-v2.1_pri vs T2T_primates_CenSat_v2.0, no overlap: chr17_mat_hsa18
- NHGRI_mPanTro3-v2.0_alt vs T2T_primates_CenSat_v2.0, no overlap: chr6_hap2_hsa7
- NHGRI_mSymSyn1-v2.0_alt vs T2T_primates_CenSat_v2.0, no overlap: chr19_hap2 (molecule 12), chr15_hap2, chr21_hap2
- NHGRI_mSymSyn1-v2.1_pri vs T2T_primates_CenSat_v2.0, no overlap: chr19_hap1 (molecule 12), chr21_hap1, chrY_hap2
<!-- END TABLES -->

### What the comparison shows

Where our call contradicts a published annotation it is labelled **contested**, never left for the reader to
find: the `agreement` column of `official_overlap.tsv`, `published_agreement=contested` with a `Note` in the
GFF, a red mark in the karyotype, and the generated list above. Two kinds: our locus does not overlap the
published interval; the authors place a centromere where we make no call.

- **CHM13, the one independent test.** Our locus overlaps the live HOR of the T2T Consortium censat v2.1
  on 24 of 24 chromosomes. Taking the largest array of each chromosome
  instead lands on the live HOR on 16 of 24 (HSat1/2/3 arrays are larger on 8).
- **Sheep, indirect.** The authors published genes in their centromeric regions, not the regions. Our
  locus overlaps the span of those genes on 16 of 25 chromosomes. chrX is the clear miss: the only large
  satellite array beside the authors' genes (0.65–7.2 Mb) is a 1.1-Mb SatII (702 bp) array at 7.22–8.35 Mb.
  SatII is on 14 of the 28 sheep chromosomes, too few to join as a co-family, so the call is a 22-kb SatI
  period-class array at 91.9 Mb, flagged `other_family_larger_on_chromosome`.
- **Apes and siamang.** CenSat `active_hor` is the longest high-confidence (SF1–SF3) alpha-satellite HOR
  array of the chromosome (the rule since CenSat v0.1), except on 15 contigs where the authors set the active
  array by hand from CENP-A localisation and the CDR ([T2TPrimates README](https://github.com/hloucks/CenSatData/blob/main/T2TPrimates/README.md),
  [`Primate_alpha_summary.sh`](https://github.com/hloucks/CenSatData/blob/main/T2TPrimates/Code/Primate_alpha_summary.sh)).
  Agreement of our locus with a size-defined label is consistency rather than validation; a contested
  chromosome is a disagreement, not proof that either annotation is wrong.
- **Macaques.** The published intervals are defined from sequence (presumed-active arrays, alpha-satellite
  spans, RepeatMasker), so agreement is consistency rather than validation.
- **Chromosomes without a call** (`outliers.tsv`): chrY of pig, goats and giraffe, and horse Chr11, which
  the T2T-horse authors describe as satellite-free. On the haplotype assemblies of goat (chr9 and chr21
  maternal, chr28 paternal) and sheep (chr7 and chr14, plus chr15 paternal) the chromosome is 2.4–20.8 Mb
  shorter than the same chromosome in the other assemblies of that individual and carries no array of the
  centromeric family: the satellite block is absent from that assembly.
- **Clades added in v1.3.** With its rules unchanged, the caller places a locus on 1,412 of 1,609 chromosomes of
  the 63 fish assemblies, 426 of 498 of the 21 added mammal assemblies, 73 of 121 of the three birds (41 of 41 in
  duck, 20 of 40 in red-legged partridge, 12 of 40 in sand martin), 27 of 34 in the turtle and 10 of 18 in the
  Gila monster. Nothing is published to compare these with. The rules were set on mammals; where the call rate is
  low, the broadest satellite covers only part of the chromosomes, and `outliers.tsv` names the rest.
- **Reproducibility.** Two runs of the caller on the same catalogues gave byte-identical files on all 17
  assemblies compared. The NCBI-download route of `reproduce.sh` gave a satellome catalogue byte-identical
  to our local run (tested on GCA_048601705.1, the smoke-test fixture).

## How a centromere is called

1. **Arrays ≥ 10 kb** from the FasTAN BED, with the arraysplitter monomer / HOR period and consensus.
2. **Sequence families** — monomer consensus searched against doubled and reverse-complement-doubled
   copies of every other consensus (rotation- and strand-aware), BLASTN, single linkage at ≥ 0.70 identity
   over ≥ 0.80 of the query.
3. **Superfamilies** — a family joins the largest family (the founder) whose base monomer period is
   equal to or 2× its own (alpha 169 / 171 / 340 join; a 32-bp subterminal satellite does not). Families
   are compared with the founder only, never member to member, so periods cannot drift along a chain.
4. **The centromeric superfamily** is the one present (one array ≥ 10 kb) on the most chromosomes —
   breadth, not size: the largest array of a chromosome is HSat1/2/3 rather than alpha on 8 of 24 CHM13
   chromosomes. Superfamilies without a consensus ≥ 30 bp (microsatellites, the telomeric repeat) are not
   eligible. Broad superfamilies co-located with it (bovid SatI + SatII) join as co-families.
5. **Per chromosome** its arrays, plus arrays of the same monomer-period class, are merged across gaps
   < 1 Mb; the cluster with the most bp is the locus. If the superfamily also caps both chromosome ends on
   most chromosomes (siamang), an interior cluster is preferred.

Every choice and the measurement behind it: [`docs/PARAMETERS.md`](docs/PARAMETERS.md).

## Published annotations used for comparison

They are compared, never used as input. Only some of them rest on evidence other than the sequence itself:

| assemblies | published set | definition | independent? |
|---|---|---|---|
| CHM13 | T2T Consortium censat v2.1, live HOR `hor_*(…L)` | CENP-A / CDR-informed | yes |
| chimpanzee, bonobo, gorilla, orangutans | T2T Primates CenSat v2.0 `active_hor` | longest high-confidence HOR array of the chromosome; set from CENP-A and CDR on 14 contigs (bonobo 9, Sumatran orangutan 3, Bornean orangutan 2) | mostly no |
| siamang | T2T Primates CenSat v2.0 `active_hor` (joined by accession via the marbl/T2T-Browser alias table) | longest high-confidence HOR array; chrY_hap2 set from CENP-A and CDR | mostly no |
| T2T-MFA8 | Supp. Table 39 presumed-active S1S2 arrays; release `centromere.bed` | sequence only | no |
| T2T-MMU8 | release `alpha_satellite.bed` (RepeatMasker) | sequence only | no |
| T2T-sheep1.0 | Supp. Table 8, genes in the authors' centromeric regions | CENP-A ChIP + methylation + satellite | yes, indirectly |
| goats, pig, horse, giraffe | — no interval file published | | |

URLs and md5 of every file: [`config/official_downloads.tsv`](config/official_downloads.tsv). The public CenSat
v2.0 files include chrY of the assemblies that carry it (gorilla and bonobo primary and paternal; chimpanzee,
orangutan and siamang primary), and those 8 are compared; ape chrY tracks from other sources are not used.
Note for siamang: the NCBI assembly report keeps pre-v2.1 chromosome labels (`chr19_hap1` is molecule 12);
join on accession or molecule, never on that name.

## Reproduce

```bash
export PYTHONNOUSERSITE=1                # else pip skips packages it finds in ~/.local (arraysplitter's binary goes missing)
conda env create -f environment.yml && conda activate t2t_centromeres   # mamba env create: 1.5 min vs > 12 min
python -m satellome --install-fastan --install-tanbed --install-rust-tools
bash workflow/check_env.sh               # versions vs config/toolchain.tsv
bash tests/smoke_test.sh work_smoke/      # ~3 min: the 17-Mb MFA8 chrY accession end to end, md5 vs tests/expected/
bash workflow/reproduce.sh work/ 5 10    # all 121: 5 assemblies at a time, 10 threads each
```

`reproduce.sh` downloads each GenBank assembly from NCBI (md5-checked), runs satellome
(`workflow/run_satellome.sh`), fetches the published annotations (md5-checked, refuses on mismatch), calls
the centromeres (`workflow/call_all.sh`) and compares the result with `annotations/` file by file
(`workflow/verify.sh`). satellome took 15–71 min wall per assembly (median 36) at 10 threads, five assemblies
at a time (logs of the 32 release runs; Ubuntu 22.04, 192 cores shared with other jobs); most of it is the
monomer decomposition (median 64 % of a run), peak memory at most 4.04 GB. The 89 assemblies added in v1.3 took
2–42 min each (median 5; most are fish genomes under 1 Gb). One, the 3.45-Gb water deer GCA_056247825.2, ran past
satellome's one-hour limit for ArraySplitter when five assemblies ran at once and was rerun alone with 30 threads
(35 min); `run_satellome.sh` now refuses a run that stopped there instead of passing on an incomplete catalogue.
Calling all 121 takes 2 h 8 min at 16 threads. Tested on Linux only.

Toolchain the annotations were produced with: satellome 1.17.2, FasTAN 0.8 (ad3002/FASTAN e4953e9),
tanbed (richarddurbin/alntools 1fb9f4f), arraysplitter 1.7.4, BLAST+ 2.12.0, Python 3.10
([`config/toolchain.tsv`](config/toolchain.tsv)).

## Changes in v1.3

- **89 assemblies added, chosen by sequence statistics.** Every GenBank vertebrate assembly without a gap or an
  unplaced sequence (NCBI Datasets, 2026-10-16), whatever its name: 63 fishes, 21 mammals, 3 birds, a turtle and a
  lizard (selection rule [above](#which-assemblies-are-t2t)). `config/assemblies.tsv` has 121 rows.
- **Completeness tier.** `src/t2t_centromeres/tier.py` writes `annotations/TIERS.tsv` (gaps, unplaced sequences,
  telomeres per assembly); `SUMMARY.tsv` and the README tables carry the tier.
- **Chromosomes of type Linkage Group** count as chromosomes (three fish assemblies report all their chromosomes
  so); a mitochondrion deposited as an unlocalized scaffold is not counted as unplaced sequence.
- **`run_satellome.sh` refuses an incomplete satellome run** (no `gaps.bed` or telomere table), which happened once
  when ArraySplitter ran past satellome's one-hour limit.
- Every file of the 32 v1.2 assemblies is byte-identical in v1.3.

## Changes in v1.2

- **32 assemblies, not 33.** GCA_048601705.1 (T2T-MFA8v1.1 chrY deposited on its own) is the same GenBank
  sequence (CP168063.1, 17,124,583 bp) as chrY of GCA_037993035.2, which the panel already contains; v1.1
  counted that chromosome twice. It stays only as the input of `tests/smoke_test.sh` (expected md5 in
  `tests/expected/`).
- **The verdict is about the locus only.** v1.1 put the largest array of our locus forward as the candidate
  active array and labelled `contested_active_array` where a published active array was another array of the
  locus. The size of an array does not identify the active array, so that label and the candidate are gone;
  `largest_array_overlap_bp` stays in `official_overlap.tsv` as a descriptive column. The 18 such v1.1 rows
  now read `agrees` (our locus overlaps the published interval).
- **Karyotype figure.** Each published interval is its own bracket (v1.1 drew one bracket from the first to
  the last interval of a set, which on CHM13 chr3 and chr4 spanned a 42-bp satellite array between pieces of
  the live HOR); other satellite arrays are a darker grey, so an array never reads as empty sequence; no
  array is marked as the largest.
- **CenSat definition corrected.** The ape `active_hor` is the longest high-confidence HOR array of the
  chromosome, with the active array set from CENP-A and CDR on 15 contigs; v1.1 cited the CenSatData root
  README, which does not state the rule, and left out the curated contigs.
- The calls themselves (`centromeres.tsv`, arrays, families, outliers, exceptions) are byte-identical to v1.1.

## Sources

Assemblies: NCBI GenBank, accessions in [`config/assemblies.tsv`](config/assemblies.tsv). Papers whose
centromere annotations are compared (DOIs):
T2T apes — Yoo et al., *Nature* 2025, 10.1038/s41586-025-08816-3;
T2T-MFA8 — Zhang et al., *Nature* 2025, 10.1038/s41586-025-08596-w;
T2T-MMU8 — Zhang et al., *Cell* 2026, 10.1016/j.cell.2026.02.018;
T2T-sheep — Luo et al., *Nature Genetics* 2025, 10.1038/s41588-024-02037-6.
Assemblies without a published centromere interval file: T2T-goat (Wu et al., *Nature Communications*
2024, 10.1038/s41467-024-54188-z), T2T-Sscrofa (bioRxiv 10.1101/2025.10.02.679976), T2T-horse1.0
(bioRxiv 10.1101/2025.09.10.675315), T2T giraffe (*Cell Genomics* 2026, 10.1016/j.xgen.2026.101279).

## Licence

Code: MIT ([`LICENSE`](LICENSE)). Our annotations: CC BY 4.0; published intervals in the GFF remain their
authors' ([`annotations/LICENSE.md`](annotations/LICENSE.md)).
