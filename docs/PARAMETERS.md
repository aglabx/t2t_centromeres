# Parameters and why they have the values they have

Every value below is a choice (or inherited) unless marked as a constraint. Where a value was changed
during development, the measured reason is given. Defaults are the `call_centromeric_arrays.py` flags.

| parameter | value | flag | why |
|---|---|---|---|
| chromosome set | NCBI assembly report rows `assembled-molecule` / `Chromosome` (MT excluded) | — | constraint: no length cut-off needed |
| minimum array | 10 kb | `--min-array` | satellome's `10kb` class. Shorter admits interstitial minisatellites; 100 kb loses fragmented FasTAN pieces (goat 816-bp arrays come as 100–300 kb pieces) |
| family identity / cover | 0.70 / 0.80 | `--ident`, `--cov` | cross-species level of our earlier rotation-aware family clustering; 0.80 identity splits alpha suprachromosomal families within one genome |
| family search | BLASTN word 7, e-value 1e-5, plus strand vs doubled + reverse-complement-doubled targets, 1 HSP per pair, 200 targets per query | `--aligner`, `--max-targets` | MMSeqs2 (the first choice) returned only self-hits for the largest active alpha arrays of CHM13 chr18, chr11, chr6 and chrX. 200 targets keeps goat (3,790 near-identical consensus) tractable; CHM13 loci identical to the unlimited search |
| short consensus | < 30 bp grouped by canonical rotation | `--short-len` | alignment is meaningless for micro/minisatellite units |
| superfamily | sequence families whose base monomer periods are k × each other, k = 1..2, ± 3 % × k | `--period-tol` | bonobo (mat) alpha split into 18 sequence families and lost the breadth contest to a 32-bp subterminal satellite; k ≤ 2 because at k = 4 the window admits 703 ≈ 4 × 171 |
| chromosome carries a family | ≥ 100 kb of it | `--chrom-min-bp` | alternatives not measured |
| centromeric superfamily | on the most chromosomes, ties by bp | — | the largest array of a chromosome is HSat1/2/3, not alpha, on 8 of 24 CHM13 chromosomes |
| co-family | breadth ≥ 0.75 × chosen; largest clusters within the merge gap on ≥ 0.80 of shared chromosomes | `--cofamily-breadth`, `--cofamily-coloc` | goat1.0 and goat2.0 (same individual) flipped between SatI and SatII by 30 vs 29 chromosomes; with co-families both carry SatI + SatII |
| period class | monomer period = k × family base period, k = 1..4, ± 3 % × k | `--period-tol` | CHM13 chr2/6/7/X: sequence hits just under 0.70 / 0.80 |
| locus merge gap | 1 Mb | `--merge-gap` | sensitivity: CHM13 identical at 250 kb / 1 Mb / 2 Mb; chimp 21 vs 23 vs 23 official-within-ours; bonobo pat 6 / 6 / 4 disjoint; siamang 6 / 3 / 3 disjoint |
| cap mode | chosen superfamily in both terminal 10 % (≥ 100 kb each) on > 0.5 of its chromosomes → interior cluster (midpoint 10–90 %) preferred | — | siamang alpha forms multi-Mb subterminal caps on 23/26 chromosomes; without cap mode every siamang locus was a cap. Off for acrocentric genomes |
| ambiguous locus | second cluster ≥ 50 % of the first | `--ambiguous-frac` | flag only |
| terminal locus | within 1 % of an end | `--end-frac` | flag only (telocentric / acrocentric) |
