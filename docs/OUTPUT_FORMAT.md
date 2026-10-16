# Output format

## `centromeres.gff3`

GFF3, 1-based closed coordinates, seqid = GenBank accession of the chromosome (`##sequence-region` lines
give every chromosome's length). Four kinds of line:

| column 2 (source) | column 3 (type) | what |
|---|---|---|
| `satellome_cen` | `centromere` | our centromeric locus, one per chromosome, `ID=cen.<chromosome>` |
| `satellome_cen` | `satellite_DNA` | one array of that locus, `Parent=cen.<chromosome>` |
| publisher (e.g. `T2T_CHM13_censat_v2.1`) | `satellite_DNA` / `centromere` / `gene` | a published interval: active/live HOR or presumed-active array → `satellite_DNA`; alpha-satellite span → `centromere`; gene in the authors' centromeric region → `gene` |
| `intersection` | `region` | the stretch where our locus and one published interval overlap |

### Attributes of our locus (`centromere`)

| attribute | meaning |
|---|---|
| `chromosome`, `molecule` | NCBI Sequence-Name and assigned molecule (they can disagree; trust `molecule`) |
| `family` | centromeric satellite superfamily id within this assembly (`S001`; `S001+S003` = with a co-family) |
| `family_bp`, `n_arrays` | satellite bp and arrays of the family inside the locus |
| `rel_pos_mid` | locus midpoint / chromosome length |
| `second_locus_bp` | bp of the next-largest cluster of the family on this chromosome |
| `colocated_families` | other superfamilies (≥ 100 kb) within 1 Mb of the locus, `S00n:bp:m<monomer period>` |
| `<class>.source` | the published set compared, per class (`active_hor`, `presumed_active_array`, `centromere_region`, `alpha_satellite_region`, `gene_in_centromeric_region`) |
| `published_agreement` | `contested` if our call contradicts any published set on this chromosome, `agrees` if it agrees with every set that has an interval here, `not_compared` otherwise |
| `Note` | `contested by <set> <kind>` for every contradicting set |
| `<class>.relation` | see below |
| `<class>.agreement` | see below |
| `<class>.overlap_bp` | bp of our locus covered by the published intervals (for genes: by the span first..last gene) |
| `<class>.official_frac_in_locus` | overlap / published bp (for genes: / published span) |
| `<class>.locus_frac_official` | overlap / our locus length |
| `<class>.largest_array_overlap_bp` | overlap of the largest array of our locus with the published intervals; descriptive only (the size of an array does not identify the active array, so no verdict uses it) |
| `official=none_published` | the assembly has no published centromere interval file |

### Relation of our locus to a published set on one chromosome

Computed against the published span on that chromosome (first start .. last end of its intervals):

| relation | meaning |
|---|---|
| `official_within_ours` | the published span lies inside our locus |
| `ours_within_official` | our locus lies inside the published span |
| `partial_overlap` | they overlap, neither contains the other |
| `disjoint` | no overlap; `gap_bp` in `official_overlap.tsv` |
| `no_official_on_chromosome` | the set has no interval on this chromosome |
| `official_not_comparable` | the published track's assembly has a different length for this chromosome |
| `no_call` | the published set has an interval, we have no locus |

### Agreement with a published set on one chromosome

Every contradiction with a published annotation is labelled; none is left to be inferred from the relation.
The verdict is about the locus only. Our locus is the centromeric satellite zone; the functional centromere
(the active HOR, the CENP-A domain) is usually a small part of it, and the size of an array does not identify
it, so no array of the locus is put forward as the active one (changed in v1.2; v1.1 also labelled
`contested_active_array` when the published active array was not the largest array of our locus).

| agreement | meaning |
|---|---|
| `agrees` | our locus overlaps the published interval |
| `contested_locus` | our locus does not overlap the published interval (`relation=disjoint`) |
| `contested_no_call` | the authors place a centromere here, we make no call |
| `not_compared` | nothing published on this chromosome, or the published track's assembly differs |

The README lists every contested chromosome (generated from `ALL_official_overlap.tsv`); the karyotype marks
each with a red dot (no overlap) or an amber dot (no call).

### Attributes of our arrays (`satellite_DNA`, source `satellome_cen`)

`family` (superfamily), `membership` (`family` | `cofamily` | `period_class_k<k>` — joined by monomer period
k × base period, not by sequence), `fastan_period`, `mono_period`, `hor_period`, `mono_consensus`,
`largest=true` on the longest array of the locus (a size fact, not a prediction of the active array).

## Tables

- `centromeres.tsv` — one row per chromosome with a call (same fields as the locus, plus the largest array and its monomer consensus).
- `centromeric_arrays.tsv` / `.bed` — every array of every locus, with HOR consensus.
- `satellite_arrays.tsv` — every array ≥ 10 kb with a usable consensus on an assembled chromosome:
  `chromosome`, `seq_name`, `start`, `end`, `length`, `family` (superfamily), `seq_family`, `membership`
  (`family` | `cofamily` | `period_class_k<k>` | `-`), `fastan_period`, `mono_period`, `hor_period`.
- `karyotype.svg` — drawn from the tables above by `plot_karyotype.py`, written as plain SVG so it is
  byte-identical between runs and listed in `MANIFEST.md5`. Panel **a**: all chromosomes to one scale (Mb
  grid). Panel **b**: the centromeric region of each chromosome, a window around our locus (and the
  published intervals that lie within 3 Mb), each to its own scale with its width printed above it; the
  nearest published interval farther away is an arrow with its distance in Mb (↑ before, ↓ after our locus).
  Arrays: blue = centromeric superfamily, teal = co-family, light blue = joined by monomer period, four
  other colours = the four largest other superfamilies, grey = the rest (darker than the chromosome body, so
  an array never reads as empty sequence). Black box (a) / bracket (b) = our locus; magenta bracket = one
  published interval (dashed for a second set; intervals closer than 1.5 px at the panel's scale are drawn
  as one bracket, and the hover text gives their number; genes in a centromeric region are drawn as the span
  first..last gene); dots above panel a, one per published set: green = overlap, red = no overlap, amber =
  published but no call, hollow = nothing published. Red name and dashed outline = no call.
- `official_overlap.tsv` — one row per chromosome × published set (all relations above, with coordinates, and the `agreement` verdict).
- `family_breadth.tsv` — every superfamily; `n_chrom_breadth` = chromosomes with ≥ `--breadth-min-bp` (10 kb) of it,
  the breadth the choice uses; `n_chrom_ge_chrom_min_bp` = chromosomes with ≥ `--chrom-min-bp` (100 kb);
  `base_period` 0 = no consensus ≥ 30 bp (not eligible); `chosen` = 1 centromeric, 2 co-family.
- `outliers.tsv` — chromosomes without an array of the centromeric family; arrays without a usable consensus.
- `exceptions.tsv` — `breadth_near_tie`, `breadth_tie_broken_by_bp`, `short_consensus_family_as_broad`
  (a superfamily without a consensus ≥ 30 bp, at least as broad as the chosen one minus one), `cofamily_joined_locus`, `cap_mode`,
  `subterminal_cap_skipped`, `cap_mode_no_interior_cluster`, `locus_at_chromosome_end`,
  `second_locus_comparable`, `locus_includes_period_class_arrays`, `other_family_larger_on_chromosome`.

## Panel tables

- `STATUS.tsv` — one row per assembly: directory (`<tier>/<Organism>__<assembly>__<accession>`, relative to
  `annotations/`), catalogue path relative to the work directory, `called` /
  `refused` / `pending`, the caller's one-line summary.
- `TIERS.tsv` — completeness of every assembly from its sequence and NCBI assembly report, never from its name
  or NCBI assembly level (`src/t2t_centromeres/tier.py`): `chromosomes` (assembled-molecule rows of type
  Chromosome), `gaps` (runs of N, satellome `gaps.bed`), `gaps_on_chromosomes`, `unplaced_sequences` (report rows
  that are not assembled molecules: unplaced and unlocalized scaffolds), `both_telomeres` / `one_telomere` /
  `no_telomere` (satellome telomere check: TTAGGG / CCCTAA in the terminal 10 kb at density ≥ 0.15),
  `completeness_tier`:

  | tier | rule |
  |---|---|
  | `real_t2t` | no gap, no unplaced sequence, telomeric repeat at both ends of every chromosome |
  | `gap_free` | no gap, no unplaced sequence, telomeric repeat missing at one or more chromosome ends |
  | `gapped` | at least one gap or one unplaced / unlocalized sequence |

  Each assembly directory sits in the folder of its tier, `annotations/<completeness_tier>/`
  (`src/t2t_centromeres/tier_folders.py`, run by `workflow/call_all.sh` after `tier.py`).

- `SUMMARY.tsv` — one row per assembly × published set: `directory` as in `STATUS.tsv`, the call (chromosomes,
  called, family period, near ties), the relation counts and verdicts against the set, and the completeness
  columns of `TIERS.tsv`
  (`completeness_tier`, `gaps`, `unplaced_sequences`, `chromosomes_both_telomeres`).
- `ALL_centromeres.tsv`, `ALL_official_overlap.tsv` — the per-assembly tables concatenated.
