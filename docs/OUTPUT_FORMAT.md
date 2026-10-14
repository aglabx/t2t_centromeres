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
| `<class>.relation` | see below |
| `<class>.overlap_bp` | bp of our locus covered by the published intervals (for genes: by the span first..last gene) |
| `<class>.official_frac_in_locus` | overlap / published bp (for genes: / published span) |
| `<class>.locus_frac_official` | overlap / our locus length |
| `<class>.largest_array_overlap_bp` | overlap of our largest array with the published intervals |
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

### Attributes of our arrays (`satellite_DNA`, source `satellome_cen`)

`family` (superfamily), `membership` (`family` | `cofamily` | `period_class_k<k>` — joined by monomer period
k × base period, not by sequence), `fastan_period`, `mono_period`, `hor_period`, `mono_consensus`,
`largest=true` on the largest array of the locus (the candidate active array).

## Tables

- `centromeres.tsv` — one row per chromosome with a call (same fields as the locus, plus the largest array and its monomer consensus).
- `centromeric_arrays.tsv` / `.bed` — every array of every locus, with HOR consensus.
- `official_overlap.tsv` — one row per chromosome × published set (all relations above, with coordinates).
- `family_breadth.tsv` — every superfamily; `chosen` = 1 centromeric, 2 co-family.
- `outliers.tsv` — chromosomes without an array of the centromeric family; arrays without a usable consensus.
- `exceptions.tsv` — `breadth_near_tie`, `breadth_tie_broken_by_bp`, `cofamily_joined_locus`, `cap_mode`,
  `subterminal_cap_skipped`, `cap_mode_no_interior_cluster`, `locus_at_chromosome_end`,
  `second_locus_comparable`, `locus_includes_period_class_arrays`, `other_family_larger_on_chromosome`.
