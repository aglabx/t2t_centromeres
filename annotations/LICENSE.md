# Annotation licence

Our calls (features with source `satellome_cen`, the `intersection` regions, the agreement verdicts, the
karyotype figures and all tables) are released under **CC BY 4.0**
(https://creativecommons.org/licenses/by/4.0/).

Features whose GFF source is a publisher are coordinates of published centromere annotations, included so
that every call can be compared with them. They remain the work of their authors; cite them (README,
Sources) and follow their terms:

| GFF source | assemblies | origin | terms |
|---|---|---|---|
| `T2T_CHM13_censat_v2.1` | CHM13 | T2T Consortium, human-pangenomics S3 | CC0, public domain (github.com/marbl/CHM13, "Data reuse and license") |
| `T2T_primates_CenSat_v2.0` | chimpanzee, bonobo, gorilla, orangutans, siamang | T2T Primates release, via github.com/hloucks/CenSatData | CC0, public domain (github.com/marbl/Primates, "Data reuse and license"); the CenSatData repository states no licence of its own |
| `T2T-MFA8_Nature2025_S39` | T2T-MFA8 | Zhang et al. 2025 *Nature*, Supplementary Table 39 | no reuse licence stated; coordinates reproduced for comparison, with citation |
| `T2T-MFA8_release_v1.1` | T2T-MFA8 | release `T2T-MFA8v1.1.centromere.bed` | no reuse licence stated; coordinates reproduced for comparison, with citation |
| `T2T-MMU8_release_v2.0` | T2T-MMU8 | release `T2T-MMU8v2.0.alpha_satellite.bed` (Zhang et al. 2026 *Cell*) | no reuse licence stated; coordinates reproduced for comparison, with citation |
| `T2T-sheep_NatGenet2025_S8` | T2T-sheep1.0 | Luo et al. 2025 *Nature Genetics*, Supplementary Table 8 | no reuse licence stated; coordinates reproduced for comparison, with citation |

Authors who prefer their coordinates not to be redistributed here: open an issue and we will remove them
from the files; `workflow/fetch_official.sh` then downloads them from the original source instead.
