# 🧬 panel-baf

BAF-/LOH- und Kopienzahlprofile aus **Hybrid-Capture-Panel-Daten (Illumina
Paired-End, Tumor-only)** – unter Nutzung der **Off-Target-Reads**. Auch bei
kleinem Panel (~150 kb) und nur ~0,1× Off-Target-Tiefe lassen sich so
LOH, Zugewinne, Verluste, Amplifikationen und der Tumorzellgehalt
abschätzen.

```text
📥 FASTQ ──▶ ✂️ fastp ──▶ 🧬 minimap2 (hg19 primär) ──▶ 🧹 markdup
                                   │
          ┌────────────────────────┼─────────────────────────┐
          ▼                        ▼                         ▼
 🔬 Allelzählung an          📈 Off-Target-Tiefe       🎯 On-Target-Tiefe
 ~7 Mio. 1000G-SNPs          (10 kb, GC + Alu)          je Gen
          └──────────▶ 🧮 BAF-Likelihood + HMM ◀──────────┘
                               │
                 📊 Segmente · Calls · Purity · Plots
```

## ✨ Auf einen Blick

- 🔬 **BAF ohne Genotypisierung:** Bei ~0,1× lassen sich SNP-Genotypen nicht
  einzeln bestimmen. Jeder SNP mit ≥ 2 Reads liefert eine Likelihood über den
  Minor-Allel-Anteil *b* heterozygoter SNPs (Hardy-Weinberg-Prior). Summiert
  über Bins/Segmente ergibt das MLE und 95 %-KI für *b*.
- 🌍 **Ancestry-Kalibrierung:** Die 1000G-Superpopulation (EUR, AFR, AMR,
  EAS, SAS oder global) wird per Likelihood gewählt, die Heterozygotie-Rate
  gefittet. Mit globalen AFs entsteht sonst genomweit scheinbare Imbalance.
- 📈 **Off-Target-logR:** Capture-Hotspots werden verworfen; Korrektur auf
  GC und **Alu-Dichte** (Off-Target-Reads sind stark Alu-getrieben).
- 🧮 **Segmentierung:** gemeinsames HMM über (logR, *b*) pro Chromosomenarm,
  dann ASCAT-artiger Fit von Tumorzellgehalt und Allelzuständen (z. B. 2+1, 1+0).
- 🏷️ **Calls:** Amplification, Gain (imbalanced/balanced), Loss/LOH, CN-LOH,
  Possible imbalance – jeweils mit Konfidenz.
- ♻️ **Wiederaufnahme:** Fertige Schritte werden erkannt (sicher auf `requeue`).
- 🖥️ **SLURM:** Auf PALMA ein Job pro Probe; lokal seriell.

## 🚀 Schnellstart (PALMA)

```bash
git clone https://github.com/ctho1/panel-baf.git
cd panel-baf

# 1. Conda-Env einmalig anlegen (minimap2/samtools kommen aus den PALMA-Modulen)
bash scripts/setup_envs.sh

# 2. Panel-BED (hg19, "chr"-Namen) bereitstellen
mkdir -p references
cp /pfad/zur/NPMS_v4_targets_hg19.bed references/panel_targets.bed

# 3. FASTQs nach input/ kopieren und starten
cp /pfad/zu/*_R1_001.fastq.gz /pfad/zu/*_R2_001.fastq.gz input/
bash run_pipeline.sh --mail christian.thomas@ukmuenster.de
```

Beim ersten Start fehlen die Referenzen → `prepare_reference.sh` läuft als
vorgeschalteter Job (~30–60 min, v. a. Downloads: hg19-Primärchromosomen +
minimap2-Index, 1000G-SNP-Panel mit 7,1 Mio. SNPs, Annotation; ca. 12 GB). Alle Probenjobs warten per
`--dependency` darauf.

Optional:

```bash
bash run_pipeline.sh --partition normal                 # Partition fest vorgeben
PANEL_BED=/cloud/.../panel.bed bash run_pipeline.sh      # BED an anderem Ort
FOCAL_GENES="EGFR;CDKN2A,CDKN2B;MDM2;CDK4" bash run_pipeline.sh
```

Alle Pfade, Module und Parameter stehen in
[`scripts/config.sh`](scripts/config.sh) und lassen sich als
Umgebungsvariablen überschreiben.

## 📥 Eingaben

| Datei | Hinweis |
| --- | --- |
| `input/<Probe>_R1_001.fastq.gz` + `_R2_001` | Illumina Paired-End, flach in `input/` |
| `references/panel_targets.bed` | Panel-Targets (hg19); **nicht im Repo** |

## 📤 Ergebnisse

Pro Probe unter `output/<Probe>/`:

| Datei | Inhalt |
| --- | --- |
| `<Probe>.genome_BAF_logR.png` | Genomweites BAF-/logR-Profil mit Calls |
| `<Probe>.focal.png` | Zoom auf `FOCAL_GENES` (Standard: EGFR, CDKN2A/B, PTEN, PDGFRA/KIT) |
| `chromosomes/<Probe>.chrN.png` | Einzelchromosomen |
| `<Probe>.segments.tsv` | Segmente: logR, *b* (MLE + KI), p-Imbalance, Call, Konfidenz, Allelzustand |
| `<Probe>.arms.tsv` | dasselbe pro Chromosomenarm |
| `<Probe>.gene_coverage.tsv` | On-Target-Tiefe (GC-korr.) + het. SNPs je Gen |
| `<Probe>.purity.json` | Tumorzellgehalt, logR-Kompression γ |
| `<Probe>.summary.json` | QC (Reads, Duplikate, On-Target-Anteil, Tiefe), Geschlecht, Ancestry, Events |
| `qc/` | fastp-Report |

BAM, Allelzählungen und Bin-Tabellen bleiben unter `tmp/<Probe>/`.

## ⚙️ Ausführung auf PALMA

Wie in `mNGS-pipeline`: `sinfo`-basierte Wahl zwischen **`normal`** und
**`requeue`**, Module `palma/2024a` (minimap2 2.29, SAMtools 1.21), alles
Weitere (fastp, bcftools, mosdepth, Python-Stack) aus der Conda-Env.
Ressourcen pro Probe: 16 CPUs, 48 GB, ≤ 2 h (typisch ~15 min für 20 Mio.
Read-Paare). Liegt `/cloud/wwu1/e_np_ngs/references/genomes/hg19.fa` vor,
wird sie statt eines Downloads genutzt.

## 🧪 Tests

```bash
conda_envs/baf/bin/python -m unittest discover -s tests -v
```

Getestet werden u. a. die BAF-Likelihood an simulierten Low-Depth-SNPs
(balanciert, Verlust bei 60 % Purity, LOH), die Hotspot-robuste
Ratenschätzung und die HMM-Segmentierung.

## ⚠️ Einschränkungen

- Kein Normal-Referenzpool: logR-Amplituden sind komprimiert (γ ≈ 0,7), und
  Zugewinne Alu-reicher Chromosomen (v. a. chr19) werden teilweise
  wegkorrigiert – dort ist die BAF die verlässlichere Evidenz.
- Off-Target-BAF erfasst LOH/Imbalance ab etwa 10 Mb; subklonale oder kleine
  Veränderungen nicht sicher. p-Werte ignorieren LD → strenge Schwellen.
- Einzelne 3–5-Mb-Fenster mit *b* ≈ 0 in sonst neutralen Regionen können
  Keimbahn-ROH sein, kein Tumor-LOH.
- On-Target-Tiefe ohne Referenz: nur große Abweichungen (Amplifikation,
  homozygote Deletion) interpretierbar.

> [!IMPORTANT]
> `input/`, `output/`, `tmp/`, `log/` und `references/` gehören nicht ins
> Repository. Ergebnisse vor diagnostischer Verwendung immer fachlich prüfen.
