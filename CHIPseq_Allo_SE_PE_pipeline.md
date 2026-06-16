# ChIP-seq Pipeline with Allo Multi-Mapping Rescue

Reusable GitHub-ready pipeline for **single-end (SE)** and **paired-end (PE)** ChIP-seq data using `SRA-tools`, `fastp`, `bowtie2`, `samtools`, `Allo`, and `MACS2`.

Project example used here:

```bash
/media/ismam/Expansion/RNA_WGBS_ChIP/D_pulex/ChIP
```

Example input/control sample:

```bash
SRR8375037
```

`SRR8375037` is a **single-end input/control ChIP-seq sample**.

---

## 1. Directory structure

Run from the ChIP project folder:

```bash
cd /media/ismam/Expansion/RNA_WGBS_ChIP/D_pulex/ChIP

mkdir -p input/raw/sra
mkdir -p input/raw
mkdir -p input/trimmed
mkdir -p input/ref
mkdir -p input/ref/bowtie2_index
mkdir -p input/aligned
mkdir -p input/allo_out
mkdir -p input/logs

mkdir -p chip/raw/sra
mkdir -p chip/raw
mkdir -p chip/trimmed
mkdir -p chip/aligned
mkdir -p chip/allo_out
mkdir -p chip/logs

mkdir -p macs2_out
```

Recommended meaning:

```text
input/ = input/control sample
chip/  = treatment ChIP sample
macs2_out/ = peak-calling output
```

---

## 2. Software environment

Create the conda environment if needed:

```bash
conda create -n chipseqTools -c conda-forge -c bioconda \
  sra-tools fastp bowtie2 samtools macs2 bedtools deeptools pip -y
```

Activate it:

```bash
conda activate chipseqTools
```

Check tools:

```bash
fastp --version
bowtie2 --version
samtools --version
macs2 --version
```

---

## 3. Install Allo

If Allo is already downloaded at `/home/ismam/allo`, install it into the active conda environment:

```bash
conda activate chipseqTools
cd /home/ismam/allo
pip install -e .
```

Test:

```bash
which allo
allo --help
```

If `allo --help` works, Allo is ready.

---

## 4. Reference genome and Bowtie2 index

Return to the project folder:

```bash
cd /media/ismam/Expansion/RNA_WGBS_ChIP/D_pulex/ChIP
```

Set variables:

```bash
GENOME=input/ref/Daphnia_pulex_gca021134715v1rs.ASM2113471v1.dna.toplevel.fa
BOWTIE2_INDEX=input/ref/bowtie2_index/Daphnia_pulex
THREADS=16
```

Build Bowtie2 index:

```bash
bowtie2-build --threads ${THREADS} ${GENOME} ${BOWTIE2_INDEX}
```

Check:

```bash
ls -lh input/ref/bowtie2_index/
```

Expected index files:

```text
Daphnia_pulex.1.bt2
Daphnia_pulex.2.bt2
Daphnia_pulex.3.bt2
Daphnia_pulex.4.bt2
Daphnia_pulex.rev.1.bt2
Daphnia_pulex.rev.2.bt2
```

---

# PART A. Single-end input/control ChIP-seq pipeline

Use this for SRA layout **SINGLE**.

Example:

```bash
SAMPLE=SRR8375037
THREADS=16
GENOME=input/ref/Daphnia_pulex_gca021134715v1rs.ASM2113471v1.dna.toplevel.fa
BOWTIE2_INDEX=input/ref/bowtie2_index/Daphnia_pulex
```

---

## A1. Download single-end data

```bash
cd /media/ismam/Expansion/RNA_WGBS_ChIP/D_pulex/ChIP

prefetch ${SAMPLE} \
  --output-directory input/raw/sra
```

Convert to FASTQ:

```bash
fasterq-dump input/raw/sra/${SAMPLE} \
  --outdir input/raw \
  --threads ${THREADS}
```

Compress:

```bash
gzip -f input/raw/${SAMPLE}.fastq
```

Check:

```bash
ls -lh input/raw/${SAMPLE}.fastq.gz
```

---

## A2. Trim single-end reads with fastp

```bash
fastp \
  -i input/raw/${SAMPLE}.fastq.gz \
  -o input/trimmed/${SAMPLE}.trimmed.fastq.gz \
  -w ${THREADS} \
  -h input/trimmed/${SAMPLE}_fastp.html \
  -j input/trimmed/${SAMPLE}_fastp.json
```

Check:

```bash
ls -lh input/trimmed/${SAMPLE}.trimmed.fastq.gz
ls -lh input/trimmed/${SAMPLE}_fastp.html
```

---

## A3. Align single-end reads with Bowtie2 for Allo

Allo needs multi-mapping locations, so use `-k 25`.

```bash
bowtie2 \
  -x ${BOWTIE2_INDEX} \
  -q input/trimmed/${SAMPLE}.trimmed.fastq.gz \
  -S input/aligned/${SAMPLE}_MMR.sam \
  -k 25 \
  -p ${THREADS} \
  2> input/logs/${SAMPLE}_bowtie2_MMR.log
```

Check:

```bash
cat input/logs/${SAMPLE}_bowtie2_MMR.log
ls -lh input/aligned/${SAMPLE}_MMR.sam
```

Example successful alignment from `SRR8375037`:

```text
76.67% overall alignment rate
```

---

## A4. Collate SAM by read name

Allo expects alignments grouped by read name.

```bash
samtools collate \
  -@ ${THREADS} \
  -o input/aligned/${SAMPLE}_MMR_collate.sam \
  input/aligned/${SAMPLE}_MMR.sam
```

Check:

```bash
ls -lh input/aligned/${SAMPLE}_MMR_collate.sam
```

---

## A5. Run Allo for single-end input/control

For input/control/background, use `--random`.

```bash
mkdir -p input/allo_out

allo input/aligned/${SAMPLE}_MMR_collate.sam \
  -seq se \
  --random \
  -o input/allo_out/${SAMPLE}_allo_random.sam \
  -p ${THREADS}
```

Check output:

```bash
ls -lh input/allo_out/${SAMPLE}_allo_random.sam
```

Check Allo assignment tags:

```bash
grep -c "ZA:" input/allo_out/${SAMPLE}_allo_random.sam
grep -c "ZZ:" input/allo_out/${SAMPLE}_allo_random.sam
```

Interpretation:

```text
ZA = multi-mapped reads assigned by Allo using unique-read evidence
ZZ = reads randomly assigned because no unique-read evidence was available
```

Example result from `SRR8375037`:

```text
ZA = 5,507,923
ZZ =   774,715
```

---

## A6. Convert Allo SAM to sorted BAM

```bash
samtools view -@ ${THREADS} -bS \
  input/allo_out/${SAMPLE}_allo_random.sam \
  > input/allo_out/${SAMPLE}_allo_random.bam

samtools sort -@ ${THREADS} \
  -o input/allo_out/${SAMPLE}_allo_random.sorted.bam \
  input/allo_out/${SAMPLE}_allo_random.bam

samtools index input/allo_out/${SAMPLE}_allo_random.sorted.bam
```

Check:

```bash
ls -lh input/allo_out/${SAMPLE}_allo_random.sorted.bam*
samtools flagstat input/allo_out/${SAMPLE}_allo_random.sorted.bam
```

---

## A7. Make primary-only BAM for MACS2

Allo output may contain secondary alignments. Before MACS2, remove secondary alignments with `-F 256`.

```bash
samtools view -@ ${THREADS} -b -F 256 \
  input/allo_out/${SAMPLE}_allo_random.sorted.bam \
  > input/allo_out/${SAMPLE}_allo_random.primary.bam

samtools sort -@ ${THREADS} \
  -o input/allo_out/${SAMPLE}_allo_random.primary.sorted.bam \
  input/allo_out/${SAMPLE}_allo_random.primary.bam

samtools index input/allo_out/${SAMPLE}_allo_random.primary.sorted.bam
```

Check:

```bash
ls -lh input/allo_out/${SAMPLE}_allo_random.primary.sorted.bam*
samtools flagstat input/allo_out/${SAMPLE}_allo_random.primary.sorted.bam
```

Remove intermediate BAMs:

```bash
rm input/allo_out/${SAMPLE}_allo_random.bam
rm input/allo_out/${SAMPLE}_allo_random.primary.bam
```

Final single-end input/control BAM:

```bash
input/allo_out/${SAMPLE}_allo_random.primary.sorted.bam
```

For `SRR8375037`, final control BAM:

```bash
input/allo_out/SRR8375037_allo_random.primary.sorted.bam
```

---

# PART B. Paired-end ChIP-seq pipeline

Use this for SRA layout **PAIRED**.

Set variables:

```bash
SAMPLE=SRRXXXXXXX
THREADS=16
GENOME=input/ref/Daphnia_pulex_gca021134715v1rs.ASM2113471v1.dna.toplevel.fa
BOWTIE2_INDEX=input/ref/bowtie2_index/Daphnia_pulex
```

---

## B1. Download paired-end data

```bash
cd /media/ismam/Expansion/RNA_WGBS_ChIP/D_pulex/ChIP

prefetch ${SAMPLE} \
  --output-directory chip/raw/sra
```

Convert to paired FASTQ:

```bash
fasterq-dump chip/raw/sra/${SAMPLE} \
  --split-files \
  --outdir chip/raw \
  --threads ${THREADS}
```

Compress:

```bash
gzip -f chip/raw/${SAMPLE}_1.fastq
gzip -f chip/raw/${SAMPLE}_2.fastq
```

Check:

```bash
ls -lh chip/raw/${SAMPLE}_1.fastq.gz
ls -lh chip/raw/${SAMPLE}_2.fastq.gz
```

---

## B2. Trim paired-end reads with fastp

```bash
fastp \
  -i chip/raw/${SAMPLE}_1.fastq.gz \
  -I chip/raw/${SAMPLE}_2.fastq.gz \
  -o chip/trimmed/${SAMPLE}_1.trimmed.fastq.gz \
  -O chip/trimmed/${SAMPLE}_2.trimmed.fastq.gz \
  --detect_adapter_for_pe \
  -w ${THREADS} \
  -h chip/trimmed/${SAMPLE}_fastp.html \
  -j chip/trimmed/${SAMPLE}_fastp.json
```

Check:

```bash
ls -lh chip/trimmed/${SAMPLE}_1.trimmed.fastq.gz
ls -lh chip/trimmed/${SAMPLE}_2.trimmed.fastq.gz
```

---

## B3. Align paired-end reads with Bowtie2 for Allo

For paired-end data, use `-1` and `-2`. Keep multiple alignments with `-k 25`.

```bash
bowtie2 \
  -x ${BOWTIE2_INDEX} \
  -1 chip/trimmed/${SAMPLE}_1.trimmed.fastq.gz \
  -2 chip/trimmed/${SAMPLE}_2.trimmed.fastq.gz \
  -S chip/aligned/${SAMPLE}_MMR.sam \
  -k 25 \
  --no-mixed \
  --no-discordant \
  -p ${THREADS} \
  2> chip/logs/${SAMPLE}_bowtie2_MMR.log
```

Check:

```bash
cat chip/logs/${SAMPLE}_bowtie2_MMR.log
ls -lh chip/aligned/${SAMPLE}_MMR.sam
```

---

## B4. Collate paired-end SAM by read name

```bash
samtools collate \
  -@ ${THREADS} \
  -o chip/aligned/${SAMPLE}_MMR_collate.sam \
  chip/aligned/${SAMPLE}_MMR.sam
```

Check:

```bash
ls -lh chip/aligned/${SAMPLE}_MMR_collate.sam
```

---

## B5. Run Allo for paired-end treatment sample

For treatment ChIP samples:

```bash
mkdir -p chip/allo_out

allo chip/aligned/${SAMPLE}_MMR_collate.sam \
  -seq pe \
  -o chip/allo_out/${SAMPLE}_allo.sam \
  -p ${THREADS}
```

For paired-end input/control samples, use `--random`:

```bash
allo chip/aligned/${SAMPLE}_MMR_collate.sam \
  -seq pe \
  --random \
  -o chip/allo_out/${SAMPLE}_allo_random.sam \
  -p ${THREADS}
```

Check tags:

```bash
grep -c "ZA:" chip/allo_out/${SAMPLE}_allo.sam
grep -c "ZZ:" chip/allo_out/${SAMPLE}_allo.sam
```

If you used the random input/control version:

```bash
grep -c "ZA:" chip/allo_out/${SAMPLE}_allo_random.sam
grep -c "ZZ:" chip/allo_out/${SAMPLE}_allo_random.sam
```

---

## B6. Convert paired-end Allo SAM to sorted BAM

For treatment sample without `--random`:

```bash
samtools view -@ ${THREADS} -bS \
  chip/allo_out/${SAMPLE}_allo.sam \
  > chip/allo_out/${SAMPLE}_allo.bam

samtools sort -@ ${THREADS} \
  -o chip/allo_out/${SAMPLE}_allo.sorted.bam \
  chip/allo_out/${SAMPLE}_allo.bam

samtools index chip/allo_out/${SAMPLE}_allo.sorted.bam
```

Make primary-only BAM:

```bash
samtools view -@ ${THREADS} -b -F 256 \
  chip/allo_out/${SAMPLE}_allo.sorted.bam \
  > chip/allo_out/${SAMPLE}_allo.primary.bam

samtools sort -@ ${THREADS} \
  -o chip/allo_out/${SAMPLE}_allo.primary.sorted.bam \
  chip/allo_out/${SAMPLE}_allo.primary.bam

samtools index chip/allo_out/${SAMPLE}_allo.primary.sorted.bam
```

Check:

```bash
ls -lh chip/allo_out/${SAMPLE}_allo.primary.sorted.bam*
samtools flagstat chip/allo_out/${SAMPLE}_allo.primary.sorted.bam
```

Remove intermediate BAMs:

```bash
rm chip/allo_out/${SAMPLE}_allo.bam
rm chip/allo_out/${SAMPLE}_allo.primary.bam
```

Final paired-end treatment BAM:

```bash
chip/allo_out/${SAMPLE}_allo.primary.sorted.bam
```

---

# PART C. MACS2 peak calling

## C1. Define treatment and control BAMs

Example:

```bash
TREAT_BAM=chip/allo_out/SRR_TREATMENT_allo.primary.sorted.bam
CONTROL_BAM=input/allo_out/SRR8375037_allo_random.primary.sorted.bam
OUTDIR=macs2_out
NAME=Dpulex_ChIP_vs_Input_Allo
GENOME_SIZE=1.8e8
THREADS=16
```

`GENOME_SIZE=1.8e8` is a starting value for *Daphnia pulex*. Adjust it if you calculate a more specific effective genome size.

---

## C2. MACS2 for single-end treatment BAM

Use this if the treatment ChIP sample is single-end:

```bash
macs2 callpeak \
  -t ${TREAT_BAM} \
  -c ${CONTROL_BAM} \
  -f BAM \
  -g ${GENOME_SIZE} \
  -n ${NAME} \
  --outdir ${OUTDIR} \
  -q 0.05 \
  --nomodel \
  --extsize 147
```

Notes:

- `-f BAM` = single-end BAM.
- `--nomodel --extsize 147` is useful when MACS2 cannot build a fragment model or when you want fixed read extension.

---

## C3. MACS2 for paired-end treatment BAM

Use this if the treatment ChIP sample is paired-end:

```bash
macs2 callpeak \
  -t ${TREAT_BAM} \
  -c ${CONTROL_BAM} \
  -f BAMPE \
  -g ${GENOME_SIZE} \
  -n ${NAME} \
  --outdir ${OUTDIR} \
  -q 0.05
```

Notes:

- `-f BAMPE` = paired-end BAM.
- For paired-end BAMPE, MACS2 uses observed fragment sizes, so `--nomodel --extsize` is usually not needed.

---

## C4. Expected MACS2 output

```bash
ls -lh macs2_out/
```

Expected files:

```text
${NAME}_peaks.narrowPeak
${NAME}_peaks.xls
${NAME}_summits.bed
${NAME}_treat_pileup.bdg
${NAME}_control_lambda.bdg
```

---

# PART D. Peak QC and filtering

## D1. Count peaks

```bash
wc -l macs2_out/${NAME}_peaks.narrowPeak
```

## D2. Sort peaks

```bash
sort -k1,1 -k2,2n \
  macs2_out/${NAME}_peaks.narrowPeak \
  > macs2_out/${NAME}_peaks.sorted.narrowPeak
```

## D3. Filter strong peaks

In narrowPeak, column 9 is `-log10(q-value)`. This keeps peaks with approximately q-value <= 0.01:

```bash
awk '$9 >= 2' \
  macs2_out/${NAME}_peaks.sorted.narrowPeak \
  > macs2_out/${NAME}_peaks.q0.01.narrowPeak
```

Count filtered peaks:

```bash
wc -l macs2_out/${NAME}_peaks.q0.01.narrowPeak
```

---

# PART E. Optional bigWig generation

## E1. Treatment bigWig

```bash
bamCoverage \
  -b ${TREAT_BAM} \
  -o macs2_out/${NAME}_treatment_RPGC.bw \
  --normalizeUsing RPGC \
  --effectiveGenomeSize 180000000 \
  --binSize 10 \
  -p ${THREADS}
```

## E2. Control bigWig

```bash
bamCoverage \
  -b ${CONTROL_BAM} \
  -o macs2_out/${NAME}_control_RPGC.bw \
  --normalizeUsing RPGC \
  --effectiveGenomeSize 180000000 \
  --binSize 10 \
  -p ${THREADS}
```

---

# PART F. Complete single-end input/control example for SRR8375037

```bash
conda activate chipseqTools

cd /media/ismam/Expansion/RNA_WGBS_ChIP/D_pulex/ChIP

SAMPLE=SRR8375037
THREADS=16
GENOME=input/ref/Daphnia_pulex_gca021134715v1rs.ASM2113471v1.dna.toplevel.fa
BOWTIE2_INDEX=input/ref/bowtie2_index/Daphnia_pulex

mkdir -p input/raw/sra input/raw input/trimmed input/ref/bowtie2_index input/aligned input/allo_out input/logs

prefetch ${SAMPLE} --output-directory input/raw/sra

fasterq-dump input/raw/sra/${SAMPLE} \
  --outdir input/raw \
  --threads ${THREADS}

gzip -f input/raw/${SAMPLE}.fastq

fastp \
  -i input/raw/${SAMPLE}.fastq.gz \
  -o input/trimmed/${SAMPLE}.trimmed.fastq.gz \
  -w ${THREADS} \
  -h input/trimmed/${SAMPLE}_fastp.html \
  -j input/trimmed/${SAMPLE}_fastp.json

bowtie2 \
  -x ${BOWTIE2_INDEX} \
  -q input/trimmed/${SAMPLE}.trimmed.fastq.gz \
  -S input/aligned/${SAMPLE}_MMR.sam \
  -k 25 \
  -p ${THREADS} \
  2> input/logs/${SAMPLE}_bowtie2_MMR.log

samtools collate \
  -@ ${THREADS} \
  -o input/aligned/${SAMPLE}_MMR_collate.sam \
  input/aligned/${SAMPLE}_MMR.sam

allo input/aligned/${SAMPLE}_MMR_collate.sam \
  -seq se \
  --random \
  -o input/allo_out/${SAMPLE}_allo_random.sam \
  -p ${THREADS}

grep -c "ZA:" input/allo_out/${SAMPLE}_allo_random.sam
grep -c "ZZ:" input/allo_out/${SAMPLE}_allo_random.sam

samtools view -@ ${THREADS} -bS \
  input/allo_out/${SAMPLE}_allo_random.sam \
  > input/allo_out/${SAMPLE}_allo_random.bam

samtools sort -@ ${THREADS} \
  -o input/allo_out/${SAMPLE}_allo_random.sorted.bam \
  input/allo_out/${SAMPLE}_allo_random.bam

samtools index input/allo_out/${SAMPLE}_allo_random.sorted.bam

samtools view -@ ${THREADS} -b -F 256 \
  input/allo_out/${SAMPLE}_allo_random.sorted.bam \
  > input/allo_out/${SAMPLE}_allo_random.primary.bam

samtools sort -@ ${THREADS} \
  -o input/allo_out/${SAMPLE}_allo_random.primary.sorted.bam \
  input/allo_out/${SAMPLE}_allo_random.primary.bam

samtools index input/allo_out/${SAMPLE}_allo_random.primary.sorted.bam

samtools flagstat input/allo_out/${SAMPLE}_allo_random.primary.sorted.bam

rm input/allo_out/${SAMPLE}_allo_random.bam
rm input/allo_out/${SAMPLE}_allo_random.primary.bam
```

Final file:

```bash
input/allo_out/SRR8375037_allo_random.primary.sorted.bam
```

---

# PART G. Complete paired-end treatment template

Replace `SRRXXXXXXX` with your treatment ChIP accession.

```bash
conda activate chipseqTools

cd /media/ismam/Expansion/RNA_WGBS_ChIP/D_pulex/ChIP

SAMPLE=SRRXXXXXXX
THREADS=16
GENOME=input/ref/Daphnia_pulex_gca021134715v1rs.ASM2113471v1.dna.toplevel.fa
BOWTIE2_INDEX=input/ref/bowtie2_index/Daphnia_pulex

mkdir -p chip/raw/sra chip/raw chip/trimmed chip/aligned chip/allo_out chip/logs

prefetch ${SAMPLE} --output-directory chip/raw/sra

fasterq-dump chip/raw/sra/${SAMPLE} \
  --split-files \
  --outdir chip/raw \
  --threads ${THREADS}

gzip -f chip/raw/${SAMPLE}_1.fastq
gzip -f chip/raw/${SAMPLE}_2.fastq

fastp \
  -i chip/raw/${SAMPLE}_1.fastq.gz \
  -I chip/raw/${SAMPLE}_2.fastq.gz \
  -o chip/trimmed/${SAMPLE}_1.trimmed.fastq.gz \
  -O chip/trimmed/${SAMPLE}_2.trimmed.fastq.gz \
  --detect_adapter_for_pe \
  -w ${THREADS} \
  -h chip/trimmed/${SAMPLE}_fastp.html \
  -j chip/trimmed/${SAMPLE}_fastp.json

bowtie2 \
  -x ${BOWTIE2_INDEX} \
  -1 chip/trimmed/${SAMPLE}_1.trimmed.fastq.gz \
  -2 chip/trimmed/${SAMPLE}_2.trimmed.fastq.gz \
  -S chip/aligned/${SAMPLE}_MMR.sam \
  -k 25 \
  --no-mixed \
  --no-discordant \
  -p ${THREADS} \
  2> chip/logs/${SAMPLE}_bowtie2_MMR.log

samtools collate \
  -@ ${THREADS} \
  -o chip/aligned/${SAMPLE}_MMR_collate.sam \
  chip/aligned/${SAMPLE}_MMR.sam

allo chip/aligned/${SAMPLE}_MMR_collate.sam \
  -seq pe \
  -o chip/allo_out/${SAMPLE}_allo.sam \
  -p ${THREADS}

grep -c "ZA:" chip/allo_out/${SAMPLE}_allo.sam
grep -c "ZZ:" chip/allo_out/${SAMPLE}_allo.sam

samtools view -@ ${THREADS} -bS \
  chip/allo_out/${SAMPLE}_allo.sam \
  > chip/allo_out/${SAMPLE}_allo.bam

samtools sort -@ ${THREADS} \
  -o chip/allo_out/${SAMPLE}_allo.sorted.bam \
  chip/allo_out/${SAMPLE}_allo.bam

samtools index chip/allo_out/${SAMPLE}_allo.sorted.bam

samtools view -@ ${THREADS} -b -F 256 \
  chip/allo_out/${SAMPLE}_allo.sorted.bam \
  > chip/allo_out/${SAMPLE}_allo.primary.bam

samtools sort -@ ${THREADS} \
  -o chip/allo_out/${SAMPLE}_allo.primary.sorted.bam \
  chip/allo_out/${SAMPLE}_allo.primary.bam

samtools index chip/allo_out/${SAMPLE}_allo.primary.sorted.bam

samtools flagstat chip/allo_out/${SAMPLE}_allo.primary.sorted.bam

rm chip/allo_out/${SAMPLE}_allo.bam
rm chip/allo_out/${SAMPLE}_allo.primary.bam
```

Final paired-end treatment file:

```bash
chip/allo_out/SRRXXXXXXX_allo.primary.sorted.bam
```

---

# PART H. Single-end versus paired-end details

## Single-end raw sequence

Single-end data has one FASTQ file:

```bash
SAMPLE.fastq.gz
```

Main differences:

```bash
fasterq-dump SAMPLE --outdir raw --threads 16
fastp -i raw/SAMPLE.fastq.gz -o trimmed/SAMPLE.trimmed.fastq.gz
bowtie2 -q trimmed/SAMPLE.trimmed.fastq.gz
allo -seq se
macs2 callpeak -f BAM
```

## Paired-end raw sequence

Paired-end data has two FASTQ files:

```bash
SAMPLE_1.fastq.gz
SAMPLE_2.fastq.gz
```

Main differences:

```bash
fasterq-dump SAMPLE --split-files --outdir raw --threads 16
fastp -i raw/SAMPLE_1.fastq.gz -I raw/SAMPLE_2.fastq.gz \
      -o trimmed/SAMPLE_1.trimmed.fastq.gz -O trimmed/SAMPLE_2.trimmed.fastq.gz
bowtie2 -1 trimmed/SAMPLE_1.trimmed.fastq.gz -2 trimmed/SAMPLE_2.trimmed.fastq.gz
allo -seq pe
macs2 callpeak -f BAMPE
```

---

# PART I. Why each major step is needed

## Why use `-k 25` in Bowtie2?

Allo rescues multi-mapped reads. Therefore, Bowtie2 must report multiple possible alignments per read:

```bash
-k 25
```

Without `-k 25`, Allo cannot evaluate many multi-mapping locations.

## Why use `samtools collate`?

Allo needs all alignments for each read grouped together. `samtools collate` groups alignments by read name:

```bash
samtools collate -o sample_collate.sam sample.sam
```

## Why use `--random` for input/control?

For background/input samples, `--random` keeps reads even when Allo cannot find unique-read support in candidate regions. These reads get a `ZZ` tag.

## Why make primary-only BAM?

MACS2 should use primary alignments. Allo/Bowtie2 output can contain secondary alignments, so remove them before MACS2:

```bash
samtools view -b -F 256 input.bam > primary.bam
```

---

# PART J. Methods text for manuscript or GitHub README

Raw ChIP-seq reads were downloaded from NCBI SRA using SRA-tools and quality-trimmed using fastp. Trimmed reads were aligned to the *Daphnia pulex* reference genome using Bowtie2 with multi-mapping retention enabled (`-k 25`). Alignments were collated by read name with samtools and processed with Allo to assign multi-mapped reads using uniquely mapped read support. For input/control samples, Allo was run with the `--random` option to retain reads in regions lacking unique-read evidence. Allo-processed alignments were converted to sorted BAM files, indexed with samtools, and filtered to retain primary alignments before peak calling. Peaks were called with MACS2 using the Allo-processed ChIP BAM as treatment and the Allo-processed input BAM as control.

---

# PART K. Troubleshooting

## `allo: command not found`

```bash
conda activate chipseqTools
cd /home/ismam/allo
pip install -e .
allo --help
```

## Bowtie2 index path error

If current folder is:

```bash
/media/ismam/Expansion/RNA_WGBS_ChIP/D_pulex/ChIP
```

use:

```bash
input/ref/bowtie2_index/Daphnia_pulex
```

If current folder is:

```bash
/media/ismam/Expansion/RNA_WGBS_ChIP/D_pulex/ChIP/input
```

use:

```bash
ref/bowtie2_index/Daphnia_pulex
```

## MACS2 format error

Use:

```bash
-f BAM      # single-end
-f BAMPE    # paired-end
```

## Disk space cleanup

Large files:

```bash
*_MMR.sam
*_MMR_collate.sam
*_allo.sam
```

Remove these only after final BAMs, BAM indexes, and MACS2 outputs are verified.

Do **not** delete:

```bash
*.primary.sorted.bam
*.primary.sorted.bam.bai
*_peaks.narrowPeak
*_summits.bed
```

---

# End of pipeline
