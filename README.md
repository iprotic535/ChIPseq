# ChIP-seq Pipeline with Allo Multi-Mapping Rescue and MACS3 Peak Calling

Reusable GitHub-ready pipeline for **single-end (SE)** and **paired-end (PE)** ChIP-seq data.

This version uses **MACS3** instead of MACS2 for peak calling.

Pipeline tools:

- `sra-tools`
- `fastp`
- `bowtie2`
- `samtools`
- `Allo`
- `MACS3`
- optional: `bedtools`, `deepTools`

Project example:

```bash
/media/ismam/Expansion/RNA_WGBS_ChIP/D_pulex/ChIP
```

Example input/control sample:

```bash
SRR8375037
```

---

## Why MACS3 instead of MACS2?

MACS3 is the current Python 3 version of the MACS peak-calling framework. The main peak-calling command is still:

```bash
macs3 callpeak
```

The logic is very similar to older `macs2 callpeak` workflows, but the executable is `macs3`.

Use MACS3 here because:

1. MACS3 is actively maintained for modern Python environments.
2. It supports standard ChIP-seq peak calling with treatment and control BAMs.
3. It supports standard input formats including `BAM` and `BAMPE`.
4. It produces standard peak files such as `narrowPeak`, `broadPeak`, `summits.bed`, bedGraph pileups, and `.xls` text summaries.
5. It is easier to install in many current conda/pip environments than older MACS2.

---

# 0. Directory structure

Start in the ChIP-seq project directory:

```bash
cd /media/ismam/Expansion/RNA_WGBS_ChIP/D_pulex/ChIP
```

Create folders:

```bash
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

mkdir -p macs3_out
mkdir -p bigwig
```

Suggested usage:

- `input/` = input/control/background sample.
- `chip/` = treatment ChIP sample.
- `macs3_out/` = MACS3 peak-calling output.
- `bigwig/` = normalized coverage tracks.

---

# 1. Install software

## 1.1 Create conda environment

```bash
conda create -n chipseqTools -c conda-forge -c bioconda \
  sra-tools fastp bowtie2 samtools macs3 bedtools deeptools pip -y
```

Activate:

```bash
conda activate chipseqTools
```

Check versions:

```bash
fasterq-dump --version
fastp --version
bowtie2 --version
samtools --version
macs3 --version
bamCoverage --version
```

If you already have the environment:

```bash
conda activate chipseqTools
conda install -c conda-forge -c bioconda macs3 -y
```

---

# 2. Install Allo

If Allo is already downloaded here:

```bash
/home/ismam/allo
```

install it inside the active conda environment:

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

If `allo --help` works, continue.

---

# 3. Reference genome setup

Go to the project folder:

```bash
cd /media/ismam/Expansion/RNA_WGBS_ChIP/D_pulex/ChIP
```

Set reference variables:

```bash
THREADS=16
GENOME=input/ref/Daphnia_pulex_gca021134715v1rs.ASM2113471v1.dna.toplevel.fa
BOWTIE2_INDEX=input/ref/bowtie2_index/Daphnia_pulex
```

Build Bowtie2 index:

```bash
bowtie2-build --threads ${THREADS} ${GENOME} ${BOWTIE2_INDEX}
```

Check index:

```bash
ls -lh input/ref/bowtie2_index/
```

Expected files:

```bash
Daphnia_pulex.1.bt2
Daphnia_pulex.2.bt2
Daphnia_pulex.3.bt2
Daphnia_pulex.4.bt2
Daphnia_pulex.rev.1.bt2
Daphnia_pulex.rev.2.bt2
```

---

# PART A: Single-end input/control ChIP-seq with Allo

Use this for **single-end** SRA data.

Example sample:

```bash
SAMPLE=SRR8375037
```

---

## A1. Download single-end SRA data

```bash
cd /media/ismam/Expansion/RNA_WGBS_ChIP/D_pulex/ChIP

conda activate chipseqTools

SAMPLE=SRR8375037
THREADS=16
```

Download:

```bash
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

## A3. Align single-end reads with Bowtie2 retaining multi-mappers

Allo needs multiple possible alignments, so use:

```bash
-k 25
```

Set index:

```bash
BOWTIE2_INDEX=input/ref/bowtie2_index/Daphnia_pulex
```

Align:

```bash
bowtie2 \
  -x ${BOWTIE2_INDEX} \
  -q input/trimmed/${SAMPLE}.trimmed.fastq.gz \
  -S input/aligned/${SAMPLE}_MMR.sam \
  -k 25 \
  -p ${THREADS} \
  2> input/logs/${SAMPLE}_bowtie2_MMR.log
```

Check alignment:

```bash
cat input/logs/${SAMPLE}_bowtie2_MMR.log
ls -lh input/aligned/${SAMPLE}_MMR.sam
```

---

## A4. Collate SAM by read name for Allo

Allo expects alignments for the same read to be grouped together.

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

## A5. Run Allo for single-end input/control data

For input/control/background samples, use:

```bash
--random
```

Run:

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

Check Allo tags:

```bash
grep -c "ZA:" input/allo_out/${SAMPLE}_allo_random.sam
grep -c "ZZ:" input/allo_out/${SAMPLE}_allo_random.sam
```

Meaning:

- `ZA` = multi-mapped reads rescued/assigned by Allo using unique-read evidence.
- `ZZ` = reads assigned randomly because no unique-read support was available.

For your `SRR8375037` run, you observed:

```text
ZA: 5,507,923
ZZ:   774,715
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

## A7. Create primary-only BAM for MACS3

Allo output may contain secondary alignments. Before MACS3, remove secondary alignments:

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
samtools flagstat input/allo_out/${SAMPLE}_allo_random.primary.sorted.bam
ls -lh input/allo_out/${SAMPLE}_allo_random.primary.sorted.bam*
```

Remove intermediate BAMs:

```bash
rm input/allo_out/${SAMPLE}_allo_random.bam
rm input/allo_out/${SAMPLE}_allo_random.primary.bam
```

Final single-end input/control BAM for MACS3:

```bash
input/allo_out/${SAMPLE}_allo_random.primary.sorted.bam
```

---

# PART B: Single-end treatment ChIP-seq with Allo

Use this if the treatment ChIP sample is also **single-end**.

Set treatment sample:

```bash
SAMPLE=SRR_TREATMENT
THREADS=16
BOWTIE2_INDEX=input/ref/bowtie2_index/Daphnia_pulex
```

Create folders:

```bash
mkdir -p chip/raw/sra chip/raw chip/trimmed chip/aligned chip/allo_out chip/logs
```

Download:

```bash
prefetch ${SAMPLE} \
  --output-directory chip/raw/sra

fasterq-dump chip/raw/sra/${SAMPLE} \
  --outdir chip/raw \
  --threads ${THREADS}

gzip -f chip/raw/${SAMPLE}.fastq
```

Trim:

```bash
fastp \
  -i chip/raw/${SAMPLE}.fastq.gz \
  -o chip/trimmed/${SAMPLE}.trimmed.fastq.gz \
  -w ${THREADS} \
  -h chip/trimmed/${SAMPLE}_fastp.html \
  -j chip/trimmed/${SAMPLE}_fastp.json
```

Align with multi-mapper retention:

```bash
bowtie2 \
  -x ${BOWTIE2_INDEX} \
  -q chip/trimmed/${SAMPLE}.trimmed.fastq.gz \
  -S chip/aligned/${SAMPLE}_MMR.sam \
  -k 25 \
  -p ${THREADS} \
  2> chip/logs/${SAMPLE}_bowtie2_MMR.log
```

Collate:

```bash
samtools collate \
  -@ ${THREADS} \
  -o chip/aligned/${SAMPLE}_MMR_collate.sam \
  chip/aligned/${SAMPLE}_MMR.sam
```

Run Allo. For treatment samples, normally run without `--random`:

```bash
allo chip/aligned/${SAMPLE}_MMR_collate.sam \
  -seq se \
  -o chip/allo_out/${SAMPLE}_allo.sam \
  -p ${THREADS}
```

Convert to sorted BAM:

```bash
samtools view -@ ${THREADS} -bS \
  chip/allo_out/${SAMPLE}_allo.sam \
  > chip/allo_out/${SAMPLE}_allo.bam

samtools sort -@ ${THREADS} \
  -o chip/allo_out/${SAMPLE}_allo.sorted.bam \
  chip/allo_out/${SAMPLE}_allo.bam

samtools index chip/allo_out/${SAMPLE}_allo.sorted.bam
```

Primary-only BAM:

```bash
samtools view -@ ${THREADS} -b -F 256 \
  chip/allo_out/${SAMPLE}_allo.sorted.bam \
  > chip/allo_out/${SAMPLE}_allo.primary.bam

samtools sort -@ ${THREADS} \
  -o chip/allo_out/${SAMPLE}_allo.primary.sorted.bam \
  chip/allo_out/${SAMPLE}_allo.primary.bam

samtools index chip/allo_out/${SAMPLE}_allo.primary.sorted.bam
```

Final single-end treatment BAM:

```bash
chip/allo_out/${SAMPLE}_allo.primary.sorted.bam
```

---

# PART C: Paired-end treatment ChIP-seq with Allo

Use this for **paired-end** treatment ChIP-seq.

Set sample:

```bash
SAMPLE=SRR_TREATMENT
THREADS=16
BOWTIE2_INDEX=input/ref/bowtie2_index/Daphnia_pulex
```

Create folders:

```bash
mkdir -p chip/raw/sra chip/raw chip/trimmed chip/aligned chip/allo_out chip/logs
```

---

## C1. Download paired-end SRA data

```bash
prefetch ${SAMPLE} \
  --output-directory chip/raw/sra

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

## C2. Trim paired-end reads with fastp

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

---

## C3. Align paired-end reads with Bowtie2 retaining multi-mappers

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

## C4. Collate paired-end SAM by read name

```bash
samtools collate \
  -@ ${THREADS} \
  -o chip/aligned/${SAMPLE}_MMR_collate.sam \
  chip/aligned/${SAMPLE}_MMR.sam
```

---

## C5. Run Allo for paired-end treatment data

```bash
allo chip/aligned/${SAMPLE}_MMR_collate.sam \
  -seq pe \
  -o chip/allo_out/${SAMPLE}_allo.sam \
  -p ${THREADS}
```

Check tags:

```bash
grep -c "ZA:" chip/allo_out/${SAMPLE}_allo.sam
grep -c "ZZ:" chip/allo_out/${SAMPLE}_allo.sam
```

---

## C6. Convert paired-end Allo SAM to sorted, primary-only BAM

```bash
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
```

Check:

```bash
samtools flagstat chip/allo_out/${SAMPLE}_allo.primary.sorted.bam
```

Final paired-end treatment BAM:

```bash
chip/allo_out/${SAMPLE}_allo.primary.sorted.bam
```

---

# PART D: MACS3 peak calling

## D1. Choose treatment and control BAMs

Example:

```bash
TREAT_BAM=chip/allo_out/SRR_TREATMENT_allo.primary.sorted.bam
CONTROL_BAM=input/allo_out/SRR8375037_allo_random.primary.sorted.bam
OUTDIR=macs3_out
NAME=Dpulex_ChIP_vs_Input_Allo_MACS3
GENOME_SIZE=1.8e8
```

For *Daphnia pulex*, update `GENOME_SIZE` if you calculate a better effective genome size.

---

## D2. MACS3 narrow peak calling for single-end ChIP-seq

Use this if the treatment BAM is **single-end**:

```bash
macs3 callpeak \
  -t ${TREAT_BAM} \
  -c ${CONTROL_BAM} \
  -f BAM \
  -g ${GENOME_SIZE} \
  -n ${NAME} \
  --outdir ${OUTDIR} \
  -q 0.05 \
  --nomodel \
  --extsize 147 \
  --keep-dup all \
  --call-summits
```

Notes:

- `-f BAM` = single-end BAM.
- `--nomodel --extsize 147` = fixed extension size; useful when model building is unstable.
- `--call-summits` gives summit positions for narrow peaks.
- `--keep-dup all` avoids MACS3 removing reads after Allo processing. You can change this if you want duplicate filtering.

---

## D3. MACS3 narrow peak calling for paired-end ChIP-seq

Use this if the treatment BAM is **paired-end**:

```bash
macs3 callpeak \
  -t ${TREAT_BAM} \
  -c ${CONTROL_BAM} \
  -f BAMPE \
  -g ${GENOME_SIZE} \
  -n ${NAME} \
  --outdir ${OUTDIR} \
  -q 0.05 \
  --keep-dup all \
  --call-summits
```

Notes:

- `-f BAMPE` = paired-end BAM.
- MACS3 uses paired-end fragment information directly.
- Do not use `--nomodel --extsize` for normal paired-end BAMPE mode.

---

## D4. MACS3 broad peak calling

For broad histone marks such as H3K27me3 or H3K36me3, use `--broad`.

Single-end broad peak example:

```bash
macs3 callpeak \
  -t ${TREAT_BAM} \
  -c ${CONTROL_BAM} \
  -f BAM \
  -g ${GENOME_SIZE} \
  -n ${NAME}_broad \
  --outdir ${OUTDIR} \
  -q 0.05 \
  --broad \
  --broad-cutoff 0.05 \
  --nomodel \
  --extsize 147 \
  --keep-dup all
```

Paired-end broad peak example:

```bash
macs3 callpeak \
  -t ${TREAT_BAM} \
  -c ${CONTROL_BAM} \
  -f BAMPE \
  -g ${GENOME_SIZE} \
  -n ${NAME}_broad \
  --outdir ${OUTDIR} \
  -q 0.05 \
  --broad \
  --broad-cutoff 0.05 \
  --keep-dup all
```

---

## D5. MACS3 output files

For narrow peaks, expected files include:

```bash
${NAME}_peaks.narrowPeak
${NAME}_peaks.xls
${NAME}_summits.bed
${NAME}_treat_pileup.bdg
${NAME}_control_lambda.bdg
```

For broad peaks, expected files include:

```bash
${NAME}_peaks.broadPeak
${NAME}_peaks.gappedPeak
${NAME}_peaks.xls
${NAME}_treat_pileup.bdg
${NAME}_control_lambda.bdg
```

Check:

```bash
ls -lh ${OUTDIR}/
```

---

# PART E: Post-peak QC

## E1. Count narrow peaks

```bash
wc -l ${OUTDIR}/${NAME}_peaks.narrowPeak
```

## E2. Sort narrow peaks

```bash
sort -k1,1 -k2,2n \
  ${OUTDIR}/${NAME}_peaks.narrowPeak \
  > ${OUTDIR}/${NAME}_peaks.sorted.narrowPeak
```

## E3. Filter strong narrow peaks

Column 9 in `narrowPeak` is `-log10(qvalue)`.

Keep peaks with q-value approximately ≤ 0.01:

```bash
awk '$9 >= 2' \
  ${OUTDIR}/${NAME}_peaks.sorted.narrowPeak \
  > ${OUTDIR}/${NAME}_peaks.q0.01.narrowPeak
```

## E4. Count broad peaks

```bash
wc -l ${OUTDIR}/${NAME}_broad_peaks.broadPeak
```

---

# PART F: Generate bigWig tracks

## F1. Treatment bigWig

```bash
bamCoverage \
  -b ${TREAT_BAM} \
  -o bigwig/${NAME}_treatment_RPGC.bw \
  --normalizeUsing RPGC \
  --effectiveGenomeSize 180000000 \
  --binSize 10 \
  -p ${THREADS}
```

## F2. Control bigWig

```bash
bamCoverage \
  -b ${CONTROL_BAM} \
  -o bigwig/${NAME}_control_RPGC.bw \
  --normalizeUsing RPGC \
  --effectiveGenomeSize 180000000 \
  --binSize 10 \
  -p ${THREADS}
```

## F3. Treatment minus control bigWig

```bash
bamCompare \
  -b1 ${TREAT_BAM} \
  -b2 ${CONTROL_BAM} \
  -o bigwig/${NAME}_treatment_vs_control_log2ratio.bw \
  --operation log2 \
  --binSize 10 \
  -p ${THREADS}
```

---

# PART G: Cleanup

After confirming final BAMs, indexes, and MACS3 outputs exist:

```bash
rm input/allo_out/*_allo_random.bam 2>/dev/null
rm input/allo_out/*_allo_random.primary.bam 2>/dev/null
rm chip/allo_out/*_allo.bam 2>/dev/null
rm chip/allo_out/*_allo.primary.bam 2>/dev/null
```

Optional, only after full confirmation:

```bash
# rm input/aligned/*_MMR.sam
# rm input/aligned/*_MMR_collate.sam
# rm input/allo_out/*_allo_random.sam
# rm chip/aligned/*_MMR.sam
# rm chip/aligned/*_MMR_collate.sam
# rm chip/allo_out/*_allo.sam
```

Do not remove:

```bash
*.primary.sorted.bam
*.primary.sorted.bam.bai
*_peaks.narrowPeak
*_peaks.broadPeak
*_summits.bed
*.bw
```

---

# PART H: Complete single-end input/control example

This is your `SRR8375037` input/control command set.

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

Final control BAM:

```bash
input/allo_out/SRR8375037_allo_random.primary.sorted.bam
```

---

# PART I: Complete MACS3 example

Set treatment and control:

```bash
TREAT_BAM=chip/allo_out/SRR_TREATMENT_allo.primary.sorted.bam
CONTROL_BAM=input/allo_out/SRR8375037_allo_random.primary.sorted.bam
OUTDIR=macs3_out
NAME=Dpulex_ChIP_vs_Input_Allo_MACS3
GENOME_SIZE=1.8e8
THREADS=16
```

For single-end treatment:

```bash
macs3 callpeak \
  -t ${TREAT_BAM} \
  -c ${CONTROL_BAM} \
  -f BAM \
  -g ${GENOME_SIZE} \
  -n ${NAME} \
  --outdir ${OUTDIR} \
  -q 0.05 \
  --nomodel \
  --extsize 147 \
  --keep-dup all \
  --call-summits
```

For paired-end treatment:

```bash
macs3 callpeak \
  -t ${TREAT_BAM} \
  -c ${CONTROL_BAM} \
  -f BAMPE \
  -g ${GENOME_SIZE} \
  -n ${NAME} \
  --outdir ${OUTDIR} \
  -q 0.05 \
  --keep-dup all \
  --call-summits
```

Check output:

```bash
ls -lh ${OUTDIR}/
wc -l ${OUTDIR}/${NAME}_peaks.narrowPeak
```

---

# PART J: Single-end versus paired-end quick reference

## Single-end raw sequence

Single-end SRA produces one FASTQ file:

```bash
SAMPLE.fastq.gz
```

Use:

```bash
fasterq-dump SAMPLE --outdir raw --threads 16
```

Fastp:

```bash
fastp -i raw/SAMPLE.fastq.gz -o trimmed/SAMPLE.trimmed.fastq.gz
```

Bowtie2:

```bash
bowtie2 -x INDEX -q trimmed/SAMPLE.trimmed.fastq.gz -S SAMPLE.sam -k 25
```

Allo:

```bash
allo SAMPLE_collate.sam -seq se -o SAMPLE_allo.sam
```

MACS3:

```bash
macs3 callpeak -t TREAT.bam -c CONTROL.bam -f BAM
```

---

## Paired-end raw sequence

Paired-end SRA produces two FASTQ files:

```bash
SAMPLE_1.fastq.gz
SAMPLE_2.fastq.gz
```

Use:

```bash
fasterq-dump SAMPLE --split-files --outdir raw --threads 16
```

Fastp:

```bash
fastp \
  -i raw/SAMPLE_1.fastq.gz \
  -I raw/SAMPLE_2.fastq.gz \
  -o trimmed/SAMPLE_1.trimmed.fastq.gz \
  -O trimmed/SAMPLE_2.trimmed.fastq.gz
```

Bowtie2:

```bash
bowtie2 \
  -x INDEX \
  -1 trimmed/SAMPLE_1.trimmed.fastq.gz \
  -2 trimmed/SAMPLE_2.trimmed.fastq.gz \
  -S SAMPLE.sam \
  -k 25 \
  --no-mixed \
  --no-discordant
```

Allo:

```bash
allo SAMPLE_collate.sam -seq pe -o SAMPLE_allo.sam
```

MACS3:

```bash
macs3 callpeak -t TREAT.bam -c CONTROL.bam -f BAMPE
```

---

# PART K: Methods text

You can use this in a manuscript or GitHub README:

> Raw ChIP-seq reads were downloaded from the NCBI Sequence Read Archive using SRA-tools and quality-trimmed with fastp. Trimmed reads were aligned to the *Daphnia pulex* reference genome using Bowtie2 while retaining up to 25 possible alignments per read (`-k 25`) to preserve multi-mapping information. Alignments were collated by read name using samtools and processed with Allo to assign multi-mapped reads based on uniquely mapped read support. For input/control samples, Allo was run with the `--random` option to retain reads in regions lacking unique-read evidence. Allo-processed SAM files were converted to sorted BAM files, indexed, and filtered to retain primary alignments before peak calling. Peaks were called using MACS3 with the Allo-processed ChIP BAM as treatment and the Allo-processed input BAM as control.

---

# PART L: Troubleshooting

## `macs3: command not found`

Install MACS3:

```bash
conda activate chipseqTools
conda install -c conda-forge -c bioconda macs3 -y
macs3 --version
```

## `allo: command not found`

Install Allo:

```bash
conda activate chipseqTools
cd /home/ismam/allo
pip install -e .
allo --help
```

## Bowtie2 index not found

If you are in:

```bash
/media/ismam/Expansion/RNA_WGBS_ChIP/D_pulex/ChIP
```

use:

```bash
input/ref/bowtie2_index/Daphnia_pulex
```

If you are in:

```bash
/media/ismam/Expansion/RNA_WGBS_ChIP/D_pulex/ChIP/input
```

use:

```bash
ref/bowtie2_index/Daphnia_pulex
```

## MACS3 format error

Single-end:

```bash
-f BAM
```

Paired-end:

```bash
-f BAMPE
```

## Too many secondary alignments

Use primary-only BAM:

```bash
samtools view -b -F 256 input.sorted.bam > input.primary.bam
```

Then sort and index:

```bash
samtools sort -o input.primary.sorted.bam input.primary.bam
samtools index input.primary.sorted.bam
```

---

# End of pipeline
