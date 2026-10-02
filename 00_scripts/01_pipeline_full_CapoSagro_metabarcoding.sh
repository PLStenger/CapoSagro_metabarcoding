#!/bin/bash
#SBATCH --job-name=CapoSagro
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=36
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=/home/plstenge/CapoSagro_metabarcoding/00_scripts/01_pipeline_full_CapoSagro_metabarcoding.err
#SBATCH --output=/home/plstenge/CapoSagro_metabarcoding/00_scripts/01_pipeline_full_CapoSagro_metabarcoding.out

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================
module load conda/4.12.0
source ~/.bashrc
conda activate rachis-qiime2-2026.7

#############################################################################
# FULL PIPELINE - CapoSagro_metabarcoding - multi-marker analysis (7 markers)
#   m1 = V4V5        (16S prokaryotes, 515F-Y / 926R, Parada et al. 2016)
#   m2 = rbcL        (rbcL_F52 / rbcL_R193)
#   m3 = mlCOI/jgHCO (mlCOIintF / jgHCO2198, Leray et al. 2013)
#   m4 = teleo       (12S, teleo_F L1848 / teleo_R H1913, Valentini et al. 2016)
#   m5 = marine fungi(ITS2, ITS7 / ITS4)
#   m6 = V3V4 18S    (Eukaryotes, 515F / Ek-NSR951)
#   m7 = ITS2        (ITS2-S2F / ITS4)
#
# Same structure as the Grand_Saint_Antoine pipeline:
#   adapters -> FastQC/MultiQC raw -> Trimmomatic -> FastQC/MultiQC cleaned
#   -> manifests -> metadata -> QIIME2 import -> cutadapt -> DADA2
#   -> taxonomy -> controls -> decontam -> filtering -> barplots/exports
#
# Changes vs GSA (needed because markers differ):
#   - Trimmomatic MINLEN per marker (teleo/rbcL amplicons < 150 bp !)
#   - cutadapt in 2 passes: (1) 5' primers + discard-untrimmed,
#     (2) 3' reverse-complement of the opposite primer (read-through on
#     short amplicons: teleo, rbcL, short ITS2)
#   - inosine (I) in jgHCO2198 converted to N for cutadapt
#   - generic DADA2 function with automatic fallback
#     (independent+consensus -> pseudo+consensus -> pseudo+none)
#   - one reference database per marker (SILVA / NCBI / MIDORI2 / UNITE / PR2)
#############################################################################

########################
# 0. General parameters
########################
# Raw data: same sequencing run / folder as Grand_Saint_Antoine
BASEDIR="/storage/groups/gdec/shared_paleo/E1739/20260916_AV241601_E1739-Ps12-Testscustom_recipe_15092026"
PROJECTDIR="/home/plstenge/CapoSagro_metabarcoding"
WORKDIR="$PROJECTDIR/20260916_AV241601_E1739-Ps12-Testscustom_recipe_15092026"
# Project-level reference DB directory (re-used across runs)
TAXODIR_GLOBAL="$PROJECTDIR/taxonomy"
# GSA taxonomy dir (contains the SILVA 515f-926r classifier already trained)
TAXODIR_GSA="/home/plstenge/Grand_Saint_Antoine/taxonomy"

THREADS=32
TRIM_THREADS=8
JAVA_MEM="60G"
USE_ILLUMINA_ADAPTER_TRIM=0
SKIP_EXISTING=1            # 1 = do not redo Trimmomatic/refs if outputs exist

# Reference database versions (RESCRIPt)
UNITE_VERSION="2025-02-19"
MIDORI2_VERSION="GenBank270_2026-02-15"
PR2_VERSION="5.1.1"
SILVA_VERSION="138.2"

mkdir -p "$PROJECTDIR/00_scripts" "$TAXODIR_GLOBAL"
mkdir -p "$WORKDIR"
mkdir -p "$WORKDIR"/{manifests,imported,trimmed,dada2,taxonomy,decontam,exports,metadata,logs,tmp}
mkdir -p "$WORKDIR"/qc/{raw,cleaned}/{m1,m2,m3,m4,m5,m6,m7}
mkdir -p "$WORKDIR"/trimmomatic
cd "$WORKDIR" || exit 1
export TMPDIR="$WORKDIR/tmp"

META="$WORKDIR/metadata/sample-metadata.tsv"

MARKERS=(m1 m2 m3 m4 m5 m6 m7)

declare -A MARKER_NAME
MARKER_NAME[m1]="V4V5"
MARKER_NAME[m2]="rbcL"
MARKER_NAME[m3]="mlCOI_jgHCO"
MARKER_NAME[m4]="teleo"
MARKER_NAME[m5]="marine_fungi_ITS2"
MARKER_NAME[m6]="V3V4_18S"
MARKER_NAME[m7]="ITS2"

# Gene-specific primers only (Illumina tails removed)
declare -A PRIMER_F
declare -A PRIMER_R
PRIMER_F[m1]="GTGYCAGCMGCCGCGGTAA"          # 515F-Y
PRIMER_R[m1]="CCGYCAATTYMTTTRAGTTT"         # 926R
PRIMER_F[m2]="GTTGGATTCAAAGCTGGTGTTA"       # rbcL_F52
PRIMER_R[m2]="CVGTCCAMACAGTWGTCCATGT"       # rbcL_R193
PRIMER_F[m3]="GGWACWGGWTGAACWGTWTAYCCYCC"   # mlCOIintF
PRIMER_R[m3]="TAIACYTCIGGRTGICCRAARAAYCA"   # jgHCO2198 (I = inosine)
PRIMER_F[m4]="ACACCGCCCGTCACTCT"            # teleo_F (L1848)
PRIMER_R[m4]="CTTCCGGTACACTTACCATG"         # teleo_R (H1913)
PRIMER_F[m5]="GTGARTCATCGAATCTTTG"          # ITS7
PRIMER_R[m5]="TCCTCCGCTTATTGATATGC"         # ITS4
PRIMER_F[m6]="GTGCCAGCMGCCGCGG"             # 515F (18S)
PRIMER_R[m6]="TTGGYRAATGCTTTCGC"            # Ek-NSR951
PRIMER_F[m7]="ATGCGATACTTGGTGTGAAT"         # ITS2-S2F
PRIMER_R[m7]="TCCTCCGCTTATTGATATGC"         # ITS4

# cutadapt does not accept inosine: I -> N
for m in "${MARKERS[@]}"; do
    PRIMER_F[$m]=$(echo "${PRIMER_F[$m]}" | tr 'Ii' 'NN')
    PRIMER_R[$m]=$(echo "${PRIMER_R[$m]}" | tr 'Ii' 'NN')
done

# IUPAC-aware reverse complement
# (pure bash: does not depend on the "rev" binary)
revcomp() {
    local seq="$1" out="" i
    for (( i=${#seq}-1; i>=0; i-- )); do out+="${seq:$i:1}"; done
    echo "$out" | tr 'ACGTRYKMBDHVNSWacgtrykmbdhvnsw' 'TGCAYRMKVHDBNSWtgcayrmkvhdbnsw'
}
declare -A PRIMER_F_RC
declare -A PRIMER_R_RC
for m in "${MARKERS[@]}"; do
    PRIMER_F_RC[$m]=$(revcomp "${PRIMER_F[$m]}")
    PRIMER_R_RC[$m]=$(revcomp "${PRIMER_R[$m]}")
done

echo "=== Primers used ==="
for m in "${MARKERS[@]}"; do
    printf "%s\t%-18s F=%s\tR=%s\tRC(R)=%s\tRC(F)=%s\n" "$m" "${MARKER_NAME[$m]}" \
        "${PRIMER_F[$m]}" "${PRIMER_R[$m]}" "${PRIMER_R_RC[$m]}" "${PRIMER_F_RC[$m]}"
done

# Trimmomatic MINLEN per marker (amplicon sizes incl. primers are very different:
# teleo ~100 bp, rbcL ~165 bp -> MINLEN 150 of GSA would remove ALL reads)
declare -A TRIM_MINLEN
TRIM_MINLEN[m1]=150
TRIM_MINLEN[m2]=80
TRIM_MINLEN[m3]=150
TRIM_MINLEN[m4]=50
TRIM_MINLEN[m5]=100
TRIM_MINLEN[m6]=150
TRIM_MINLEN[m7]=100
TRIM_LEADING=30
TRIM_TRAILING=30
TRIM_SLIDINGWINDOW="26:30"

ILLUMINA_FWD="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ILLUMINA_REV="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

#####################################################
# 0bis. Sample directory list (CapoSagro project)
#   <site>-a1..a3   : 1069 1070 1177 1178
#   <site>-ab1..ab3 : 1170 1171 1172 1173 1174 1175 1176
#   ntc1..ntc6 + neg-pcr  (controls)
#   x 7 markers (-m1 ... -m7)  => 40 x 7 = 280 libraries
#####################################################
SITES_A=(1069 1070 1177 1178)
SITES_AB=(1170 1171 1172 1173 1174 1175 1176)

CAPO_DIRS=()
for m in "${MARKERS[@]}"; do
    for s in "${SITES_A[@]}";  do for r in 1 2 3; do CAPO_DIRS+=("${s}-a${r}-${m}");  done; done
    for s in "${SITES_AB[@]}"; do for r in 1 2 3; do CAPO_DIRS+=("${s}-ab${r}-${m}"); done; done
    for n in 1 2 3 4 5 6; do CAPO_DIRS+=("ntc${n}-${m}"); done
    CAPO_DIRS+=("neg-pcr-${m}")
done

printf "%s\n" "${CAPO_DIRS[@]}" > "$WORKDIR/metadata/expected_samples.txt"
echo "Expected libraries: ${#CAPO_DIRS[@]} (list: $WORKDIR/metadata/expected_samples.txt)"

#####################################################
# 1. Illumina universal adapter FASTA (Trimmomatic)
#####################################################
ADAPTERFILE="$WORKDIR/trimmomatic/illumina_universal_adapters.fa"
echo "=== Writing Illumina universal adapter file: $ADAPTERFILE ==="
cat > "$ADAPTERFILE" << 'ADAPTER_EOF'
>PrefixPE/1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>PrefixPE/2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>TruSeq2_SE
AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
>TruSeq2_PE_fwd
AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
>TruSeq2_PE_rev
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
>TruSeq3_IndexedAdapter
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>TruSeq3_UniversalAdapter
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>Nextera_Trans1
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
>Nextera_Trans2
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>NexteraPE-PE/1
CTGTCTCTTATACACATCT
>NexteraPE-PE/2
CTGTCTCTTATACACATCT
>Illumina_Single_End_Adapter_1
GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
>Illumina_Single_End_Adapter_2
GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
>Illumina_Paired_End_Adapter_1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>Illumina_Paired_End_Adapter_2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>Illumina_Multiplexing_Adapter_1
GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>Illumina_Multiplexing_Adapter_2
GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>Illumina_Multiplexing_Index_Sequencing_Primer
GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>Illumina_Multiplexing_Read2_Sequencing_Primer
GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>Illumina_DpnII_expression_Adapter_1
GATCGGAAGAGCACACGTCTGAACTCCAGTCA
>Illumina_DpnII_expression_Adapter_2
CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT
>Illumina_DpnII_Gex_Sequencing_Primer
CGACAGGTTCAGAGTTCTACAGTCCGACGATC
>Illumina_NlaIII_expression_Adapter_1
GATCGGAAGAGCACACGTCTGAACTCCAGTCA
>Illumina_NlaIII_expression_Adapter_2
CATGCAGAAGACGGCATACGAGCTCTTCCGATCT
>Illumina_NlaIII_Gex_Sequencing_Primer
CCGACTATGCCGTCTGTTCCGAAGGTCCGACGATC
>Illumina_Small_RNA_Adapter_1
ATCTCGTATGCCGTCTTCTGCTTG
>Illumina_Small_RNA_RT_Primer
CAAGCAGAAGACGGCATACGA
>Illumina_Small_RNA_PCR_Primer_1
CAAGCAGAAGACGGCATACGA
>Illumina_Small_RNA_PCR_Primer_2
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
ADAPTER_EOF
[[ -s "$ADAPTERFILE" ]] || { echo "ERROR: adapter file was not created correctly." >&2; exit 1; }
echo "Adapter file ready with $(grep -c '^>' "$ADAPTERFILE") sequences."

#####################################################
# 2. Resolve raw R1/R2 fastq paths for every sample
#    (exact name first, then case-insensitive search)
#####################################################
declare -A RAW_R1
declare -A RAW_R2
MISSING_LOG="$WORKDIR/logs/missing_samples.txt"
: > "$MISSING_LOG"

echo "=== Resolving raw fastq paths for all samples ==="
for sampledir in "${CAPO_DIRS[@]}"; do
    d="$BASEDIR/$sampledir"
    if [[ ! -d "$d" ]]; then
        d=$(find "$BASEDIR" -mindepth 1 -maxdepth 1 -type d -iname "$sampledir" | head -n1)
    fi
    if [[ -z "$d" || ! -d "$d" ]]; then
        echo "WARNING: directory missing: $BASEDIR/$sampledir" >&2
        echo -e "$sampledir\tdirectory_missing" >> "$MISSING_LOG"
        continue
    fi
    r1=$(find "$d" -maxdepth 1 -type f -iname "*_R1*.fastq.gz" | sort | head -n1)
    r2=$(find "$d" -maxdepth 1 -type f -iname "*_R2*.fastq.gz" | sort | head -n1)
    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "WARNING: fastq files missing in $d" >&2
        echo -e "$sampledir\tfastq_missing" >> "$MISSING_LOG"
        continue
    fi
    RAW_R1[$sampledir]="$r1"
    RAW_R2[$sampledir]="$r2"
done

echo "=== Samples found per marker ==="
for marker in "${MARKERS[@]}"; do
    n=0
    for sampledir in "${CAPO_DIRS[@]}"; do
        [[ "$sampledir" == *-"$marker" && -n "${RAW_R1[$sampledir]:-}" ]] && n=$((n+1))
    done
    echo "$marker (${MARKER_NAME[$marker]}): $n / 40"
done
echo "Missing samples (if any): $MISSING_LOG"
[[ ${#RAW_R1[@]} -gt 0 ]] || { echo "ERROR: no fastq found in $BASEDIR for CapoSagro samples." >&2; exit 1; }

#####################################################
# 3. FastQC + MultiQC on RAW data, PER MARKER
#####################################################
echo "=== FastQC on raw reads ==="
for marker in "${MARKERS[@]}"; do
    outdir="$WORKDIR/qc/raw/${marker}"
    mkdir -p "$outdir"
    for sampledir in "${CAPO_DIRS[@]}"; do
        [[ "$sampledir" == *-"$marker" ]] || continue
        [[ -n "${RAW_R1[$sampledir]:-}" && -n "${RAW_R2[$sampledir]:-}" ]] || continue
        fastqc -t "$TRIM_THREADS" -o "$outdir" "${RAW_R1[$sampledir]}" "${RAW_R2[$sampledir]}" \
            2>&1 | tee -a "$WORKDIR/logs/fastqc_raw_${marker}.log"
    done
done

echo "=== MultiQC on raw reads, one report per marker ==="
for marker in "${MARKERS[@]}"; do
    outdir="$WORKDIR/qc/raw/${marker}"
    multiqc "$outdir" -o "$outdir" -n "multiqc_raw_${marker}_${MARKER_NAME[$marker]}" -f \
        2>&1 | tee -a "$WORKDIR/logs/multiqc_raw_${marker}.log"
done

#####################################################
# 4. Trimmomatic cleaning (MINLEN adapted per marker)
#####################################################
echo "=== Trimmomatic cleaning ==="
declare -A CLEAN_R1
declare -A CLEAN_R2

for sampledir in "${CAPO_DIRS[@]}"; do
    [[ -n "${RAW_R1[$sampledir]:-}" && -n "${RAW_R2[$sampledir]:-}" ]] || continue
    marker="${sampledir##*-}"
    minlen="${TRIM_MINLEN[$marker]}"

    out_paired_r1="$WORKDIR/trimmomatic/${sampledir}_R1.paired.fastq.gz"
    out_single_r1="$WORKDIR/trimmomatic/${sampledir}_R1.single.fastq.gz"
    out_paired_r2="$WORKDIR/trimmomatic/${sampledir}_R2.paired.fastq.gz"
    out_single_r2="$WORKDIR/trimmomatic/${sampledir}_R2.single.fastq.gz"

    if [[ "$SKIP_EXISTING" -eq 1 && -s "$out_paired_r1" && -s "$out_paired_r2" ]]; then
        echo "Trimmomatic output exists for $sampledir, skipped."
    else
        trimmomatic PE -Xmx"$JAVA_MEM" -threads "$TRIM_THREADS" -phred33 \
            "${RAW_R1[$sampledir]}" "${RAW_R2[$sampledir]}" \
            "$out_paired_r1" "$out_single_r1" \
            "$out_paired_r2" "$out_single_r2" \
            ILLUMINACLIP:"$ADAPTERFILE":2:30:10 \
            LEADING:"$TRIM_LEADING" TRAILING:"$TRIM_TRAILING" \
            SLIDINGWINDOW:"$TRIM_SLIDINGWINDOW" MINLEN:"$minlen" \
            2>&1 | tee "$WORKDIR/logs/trimmomatic_${sampledir}.log"
    fi

    if [[ -s "$out_paired_r1" && -s "$out_paired_r2" ]]; then
        CLEAN_R1[$sampledir]="$out_paired_r1"
        CLEAN_R2[$sampledir]="$out_paired_r2"
    else
        echo "ERROR: Trimmomatic did not produce paired output for $sampledir" >&2
    fi
done

#####################################################
# 5. FastQC + MultiQC on CLEANED data, PER MARKER
#####################################################
echo "=== FastQC on Trimmomatic-cleaned reads ==="
for marker in "${MARKERS[@]}"; do
    outdir="$WORKDIR/qc/cleaned/${marker}"
    mkdir -p "$outdir"
    for sampledir in "${CAPO_DIRS[@]}"; do
        [[ "$sampledir" == *-"$marker" ]] || continue
        [[ -n "${CLEAN_R1[$sampledir]:-}" && -n "${CLEAN_R2[$sampledir]:-}" ]] || continue
        fastqc -t "$TRIM_THREADS" -o "$outdir" "${CLEAN_R1[$sampledir]}" "${CLEAN_R2[$sampledir]}" \
            2>&1 | tee -a "$WORKDIR/logs/fastqc_cleaned_${marker}.log"
    done
done

echo "=== MultiQC on cleaned reads, one report per marker ==="
for marker in "${MARKERS[@]}"; do
    outdir="$WORKDIR/qc/cleaned/${marker}"
    multiqc "$outdir" -o "$outdir" -n "multiqc_cleaned_${marker}_${MARKER_NAME[$marker]}" -f \
        2>&1 | tee -a "$WORKDIR/logs/multiqc_cleaned_${marker}.log"
done

#####################################################
# 6. QIIME2 manifests from CLEANED fastq
#    (only non-empty gz files are kept: an empty
#    control would make the import fail)
#####################################################
echo "=== Building manifests from cleaned reads ==="
for marker in "${MARKERS[@]}"; do
    manifest="$WORKDIR/manifests/manifest_${marker}.tsv"
    printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$manifest"
    for sampledir in "${CAPO_DIRS[@]}"; do
        [[ "$sampledir" == *-"$marker" ]] || continue
        [[ -n "${CLEAN_R1[$sampledir]:-}" && -n "${CLEAN_R2[$sampledir]:-}" ]] || continue
        nreads=$(zcat "${CLEAN_R1[$sampledir]}" | head -n 4 | wc -l)
        if [[ "$nreads" -lt 4 ]]; then
            echo "WARNING: $sampledir has 0 reads after Trimmomatic, excluded from manifest." >&2
            echo -e "$sampledir\tzero_reads_after_trimmomatic" >> "$MISSING_LOG"
            continue
        fi
        printf "%s\t%s\t%s\n" "$sampledir" "${CLEAN_R1[$sampledir]}" "${CLEAN_R2[$sampledir]}" >> "$manifest"
    done
    echo "Manifest ${marker} (${MARKER_NAME[$marker]}): $manifest  ($(($(wc -l < "$manifest")-1)) samples)"
    column -t -s $'\t' "$manifest" | head
    echo
done

#####################################################
# 7. Metadata (shared across markers)
#    sample-id | site | replicate | sample_type | marker
#    | marker_name | sample-or-control | control_type
#####################################################
printf "sample-id\tsite\treplicate\tsample_type\tmarker\tmarker_name\tsample-or-control\tcontrol_type\n" > "$META"

for sampledir in "${CAPO_DIRS[@]}"; do
    marker="${sampledir##*-}"          # m1
    rest="${sampledir%-*}"             # 1069-a1 | 1170-ab2 | ntc3 | neg-pcr
    case "$rest" in
        ntc*)
            site="$rest"; replicate="${rest#ntc}"; stype="ntc"
            soc="control"; ctype="NTC" ;;
        neg-pcr*)
            site="neg-pcr"; replicate="1"; stype="neg-pcr"
            soc="control"; ctype="PCR_negative" ;;
        *)
            site="${rest%%-*}"         # 1069
            rep="${rest#*-}"           # a1 / ab2
            replicate="${rep: -1}"     # 1
            stype="${rep%?}"           # a / ab
            soc="sample"; ctype="none" ;;
    esac
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$sampledir" "$site" "$replicate" "$stype" \
        "$marker" "${MARKER_NAME[$marker]}" "$soc" "$ctype" >> "$META"
done
echo "Metadata written: $META"

#####################################################
# 8. QIIME2 import + primer removal (cutadapt)
#   pass 1: 5' primers, discard reads without primer
#   pass 2: 3' read-through (RC of opposite primer),
#           reads kept even if not found
#####################################################
for marker in "${MARKERS[@]}"; do
    manifest="$WORKDIR/manifests/manifest_${marker}.tsv"
    demux="$WORKDIR/imported/demux_${marker}.qza"
    trimmed5="$WORKDIR/trimmed/trimmed5p_${marker}.qza"
    trimmed="$WORKDIR/trimmed/trimmed_${marker}.qza"
    trimmed_qzv="$WORKDIR/trimmed/trimmed_${marker}.qzv"

    if [[ $(wc -l < "$manifest") -lt 2 ]]; then
        echo "WARNING: manifest for $marker is empty, marker skipped." >&2
        continue
    fi

    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$manifest" \
        --input-format PairedEndFastqManifestPhred33V2 \
        --output-path "$demux" \
        2>&1 | tee "$WORKDIR/logs/import_${marker}.log"
    [[ -f "$demux" ]] || { echo "ERROR: import failed for $marker" >&2; continue; }

    ADAPT_OPTS=()
    if [[ "$USE_ILLUMINA_ADAPTER_TRIM" -eq 1 ]]; then
        ADAPT_OPTS=(--p-adapter-f "$ILLUMINA_FWD" --p-adapter-r "$ILLUMINA_REV")
    fi

    # pass 1 : 5' primers
    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences "$demux" \
        --p-cores "$THREADS" \
        --p-front-f "${PRIMER_F[$marker]}" \
        --p-front-r "${PRIMER_R[$marker]}" \
        "${ADAPT_OPTS[@]}" \
        --p-match-read-wildcards \
        --p-match-adapter-wildcards \
        --p-discard-untrimmed \
        --p-no-indels \
        --o-trimmed-sequences "$trimmed5" \
        --o-stats "$WORKDIR/trimmed/cutadapt5p_${marker}_stats.qza" \
        --verbose \
        2>&1 | tee "$WORKDIR/logs/cutadapt5p_${marker}.log"
    [[ -f "$trimmed5" ]] || { echo "ERROR: cutadapt pass 1 did not create $trimmed5" >&2; continue; }

    # pass 2 : 3' read-through of the opposite primer
    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences "$trimmed5" \
        --p-cores "$THREADS" \
        --p-adapter-f "${PRIMER_R_RC[$marker]}" \
        --p-adapter-r "${PRIMER_F_RC[$marker]}" \
        --p-match-read-wildcards \
        --p-match-adapter-wildcards \
        --p-overlap 10 \
        --p-minimum-length 30 \
        --o-trimmed-sequences "$trimmed" \
        --o-stats "$WORKDIR/trimmed/cutadapt3p_${marker}_stats.qza" \
        --verbose \
        2>&1 | tee "$WORKDIR/logs/cutadapt3p_${marker}.log"

    if [[ ! -f "$trimmed" ]]; then
        echo "WARNING: cutadapt pass 2 failed for $marker, using pass-1 output." >&2
        cp "$trimmed5" "$trimmed"
    fi

    qiime demux summarize \
        --i-data "$trimmed" \
        --o-visualization "$trimmed_qzv"
done

#####################################################
# 9. DADA2 denoise-paired, per marker
#    trunc = 0 everywhere except V4V5 (same as GSA m3)
#    because ITS / teleo / rbcL have variable or short
#    lengths. Automatic fallback strategy (GSA m2 fix).
#####################################################
declare -A TRUNC_F TRUNC_R MAXEE_F MAXEE_R
TRUNC_F[m1]=220; TRUNC_R[m1]=180; MAXEE_F[m1]=5.0;  MAXEE_R[m1]=5.0    # V4V5
TRUNC_F[m2]=0;   TRUNC_R[m2]=0;   MAXEE_F[m2]=5.0;  MAXEE_R[m2]=5.0    # rbcL
TRUNC_F[m3]=0;   TRUNC_R[m3]=0;   MAXEE_F[m3]=5.0;  MAXEE_R[m3]=5.0    # COI
TRUNC_F[m4]=0;   TRUNC_R[m4]=0;   MAXEE_F[m4]=5.0;  MAXEE_R[m4]=5.0    # teleo
TRUNC_F[m5]=0;   TRUNC_R[m5]=0;   MAXEE_F[m5]=5.0;  MAXEE_R[m5]=5.0    # marine fungi ITS2
TRUNC_F[m6]=0;   TRUNC_R[m6]=0;   MAXEE_F[m6]=5.0;  MAXEE_R[m6]=5.0    # V3V4 18S
TRUNC_F[m7]=0;   TRUNC_R[m7]=0;   MAXEE_F[m7]=5.0;  MAXEE_R[m7]=5.0    # ITS2

run_dada2() {
    local marker="$1" pooling="$2" chimera="$3" tag="$4"
    local trimmed="$WORKDIR/trimmed/trimmed_${marker}.qza"
    rm -f "$WORKDIR/dada2/table_${marker}.qza" "$WORKDIR/dada2/rep-seqs_${marker}.qza" \
          "$WORKDIR/dada2/stats_${marker}.qza" "$WORKDIR/dada2/base-transition-stats_${marker}.qza"
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs "$trimmed" \
        --p-trunc-len-f "${TRUNC_F[$marker]}" \
        --p-trunc-len-r "${TRUNC_R[$marker]}" \
        --p-max-ee-f "${MAXEE_F[$marker]}" \
        --p-max-ee-r "${MAXEE_R[$marker]}" \
        --p-pooling-method "$pooling" \
        --p-chimera-method "$chimera" \
        --p-n-threads "$THREADS" \
        --o-table "$WORKDIR/dada2/table_${marker}.qza" \
        --o-representative-sequences "$WORKDIR/dada2/rep-seqs_${marker}.qza" \
        --o-denoising-stats "$WORKDIR/dada2/stats_${marker}.qza" \
        --o-base-transition-stats "$WORKDIR/dada2/base-transition-stats_${marker}.qza" \
        --verbose \
        2>&1 | tee "$WORKDIR/logs/dada2_${marker}_${tag}.log"
    [[ -f "$WORKDIR/dada2/table_${marker}.qza" ]]
}

for marker in "${MARKERS[@]}"; do
    trimmed="$WORKDIR/trimmed/trimmed_${marker}.qza"
    [[ -f "$trimmed" ]] || { echo "ERROR: missing file $trimmed" >&2; continue; }

    echo "=== DADA2 $marker (${MARKER_NAME[$marker]}) ==="
    if run_dada2 "$marker" independent consensus "indep_consensus"; then
        echo "$marker: DADA2 OK (independent + consensus)"
    elif run_dada2 "$marker" pseudo consensus "pseudo_consensus"; then
        echo "$marker: DADA2 OK (pseudo + consensus)"
    elif run_dada2 "$marker" pseudo none "pseudo_nochimera"; then
        echo "$marker: DADA2 OK (pseudo, NO chimera removal)"
    else
        echo "ERROR: $marker DADA2 failed with all strategies." >&2
        continue
    fi

    qiime metadata tabulate \
        --m-input-file "$WORKDIR/dada2/stats_${marker}.qza" \
        --o-visualization "$WORKDIR/dada2/stats_${marker}.qzv" || true

    qiime feature-table summarize \
        --i-table "$WORKDIR/dada2/table_${marker}.qza" \
        --m-metadata-file "$META" \
        --o-feature-frequencies "$WORKDIR/dada2/feature-freq_${marker}.qza" \
        --o-sample-frequencies "$WORKDIR/dada2/sample-freq_${marker}.qza" \
        --o-summary "$WORKDIR/dada2/table_${marker}.qzv" || true

    qiime feature-table tabulate-seqs \
        --i-data "$WORKDIR/dada2/rep-seqs_${marker}.qza" \
        --o-visualization "$WORKDIR/dada2/rep-seqs_${marker}.qzv" || true
done

#####################################################
# 10. Taxonomy
#   m1 V4V5     : SILVA 138.2 515f-926r sklearn classifier
#                 (re-used from GSA, or built if missing)
#   m2 rbcL     : NCBI (RESCRIPt) + vsearch consensus
#   m3 COI      : MIDORI2 CO1 + vsearch consensus
#   m4 teleo    : MIDORI2 srRNA (12S) + vsearch consensus
#   m5 fungi    : UNITE (all eukaryotes, 99%) + vsearch
#   m6 18S V3V4 : PR2 + vsearch consensus
#   m7 ITS2     : UNITE (all eukaryotes, 99%) + vsearch
#####################################################
mkdir -p "$TAXODIR_GLOBAL"

# --- m1: SILVA classifier ---
CLF_NAME="silva-138.2-ssu-nr99-515f-926r-classifier.qza"
M1_CLASSIFIER=""
for c in "$WORKDIR/taxonomy/$CLF_NAME" "$TAXODIR_GLOBAL/$CLF_NAME" "$TAXODIR_GSA/$CLF_NAME"; do
    [[ -f "$c" ]] && { M1_CLASSIFIER="$c"; break; }
done

if [[ -z "$M1_CLASSIFIER" ]]; then
    echo "SILVA 515f-926r classifier not found -> building it (long step)..."
    S="$TAXODIR_GLOBAL/silva"; mkdir -p "$S"
    qiime rescript get-silva-data \
        --p-version "$SILVA_VERSION" --p-target SSURef_NR99 \
        --o-silva-sequences "$S/silva-rna-seqs.qza" \
        --o-silva-taxonomy "$S/silva-tax.qza" --verbose
    qiime rescript reverse-transcribe --i-rna-sequences "$S/silva-rna-seqs.qza" --o-dna-sequences "$S/silva-seqs.qza"
    qiime rescript cull-seqs --i-sequences "$S/silva-seqs.qza" --p-n-jobs "$THREADS" --o-clean-sequences "$S/silva-seqs-cleaned.qza"
    qiime feature-classifier extract-reads \
        --i-sequences "$S/silva-seqs-cleaned.qza" \
        --p-f-primer "${PRIMER_F[m1]}" --p-r-primer "${PRIMER_R[m1]}" \
        --p-n-jobs "$THREADS" --p-read-orientation forward \
        --o-reads "$S/silva-seqs-515f-926r.qza"
    qiime rescript dereplicate \
        --i-sequences "$S/silva-seqs-515f-926r.qza" --i-taxa "$S/silva-tax.qza" \
        --p-mode uniq --p-threads "$THREADS" \
        --o-dereplicated-sequences "$S/silva-seqs-515f-926r-derep.qza" \
        --o-dereplicated-taxa "$S/silva-tax-515f-926r-derep.qza"
    qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads "$S/silva-seqs-515f-926r-derep.qza" \
        --i-reference-taxonomy "$S/silva-tax-515f-926r-derep.qza" \
        --o-classifier "$TAXODIR_GLOBAL/$CLF_NAME"
    [[ -f "$TAXODIR_GLOBAL/$CLF_NAME" ]] && M1_CLASSIFIER="$TAXODIR_GLOBAL/$CLF_NAME"
fi

if [[ -n "$M1_CLASSIFIER" && -f "$WORKDIR/dada2/rep-seqs_m1.qza" ]]; then
    qiime feature-classifier classify-sklearn \
        --i-classifier "$M1_CLASSIFIER" \
        --i-reads "$WORKDIR/dada2/rep-seqs_m1.qza" \
        --p-n-jobs "$THREADS" \
        --o-classification "$WORKDIR/taxonomy/taxonomy_m1.qza" \
        --verbose \
        2>&1 | tee "$WORKDIR/logs/taxonomy_m1.log"
else
    echo "WARNING: m1 taxonomy skipped (classifier: ${M1_CLASSIFIER:-NONE})" >&2
fi

# --- reference DB builders (cached in $TAXODIR_GLOBAL) ---
declare -A REF_SEQS REF_TAX

build_unite() {
    local s="$TAXODIR_GLOBAL/unite_${UNITE_VERSION}_euk99_seqs.qza"
    local t="$TAXODIR_GLOBAL/unite_${UNITE_VERSION}_euk99_tax.qza"
    if [[ ! -f "$s" || ! -f "$t" || "$SKIP_EXISTING" -eq 0 ]]; then
        qiime rescript get-unite-data \
            --p-version "$UNITE_VERSION" --p-taxon-group eukaryotes --p-cluster-id 99 \
            --o-taxonomy "$t" --o-sequences "$s" --verbose \
            2>&1 | tee "$WORKDIR/logs/ref_unite.log"
    fi
    UNITE_SEQS="$s"; UNITE_TAX="$t"
}

build_midori2() {   # $1 = gene (CO1 / srRNA) ; $2 = marker
    local gene="$1" marker="$2"
    local s="$TAXODIR_GLOBAL/midori2_${MIDORI2_VERSION}_${gene}_seqs.qza"
    local t="$TAXODIR_GLOBAL/midori2_${MIDORI2_VERSION}_${gene}_tax.qza"
    if [[ ! -f "$s" || ! -f "$t" || "$SKIP_EXISTING" -eq 0 ]]; then
        local d="$TAXODIR_GLOBAL/midori2_${gene}_dl"
        rm -rf "$d"; mkdir -p "$d"
        qiime rescript get-midori2-data \
            --p-mito-gene "$gene" --p-version "$MIDORI2_VERSION" --p-ref-seq-type uniq \
            --o-midori2-sequences "$d/seqs" --o-midori2-taxonomy "$d/tax" --verbose \
            2>&1 | tee "$WORKDIR/logs/ref_midori2_${gene}.log"
        # Collections: take the .qza inside each output directory
        local fs ft
        fs=$(find "$d/seqs" -name "*.qza" | head -n1)
        ft=$(find "$d/tax"  -name "*.qza" | head -n1)
        if [[ -n "$fs" && -n "$ft" ]]; then
            qiime rescript cull-seqs --i-sequences "$fs" --p-num-degenerates 5 \
                --p-homopolymer-length 12 --p-n-jobs "$THREADS" --o-clean-sequences "$s"
            cp "$ft" "$t"
        fi
    fi
    REF_SEQS[$marker]="$s"; REF_TAX[$marker]="$t"
}

build_pr2() {
    local s="$TAXODIR_GLOBAL/pr2_${PR2_VERSION}_seqs.qza"
    local t="$TAXODIR_GLOBAL/pr2_${PR2_VERSION}_tax.qza"
    if [[ ! -f "$s" || ! -f "$t" || "$SKIP_EXISTING" -eq 0 ]]; then
        qiime rescript get-pr2-data \
            --p-version "$PR2_VERSION" \
            --o-pr2-sequences "$s" --o-pr2-taxonomy "$t" --verbose \
            2>&1 | tee "$WORKDIR/logs/ref_pr2.log"
    fi
    REF_SEQS[m6]="$s"; REF_TAX[m6]="$t"
}

build_ncbi_rbcl() {
    local s="$TAXODIR_GLOBAL/ncbi_rbcL_seqs.qza"
    local t="$TAXODIR_GLOBAL/ncbi_rbcL_tax.qza"
    if [[ ! -f "$s" || ! -f "$t" || "$SKIP_EXISTING" -eq 0 ]]; then
        local q='rbcL[Gene] AND 150:3000[SLEN] AND (Bacillariophyta[Organism] OR Phaeophyceae[Organism] OR Rhodophyta[Organism] OR Chlorophyta[Organism] OR Streptophyta[Organism] OR Haptophyta[Organism] OR Cryptophyceae[Organism] OR Dinophyceae[Organism] OR Ochrophyta[Organism])'
        qiime rescript get-ncbi-data \
            --p-query "$q" --p-n-jobs 3 \
            --o-sequences "$TAXODIR_GLOBAL/ncbi_rbcL_raw_seqs.qza" \
            --o-taxonomy "$t" --verbose \
            2>&1 | tee "$WORKDIR/logs/ref_ncbi_rbcL.log"
        qiime rescript cull-seqs --i-sequences "$TAXODIR_GLOBAL/ncbi_rbcL_raw_seqs.qza" \
            --p-num-degenerates 5 --p-homopolymer-length 12 --p-n-jobs "$THREADS" \
            --o-clean-sequences "$TAXODIR_GLOBAL/ncbi_rbcL_culled_seqs.qza"
        qiime rescript filter-seqs-length \
            --i-sequences "$TAXODIR_GLOBAL/ncbi_rbcL_culled_seqs.qza" \
            --p-global-min 150 --p-global-max 3000 --p-threads "$THREADS" \
            --o-filtered-seqs "$s" \
            --o-discarded-seqs "$TAXODIR_GLOBAL/ncbi_rbcL_discarded_seqs.qza"
    fi
    REF_SEQS[m2]="$s"; REF_TAX[m2]="$t"
}

[[ -f "$WORKDIR/dada2/rep-seqs_m2.qza" ]] && build_ncbi_rbcl
[[ -f "$WORKDIR/dada2/rep-seqs_m3.qza" ]] && build_midori2 CO1 m3
[[ -f "$WORKDIR/dada2/rep-seqs_m4.qza" ]] && build_midori2 srRNA m4
if [[ -f "$WORKDIR/dada2/rep-seqs_m5.qza" || -f "$WORKDIR/dada2/rep-seqs_m7.qza" ]]; then
    build_unite
    REF_SEQS[m5]="$UNITE_SEQS"; REF_TAX[m5]="$UNITE_TAX"
    REF_SEQS[m7]="$UNITE_SEQS"; REF_TAX[m7]="$UNITE_TAX"
fi
[[ -f "$WORKDIR/dada2/rep-seqs_m6.qza" ]] && build_pr2

# vsearch identity thresholds per marker (adjust after first inspection)
declare -A PERC_ID
PERC_ID[m2]=0.90   # rbcL
PERC_ID[m3]=0.90   # COI
PERC_ID[m4]=0.97   # teleo 12S (fish, near-species level)
PERC_ID[m5]=0.90   # fungi ITS2
PERC_ID[m6]=0.90   # 18S V3V4
PERC_ID[m7]=0.90   # ITS2

for marker in m2 m3 m4 m5 m6 m7; do
    REP="$WORKDIR/dada2/rep-seqs_${marker}.qza"
    RS="${REF_SEQS[$marker]:-}"; RT="${REF_TAX[$marker]:-}"
    if [[ -f "$REP" && -n "$RS" && -f "$RS" && -n "$RT" && -f "$RT" ]]; then
        qiime feature-classifier classify-consensus-vsearch \
            --i-query "$REP" \
            --i-reference-reads "$RS" \
            --i-reference-taxonomy "$RT" \
            --p-perc-identity "${PERC_ID[$marker]}" \
            --p-min-consensus 0.51 \
            --p-top-hits-only \
            --p-maxaccepts 10 \
            --p-threads "$THREADS" \
            --o-classification "$WORKDIR/taxonomy/taxonomy_${marker}.qza" \
            --o-search-results "$WORKDIR/taxonomy/search_${marker}.qza" \
            --verbose \
            2>&1 | tee "$WORKDIR/logs/taxonomy_${marker}.log"
    else
        echo "WARNING: ${marker} (${MARKER_NAME[$marker]}) classification skipped (rep-seqs or reference DB unavailable)." >&2
    fi
done

for marker in "${MARKERS[@]}"; do
    TAXO="$WORKDIR/taxonomy/taxonomy_${marker}.qza"
    [[ -f "$TAXO" ]] && qiime metadata tabulate --m-input-file "$TAXO" \
        --o-visualization "$WORKDIR/taxonomy/taxonomy_${marker}.qzv" || true
done

#####################################################
# 11. Inspect negative controls BEFORE decontamination
#####################################################
for marker in "${MARKERS[@]}"; do
    TABLE="$WORKDIR/dada2/table_${marker}.qza"
    TAXO="$WORKDIR/taxonomy/taxonomy_${marker}.qza"
    [[ -f "$TABLE" && -f "$TAXO" ]] || continue

    qiime feature-table filter-samples \
        --i-table "$TABLE" \
        --m-metadata-file "$META" \
        --p-where "[marker]='${marker}' AND [sample-or-control]='control'" \
        --o-filtered-table "$WORKDIR/decontam/controls_table_${marker}.qza" \
        2>&1 | tee -a "$WORKDIR/logs/controls_${marker}.log" || true

    if [[ -f "$WORKDIR/decontam/controls_table_${marker}.qza" ]]; then
        qiime taxa barplot \
            --i-table "$WORKDIR/decontam/controls_table_${marker}.qza" \
            --i-taxonomy "$TAXO" \
            --m-metadata-file "$META" \
            --o-visualization "$WORKDIR/decontam/controls_barplot_${marker}.qzv" \
            2>&1 | tee -a "$WORKDIR/logs/controls_${marker}.log" || true
    fi
done

#####################################################
# 12. Decontam identification (prevalence method)
#####################################################
for marker in "${MARKERS[@]}"; do
    TABLE="$WORKDIR/dada2/table_${marker}.qza"
    REP="$WORKDIR/dada2/rep-seqs_${marker}.qza"
    SCORE="$WORKDIR/decontam/decontam-scores_${marker}.qza"
    SCORE_VIZ="$WORKDIR/decontam/decontam-scoreviz_${marker}.qzv"
    [[ -f "$TABLE" ]] || continue

    qiime quality-control decontam-identify \
        --i-table "$TABLE" \
        --m-metadata-file "$META" \
        --p-method prevalence \
        --p-prev-control-column sample-or-control \
        --p-prev-control-indicator control \
        --o-decontam-scores "$SCORE" \
        2>&1 | tee -a "$WORKDIR/logs/decontam_${marker}.log" || true

    if [[ -f "$SCORE" ]]; then
        if [[ -f "$REP" ]]; then
            qiime quality-control decontam-score-viz \
                --i-decontam-scores "$SCORE" --i-table "$TABLE" --i-rep-seqs "$REP" \
                --p-threshold 0.1 --o-visualization "$SCORE_VIZ" \
                2>&1 | tee -a "$WORKDIR/logs/decontam_${marker}.log" || true
        else
            qiime quality-control decontam-score-viz \
                --i-decontam-scores "$SCORE" --i-table "$TABLE" \
                --p-threshold 0.1 --o-visualization "$SCORE_VIZ" \
                2>&1 | tee -a "$WORKDIR/logs/decontam_${marker}.log" || true
        fi
    else
        echo "$marker: decontam-identify did not produce scores (not enough controls or table too small)."
    fi
done

#####################################################
# 13. Remove contaminant features and control samples
#     (manual filtering on decontam p-score, as in GSA)
#####################################################
for marker in "${MARKERS[@]}"; do
    TABLE="$WORKDIR/dada2/table_${marker}.qza"
    REP="$WORKDIR/dada2/rep-seqs_${marker}.qza"
    SCORE="$WORKDIR/decontam/decontam-scores_${marker}.qza"
    TABLE_DC="$WORKDIR/decontam/table_${marker}_decontam.qza"
    REP_DC="$WORKDIR/decontam/rep-seqs_${marker}_decontam.qza"
    TABLE_FINAL="$WORKDIR/decontam/table_${marker}_final.qza"
    REP_FINAL="$WORKDIR/decontam/rep-seqs_${marker}_final.qza"
    [[ -f "$TABLE" && -f "$REP" ]] || continue

    if [[ -f "$SCORE" ]]; then
        qiime feature-table filter-features \
            --i-table "$TABLE" \
            --m-metadata-file "$SCORE" \
            --p-where '[p] > 0.1 OR [p] IS NULL' \
            --o-filtered-table "$TABLE_DC" \
            2>&1 | tee -a "$WORKDIR/logs/filter_${marker}.log" || true
        if [[ -f "$TABLE_DC" ]]; then
            qiime feature-table filter-seqs \
                --i-data "$REP" --i-table "$TABLE_DC" --o-filtered-data "$REP_DC" \
                2>&1 | tee -a "$WORKDIR/logs/filter_${marker}.log" || true
        fi
    else
        echo "$marker: no decontam score available, using raw table/rep-seqs before sample filtering."
        cp "$TABLE" "$TABLE_DC"
        cp "$REP" "$REP_DC"
    fi

    if [[ -f "$TABLE_DC" ]]; then
        qiime feature-table filter-samples \
            --i-table "$TABLE_DC" \
            --m-metadata-file "$META" \
            --p-where "[sample-or-control]='sample' AND [marker]='${marker}'" \
            --o-filtered-table "$TABLE_FINAL" \
            2>&1 | tee -a "$WORKDIR/logs/filter_${marker}.log" || true
    fi

    if [[ -f "$TABLE_FINAL" && -f "$REP_DC" ]]; then
        qiime feature-table filter-seqs \
            --i-data "$REP_DC" --i-table "$TABLE_FINAL" --o-filtered-data "$REP_FINAL" \
            2>&1 | tee -a "$WORKDIR/logs/filter_${marker}.log" || true
    fi
done

#####################################################
# 14. Final taxonomy barplots + exports
#####################################################
for marker in "${MARKERS[@]}"; do
    TAXO="$WORKDIR/taxonomy/taxonomy_${marker}.qza"
    TABLE_FINAL="$WORKDIR/decontam/table_${marker}_final.qza"
    REP_FINAL="$WORKDIR/decontam/rep-seqs_${marker}_final.qza"
    TABLE_RAW="$WORKDIR/dada2/table_${marker}.qza"
    REP_RAW="$WORKDIR/dada2/rep-seqs_${marker}.qza"

    EXPORT_DIR="$WORKDIR/exports/${marker}_${MARKER_NAME[$marker]}_final"
    mkdir -p "$EXPORT_DIR"

    TABLE_TO_USE="$TABLE_RAW"; REP_TO_USE="$REP_RAW"
    if [[ -f "$TABLE_FINAL" && -f "$REP_FINAL" ]]; then
        TABLE_TO_USE="$TABLE_FINAL"; REP_TO_USE="$REP_FINAL"
    fi
    [[ -f "$TABLE_TO_USE" ]] || { echo "$marker: no table available, export skipped."; continue; }
    [[ -f "$REP_TO_USE" ]]   || { echo "$marker: no rep-seqs available, export skipped."; continue; }

    qiime feature-table summarize \
        --i-table "$TABLE_TO_USE" --m-metadata-file "$META" \
        --o-feature-frequencies "$EXPORT_DIR/feature-frequencies_${marker}.qza" \
        --o-sample-frequencies "$EXPORT_DIR/sample-frequencies_${marker}.qza" \
        --o-summary "$EXPORT_DIR/table_${marker}.qzv" || true

    qiime feature-table tabulate-seqs \
        --i-data "$REP_TO_USE" --o-visualization "$EXPORT_DIR/rep-seqs_${marker}.qzv" || true

    if [[ -f "$TAXO" ]]; then
        qiime taxa barplot \
            --i-table "$TABLE_TO_USE" --i-taxonomy "$TAXO" --m-metadata-file "$META" \
            --o-visualization "$EXPORT_DIR/taxa-barplot_${marker}.qzv" || true
    fi

    qiime tools export --input-path "$TABLE_TO_USE" --output-path "$EXPORT_DIR/table_export"
    qiime tools export --input-path "$REP_TO_USE"   --output-path "$EXPORT_DIR/repseq_export"
    [[ -f "$TAXO" ]] && qiime tools export --input-path "$TAXO" --output-path "$EXPORT_DIR/taxonomy_export"

    if [[ -f "$EXPORT_DIR/table_export/feature-table.biom" && -f "$EXPORT_DIR/taxonomy_export/taxonomy.tsv" ]]; then
        biom add-metadata \
            -i "$EXPORT_DIR/table_export/feature-table.biom" \
            -o "$EXPORT_DIR/feature-table-with-tax.biom" \
            --observation-metadata-fp "$EXPORT_DIR/taxonomy_export/taxonomy.tsv" \
            --observation-header "OTUID,taxonomy,confidence" \
            --sc-separated taxonomy
        biom convert \
            -i "$EXPORT_DIR/feature-table-with-tax.biom" \
            -o "$EXPORT_DIR/ASV_table_${marker}_${MARKER_NAME[$marker]}_taxonomy.tsv" \
            --to-tsv --header-key taxonomy
    fi
done

echo "=== Full CapoSagro pipeline finished ==="
