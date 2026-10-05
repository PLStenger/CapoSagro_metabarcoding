#!/bin/bash
#SBATCH --job-name=CapoSagro_taxo
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=36
#SBATCH --mem=500G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=/home/plstenge/CapoSagro_metabarcoding/00_scripts/02_taxonomy_m2_m7_CapoSagro.err
#SBATCH --output=/home/plstenge/CapoSagro_metabarcoding/00_scripts/02_taxonomy_m2_m7_CapoSagro.out

#############################################################################
# CapoSagro_metabarcoding - TAXONOMY m2 -> m7 (resume after 01_pipeline_full)
#
#   m2 rbcL mini (F52/R193)   : NCBI rbcL (RESCRIPt get-ncbi-data, 5 groups)
#   m3 mlCOI / jgHCO2198      : MIDORI2 CO1   (uniq)
#   m4 teleo (12S)            : MIDORI2 srRNA (uniq)
#   m5 marine fungi ITS7/ITS4 : UNITE 2025-02-19, all eukaryotes, 99 %
#   m6 18S V3V4 515F/NSR951   : PR2 5.1.1
#   m7 ITS2 S2F/ITS4          : UNITE 2025-02-19, all eukaryotes, 99 %
#   Classification: classify-consensus-vsearch (no training needed)
#
# WHY THE FIRST RUN FAILED: all RESCRIPt "get-*" actions need INTERNET.
# On most clusters the compute nodes (here partition smp) have no internet
# -> no .qza created -> "does not exist" for cull-seqs / find.
#
# USAGE (2 steps):
#  1) DOWNLOAD on a node WITH internet (login node), in screen/tmux or nohup:
#       conda activate rachis-qiime2-2026.7
#       nohup bash 02_taxonomy_m2_m7_CapoSagro.sh download > 02_download.log 2>&1 &
#     (proxy if needed: export https_proxy=http://proxy:port http_proxy=...)
#  2) CLASSIFICATION on the cluster:
#       sbatch 02_taxonomy_m2_m7_CapoSagro.sh classify
#  ("sbatch 02_taxonomy_m2_m7_CapoSagro.sh" without argument = "all":
#   downloads what is missing if internet is available, then classifies)
#############################################################################

MODE="${1:-all}"     # download | classify | all

# ==============================================================================
# ENVIRONMENT
# ==============================================================================
if [[ -z "${CONDA_DEFAULT_ENV:-}" || "${CONDA_DEFAULT_ENV}" != "rachis-qiime2-2026.7" ]]; then
    module load conda/4.12.0 2>/dev/null
    source ~/.bashrc
    conda activate rachis-qiime2-2026.7
fi

########################
# 0. Parameters
########################
PROJECTDIR="/home/plstenge/CapoSagro_metabarcoding"
WORKDIR="$PROJECTDIR/20260916_AV241601_E1739-Ps12-Testscustom_recipe_15092026"
TAXODIR_GLOBAL="$PROJECTDIR/taxonomy"
RAWDB="$TAXODIR_GLOBAL/raw_downloads"      # raw RESCRIPt downloads
REFDB="$TAXODIR_GLOBAL/ref_ready"          # cleaned DBs used for classification
META="$WORKDIR/metadata/sample-metadata.tsv"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
[[ "$THREADS" -gt 32 ]] && THREADS=32
SKIP_EXISTING=1
MAX_TRIES=3

UNITE_VERSION="2025-02-19"
MIDORI2_VERSION="GenBank270_2026-02-15"
PR2_VERSION="5.1.1"

MARKERS_TAXO=(m2 m3 m4 m5 m6 m7)
declare -A MARKER_NAME=(
  [m2]="rbcL" [m3]="mlCOI_jgHCO" [m4]="teleo"
  [m5]="marine_fungi_ITS2" [m6]="V3V4_18S" [m7]="ITS2" )

# vsearch parameters per marker (adjust after first inspection)
declare -A PERC_ID=( [m2]=0.90 [m3]=0.85 [m4]=0.95 [m5]=0.85 [m6]=0.85 [m7]=0.85 )
declare -A QUERY_COV=( [m2]=0.80 [m3]=0.80 [m4]=0.80 [m5]=0.70 [m6]=0.80 [m7]=0.70 )
MIN_CONSENSUS=0.51
MAXACCEPTS=10

# NCBI rbcL: split by group (smaller queries = fewer NCBI timeouts).
# rbcL mini F52/R193 (Little et al. 2014): land plants, seagrasses, macroalgae.
# Remove "Streptophyta" here if the download is too long (largest group).
RBCL_GROUPS=(Streptophyta Chlorophyta Rhodophyta Phaeophyceae Bacillariophyta)

mkdir -p "$PROJECTDIR/00_scripts" "$TAXODIR_GLOBAL" "$RAWDB" "$REFDB" \
         "$WORKDIR"/{taxonomy,logs,exports,decontam}

# IMPORTANT : temporaires courts et locaux au nœud.
# Ne pas employer $WORKDIR/tmp : chemin trop long + FS partagé/NFS.
export TMPDIR="/tmp/${USER}_q2"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"

echo "TMPDIR=${TMPDIR}"
python - <<'PY'
import tempfile
print("Python temporary directory:", tempfile.gettempdir())
PY

LOG="$WORKDIR/logs"
cd "$WORKDIR" || exit 1

echo "=== MODE: $MODE | THREADS: $THREADS | host: $(hostname) ==="

#####################################################
# Helpers
#####################################################
has_internet() {
    if command -v curl >/dev/null 2>&1; then
        curl -s --max-time 20 -o /dev/null "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi"
    else
        wget -q -T 20 -O /dev/null "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi"
    fi
}

ok() { [[ -s "$1" ]]; }     # file exists and is not empty

#####################################################
# 1. DOWNLOAD functions (need internet)
#####################################################
dl_unite() {
    local s="$RAWDB/unite_${UNITE_VERSION}_euk99_seqs.qza"
    local t="$RAWDB/unite_${UNITE_VERSION}_euk99_tax.qza"
    if [[ "$SKIP_EXISTING" -eq 1 ]] && ok "$s" && ok "$t"; then echo "UNITE already downloaded."; return 0; fi
    for try in $(seq 1 $MAX_TRIES); do
        echo "--- UNITE download, attempt $try ---"
        qiime rescript get-unite-data \
            --p-version "$UNITE_VERSION" --p-taxon-group eukaryotes --p-cluster-id 99 \
            --o-sequences "$s" --o-taxonomy "$t" --verbose \
            2>&1 | tee "$LOG/dl_unite_try${try}.log"
        ok "$s" && ok "$t" && return 0
        sleep 60
    done
    echo "ERROR: UNITE download failed (see $LOG/dl_unite_try*.log)" >&2; return 1
}

dl_pr2() {
    local s="$RAWDB/pr2_${PR2_VERSION}_seqs.qza"
    local t="$RAWDB/pr2_${PR2_VERSION}_tax.qza"
    if [[ "$SKIP_EXISTING" -eq 1 ]] && ok "$s" && ok "$t"; then echo "PR2 already downloaded."; return 0; fi
    for try in $(seq 1 $MAX_TRIES); do
        echo "--- PR2 download, attempt $try ---"
        qiime rescript get-pr2-data \
            --p-version "$PR2_VERSION" \
            --o-pr2-sequences "$s" --o-pr2-taxonomy "$t" --verbose \
            2>&1 | tee "$LOG/dl_pr2_try${try}.log"
        ok "$s" && ok "$t" && return 0
        sleep 60
    done
    echo "ERROR: PR2 download failed (see $LOG/dl_pr2_try*.log)" >&2; return 1
}

# MIDORI2 outputs are COLLECTIONS (directories containing .qza files)
dl_midori2() {
    local gene="$1"
    local s="$RAWDB/midori2_${MIDORI2_VERSION}_${gene}_seqs.qza"
    local t="$RAWDB/midori2_${MIDORI2_VERSION}_${gene}_tax.qza"
    if [[ "$SKIP_EXISTING" -eq 1 ]] && ok "$s" && ok "$t"; then echo "MIDORI2 $gene already downloaded."; return 0; fi
    local d="$RAWDB/midori2_${gene}_collection"
    for try in $(seq 1 $MAX_TRIES); do
        echo "--- MIDORI2 $gene download, attempt $try ---"
        rm -rf "$d"; mkdir -p "$d"
        qiime rescript get-midori2-data \
            --p-mito-gene "$gene" \
            --p-version "$MIDORI2_VERSION" \
            --p-ref-seq-type uniq \
            --o-midori2-sequences "$d/seqs" \
            --o-midori2-taxonomy "$d/tax" \
            --verbose \
            2>&1 | tee "$LOG/dl_midori2_${gene}_try${try}.log"
        local fs="" ft=""
        [[ -d "$d/seqs" ]] && fs=$(find "$d/seqs" -name "*.qza" | head -n1)
        [[ -d "$d/tax"  ]] && ft=$(find "$d/tax"  -name "*.qza" | head -n1)
        if [[ -n "$fs" && -n "$ft" ]]; then
            cp "$fs" "$s"; cp "$ft" "$t"
            return 0
        fi
        sleep 60
    done
    echo "ERROR: MIDORI2 $gene download failed (see $LOG/dl_midori2_${gene}_try*.log)" >&2; return 1
}

dl_ncbi_rbcl() {
    local all_ok=0
    for grp in "${RBCL_GROUPS[@]}"; do
        local s="$RAWDB/ncbi_rbcL_${grp}_seqs.qza"
        local t="$RAWDB/ncbi_rbcL_${grp}_tax.qza"
        if [[ "$SKIP_EXISTING" -eq 1 ]] && ok "$s" && ok "$t"; then echo "NCBI rbcL $grp already downloaded."; continue; fi
        local q="rbcL[Gene] AND ${grp}[Organism] AND 100:2000[SLEN]"
        local done_grp=0
        for try in $(seq 1 $MAX_TRIES); do
            echo "--- NCBI rbcL $grp, attempt $try : $q ---"
            qiime rescript get-ncbi-data \
                --p-query "$q" \
                --p-n-jobs 2 \
                --p-logging-level INFO \
                --o-sequences "$s" --o-taxonomy "$t" --verbose \
                2>&1 | tee "$LOG/dl_ncbi_rbcL_${grp}_try${try}.log"
            if ok "$s" && ok "$t"; then done_grp=1; break; fi
            sleep 120
        done
        [[ "$done_grp" -eq 1 ]] || { echo "ERROR: NCBI rbcL $grp failed." >&2; all_ok=1; }
    done
    return $all_ok
}

run_downloads() {
    if ! has_internet; then
        echo "ERROR: no internet access on $(hostname)." >&2
        echo "  -> run: bash $0 download   on the LOGIN node (or set https_proxy)." >&2
        return 1
    fi
    echo "Internet OK on $(hostname)."
    dl_unite
    dl_pr2
    dl_midori2 CO1
    dl_midori2 srRNA
    dl_ncbi_rbcl
    echo "=== Downloads finished. Content of $RAWDB : ==="
    ls -lh "$RAWDB"/*.qza 2>/dev/null
}

#####################################################
# 2. PREPARE reference DBs (no internet needed)
#####################################################
declare -A REF_SEQS REF_TAX

prep_simple() {   # $1 marker key name ; $2 raw seqs ; $3 raw tax ; $4 out prefix ; $5 cull(1/0)
    local s_raw="$2" t_raw="$3" pref="$4" cull="$5"
    local s_out="$REFDB/${pref}_seqs.qza" t_out="$REFDB/${pref}_tax.qza"
    if [[ "$SKIP_EXISTING" -eq 1 ]] && ok "$s_out" && ok "$t_out"; then
        echo "$pref ready."; return 0
    fi
    ok "$s_raw" && ok "$t_raw" || { echo "WARNING: raw DB missing for $pref ($s_raw)" >&2; return 1; }
    if [[ "$cull" -eq 1 ]]; then
        qiime rescript cull-seqs \
            --i-sequences "$s_raw" \
            --p-num-degenerates 5 --p-homopolymer-length 12 \
            --p-n-jobs "$THREADS" \
            --o-clean-sequences "$s_out" \
            2>&1 | tee "$LOG/prep_${pref}.log"
    else
        cp "$s_raw" "$s_out"
    fi
    cp "$t_raw" "$t_out"
    ok "$s_out" && ok "$t_out"
}

prep_rbcl() {
    local s_out="$REFDB/ncbi_rbcL_seqs.qza" t_out="$REFDB/ncbi_rbcL_tax.qza"
    if [[ "$SKIP_EXISTING" -eq 1 ]] && ok "$s_out" && ok "$t_out"; then echo "rbcL ready."; return 0; fi
    local seqs=() taxs=()
    for grp in "${RBCL_GROUPS[@]}"; do
        local s="$RAWDB/ncbi_rbcL_${grp}_seqs.qza" t="$RAWDB/ncbi_rbcL_${grp}_tax.qza"
        if ok "$s" && ok "$t"; then seqs+=("$s"); taxs+=("$t")
        else echo "WARNING: rbcL group $grp missing, ignored." >&2; fi
    done
    [[ ${#seqs[@]} -gt 0 ]] || { echo "WARNING: no NCBI rbcL group available." >&2; return 1; }

    local merged_s="$REFDB/tmp_rbcL_merged_seqs.qza" merged_t="$REFDB/tmp_rbcL_merged_tax.qza"
    if [[ ${#seqs[@]} -eq 1 ]]; then
        cp "${seqs[0]}" "$merged_s"; cp "${taxs[0]}" "$merged_t"
    else
        qiime feature-table merge-seqs --i-data "${seqs[@]}" --o-merged-data "$merged_s" \
            2>&1 | tee "$LOG/prep_rbcL.log"
        qiime feature-table merge-taxa --i-data "${taxs[@]}" --o-merged-data "$merged_t" \
            2>&1 | tee -a "$LOG/prep_rbcL.log"
    fi
    qiime rescript cull-seqs \
        --i-sequences "$merged_s" \
        --p-num-degenerates 5 --p-homopolymer-length 12 --p-n-jobs "$THREADS" \
        --o-clean-sequences "$REFDB/tmp_rbcL_culled_seqs.qza" \
        2>&1 | tee -a "$LOG/prep_rbcL.log"
    qiime rescript dereplicate \
        --i-sequences "$REFDB/tmp_rbcL_culled_seqs.qza" \
        --i-taxa "$merged_t" \
        --p-mode uniq --p-threads "$THREADS" \
        --o-dereplicated-sequences "$s_out" \
        --o-dereplicated-taxa "$t_out" \
        2>&1 | tee -a "$LOG/prep_rbcL.log"
    rm -f "$REFDB"/tmp_rbcL_*.qza
    ok "$s_out" && ok "$t_out"
}

prepare_refs() {
    echo "=== Preparing reference databases ==="
    prep_rbcl \
        && { REF_SEQS[m2]="$REFDB/ncbi_rbcL_seqs.qza"; REF_TAX[m2]="$REFDB/ncbi_rbcL_tax.qza"; }
    prep_simple m3 "$RAWDB/midori2_${MIDORI2_VERSION}_CO1_seqs.qza" \
                   "$RAWDB/midori2_${MIDORI2_VERSION}_CO1_tax.qza" "midori2_CO1" 1 \
        && { REF_SEQS[m3]="$REFDB/midori2_CO1_seqs.qza"; REF_TAX[m3]="$REFDB/midori2_CO1_tax.qza"; }
    prep_simple m4 "$RAWDB/midori2_${MIDORI2_VERSION}_srRNA_seqs.qza" \
                   "$RAWDB/midori2_${MIDORI2_VERSION}_srRNA_tax.qza" "midori2_srRNA" 1 \
        && { REF_SEQS[m4]="$REFDB/midori2_srRNA_seqs.qza"; REF_TAX[m4]="$REFDB/midori2_srRNA_tax.qza"; }
    prep_simple m5 "$RAWDB/unite_${UNITE_VERSION}_euk99_seqs.qza" \
                   "$RAWDB/unite_${UNITE_VERSION}_euk99_tax.qza" "unite_euk99" 0 \
        && { REF_SEQS[m5]="$REFDB/unite_euk99_seqs.qza"; REF_TAX[m5]="$REFDB/unite_euk99_tax.qza";
             REF_SEQS[m7]="$REFDB/unite_euk99_seqs.qza"; REF_TAX[m7]="$REFDB/unite_euk99_tax.qza"; }
    prep_simple m6 "$RAWDB/pr2_${PR2_VERSION}_seqs.qza" \
                   "$RAWDB/pr2_${PR2_VERSION}_tax.qza" "pr2" 0 \
        && { REF_SEQS[m6]="$REFDB/pr2_seqs.qza"; REF_TAX[m6]="$REFDB/pr2_tax.qza"; }

    echo "=== Reference DBs available ==="
    for m in "${MARKERS_TAXO[@]}"; do
        printf "%s %-18s seqs=%s\n" "$m" "${MARKER_NAME[$m]}" "${REF_SEQS[$m]:-MISSING}"
    done
}

#####################################################
# 3. CLASSIFICATION + controls + exports
#####################################################
classify_marker() {
    local marker="$1"
    local REP="$WORKDIR/dada2/rep-seqs_${marker}.qza"
    local TAXO="$WORKDIR/taxonomy/taxonomy_${marker}.qza"
    local SEARCH="$WORKDIR/taxonomy/search_${marker}.qza"
    local RS="${REF_SEQS[$marker]:-}" RT="${REF_TAX[$marker]:-}"

    ok "$REP" || { echo "WARNING: $marker no rep-seqs ($REP)" >&2; return 1; }
    [[ -n "$RS" && -n "$RT" ]] && ok "$RS" && ok "$RT" \
        || { echo "WARNING: $marker reference DB unavailable" >&2; return 1; }

    if [[ "$SKIP_EXISTING" -eq 1 ]] && ok "$TAXO"; then
        echo "$marker taxonomy already exists, skipped."
    else
        echo "=== vsearch classification $marker (${MARKER_NAME[$marker]}) id=${PERC_ID[$marker]} ==="
        qiime feature-classifier classify-consensus-vsearch \
            --i-query "$REP" \
            --i-reference-reads "$RS" \
            --i-reference-taxonomy "$RT" \
            --p-perc-identity "${PERC_ID[$marker]}" \
            --p-query-cov "${QUERY_COV[$marker]}" \
            --p-min-consensus "$MIN_CONSENSUS" \
            --p-maxaccepts "$MAXACCEPTS" \
            --p-top-hits-only \
            --p-strand both \
            --p-threads "$THREADS" \
            --o-classification "$TAXO" \
            --o-search-results "$SEARCH" \
            --verbose \
            2>&1 | tee "$LOG/taxonomy_${marker}.log"
    fi
    ok "$TAXO" || { echo "ERROR: classification failed for $marker" >&2; return 1; }

    qiime metadata tabulate --m-input-file "$TAXO" \
        --o-visualization "$WORKDIR/taxonomy/taxonomy_${marker}.qzv" || true

    # vsearch hits (blast6: % identity per ASV) -> useful to filter species-level calls
    if ok "$SEARCH"; then
        rm -rf "$WORKDIR/taxonomy/search_${marker}_export"
        qiime tools export --input-path "$SEARCH" \
            --output-path "$WORKDIR/taxonomy/search_${marker}_export" || true
    fi
    return 0
}

controls_barplot() {
    local marker="$1"
    local TABLE="$WORKDIR/dada2/table_${marker}.qza"
    local TAXO="$WORKDIR/taxonomy/taxonomy_${marker}.qza"
    local CTRL="$WORKDIR/decontam/controls_table_${marker}.qza"
    ok "$TABLE" && ok "$TAXO" || return 0
    if ! ok "$CTRL"; then
        qiime feature-table filter-samples \
            --i-table "$TABLE" --m-metadata-file "$META" \
            --p-where "[marker]='${marker}' AND [sample-or-control]='control'" \
            --o-filtered-table "$CTRL" \
            2>&1 | tee -a "$LOG/controls_${marker}.log" || true
    fi
    ok "$CTRL" && qiime taxa barplot \
        --i-table "$CTRL" --i-taxonomy "$TAXO" --m-metadata-file "$META" \
        --o-visualization "$WORKDIR/decontam/controls_barplot_${marker}.qzv" \
        2>&1 | tee -a "$LOG/controls_${marker}.log" || true
}

export_marker() {
    local marker="$1"
    local TAXO="$WORKDIR/taxonomy/taxonomy_${marker}.qza"
    local TABLE_FINAL="$WORKDIR/decontam/table_${marker}_final.qza"
    local REP_FINAL="$WORKDIR/decontam/rep-seqs_${marker}_final.qza"
    local TABLE_RAW="$WORKDIR/dada2/table_${marker}.qza"
    local REP_RAW="$WORKDIR/dada2/rep-seqs_${marker}.qza"
    local EXPORT_DIR="$WORKDIR/exports/${marker}_${MARKER_NAME[$marker]}_final"
    mkdir -p "$EXPORT_DIR"
    ok "$TAXO" || return 1

    local TABLE_TO_USE="$TABLE_RAW" REP_TO_USE="$REP_RAW"
    if ok "$TABLE_FINAL" && ok "$REP_FINAL"; then
        TABLE_TO_USE="$TABLE_FINAL"; REP_TO_USE="$REP_FINAL"
    fi
    ok "$TABLE_TO_USE" || { echo "$marker: no table, export skipped."; return 1; }
    echo "$marker: export from $(basename "$TABLE_TO_USE")"

    qiime taxa barplot \
        --i-table "$TABLE_TO_USE" --i-taxonomy "$TAXO" --m-metadata-file "$META" \
        --o-visualization "$EXPORT_DIR/taxa-barplot_${marker}.qzv" || true

    rm -rf "$EXPORT_DIR/table_export" "$EXPORT_DIR/repseq_export" "$EXPORT_DIR/taxonomy_export"
    qiime tools export --input-path "$TABLE_TO_USE" --output-path "$EXPORT_DIR/table_export"
    qiime tools export --input-path "$REP_TO_USE"   --output-path "$EXPORT_DIR/repseq_export"
    qiime tools export --input-path "$TAXO"         --output-path "$EXPORT_DIR/taxonomy_export"

    if ok "$EXPORT_DIR/table_export/feature-table.biom" && ok "$EXPORT_DIR/taxonomy_export/taxonomy.tsv"; then
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

    # copy vsearch hits next to the ASV table
    local hits
    hits=$(find "$WORKDIR/taxonomy/search_${marker}_export" -name "*.tsv" 2>/dev/null | head -n1)
    [[ -n "$hits" ]] && cp "$hits" "$EXPORT_DIR/vsearch_hits_${marker}_blast6.tsv"

    # quick assignment summary
    local tx="$EXPORT_DIR/taxonomy_export/taxonomy.tsv"
    local n_tot n_un
    n_tot=$(tail -n +2 "$tx" | wc -l)
    n_un=$(tail -n +2 "$tx" | awk -F'\t' '$2 ~ /^Unassigned/' | wc -l)
    printf "%s\t%s\t%s\t%s\n" "$marker" "${MARKER_NAME[$marker]}" "$n_tot" "$n_un" >> "$SUMMARY"
}

#####################################################
# MAIN
#####################################################
case "$MODE" in
    download)
        run_downloads
        exit $?
        ;;
    classify|all)
        if [[ "$MODE" == "all" ]]; then
            if has_internet; then run_downloads
            else echo "No internet on $(hostname): using already downloaded DBs in $RAWDB"; fi
        fi
        ;;
    *)
        echo "Usage: $0 [download|classify|all]" >&2; exit 1 ;;
esac

ok "$META" || { echo "ERROR: metadata missing: $META (run 01_pipeline first)" >&2; exit 1; }

prepare_refs

SUMMARY="$WORKDIR/taxonomy/taxonomy_summary_m2_m7.tsv"
printf "marker\tmarker_name\tn_ASV\tn_Unassigned\n" > "$SUMMARY"

for marker in "${MARKERS_TAXO[@]}"; do
    classify_marker "$marker" || continue
    controls_barplot "$marker"
    export_marker "$marker"
done

echo "=== Assignment summary ==="
column -t -s $'\t' "$SUMMARY"
echo "=== Taxonomy m2-m7 finished ==="
