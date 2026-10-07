# bash

module load conda/4.12.0
source ~/.bashrc
conda activate rachis-qiime2-2026.7

WORKDIR="/home/plstenge/CapoSagro_metabarcoding/20260916_AV241601_E1739-Ps12-Testscustom_recipe_15092026"
METADATA="${WORKDIR}/metadata/sample-metadata.tsv"
OUTDIR="${WORKDIR}/rarefaction"
mkdir -p "$OUTDIR"

declare -A MIN_READS
MIN_READS[m1]=47000
MIN_READS[m2]=16000 # suppression de 1069-a1-m2
MIN_READS[m3]=17000 # suppression de 1069-a2-m3, 1070-a3-m3, 1070-a2-m3, 1170-ab3-m3
MIN_READS[m4]=19900 # suppression de 1069-a1-m4
MIN_READS[m5]=96000 # suppression de 1070-a3-m5 et 1070-a2-m5
MIN_READS[m6]=26000 # suppression de 1069-a3-m6, 1069-a2-m6 et 1069-a1-m6
MIN_READS[m7]=20000 # suppression de 1069-a1-m7 et 1069-a2-m7 et 1070-a2-m7 et 1070-a3-m7 et 1069-a3-m7 et 1177-a3-m7

for marker in m1 m2 m3 m4 m5 m6 m7; do
    TABLE="${WORKDIR}/decontam/table_${marker}_final.qza"
    DEPTH="${MIN_READS[$marker]}"
    FILTERED="${OUTDIR}/table_${marker}_min${DEPTH}.qza"
    QZV="${OUTDIR}/alpha-rarefaction_${marker}_min${DEPTH}.qzv"

    if [[ ! -f "$TABLE" ]]; then
        echo "WARNING: table absente : ${TABLE}" >&2
        continue
    fi

    echo "=== ${marker}: exclusion des échantillons < ${DEPTH} lectures ==="

    qiime feature-table filter-samples \
      --i-table "$TABLE" \
      --p-min-frequency "$DEPTH" \
      --o-filtered-table "$FILTERED"

    echo "=== ${marker}: alpha-rarefaction jusqu'à ${DEPTH} ==="

    qiime diversity alpha-rarefaction \
      --i-table "$FILTERED" \
      --p-min-depth 1 \
      --p-max-depth "$DEPTH" \
      --m-metadata-file "$METADATA" \
      --o-visualization "$QZV"
done

echo "=== Tables filtrées et alpha-rarefactions terminées ==="


module load conda/4.12.0
source ~/.bashrc
conda activate rachis-qiime2-2026.7

set -euo pipefail

WORKDIR="/home/plstenge/CapoSagro_metabarcoding/20260916_AV241601_E1739-Ps12-Testscustom_recipe_15092026"
OUTDIR="${WORKDIR}/rarefied_ASV_tables"
MARKERS=(m1 m2 m3 m4 m5 m6 m7)

declare -A MARKERNAME
MARKERNAME[m1]="V4V5"
MARKERNAME[m2]="rbcL"
MARKERNAME[m3]="mlCOI_jgHCO"
MARKERNAME[m4]="teleo"
MARKERNAME[m5]="marine_fungi_ITS2"
MARKERNAME[m6]="V3V4_18S"
MARKERNAME[m7]="ITS2"

declare -A DEPTH
DEPTH[m1]=47000
DEPTH[m2]=16000
DEPTH[m3]=17000
DEPTH[m4]=19900
DEPTH[m5]=96000
DEPTH[m6]=26000
DEPTH[m7]=20000

mkdir -p "$OUTDIR"
printf "marker\tmarker_name\trarefaction_depth\tsamples_before\tsamples_retained\tsamples_excluded\n" \
    > "${OUTDIR}/rarefaction_summary.tsv"

for marker in "${MARKERS[@]}"; do
    TABLE="${WORKDIR}/decontam/table_${marker}_final.qza"
    TAXONOMY="${WORKDIR}/taxonomy/taxonomy_${marker}.qza"
    D="${DEPTH[$marker]}"
    LABEL="${MARKERNAME[$marker]}"

    MARKEROUT="${OUTDIR}/${marker}_${LABEL}_depth${D}"
    RAREFIED="${MARKEROUT}/table_${marker}_rarefied_${D}.qza"
    TABLE_EXPORT="${MARKEROUT}/table_export"
    TAX_EXPORT="${MARKEROUT}/taxonomy_export"
    TSV_RAW="${MARKEROUT}/ASV_table_${marker}_${LABEL}_rarefied_${D}.tsv"
    TSV_TAX="${MARKEROUT}/ASV_table_${marker}_${LABEL}_rarefied_${D}_taxonomy.tsv"
    BEFORE="${MARKEROUT}/sample_depths_before.tsv"
    AFTER="${MARKEROUT}/sample_depths_rarefied.tsv"
    EXCLUDED="${MARKEROUT}/excluded_samples_below_${D}.tsv"

    if [[ ! -f "$TABLE" ]]; then
        echo "WARNING: table absente : ${TABLE}" >&2
        continue
    fi

    mkdir -p "$MARKEROUT"
    rm -rf "$TABLE_EXPORT" "$TAX_EXPORT"

    echo "======================================================"
    echo "Marqueur : ${marker} (${LABEL})"
    echo "Profondeur de rarefaction : ${D}"
    echo "======================================================"

    # Profondeurs avant rarefaction et liste explicite des échantillons exclus.
    qiime tools export \
        --input-path "$TABLE" \
        --output-path "${MARKEROUT}/before_export"

    python - "${MARKEROUT}/before_export/feature-table.biom" "$BEFORE" <<'PY'
import sys
import biom

biom_file, output_file = sys.argv[1], sys.argv[2]
table = biom.load_table(biom_file)

with open(output_file, "w") as out:
    out.write("sample-id\tdepth_before_rarefaction\n")
    for sample_id, depth in zip(table.ids(axis="sample"), table.sum(axis="sample")):
        out.write(f"{sample_id}\t{int(depth)}\n")
PY

    awk -F'\t' -v depth="$D" \
        'BEGIN {OFS="\t"} NR == 1 {print $0, "rarefaction_depth"; next}
         $2 < depth {print $0, depth}' \
        "$BEFORE" > "$EXCLUDED"

    # Création de la table rarefiée.
    qiime feature-table rarefy \
        --i-table "$TABLE" \
        --p-sampling-depth "$D" \
        --o-rarefied-table "$RAREFIED"

    # Export BIOM puis conversion en table TSV.
    qiime tools export \
        --input-path "$RAREFIED" \
        --output-path "$TABLE_EXPORT"

    biom convert \
        -i "${TABLE_EXPORT}/feature-table.biom" \
        -o "$TSV_RAW" \
        --to-tsv

    # Contrôle : toutes les bibliothèques conservées doivent totaliser D lectures.
    python - "${TABLE_EXPORT}/feature-table.biom" "$AFTER" <<'PY'
import sys
import biom

biom_file, output_file = sys.argv[1], sys.argv[2]
table = biom.load_table(biom_file)

with open(output_file, "w") as out:
    out.write("sample-id\tdepth_after_rarefaction\n")
    for sample_id, depth in zip(table.ids(axis="sample"), table.sum(axis="sample")):
        out.write(f"{sample_id}\t{int(depth)}\n")
PY

    # Ajout de la taxonomie dans une seconde table TSV, uniquement si disponible.
    if [[ -f "$TAXONOMY" ]]; then
        qiime tools export \
            --input-path "$TAXONOMY" \
            --output-path "$TAX_EXPORT"

        biom add-metadata \
            -i "${TABLE_EXPORT}/feature-table.biom" \
            -o "${MARKEROUT}/feature-table-with-taxonomy.biom" \
            --observation-metadata-fp "${TAX_EXPORT}/taxonomy.tsv" \
            --observation-header "Feature ID,taxon,confidence" \
            --sc-separated "taxon"

        biom convert \
            -i "${MARKEROUT}/feature-table-with-taxonomy.biom" \
            -o "$TSV_TAX" \
            --to-tsv \
            --header-key taxon
    else
        echo "WARNING: taxonomie absente pour ${marker} : ${TAXONOMY}" >&2
    fi

    N_BEFORE=$(( $(wc -l < "$BEFORE") - 1 ))
    N_AFTER=$(( $(wc -l < "$AFTER") - 1 ))
    N_EXCLUDED=$(( N_BEFORE - N_AFTER ))

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$marker" "$LABEL" "$D" "$N_BEFORE" "$N_AFTER" "$N_EXCLUDED" \
        >> "${OUTDIR}/rarefaction_summary.tsv"

    rm -rf "${MARKEROUT}/before_export"

    echo "OK : ${TSV_RAW}"
    [[ -f "$TSV_TAX" ]] && echo "OK : ${TSV_TAX}"
done

echo
echo "Terminé."
echo "Résumé : ${OUTDIR}/rarefaction_summary.tsv"



