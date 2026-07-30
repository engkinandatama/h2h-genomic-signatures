#!/usr/bin/env bash
# Report how far the pipeline has got, without touching the running job.
#
# Counts real output files rather than parsing the Snakemake log, so it works
# whether the run is attached to a terminal, detached, or inside sbatch.
#
# Usage:  bash scripts/progress.sh
set -uo pipefail

R="${1:-results/Zoonotic_Convergent_Signatures}"

if [ ! -d "$R" ]; then
    echo "Belum ada direktori hasil di: $R"
    echo "Jalankan dari root repo, atau berikan path sebagai argumen."
    exit 1
fi

# Number of virus/protein datasets and per-virus jobs, read from the config so
# the expected totals stay correct when the design changes.
read -r N_DATASETS N_VIRUSES < <(python3 - <<'PY'
import yaml
c = yaml.safe_load(open("config/config.yaml"))
d = sum(len(s["proteins"]) for s in c["virus_groups"].values())
# A virus can belong to more than one group (Nipah_Reservoir serves both the
# NiV-B and the Malaysia contrast), so count unique virus/protein pairs rather
# than summing per group, which double-counted the shared reservoir.
pairs = {(x, p) for g in c["virus_groups"].values()
         for x in g["viruses"] if x in c["viruses"]
         for p in c["viruses"][x]["proteins"]}
v = len(pairs)
print(d, v)
PY
)

bar() {  # bar <done> <total>
    local d=$1 t=$2 w=28
    [ "$t" -eq 0 ] && t=1
    local f=$(( d * w / t ))
    printf "["
    printf "%${f}s" | tr ' ' '#'
    printf "%$(( w - f ))s" | tr ' ' '.'
    printf "] %3d/%-3d" "$d" "$t"
}

# Count matching files. `exclude` drops paths containing that substring, needed
# because *_fel.json also matches *_contrast_fel.json.
count() {
    local pattern=$1 exclude=${2:-}
    if [ -n "$exclude" ]; then
        ls $pattern 2>/dev/null | grep -v "$exclude" | wc -l
    else
        ls $pattern 2>/dev/null | wc -l
    fi
}

echo "=============================================================="
echo " Progres pipeline   $(date '+%Y-%m-%d %H:%M:%S')"
echo " $R"
echo "=============================================================="
echo
printf " %-26s " "1. unduh & ekstrak";  bar "$(count "$R/01_raw_fasta/*/*_filtered.fasta")" "$N_VIRUSES"; echo
printf " %-26s " "2. alignment";        bar "$(count "$R/02_aligned/*/*_codon_aligned.fasta")" "$N_DATASETS"; echo
printf " %-26s " "3. pohon filogeni";   bar "$(count "$R/03_trees/*/*_labeled.treefile")" "$N_DATASETS"; echo
echo
echo " 4. analisis seleksi (per metode):"
for m in meme fel fubar slac contrast_fel prime busted absrel relax gard; do
    if [ "$m" = "fel" ]; then
        n=$(count "$R/04_selection/*/*_fel.json" "contrast_fel")
    else
        n=$(count "$R/04_selection/*/*_${m}.json")
    fi
    printf "    %-23s " "$m"; bar "$n" "$N_DATASETS"; echo
done
echo
printf " %-26s " "5. tabel gabungan";   bar "$(count "$R/04_selection/*/*_all_sites.tsv")" "$N_DATASETS"; echo
printf " %-26s " "6. struktur 3D";      bar "$(count "$R/05_structure/*/*_3d_view.html")" "$N_DATASETS"; echo
printf " %-26s " "7. Manhattan plot";   bar "$(count "$R/07_figures/*/*_manhattan.png")" "$N_DATASETS"; echo
printf " %-26s " "8. konvergensi";      bar "$(count "$R/06_statistics/convergence_summary.json")" 1; echo
echo
echo "--------------------------------------------------------------"
echo " Sedang berjalan sekarang:"
if ! ps -eo etime=,args= 2>/dev/null | grep -E "hyphy (gard|prime|meme|fel|fubar|slac|absrel|busted|relax|contrast-fel)|iqtree|mafft" | grep -v grep | \
     awk '{printf "   %-12s %s %s\n", $1, $2, $3}' | head -12; then
    echo "   (tidak ada proses tool yang terdeteksi)"
fi
echo
echo "--------------------------------------------------------------"
echo " Job termahal yang sudah selesai:"
python3 - "$R" <<'PY'
import sys, os, glob, csv
R = sys.argv[1]
rows = []
for f in glob.glob(os.path.join(R, "benchmarks", "*", "*.tsv")):
    try:
        r = list(csv.DictReader(open(f), delimiter="\t"))
        if r:
            rows.append((float(r[0]["s"]),
                         os.path.basename(os.path.dirname(f)),
                         os.path.basename(f)[:-4]))
    except Exception:
        pass
rows.sort(reverse=True)
total = sum(r[0] for r in rows)
print(f"   {len(rows)} job selesai, total {total/3600:.1f} jam CPU-wall")
for s, rule, name in rows[:5]:
    print(f"   {s/3600:6.2f} jam  {rule}/{name}")
PY
echo
echo " Kegagalan sejauh ini:"
n_failed=$(grep -l '"status": "FAILED"' "$R"/04_selection/*/*.json 2>/dev/null | wc -l)
n_na=$(grep -l '"status": "NOT_APPLICABLE"' "$R"/04_selection/*/*.json 2>/dev/null | wc -l)
echo "   FAILED         : $n_failed   (harus 0)"
echo "   NOT_APPLICABLE : $n_na   (normal: Contrast-FEL & RELAX untuk Ebola/Sudan)"
echo "=============================================================="
