# Runbook: menjalankan ulang pipeline di HPC

Branch: `fix/pipeline-correctness`
Perkiraan waktu: ~1 jam membangun environment, ~4–5 jam menjalankan pipeline.

Jalankan blok per blok. Setiap tahap punya pemeriksaan; **jangan lanjut kalau
pemeriksaannya gagal.**

---

## 0. Sebelum mulai — yang perlu disiapkan

Satu hal tidak ada di repo dan harus kamu bawa manual: folder data GISAID.
Sekuens GISAID tidak boleh diredistribusikan, jadi sengaja tidak di-commit.

Dari laptop, salin ke HPC:

Folder itu harus mendarat **di dalam repo**, bukan di home, karena
`config/config.yaml` menunjuknya dengan path relatif (`.dev/gisaid_epiniv/...`)
dan Snakemake dijalankan dari root repo.

Lakukan langkah 1 (clone) lebih dulu, lalu dari laptop:

```bash
# jalankan di laptop, bukan di HPC
# ganti <repo-di-hpc> dengan path repo hasil clone, mis:
#   /datadrive/drive_a/engkinandatama/h2h
scp -r ~/projects/h2h-genomic-signatures/.dev/gisaid_epiniv \
    <user>@<hpc-host>:<repo-di-hpc>/.dev/
```

Kalau `.dev` belum ada di HPC, buat dulu di sana: `mkdir -p <repo-di-hpc>/.dev`.

Kalau kamu memilih menjalankan tanpa GISAID, lewati langkah ini dan kosongkan
blok `gisaid_sources` di config nanti. Pipeline tetap jalan, hanya kehilangan
sekitar 6–8 sekuens India dan babi Malaysia.

---

## 1. Clone dan siapkan environment

```bash
cd ~
git clone git@github.com:engkinandatama/h2h-genomic-signatures.git h2h
cd h2h
git checkout fix/pipeline-correctness

# konfirmasi branch dan commit terakhir
git log --oneline -3
```

Buat environment untuk Snakemake sendiri (bukan environment tool-nya, itu dibuat
otomatis oleh `--use-conda`):

```bash
conda create -y -n snakemake -c conda-forge -c bioconda snakemake mamba
conda activate snakemake
snakemake --version   # harus keluar angka, bukan error
```

---

## 2. Pasang data GISAID dan kunci API

```bash
# verifikasi kedua file yang dibaca config benar-benar ada
python3 -c "
import yaml, os
c = yaml.safe_load(open('config/config.yaml'))
for tax, src in c.get('gisaid_sources', {}).items():
    for k, v in src.items():
        print(f'  {k:9s}', 'ADA   ' if os.path.exists(v) else 'HILANG', v)
"

export NCBI_API_KEY="<kunci-anda>"
echo "${NCBI_API_KEY:0:6}..."   # cek terisi, jangan cetak penuh
```

`NCBI_API_KEY` dibaca dari environment, bukan dari file `.env`. Tanpa kunci ini
laju permintaan NCBI turun dari 10 ke 3 per detik dan unduhan besar bisa gagal.

Supaya tidak hilang saat sesi terputus, simpan di `~/.bashrc` atau di skrip
sbatch-mu.

---

## 3. Periksa konfigurasi sebelum menjalankan apa pun

```bash
python3 -c "
import yaml
c = yaml.safe_load(open('config/config.yaml'))
print('viruses      :', len(c['viruses']))
print('virus_groups :', len(c['virus_groups']))
print('min_methods  :', c['params']['hyphy']['min_methods'])
print('gisaid       :', list(c.get('gisaid_sources', {}).keys()) or 'tidak dipakai')
for g, s in c['virus_groups'].items():
    cats = [c['viruses'][v]['category'] for v in s['viruses'] if v in c['viruses']]
    print(f'  {g:20s} {cats}')
"
```

Yang harus terlihat: 13 viruses, 8 virus_groups, `min_methods: 2`, dan Ebola serta
Sudan hanya punya `['H2H']` tanpa Reservoir. Kalau kamu tidak memakai GISAID,
kosongkan blok `gisaid_sources` di `config/config.yaml`.

---

## 4. Dry run

```bash
snakemake -n --quiet -s workflow/Snakefile | tail -20
```

Harus berakhir dengan `total` sekitar **520 job** dan tanpa `Error`. Kalau ada
`KeyError` atau `MissingInputException`, berhenti di sini dan kirim keluarannya.

---

## 5. Bangun environment conda lebih dulu

Pisahkan tahap ini supaya kegagalan solver tidak bercampur dengan kegagalan
analisis:

```bash
snakemake --use-conda --conda-create-envs-only -c 8 -s workflow/Snakefile
```

Sekitar 30–60 menit. Setelah selesai, kunci versinya supaya run berikutnya
identik:

```bash
for env in $(ls .snakemake/conda/*.yaml 2>/dev/null); do :; done
conda env list | head
# opsional tapi disarankan:
# conda list --explicit -p <prefix-env> > workflow/envs/<nama>.lock
```

---

## 6. Jalankan pipeline

```bash
snakemake --use-conda \
          --cores 64 \
          --resources hyphy_jobs=4 iqtree_jobs=2 \
          --rerun-incomplete \
          --keep-going \
          -s workflow/Snakefile \
          2>&1 | tee run_$(date +%Y%m%d_%H%M).log
```

`--cores` hanya mengatur jumlah slot job. IQ-TREE dan HyPhy membuat thread sendiri,
jadi `--resources` yang menjaga total tetap di dalam 64: 2×16 + 4×8 = 64.

`--keep-going` membuat satu dataset yang gagal tidak menghentikan yang lain. Itu
disengaja, dan karena itulah langkah 7 wajib.

Untuk SLURM, bungkus dengan sbatch:

```bash
sbatch --job-name=h2h --cpus-per-task=64 --mem=200G --time=12:00:00 \
       --wrap="cd ~/h2h && conda run -n snakemake snakemake --use-conda \
               --cores 64 --resources hyphy_jobs=4 iqtree_jobs=2 \
               --rerun-incomplete --keep-going -s workflow/Snakefile"
```

---

## 7. WAJIB — periksa kegagalan senyap sebelum melihat hasil

Pipeline sengaja tidak berhenti saat sebuah tool gagal; ia menulis penanda.
Kalau langkah ini dilewati, kegagalan bisa terbaca sebagai hasil negatif yang
bersih — persis yang terjadi pada run sebelumnya.

```bash
R=results/Zoonotic_Convergent_Signatures

echo "== JSON yang gagal =="
grep -l '"status": "FAILED"' $R/04_selection/*/*.json 2>/dev/null || echo "  tidak ada"

echo "== GARD yang tidak OK =="
awk -F'\t' 'FNR==2 && $3!="OK" {print FILENAME": "$3}' \
    $R/04_selection/*/*_gard_summary.txt 2>/dev/null || echo "  semua OK"

echo "== file hasil yang kosong =="
find $R -name "*.json" -size -10c 2>/dev/null | head
find $R -name "*_mapped.pdb" -size -1k 2>/dev/null | head

echo "== struktur yang tidak tersedia =="
grep -l "NO STRUCTURE AVAILABLE" $R/05_structure/*/*_mapped.pdb 2>/dev/null || echo "  tidak ada"

echo "== error di log =="
grep -rl "^Error\|Traceback" $R/logs/ 2>/dev/null | head
```

Keluaran kosong di semua bagian berarti tidak ada tahap yang gagal diam-diam.
Kalau ada isinya, catat mana yang gagal — hasil untuk dataset itu tidak boleh
dipakai sebelum penyebabnya jelas.

---

## 8. Periksa ukuran dataset yang benar-benar masuk

```bash
R=results/Zoonotic_Convergent_Signatures
for d in $R/01_raw_fasta/*/; do
  v=$(basename "$d")
  for f in "$d"*_filtered.fasta; do
    [ -e "$f" ] || continue
    printf "  %-28s %-12s %3s sekuens\n" "$v" \
      "$(basename "$f" _filtered.fasta)" "$(grep -c '^>' "$f")"
  done
done
```

Perhatikan kelompok dengan n di bawah 10 — kontras berbasis cabang pada ukuran itu
tidak punya daya, dan harus ditulis sebagai keterbatasan, bukan dilaporkan sebagai
hasil nol.

---

## 9. Baca hasilnya

```bash
R=results/Zoonotic_Convergent_Signatures
cat $R/06_statistics/convergence_summary.json
```

Empat angka yang menentukan arah manuskrip:

- `n_differential_raw_p` — jumlah situs pada ambang p mentah
- `n_differential_after_fdr` — berapa yang bertahan setelah koreksi
- `permutation_tests[].p_value` — apakah tumpang tindih bin lebih dari kebetulan
- `recombination.datasets_failed_or_missing` — harus 0

Lalu:

```bash
column -t -s$'\t' $R/06_statistics/hotspot_bins.tsv | head
column -t -s$'\t' $R/06_statistics/hotspot_permutation_test.tsv
column -t -s$'\t' $R/06_statistics/differential_site_counts.tsv
head -3 $R/05_structure/*/*_coordinate_map.tsv
```

`_coordinate_map.tsv` adalah yang menghubungkan nomor di gambar dengan nomor di
tabel. Kolom `alignment_column` dan `uniprot_residue` harus dipakai konsisten:
tabel manuskrip memakai nomor UniProt.

---

## 10. Bawa hasilnya kembali

```bash
# di HPC — arsipkan yang ringan saja, bukan seluruh direktori
tar czf h2h_results_$(date +%Y%m%d).tar.gz \
    results/Zoonotic_Convergent_Signatures/06_statistics \
    results/Zoonotic_Convergent_Signatures/07_figures \
    results/Zoonotic_Convergent_Signatures/04_selection/*/*_all_sites.tsv \
    results/Zoonotic_Convergent_Signatures/05_structure/*/*_coordinate_map.tsv \
    results/Zoonotic_Convergent_Signatures/logs

# dari laptop
scp <user>@<hpc-host>:~/h2h/h2h_results_*.tar.gz ~/Downloads/
```

Direktori penuh berukuran ratusan MB sampai beberapa GB; yang di atas cukup untuk
menilai hasil dan menyusun revisi.

---

## Kalau ada yang gagal

Jalankan ulang hanya bagian yang gagal, Snakemake akan melewati yang sudah jadi:

```bash
snakemake --use-conda --cores 64 --resources hyphy_jobs=4 iqtree_jobs=2 \
          --rerun-incomplete -s workflow/Snakefile
```

Untuk memaksa satu dataset dihitung ulang, hapus outputnya lalu jalankan lagi:

```bash
rm results/Zoonotic_Convergent_Signatures/04_selection/<grup>/<protein>_*.json
```

Untuk melihat alasan sebuah job akan dijalankan:

```bash
snakemake -n -r -s workflow/Snakefile --until <nama_rule> | head -40
```

---

## 11. Melanjutkan run yang terputus

Snakemake melacak output yang sudah jadi, jadi menjalankan ulang **tidak** memulai
dari nol — yang sudah selesai dilewati. Tapi run yang mati mendadak (server
restart, SIGKILL, koneksi putus) meninggalkan dua masalah yang harus dibereskan
lebih dulu.

### 11a. Lepas lock

Snakemake mengunci direktori kerja selama berjalan. Kalau prosesnya mati tanpa
sempat melepasnya, run berikutnya menolak start dengan pesan
`Directory cannot be locked`.

```bash
cd /datadrive/drive_a/engkinandatama/h2h
snakemake --unlock -s workflow/Snakefile
```

### 11b. Buang output yang terpotong

Proses yang dibunuh di tengah penulisan bisa meninggalkan JSON separuh. File itu
ada dan ukurannya wajar, jadi Snakemake menganggapnya selesai, tetapi parser akan
gagal membacanya. Periksa dan hapus:

```bash
python3 - <<'PY'
import json, glob, os
R = "results/Zoonotic_Convergent_Signatures"
bad = []
for f in glob.glob(f"{R}/04_selection/*/*.json"):
    try:
        json.load(open(f))
    except Exception:
        bad.append(f)
print(f"JSON tidak bisa diparse: {len(bad)}")
for b in bad:
    print("  ", b, os.path.getsize(b), "byte")
    os.remove(b)
print("Sudah dihapus; Snakemake akan menghitungnya ulang.")
PY
```

Periksa juga alignment dan tree yang mungkin terpotong:

```bash
R=results/Zoonotic_Convergent_Signatures
find $R/02_aligned $R/03_trees -type f -size -100c 2>/dev/null   # curigai yang kosong
```

### 11c. Lanjutkan

```bash
snakemake --use-conda --cores 64 \
          --resources hyphy_jobs=4 iqtree_jobs=2 \
          --rerun-incomplete --keep-going \
          -s workflow/Snakefile \
          2>&1 | tee -a run_resume_$(date +%Y%m%d_%H%M).log
```

`--rerun-incomplete` mengulang job yang ditandai Snakemake sebagai belum tuntas.
Cek dulu berapa yang tersisa sebelum menjalankan:

```bash
snakemake -n -s workflow/Snakefile | tail -5
```

---

## 12. Supaya tidak terputus lagi

Sesi interaktif mati bersama koneksi SSH dan bersama server. Pakai salah satu:

**tmux** — paling sederhana, bisa dilepas lalu disambung lagi:

```bash
tmux new -s h2h
# jalankan snakemake di dalamnya, lalu lepas dengan Ctrl-b lalu d
# menyambung kembali:
tmux attach -t h2h
```

**nohup** — tetap jalan setelah logout, tapi tidak bisa dilihat langsung:

```bash
nohup snakemake --use-conda --cores 64 \
      --resources hyphy_jobs=4 iqtree_jobs=2 \
      --rerun-incomplete --keep-going -s workflow/Snakefile \
      > run_$(date +%Y%m%d_%H%M).log 2>&1 &
echo $! > snakemake.pid
```

Pantau dengan `bash scripts/progress.sh` atau `tail -f run_*.log`.

**sbatch** — kalau ada SLURM, ini yang paling tahan restart karena scheduler akan
menjadwalkan ulang. Lihat langkah 6.

Apa pun pilihannya, pastikan `NCBI_API_KEY` ikut terbawa. Untuk nohup dan tmux,
`export` di shell yang sama sudah cukup; untuk sbatch, tulis di dalam skripnya.
