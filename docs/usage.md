# Pipeline Usage Guide

Dokumentasi ini menjelaskan cara mengeksekusi pipeline Snakemake di berbagai environment (Lokal, HPC Single Node, dan HPC dengan SLURM).

## 1. Persiapan (Prerequisites)
Pastikan Anda telah menginstal `conda` (disarankan menggunakan `mamba` atau `micromamba` untuk instalasi yang lebih cepat) dan `snakemake`.

```bash
# Membuat environment khusus untuk menjalankan snakemake
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
```

### Kunci API NCBI

`rule download_genomes` membaca `NCBI_API_KEY` dari environment, bukan dari file
`.env`. Tanpa kunci ini rate limit NCBI turun dari 10 ke 3 request per detik dan
unduhan bisa gagal pada dataset besar:

```bash
export NCBI_API_KEY="<kunci-anda>"
```

### Sumber data GISAID (opsional)

GISAID tidak punya API bulk terbuka. Unduh manual dari antarmuka web dengan opsi
**compl CDS** dan **high cov.** dicentang, lalu catat path FASTA dan metadata di
`gisaid_sources` pada `config/config.yaml`. Sekuens GISAID **tidak boleh**
diredistribusikan; bagikan daftar accession `EPI_ISL_...` beserta acknowledgement
laboratorium, bukan alignment-nya.

## 2. Konfigurasi (Config)
Semua parameter dieksekusi berdasarkan satu file: `config/config.yaml`.
Sebelum menjalankan pipeline, pastikan Anda:
1. Menyetel `project_name` (untuk memisahkan output folder antar eksperimen).
2. Menyetel `max_sequences_per_group` (sesuaikan dengan kekuatan komputasi Anda).
3. Mengecek alokasi `resources` (threads dan mem_mb) di bagian bawah konfigurasi.

---

## 3. Menjalankan di Lokal / HPC Tanpa SLURM (Single Node)
Gunakan mode ini jika Anda menjalankan pipeline di laptop, PC, atau telah mem- *booking* satu node penuh di HPC.

```bash
snakemake --use-conda -c 64 --resources hyphy_jobs=4 iqtree_jobs=2
```
**Penjelasan:**
*   `--use-conda`: Memerintahkan Snakemake untuk otomatis membuat *environment* independen per- *rule* (misal: env khusus HyPhy, env khusus IQ-TREE) di *background*.
*   `-c 64`: Mengizinkan Snakemake menggunakan maksimal 64 *cores*. Snakemake akan otomatis membagikan *core* tersebut ke aplikasi di bawahnya berdasarkan nilai `threads` yang ada di `config.yaml` tanpa menyebabkan *oversubscription* / *crash*.

---

## 4. Menjalankan di HPC dengan SLURM (Multi Node)
Jika kluster HPC Anda menggunakan penjadwal pekerjaan SLURM, Snakemake dapat bertindak sebagai *master* yang secara otomatis men-*submit* skrip `sbatch` untuk setiap tahapan.

**Langkah 1:** Instal plugin executor SLURM untuk Snakemake:
```bash
pip install snakemake-executor-plugin-slurm
```

**Langkah 2:** Eksekusi pipeline dengan argumen SLURM:
```bash
snakemake --use-conda --executor slurm --jobs 100
```
**Penjelasan:**
*   `--executor slurm`: Memberi tahu Snakemake untuk melempar pekerjaan ke SLURM, bukan dijalankan di node yang sedang aktif.
*   `--jobs 100`: Membatasi agar Snakemake tidak men-*submit* lebih dari 100 antrean SLURM sekaligus.
*   *Magic*: Snakemake akan membaca variabel `threads: 16` dan `mem_mb: 16000` dari `config.yaml` dan secara transparan mengonversinya menjadi *SBATCH directive* (`#SBATCH --cpus-per-task=16 --mem=16000M`) saat mengirim *job* ke antrean SLURM.

---

## 5. Output & Benchmarking
Seluruh hasil akan otomatis masuk ke folder `results/<project_name>/`.
Jika ada *error*, buka file `.log` di dalam `results/<project_name>/logs/`.
Untuk penggunaan RAM dan waktu tiap program, periksa `.tsv` di `results/<project_name>/benchmarks/`.

### Menandai kegagalan yang tidak menghentikan pipeline

Beberapa rule menulis sentinel `{"status": "FAILED"}` ketika tool-nya gagal, agar
DAG tetap berjalan tanpa menyamarkan kegagalan sebagai hasil negatif. Setelah run
selesai, periksa:

```bash
grep -l '"status": "FAILED"' results/<project_name>/04_selection/*/*.json
awk 'FNR==2 && $3!="OK"' results/<project_name>/04_selection/*/*_gard_summary.txt
```

Keluaran kosong berarti tidak ada tahap yang gagal diam-diam.
