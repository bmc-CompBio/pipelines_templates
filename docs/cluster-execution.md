---
name: cluster-execution
description: Transfer a local Snakemake project to the LMU bioinformatics cluster (<cluster-login-host>), execute it via SLURM, retrieve output directories inferred from rule all, and clean up on success. Use when running computationally intensive Snakemake workflows that need cluster resources.
---

# Cluster Execution

> **Sanitized copy.** User name, login host and personal paths are replaced by placeholders: `<user>` (cluster account), `<cluster-login-host>` (login node), `<gpu-node>`, `<path-to>`. Replace them with your own values before using this file as a Claude Code skill (see the repository README).


Transfers a local Snakemake project to the LMU cluster, runs it under SLURM via snakemake-executor-plugin-cluster-generic, retrieves the output directories defined by `rule all`, and removes the project from the cluster on success. On failure, leaves everything intact for debugging.

## Cluster Details

| Parameter | Value |
|-----------|-------|
| Host | `<user>@<cluster-login-host>` |
| Base work dir | `/work/project/<user>/work` |
| Conda/mamba root | `/work/project/<user>/miniforge3/` (current installation) |
| Conda init script | `/work/project/<user>/miniforge3/etc/profile.d/conda.sh` |
| OS | CentOS 7 |

## Queue Reference

| Queue | Cores | RAM | Notes |
|-------|-------|-----|-------|
| slim16 | 32 | 256 GB | Default |
| slim18 | 36 | 256 GB | |
| fat | 48 | 1.5 TB | Memory-intensive jobs |
| gpu | 64 | 128 GB | Node `<gpu-node>`, 2× NVIDIA RTX A5000 (24GB each), driver 525.60.13 → max CUDA 12.0. GRES must be requested explicitly (see Hard Constraints and GPU ML Frameworks section below) or the job gets no GPU. |

## Environment Design Philosophy

**Default: one unified conda environment per pipeline.**

- Create ONE environment containing all tools + `snakemake-minimal` + `snakemake-executor-plugin-cluster-generic`.
- Activate it, then run `snakemake --profile profile/` — nothing else.
- Do **not** use `--use-conda` or per-rule environments by default.
- Only use `--use-conda` when there are genuine dependency conflicts between tools that cannot coexist in one environment.

**Environment yaml rules:**
- Channels: `conda-forge` and `bioconda` only — never `defaults`.
- Always pin `python=3.11`.
- Prefer `>=` version specs; avoid exact pins unless necessary (e.g., glibc-constrained tools like sra-tools).
- `channel_priority: strict` must be set in `~/.condarc` on the cluster.
- All environments live at prefix `/work/project/<user>/conda/envs/<name>` (shared storage).

## Hard Constraints

- **Never upgrade or modify the system conda/mamba installation.** If a version incompatibility is detected (e.g., conda too old for the requested snakemake version), stop and report — do not attempt to upgrade conda, mamba, or any base-environment package.
- **Never use `fasterq-dump`** on the cluster. Use `prefetch` + `fastq-dump` in separate rules instead.
- **Web access is only available on the master node** — compute nodes have no internet. `prefetch`, `wget`, `curl`, and other download tools **must** run as `localrule`. `fastq-dump` reads local files already fetched by `prefetch` and has no network requirement — it can and should run as a normal compute job. Keep `prefetch` and `fastq-dump` in **separate rules**.
- The cluster runs an older SLURM — do not use SLURM flags introduced after ~2018.
- CentOS 7 uses older glibc — keep binaries compatible.
- **Dorado is not installed natively** — it runs via apptainer container. See Dorado section below.
- **GPU jobs need explicit `--gres=gpu:N` in the sbatch submit command**, or they land on `<gpu-node>` with zero GPUs visible. `GresTypes=gpu` and `SelectTypeParameters=CR_CORE_MEMORY` are configured cluster-wide — GRES is never granted implicitly. Add a `gres` resource key to every rule that needs a GPU (`gres: "gpu:1"`) and a `gres=gpu:0` default-resource for everything else, then reference `{resources.gres}` in the profile's `cluster-generic-submit-cmd` (see Step 6 and the GPU ML Frameworks section below — the profile template shown there already includes this).

---

## Workflow

### Step 1 — Resolve project directory

Use the path passed as argument. If none given, use the current working directory.

Verify a `Snakefile` exists: check `<project_dir>/Snakefile` first, then `<project_dir>/workflow/Snakefile`. If neither exists, stop and report.

Derive `<project_name>` from the directory basename.

**The cluster working directory `/work/project/<user>/work/<project_name>/` is the canonical home for the pipeline.** All pipeline files (Snakefile, config, samples, envs, scripts) live there, and all relative paths inside those files are interpreted relative to it. Raw data files (FASTQs, BAMs, etc.) that are needed by the pipeline must also be accessible from that directory — either co-located (e.g., `fastq/`) or on a shared filesystem the compute nodes can read. If raw data exists only on the local Mac (not on an NFS share visible to the cluster), it must be rsynced to the cluster project directory in Step 5 before the pipeline can run.

### Step 2 — Parse `rule all` for output file patterns

Read the Snakefile. Locate the `rule all:` block and its `input:` section. Extract every file path string listed there and convert each to a shell glob pattern:

- For literal paths (e.g. `"results/multiqc/multiqc_report.html"`): use as-is.
- For `expand("path/{wc}/file.ext", wc=LIST)` calls: replace every `{wildcard}` with `*`, yielding e.g. `results/*/bigwig/*.bw`.
- For paths with inline wildcards (e.g. `"results/{condition}.bed"`): replace `{...}` with `*`.

Deduplicate. These glob patterns are `<output_patterns>`.

Example: rule all inputs
```python
expand("results/{sample}/bigwig/{sample}.unique.bw", sample=SAMPLES),
expand("results/{sample}/peaks/{sample}.narrowPeak", sample=DRIP_SAMPLES),
"results/consensus/{condition}.consensus_rloops.bed",
"results/multiqc/multiqc_report.html",
```
→ `<output_patterns>`:
```
results/*/bigwig/*.unique.bw
results/*/peaks/*.narrowPeak
results/consensus/*.bed
results/multiqc/multiqc_report.html
```

**Never collapse to a top-level directory.** Syncing `results/` would pull back every intermediate file (BAMs, FASTQs, trimmed reads) that is not a declared output. The glob patterns are the retrieval scope — nothing else is ever copied back.

### Step 3 — Evaluate and assign per-rule resources

Read the full Snakefile and all included `.smk` files. For each rule (excluding `rule all`), determine appropriate SLURM resources.

#### What to read per rule

- Rule name
- `shell:` or `run:` block — identify the tool(s) being called
- Existing `threads:` directive (respect it; do not lower it)
- Existing `resources:` block (any already-specified values take priority; only fill in what is missing)

#### Resource assignment heuristics

**Partition and memory — use the lowest tier that fits:**

| Tool category | Examples | Partition | mem_mb |
|---------------|----------|-----------|--------|
| De novo assembly (large genome) | SPAdes, Flye, Hifiasm, Canu, wtdbg2 | fat | 512000–1000000 |
| De novo assembly (small/microbial) | SPAdes -k, Unicycler | slim16 | 32000–64000 |
| Genome alignment index build | STAR --genomeGenerate, BWA index (large) | slim16 | 32000–64000 |
| Short-read alignment | BWA mem, Bowtie2, HISAT2 | slim16 | 8000–16000 |
| Long-read alignment | Minimap2 | slim16 | 16000–32000 |
| RNA-seq quantification | STAR 2-pass, Salmon, kallisto | slim16 | 16000–32000 |
| Variant calling | GATK HaplotypeCaller, DeepVariant, bcftools call | slim16 | 8000–16000 |
| BAM processing | samtools sort, index, markdup, view | slim16 | 4000–8000 |
| QC / trimming | FastQC, fastp, Trimmomatic, cutadapt | slim16 | 2000–4000 |
| SRA prefetch | prefetch | **localrule** (master node only — needs internet) | 2000 |
| SRA dump | fastq-dump | slim16 | 2000 |
| Text / tabular processing | awk, sort, bedtools, grep | slim16 | 1000–2000 |
| R / Python scripts (generic) | Rscript, python | slim16 | 4000–16000 (scale with expected data size) |
| MultiQC, reports | multiqc | slim16 | 2000–4000 |
| GPU tools | DeepVariant GPU, Dorado, Guppy | gpu | 16000 |
| Dorado (apptainer) | basecalling, demux | gpu | 16000 |

**Threads:**

| Tool category | Threads |
|---------------|---------|
| Assemblers | 16–32 |
| Aligners (BWA, STAR, HISAT2, Minimap2) | 8–16 |
| samtools sort / markdup | 4–8 |
| Variant callers | 4–8 |
| fastp, Trimmomatic | 4–8 |
| FastQC | 2–4 |
| SRA prefetch, text processing | 1–2 |
| SRA fastq-dump | 1–2 |
| R / Python scripts | 1–4 |
| MultiQC | 1 |

**Time (minutes):**

| Tool category | Time |
|---------------|------|
| De novo assembly (large) | 1440–4320 |
| De novo assembly (small) | 120–480 |
| Alignment (per sample) | 60–180 |
| Index build | 60–120 |
| SRA prefetch | 120–480 |
| SRA fastq-dump | 30–120 |
| QC / trimming | 30–60 |
| BAM processing | 30–60 |
| Variant calling (per sample) | 60–240 |
| R / Python scripts | 30–120 |
| Text processing | 10–30 |

#### Localrules — web-access rules must not be submitted to SLURM

Compute nodes have no internet access. Identify every rule whose shell block calls a tool that requires network access: `prefetch`, `wget`, `curl`, `rsync` to a remote host, or any other download utility. These rules **must** appear in a `localrule` (or `localrules:`) declaration in the Snakefile so Snakemake runs them on the master node.

`fastq-dump` reads SRA files already fetched by `prefetch` and does **not** need internet — it must be a separate rule and runs as a normal compute job on SLURM. If a Snakefile combines `prefetch` and `fastq-dump` in a single rule, split them before transferring.

When preparing the job (Step 3), check whether download rules are already declared as localrules. If not, report this to the user and add the missing `localrules:` line to the Snakefile before transferring.

#### Rules that already have complete resources

If a rule already declares both `threads:` and a `resources:` block covering `mem_mb`, `partition`, and `time`, skip it entirely — do not override.

#### Output

Produce a resource table for confirmation:

```
Rule              | partition | mem_mb  | threads | time (min)
------------------|-----------|---------|---------|----------
align_reads       | slim16    | 16000   | 8       | 120
assemble_genome   | fat       | 800000  | 32      | 2880
download_sra      | slim16    | 2000    | 1       | 240
...
```

Show this table to the user and explain any non-obvious choices (e.g., why fat was chosen). Wait briefly for the user to object before proceeding. If the user corrects any allocation, apply the correction.

Store the final allocations as `<rule_resources>` — they will be written into the profile config as `set-resources` and `set-threads` entries in Step 5.

### Step 4 — Verify or create the pipeline environment file

Look for, in order: `environment.yml`, `environment.yaml`, `env.yml`, `env.yaml` in the project root.

The file must satisfy:
- Channels: `conda-forge` and `bioconda` only — no `defaults`
- `python=3.11` pinned
- `snakemake-minimal=9.19.0` pinned to this exact version
- `conda>=24.7.1` included in the runner env (snakemake 9.x requires this — installing it into the runner env satisfies the check without touching the system conda). **Exception:** it cannot be solved together with `sra-tools=2.11.0` (old hdf5/libxml2 pins). Pipelines that need `sra-tools` omit it; Snakemake 9.19.0 runs fine without `--use-conda`. For the same reason `samtools` must stay `<1.17` in those environments.
- `snakemake-executor-plugin-cluster-generic` listed (installed by mamba, not pip)
- All pipeline tools present
- `sra-tools=2.11.0` pinned exactly if sra-tools is needed (glibc 2.17 on CentOS 7 cannot run ≥3.0)

Expected structure:
```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - snakemake-minimal=9.19.0
  - conda>=24.7.1
  - snakemake-executor-plugin-cluster-generic
  - <all pipeline tools, with >= version specs where possible>
```

If no file is found, stop and ask the user. If one is found but violates the rules above (wrong channels, missing snakemake components, unpinned python), fix it before proceeding.

**Known structural conflicts requiring `--use-conda`:** picard + macs2 cannot coexist in one environment. Picard's bioconda recipe bundles R scripts and pulls in `r-base`, which creates a `krb5` version conflict with macs2's numpy dependency chain. This is a bioconda packaging issue and cannot be resolved by loosening pins. When a pipeline uses both, `--use-conda` with separate per-rule envs is required — see Step 7b.

### Step 5 — Transfer the project to the cluster

Transfer pipeline files (Snakefile, config, envs, scripts, samples):

```bash
rsync -avz --progress \
  --exclude='.snakemake' \
  --exclude='__pycache__' \
  --exclude='*.egg-info' \
  --exclude='.git' \
  <project_dir>/ \
  <user>@<cluster-login-host>:/work/project/<user>/work/<project_name>/
```

Do **not** exclude output directories — partial prior results may be present and rerunnable.

**Raw data transfer:** If the pipeline reads locally-generated data (FASTQs, BAMs, etc.) that are not downloaded from public archives (SRA, ENA), those files must also be copied to the cluster. Check the samples.tsv or config for data paths. Transfer raw data to a subdirectory of the project (e.g., `fastq/`) so all paths resolve relative to the cluster project root:

```bash
rsync -avz --progress \
  <local_raw_data_dir>/ \
  <user>@<cluster-login-host>:/work/project/<user>/work/<project_name>/fastq/
```

Update any absolute paths in samples.tsv / config.yaml to the cluster paths before or after transfer.

### Step 6 — Write the SLURM profile to the cluster

Create the profile directory:

```bash
ssh <user>@<cluster-login-host> "mkdir -p /work/project/<user>/work/<project_name>/profile"
```

Write `/work/project/<user>/work/<project_name>/profile/config.yaml` via SSH heredoc. The file contains the base profile plus `set-threads` and `set-resources` sections populated from `<rule_resources>`:

```yaml
executor: cluster-generic
cluster-generic-submit-cmd: >-
  mkdir -p logs/{rule} &&
  sbatch
  --partition={resources.partition}
  --qos={resources.qos}
  --cpus-per-task={threads}
  --mem={resources.mem_mb}
  --gres={resources.gres}
  --job-name=smk-{rule}-{wildcards}
  --output=logs/{rule}/{rule}-{wildcards}-%j.out
default-resources:
  - partition=slim16
  - qos=normal
  - mem_mb=1000
  - time=60
  - gres=gpu:0
set-threads:
  - align_reads=8          # example — replace with actual rules and values
  - assemble_genome=32
set-resources:
  - align_reads:mem_mb=16000
  - align_reads:time=120
  - assemble_genome:mem_mb=800000
  - assemble_genome:partition=fat
  - assemble_genome:time=2880
  - train_model:partition=gpu
  - train_model:gres=gpu:1
  # ... one entry per rule per resource that differs from default
restart-times: 3
max-jobs-per-second: 10
max-status-checks-per-second: 1
local-cores: 1
latency-wait: 60
jobs: 200
keep-going: True
rerun-incomplete: True
printshellcmds: True
scheduler: greedy
```

The `--gres={resources.gres}` flag and the `gres=gpu:0` default are **always included**, even for pipelines with no GPU rules — `gpu:0` is a valid, harmless no-op for CPU jobs, and it means the template never needs to change between GPU and non-GPU pipelines. Only rules that actually need a GPU get a `gres=gpu:1` (or `gpu:2`, etc.) override, exactly like any other per-rule resource. Alternatively, a rule can declare `gres` directly in its own `resources:` block in the Snakefile instead of via `set-resources` — same effect, and it keeps the requirement visible next to the rule itself (this is what was done for the ML training rule below).

Rules that already declare their own `resources:` in the Snakefile take priority over `set-resources` entries at runtime — do not duplicate them unless the Snakefile values are clearly wrong.

### Step 7 — Create a fresh conda environment on the cluster

The environment prefix is `/work/project/<user>/conda/envs/smk-<project_name>` (shared storage, required).

First, discover the mamba binary path (the installation location may change):

```bash
MAMBA=$(ssh <user>@<cluster-login-host> "which mamba 2>/dev/null || echo /work/project/<user>/miniforge3/bin/mamba")
```

Then create the environment:

```bash
ssh <user>@<cluster-login-host> bash -s << 'ENVCREATE'
  MAMBA=$(which mamba 2>/dev/null || echo /work/project/<user>/miniforge3/bin/mamba)
  ENV_PREFIX=/work/project/<user>/conda/envs/smk-<project_name>
  "$MAMBA" env remove --prefix "$ENV_PREFIX" -y 2>/dev/null || true
  "$MAMBA" env create \
    --prefix "$ENV_PREFIX" \
    -f /work/project/<user>/work/<project_name>/<env_file> \
    -y
ENVCREATE
```

The environment yaml includes `snakemake-minimal` and `snakemake-executor-plugin-cluster-generic` — no pip install needed. Do not use `pip install` to add anything; it bypasses the solver and can break the environment.

If env creation fails, stop and report — do not proceed to execution.

### Step 7b — Pre-build per-rule conda environments with mamba [EXCEPTION ONLY]

**Only execute this step if `--use-conda` is explicitly required** (genuine dependency conflicts between tools, e.g. picard + macs2). The default single-env model does not use `--use-conda` and skips this step entirely.

**⚠️ $CONDA_PREFIX shadowing:** activating a conda env overwrites `$CONDA_PREFIX` to point at that env's directory. Never use `$CONDA_PREFIX` as the argument to `--conda-prefix`. Always use a separately-named variable (`SMK_CONDA_PREFIX`).

Snakemake's internal env creation hangs indefinitely on this cluster. Pre-build all per-rule environments with mamba before starting Snakemake.

**Stub location — snakemake 9.x:** writes yaml stubs to `$SMK_CONDA_PREFIX/{hash}.yaml` (not `.snakemake/conda/` as in 8.x). The script below uses `$SMK_CONDA_PREFIX` as the watch location, which is correct for 9.x. The system mamba (miniforge3) is used for env creation; the newer `conda` inside the runner env satisfies snakemake's version check but is not used for env builds.

```bash
ssh <user>@<cluster-login-host> bash -s << 'PREBUILD'
  MAMBA=$(which mamba 2>/dev/null || echo /work/project/<user>/miniforge3/bin/mamba)
  SMK_CONDA_PREFIX=/work/project/<user>/conda/envs   # NOT $CONDA_PREFIX
  cd /work/project/<user>/work/<project_name>
  source /work/project/<user>/miniforge3/etc/profile.d/conda.sh
  conda activate /work/project/<user>/conda/envs/smk-<project_name>

  # Generate hash yaml stubs — snakemake 9.x writes to $SMK_CONDA_PREFIX/{hash}.yaml
  snakemake --use-conda --conda-prefix "$SMK_CONDA_PREFIX" --conda-create-envs-only &
  SMK_PID=$!
  for i in $(seq 1 30); do
    sleep 2
    NEW=$(for f in "$SMK_CONDA_PREFIX"/*.yaml; do
      [ -f "$f" ] || continue
      prefix="${f%.yaml}"
      [ -d "${prefix}/conda-meta" ] || echo "$f"
    done)
    [ -n "$NEW" ] && kill $SMK_PID && break
  done
  kill $SMK_PID 2>/dev/null

  # Build each env that doesn't already have conda-meta/ populated
  for yaml in "$SMK_CONDA_PREFIX"/*.yaml; do
    [ -f "$yaml" ] || continue
    prefix="${yaml%.yaml}"
    [ -d "${prefix}/conda-meta" ] && continue
    "$MAMBA" env create --file "$yaml" --prefix "$prefix" -y
  done
PREBUILD
```

**Symlink shortcut:** if an env was already built under a different name/path, a symlink from the expected `{hash}_/` directory to the actual env is valid. Snakemake detects env readiness via `{hash}_/conda-meta/` presence, not by any script marker.

If any env fails to build, stop and report — do not proceed to execution.

### Step 8 — Execute Snakemake

Run via SSH, streaming output to the local terminal. The local terminal is kept open for the duration of the run.

```bash
ssh <user>@<cluster-login-host> bash -s << 'SNAKEMAKE'
  source /work/project/<user>/miniforge3/etc/profile.d/conda.sh
  conda activate /work/project/<user>/conda/envs/smk-<project_name>
  cd /work/project/<user>/work/<project_name>
  snakemake --profile profile/
SNAKEMAKE
```

Capture the exit code.

### Step 9a — Success: retrieve output and clean up

If exit code is **0**:

1. Retrieve **only** the files matching `<output_patterns>`. Expand each pattern on the cluster, then rsync each matched file individually — preserving its relative path under `<project_dir>`.

```bash
REMOTE="<user>@<cluster-login-host>"
RBASE="/work/project/<user>/work/<project_name>"
LBASE="<project_dir>"

# For each output pattern from rule all (e.g. "results/*/bigwig/*.bw"):
for pattern in \
  "results/*/bigwig/*.bw" \
  "results/*/peaks/*.narrowPeak" \
  "results/consensus/*.bed" \
  "results/multiqc/multiqc_report.html"; do

  # Expand the glob on the cluster
  ssh "$REMOTE" "ls ${RBASE}/${pattern} 2>/dev/null" | while read -r remote_file; do
    rel="${remote_file#${RBASE}/}"           # strip remote base → relative path
    local_file="${LBASE}/${rel}"
    mkdir -p "$(dirname "$local_file")"
    rsync -avz --progress "${REMOTE}:${remote_file}" "${local_file}"
  done
done
```

This pattern is general: it works for any `rule all` regardless of directory structure, never transfers files not listed as outputs, and handles per-sample wildcards automatically without knowing sample names in advance.

2. Remove the project from the cluster:

```bash
ssh <user>@<cluster-login-host> \
  "rm -rf /work/project/<user>/work/<project_name>"
```

3. Remove the conda environment:

```bash
ssh <user>@<cluster-login-host> \
  "/work/project/<user>/miniforge3/bin/mamba env remove --prefix /work/project/<user>/conda/envs/smk-<project_name> -y"
```

Report success and which directories were retrieved.

### Step 9b — Failure: stop and preserve

If exit code is **non-zero**:

- Do **not** rsync anything back.
- Do **not** remove anything from the cluster.
- Report the exit code and the preserved cluster path:
  `<user>@<cluster-login-host>:/work/project/<user>/work/<project_name>/`
- Inform the user that the cluster state is preserved and they can ask for debugging help. Snakemake logs are in `logs/` and `.snakemake/log/` within the project directory on the cluster.

---

## Debugging

**Debug outside Snakemake first.** Activate the pipeline environment and run the failing tool manually with the exact same inputs and parameters. This isolates whether the problem is in the tool, the data, or Snakemake/SLURM plumbing.

```bash
ssh <user>@<cluster-login-host>
source /work/project/<user>/miniforge3/etc/profile.d/conda.sh
conda activate /work/project/<user>/conda/envs/smk-<project_name>
cd /work/project/<user>/work/<project_name>
# run the failing command directly
```

Once the tool works manually, then inspect Snakemake/SLURM logs:

- `.snakemake/log/` — main snakemake log
- `logs/<rule>/` — per-rule SLURM stdout/stderr
- `squeue -u <user>` — check for stuck jobs
- `sacct -u <user> --format=JobID,JobName,State,ExitCode,Elapsed` — recent job history

To check a specific SLURM job: `scontrol show job <jobid>`

## Common Issues

- **conda activate fails in non-interactive SSH**: always source the conda init script first (`source .../etc/profile.d/conda.sh`) before calling `conda activate`
- **mamba binary path may change**: always discover mamba at runtime with `$(which mamba 2>/dev/null || echo /work/project/<user>/miniforge3/bin/mamba)` rather than hardcoding the path
- **never modify the system conda/mamba installation**: if a version conflict is found, stop and report — do not upgrade conda, mamba, or any base-env package
- **fasterq-dump**: replace with `prefetch <accession>` (localrule) + `fastq-dump --split-files <accession>` (separate compute rule)
- **No internet on compute nodes**: `prefetch`, `wget`, `curl`, and other download tools will fail on compute nodes — declare them as `localrules`. `fastq-dump` does not need internet (it reads local SRA files) and runs fine as a compute job; do not make it a localrule. Keep `prefetch` and `fastq-dump` in separate rules. Both rules need a `conda:` directive even though `fastq-dump` is a compute job.
- **picard + macs2 cannot share one environment**: bioconda's picard recipe bundles R scripts that pull in `r-base`, creating a `krb5` conflict with macs2's numpy chain. This is structural, not fixable by loosening pins — use `--use-conda` with separate per-rule envs.
- **Old SLURM**: avoid `--parsable`, `--array` with step syntax, or `srun` flags added after SLURM 17.x
- **OOM kill**: if a job is killed by the OOM killer (exit code 137), increase `mem_mb` for that rule — move to fat if needed
- **set-resources ignored**: Snakemake only applies `set-resources` if the rule does not already declare a matching resource in its `resources:` block — check the Snakefile if overrides seem to have no effect
- **sra-tools must be pinned to 2.11.0**: pin as `sra-tools=2.11.0` in the environment yaml. sra-tools ≥3.0 requires glibc 2.25+; CentOS 7 ships glibc 2.17.
- **never use the `defaults` channel**: always use `conda-forge` + `bioconda` only. `defaults` causes solver conflicts with bioconda packages and produces non-reproducible environments.
- **channel_priority must be strict**: if env creation produces unexpected package versions or solver conflicts, check that `~/.condarc` on the cluster contains `channel_priority: strict`.
- **do not pip install into conda envs**: use the conda yaml for everything. pip bypasses the solver and causes silent dependency breakage.
- **--conda-frontend mamba is deprecated**: Snakemake 8+ removed this setting — it prints a warning and has no effect. Do not include it.
- **do not use --use-conda by default**: the single unified env model is simpler, more debuggable, and avoids the cluster-specific hang in Snakemake's internal env creation. Only switch to `--use-conda` when tools have genuine conflicts that cannot be resolved in one environment.
- **Dorado is not a conda package — use apptainer**: `apptainer` is not on the default PATH; load it first with `module load apptainer/1.2.5`. The container image is at `<path-to>/dorado.sif`. Run dorado as: `module load apptainer/1.2.5 && apptainer exec --nv <path-to>/dorado.sif dorado <args>`. The `--nv` flag is required to pass GPU access into the container. Do not use `singularity` directly (it exists at `/usr/bin/singularity` but apptainer via module is the correct path). Dorado rules must target the `gpu` partition.
- **conda-forge TensorFlow/PyTorch GPU builds need `CONDA_OVERRIDE_CUDA` set when solving from the login node**: the login node has no NVIDIA driver, so conda's `__cuda` virtual package (used to gate GPU package variants) is undetectable there and any GPU build gets rejected as "requires __cuda, which is missing on the system" — even though the *target* GPU node has a perfectly good driver. Fix: `export CONDA_OVERRIDE_CUDA=12.0` (match the gpu node's driver ceiling, see GPU ML Frameworks section) before `mamba env create`/`mamba install`.
- **match the conda CUDA build to the driver, not just to what's newest**: `<gpu-node>`'s driver (525.60.13) tops out at CUDA 12.0. conda-forge TensorFlow builds are tagged by CUDA version in the build string (e.g. `cuda120py311...`, `cuda126py311...`) — installing a `cuda126`/`cuda128` build on this node fails or silently falls back to CPU. Always pick the build whose CUDA tag is `<=` the driver's max (check via `srun --partition=gpu --gres=gpu:1 --time=00:02:00 nvidia-smi`), not simply the latest available.
- **snakemake-minimal and TensorFlow can conflict on `wrapt`**: `snakemake-minimal>=9` pulls in `snakemake-interface-storage-plugins`, which requires `wrapt>=1.15`; TensorFlow `<=2.15` caps `wrapt<1.15` — genuinely unsatisfiable together in one env. TensorFlow `>=2.16` relaxed this cap. If a GPU-TensorFlow env fails to solve with a `wrapt` conflict in the error tree, bump TensorFlow to `2.16`/`2.17`/`2.18.x` (pick the newest version that still has a CUDA build matching the driver — see above) rather than downgrading snakemake.

---

## Dorado (GPU basecalling via apptainer)

Dorado is available as an apptainer/Singularity container, not as a conda or system package.

**Container path:** `<path-to>/dorado.sif`  
**Version:** 1.0.0+4a76f8aa (as of 2026-07)

**Required setup in every shell block that calls dorado:**

```bash
module load apptainer/1.2.5
export APPTAINERENV_CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-}"
apptainer exec --nv <path-to>/dorado.sif dorado <subcommand> <args>
```

**Snakemake rule template:**

```python
rule dorado_basecall:
    input:  pod5 = "pod5/{sample}/"
    output: bam  = "basecalls/{sample}.bam"
    resources:
        partition = "gpu",
        gres      = "gpu:1",   # required — see Hard Constraints; omitting this gets zero GPUs
        mem_mb    = 32000,
        time      = 480
    threads: 4
    shell:
        """
        module load apptainer/1.2.5
        export APPTAINERENV_CUDA_VISIBLE_DEVICES="${{CUDA_VISIBLE_DEVICES:-}}"
        apptainer exec --nv <path-to>/dorado.sif \
            dorado basecaller \
            --device cuda:all \
            sup,5mC_5hmC \
            {input.pod5} > {output.bam}
        """
```

---

## GPU ML Frameworks (TensorFlow/PyTorch via conda)

Unlike Dorado, TensorFlow and PyTorch install as ordinary conda-forge packages — no container needed. Getting real GPU acceleration (not a silent CPU fallback) requires getting three things right together: the CUDA build tag, solving from a GPU-visible context, and any transitive dependency conflicts with Snakemake itself. Validated end-to-end on `<gpu-node>` (2× RTX A5000, driver 525.60.13, CUDA 12.0 max) with `tensorflow=2.17.0=cuda120*`.

**1. Check the actual driver ceiling before picking a package version** — don't assume the newest build is installable:

```bash
srun --partition=gpu --gres=gpu:1 --time=00:02:00 nvidia-smi   # shows "CUDA Version: 12.0" etc.
```

**2. Find a conda-forge build matching that ceiling, not the newest one:**

```bash
mamba search -c conda-forge 'tensorflow=2.17.0' | grep cuda120
```

conda-forge tags GPU builds by CUDA version in the build string (`cuda120py311...`, `cuda126py311...`, ...). A `cuda126`/`cuda128` build will fail or silently fall back to CPU on a driver that only supports CUDA 12.0.

**3. Set `CONDA_OVERRIDE_CUDA` before solving** — the login node (where you run `mamba env create`) has no GPU, so conda can't auto-detect the `__cuda` virtual package, and any GPU build gets rejected as "requires `__cuda`, which is missing on the system" even though the target GPU node is fine:

```bash
export CONDA_OVERRIDE_CUDA=12.0   # match step 1's driver ceiling
mamba env create --prefix /work/project/<user>/conda/envs/<env> -f environment.yml -y
```

**4. Watch for a `wrapt` conflict with `snakemake-minimal`:** `snakemake-minimal>=9` requires `wrapt>=1.15` (via `snakemake-interface-storage-plugins`); TensorFlow `<=2.15` caps `wrapt<1.15`. These are genuinely unsatisfiable in one environment. TensorFlow `2.16+` relaxed the cap — if the solver reports a `wrapt` conflict, bump the TensorFlow version (staying within the CUDA build ceiling from step 1/2) rather than downgrading snakemake.

**5. Environment yaml example** (validated combination):

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - snakemake-minimal=9.19.0
  - conda>=24.7.1
  - snakemake-executor-plugin-cluster-generic
  - tensorflow=2.17.0=cuda120*
  - "numpy<2"          # TF 2.17 predates numpy 2.x support
```

**6. Verify GPU is actually visible before trusting any training run** — `tf.config.list_physical_devices('GPU')` returning non-empty is necessary but not sufficient; run a real op on it:

```bash
srun --partition=gpu --gres=gpu:1 --time=00:05:00 bash -c '
  source /work/project/<user>/miniforge3/etc/profile.d/conda.sh
  conda activate /work/project/<user>/conda/envs/<env>
  python3 -c "
import tensorflow as tf
print(tf.config.list_physical_devices())
with tf.device(\"/GPU:0\"):
    a = tf.random.normal([2000,2000]); b = tf.random.normal([2000,2000])
    print((tf.matmul(a,b)).device)
"'
```

A successful run prints the GPU model (e.g. `NVIDIA RTX A5000`) and shows the matmul's device as `/job:localhost/replica:0/task:0/device:GPU:0`.

**7. Snakemake rule template** (same GRES requirement as any other GPU rule — see Hard Constraints):

```python
rule train_model:
    input:  data = "data/training_set.npz"
    output: model = "models/model.keras"
    resources:
        partition = "gpu",
        gres      = "gpu:1",
        mem_mb    = 32000,
        time      = 1440
    threads: 4
    shell:
        "MPLBACKEND=Agg python -u train.py"   # Agg avoids a headless plt.show() hang — see below
```

**8. Headless matplotlib**: any script that calls `plt.show()` (directly or via a library) will hang indefinitely on a compute node with no display — `NSApplication run`/GTK event loops block forever waiting for a window that can never appear. Always set `MPLBACKEND=Agg` (env var, as above) or `matplotlib.use("Agg")` before the first `pyplot` import for any script run non-interactively. Under `Agg`, `plt.show()` is a no-op — replace it with `plt.savefig(...)` (or monkey-patch `plt.show` to do so) if the plots themselves are needed as output.

**Notes:**
- `${{CUDA_VISIBLE_DEVICES:-}}` — double braces escape Snakemake's string interpolation; the `:-` default handles unset variable on nodes where SLURM does not export it.
- Barcode data (adapter sequences including bc_bottom_flank) is compiled into the dorado binary, not stored as external files — cannot be inspected via `find`.

---

## Additions from the pipeline templates (SRA ChIP-seq / ATAC)

- **Wall time is not passed to sbatch by the profile above.** The `time` resource is set per rule but `cluster-generic-submit-cmd` has no `--time`. The templates in this repository add `--time={resources.time}` (minutes; accepted by the old SLURM). Without it, jobs get the partition default.
- **`sra-tools=2.11.0` environments:** drop `conda>=24.7.1` and use `samtools>=1.15,<1.17` (see Step 4).
- **Replace picard by samtools** (`fixmate -m | markdup -r`) to keep a single environment: picard + macs2 would otherwise need `--use-conda`.
- **Quote version specs in shell tests** (`'samtools>=1.15'`); an unquoted `>=` creates stray files named `=1.15` in the current directory.
- **Env creation takes ~25 minutes** for an `sra-tools` + R/Bioconductor + macs2 + multiqc environment; run it in the background.
