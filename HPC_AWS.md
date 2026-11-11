# Running on HPC (SLURM) and AWS

Scripts and settings below were validated locally (see "What was verified").

Two images: `wnn/` (Seurat WNN) then `trajectory/` (Monocle3).

## HPC (SLURM + Apptainer)

Apptainer runs the same images (no Docker needed). Build the image(s) as
`.sif` once on a machine with Docker, copy to the cluster, then submit:

```bash
docker build -t IMAGE .        # see run_analysis.sh for image names
apptainer build IMAGE.sif docker-daemon://IMAGE:latest
scp IMAGE.sif cluster:/path/to/repo/
# on a login node (internet): fetch data
RUNTIME=apptainer ./run_analysis.sh prefetch
mkdir -p logs && sbatch hpc/submit_slurm.sbatch
```

## AWS (single EC2 instance + S3)

One-shot instance: boots Amazon Linux 2023, installs Docker, clones this repo,
runs the analysis, syncs outputs to S3, then **terminates itself** (a
`shutdown -h +MAX_HOURS` dead-man switch also guards against hangs).

```bash
export S3_BUCKET=my-bucket          # existing bucket you can write to
export INSTANCE_PROFILE=my-ec2-role # IAM instance profile allowing s3:PutObject
./aws/launch_ec2.sh
```

Default instance `m6i.2xlarge` (60 GB gp3). Outputs land in
`s3://$S3_BUCKET/cite-seq-wnn-trajectory/<timestamp>/`, with `analysis.log`. Remember to
check cost for your instance type and region; nothing is created until you run
the script, and the instance stops billing when it terminates.

## What was verified

- `run_analysis.sh` passes `bash -n`; `DRY_RUN=1 ./run_analysis.sh` prints the container commands without running them.
