#!/bin/bash
set -euo pipefail

# Output directory
WORKDIR=`pwd`

# Path to snakemake config file
CONFIG="${WORKDIR}/config.yaml"

# Path to snakemake cluster configuration file
CLUSTER_CONFIG="cluster.json"

# Path to write the final HTML snakemake report
REPORT="${WORKDIR}/auto-asm.html"

# snakemake cluster job submission code
CLUSTER='qsub'
CLUSTER+=' -cwd'
CLUSTER+=' -terse'
CLUSTER+=' -pe parallel {cluster.cores}'
CLUSTER+=' -N {cluster.jobname}'
CLUSTER+=' -l h_vmem={cluster.mem_mb}M'
CLUSTER+=' -l mem_free={cluster.mem_mb}M'
CLUSTER+=' -l h_rt={cluster.walltime}'
CLUSTER+=' -o '$WORKDIR'/{cluster.error}'
CLUSTER+=' -e '$WORKDIR'/{cluster.output}'
CLUSTER+=' {cluster.other_options}'

# Maximum number of concurrent snakemake jobs (local + cluster)
MAX_JOBS=12

# Number of cores on submit node, on which local rules are executed
LOCAL_CORES=4

# Resource limits for auto-asm
RESOURCES='cores=128'
RESOURCES+=' mem_mb=1000000'

# Number of seconds to wait for output files upon rule completion
LATENCY_WAIT=60

# Run snakemake
if snakemake \
  --use-conda \
  --timestamp \
  --configfile $CONFIG \
  --cluster-config $CLUSTER_CONFIG \
  --cluster "$CLUSTER" \
  --jobs $MAX_JOBS \
  --local-cores $LOCAL_CORES \
  --resources $RESOURCES \
  --directory $WORKDIR \
  --latency-wait $LATENCY_WAIT \
  $*
then
  # Generate the report if snakemake finished without error
  snakemake \
    --timestamp \
    --configfile $CONFIG \
    --directory $WORKDIR \
    --report $REPORT \
    $*
fi
© 2019 GitHub, Inc.
Terms
Privacy
Security
Status
Help
Contact GitHub
Pricing
API
Training
Blog
About
