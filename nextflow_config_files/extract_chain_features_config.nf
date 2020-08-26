// SLURM config for chain features extraction jobs
// relatively lightweighted jobs
process.executor = 'slurm'
process.queue = 'batch'
process.memory = '10G'
process.time = '1h'
process.cpus = 1
