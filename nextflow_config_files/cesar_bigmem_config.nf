// SLURM config for bigmem CESAR jobs
// these jobs require huge amounts of memory and
// take too long
process.executor = 'slurm'
process.queue = 'bigmem'
process.memory = '500G'
process.time = '24h'
process.cpus = 1
