// SLURM config file for CESAR jobs
// since CESAR have various memory requirements, this
// is just a template, TOGA will fill this itself
// depending on the CESAR job bucket
process.executor = 'slurm'
process.queue = 'batch'
process.time = '24h'  // mostly 8h is enough, just for robustness
process.memory = "${_MEMORY_}G"  // to be replaced
process.cpus = 1  // CESAR utilizes a single core only
executor.queueSize = 1000  // nextflow default is 100 - too few
