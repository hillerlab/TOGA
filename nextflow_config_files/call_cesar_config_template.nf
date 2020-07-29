process.executor = 'slurm'
process.queue = 'batch'
process.time = '24h'
process.memory = "${_MEMORY_}G"  // to be replaced
process.cpus = 1
executor.queueSize = 1000
