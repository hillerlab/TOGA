#!/usr/bin/env nextflow
// Author: Bogdan Kirilenko, 2002
// Nextflow procedure to execute chain feature extraction jobs
// Joblist contains a file where each line is a separate command
// We just call these lines in parallel

// params section: basically command line arguments
params.joblist = 'NONE'  // file containing jobs

// create channel lines -> we need to execute lines in parallel
joblist = file(params.joblist)
lines = Channel.from(joblist.readLines())

process execute_jobs {
    
    // allow each process to fail 3 times
    errorStrategy 'retry'
    maxRetries 3

    input:
    val line from lines

    "${line}"  // one line represents an independent command
}
