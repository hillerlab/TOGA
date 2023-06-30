#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Author: Bogdan Kirilenko, 2023
// Nextflow procedure to execute chain feature extraction jobs
// Joblist contains a file where each line is a separate command
// We just call these lines in parallel

// params section: basically command line arguments
params.joblist = 'NONE'  // file containing jobs

// if still default -> nothing assigned: show usage message and quit
if (params.joblist == "NONE"){
    println("Usage: nextflow execute_joblist.nf  --joblist [joblist file] -c [config file]")
    System.exit(2);
}

// create channel lines -> we need to execute lines in parallel
lines = Channel.fromPath(params.joblist).splitText()

process execute_jobs {

    // allow each process to fail 3 times
    errorStrategy 'retry'
    maxRetries 3

    input:
    val line

    // one line represents an independent command
    script:
    """
    ${line}
    """
}

workflow {
    execute_jobs(lines)
}
