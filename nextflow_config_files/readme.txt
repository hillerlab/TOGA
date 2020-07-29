Nextflow config directory must contain the following files:
1) extract_chain_features_config.nf
Parametes for extracting chain features.
Quite fast jobs requiring minimal memory.
One hour should be enough.
2) call_cesar_config_template.nf
Longer jobs (for example 24h)
Make sure you have the following line in this file:
process.memory = ${_MEMORY_}
Since CESAR part might split to different joblist depending on the memory requirements.
If you still like to assign a constant value here, like this:
process.memory = '50G'
then do not use --cesar_buckets parameter.
