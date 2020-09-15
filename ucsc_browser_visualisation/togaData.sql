CREATE TABLE TOGAData (
	transcript varchar(50) not null,  -- unique projection ID: ${ref_transcript_ID}.${chain_ID}
    ref_transcript_ID varchar(50) not null,  -- transcript ID in reference
    ref_region varchar(100) not null,  -- transcript region in reference
    query_region varchar(100) not null,  -- projection region in query
    chain_score float not null,  -- chain orthology score
    chain_synteny int unsigned not null,  -- chain synteny
    chain_flank float not null,  -- flank coverage
    chain_gl_cds_fract float not null,  -- global CDS fraction
    chain_loc_cds_fract float not null,  -- local CDS fraction
    chain_exon_cov float not null,  -- local CDS coverage
    chain_intron_cov float not null,  -- local intron coverage
    status varchar(24) not null,  -- projection GLP status: loss, intact, etc
    perc_intact_ign_M float not null,  -- %intact ignoring Missing
    perc_intact_int_M float not null,  -- %intact considering missing seq intact
    intact_codon_prop float not null,  -- % of intact codons
    ouf_prop float not null,  -- out of chain proportion
    mid_intact tinyint unsigned not null,  -- middle 80% intact? 1 - True 0 - False, else - undefined
    mid_pres tinyint unsigned not null,  -- middle 80% present? 1 - True 0 - False, else - undefined
    prot_sequence longblob not null,  -- protein sequence
    svg_plot longblob not null,  -- svg string
    PRIMARY KEY(transcript)
);
