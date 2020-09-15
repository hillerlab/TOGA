CREATE TABLE TOGAInactMut (
	transcript varchar(50) not null,  -- unique projection ID: ${ref_transcript_ID}.${chain_ID}
    exon_num int unsigned not null,  -- exon number
    position int unsigned not null,  -- possition where mutation happened
    mut_class varchar(15) not null,  -- mutation class such as FS deletion
    mutation varchar(20) not null,  -- what exactly happened
    is_inact tinyint unsigned not null,  -- is this mutation inactivating, yes 1 or not 0
    mut_id varchar(10) not null  -- mut identifier
);
