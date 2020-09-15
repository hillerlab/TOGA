CREATE TABLE TOGANucl (
	transcript varchar(50) not null,  -- unique projection ID: ${ref_transcript_ID}.${chain_ID}
    exon_num int unsigned not null,  -- exon number
    exon_region varchar(100) not null,  -- region where exon was detected
    pid float not null,  -- nucleotide %id
    blosum float not null,  -- normalized blosum score
    gaps tinyint unsigned not null,  -- are there any asm gaps near? 1 - yes 0 - no
    ali_class varchar(4) not null,  -- alignemnt class: A, B, C, A+
    exp_region varchar(50) not null,  -- where exon was expected
    in_exp_region tinyint unsigned not null,  -- detected in expected region or not 1 yes 0 no
    alignment longblob not null  -- exon sequence in query
);
