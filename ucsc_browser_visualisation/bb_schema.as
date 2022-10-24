table togaBigBed
"TOGA predicted gene model"
(
    string chrom;       "Reference sequence chromosome or scaffold"
    uint   chromStart;  "Start position in chromosome"
    uint   chromEnd;    "End position in chromosome"
    string name;        "Name or ID of item, ideally both human readable and unique"
    uint score;         "Score (0-1000)"
    char[1] strand;     "+ or - for strand"
    uint thickStart;    "Start of where display should be thick (start codon)"
    uint thickEnd;      "End of where display should be thick (stop codon)"
    uint itemRgb;       "RGB value (use R,G,B string in input file)"
    int blockCount;     "Number of blocks"
    int[blockCount] blockSizes; "Comma separated list of block sizes"
    int[blockCount] chromStarts; "Start positions relative to chromStart"

    string ref_trans_id;  "Reference transcript ID"
    string ref_region; "Transcript region in the reference"
    string query_region; "Region in the query"
    float chain_score; "Chain orthology probability score"
    
    float chain_synteny; "Chain synteny log10 value"
    float chain_flank; "Chain flank feature"
    float chain_gl_cds_fract; "Chain global CDS fraction value"
    float chain_loc_cds_fract; "Chain local CDS fraction value"
    float chain_exon_cov; "Chain exon coverage value"
    
    float chain_intron_cov;  "Chain intron coverage value"
    string status;  "Gene loss classification"
    float perc_intact_ign_M;  "% intact ignoring missing"
    float perc_intact_int_M;  "% intact considering missing as intact"
    float intact_codon_prop;  "% intact codons"
    
    float ouf_prop;  "% out of chain"
    string mid_intact;  "Is middle 80% intact"
    string mid_pres;  "Is middle 80% fully present"
    lstring prot_alignment; "HTML-formatted protein alignment"
    lstring svg_line;  "SVG inactivating mutations visualization"
    
    lstring ref_link;  "Reference transcript link"
    lstring inact_mut_html_table;  "HTML-formatted inactivating mutations table"
    lstring exon_ali_html; "HTML-formatted exon alignment"
)
