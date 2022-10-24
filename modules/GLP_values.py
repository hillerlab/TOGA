"""Constants to be used in the gene loss detection modules."""

# mut classes:
MISS_EXON = "Missing exon"
DEL_EXON = "Deleted exon"
DEL_MISS = {MISS_EXON, DEL_EXON}
COMPENSATION = "COMPENSATION"
SSM = "SSM"
# (ag)acceptor-EXON-donor(gt)
SSM_D = "SSMD"  # Donor, right, GT,GC
SSM_A = "SSMA"  # Acceptor, left, AG 

START_MISSING = "START_MISSING"
ATG = "ATG"
FS_DEL = "FS_DEL"
FS_INS = "FS_INS"
BIG_DEL = "BIG_DEL"
BIG_INS = "BIG_INS"
STOP = "STOP"

# colors
BLUE = "0,0,200"
LIGHT_BLUE = "0,200,255"
LIGHT_RED = "255,50,50"
SALMON = "255,160,120"
GREY = "130,130,130"
BROWN = "159,129,112"
BLACK = "10,10,10"

# mutations-related constants
STOPS = {"TAG", "TAA", "TGA"}
D_M = {"D", "M"}
LEFT_SPLICE_CORR = ("ag",)  # acceptor
RIGHT_SPLICE_CORR = (
    "gt",
    "gc",
)  # donor
LEFT_SSID = 0
RIGHT_SSID = 1
ACCEPTOR = 0
DONOR = 1

BIG_INDEL_SIZE = 50
SAFE_EXON_DEL_SIZE = 40  # actually 39
FIRST_LAST_DEL_SIZE = 20
BIG_EXON_THR = BIG_INDEL_SIZE * 5
