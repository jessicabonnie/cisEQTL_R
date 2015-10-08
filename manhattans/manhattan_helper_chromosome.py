import manhattan_helper
import sys
import os
import getopt

tablefolder = None#'/home/jkb4y/ubs/work/data/Achilleas/eQTLs_Feb2013_pcaCorrected'
exon = False
cis = True
chromosome = None
cell_type = None
try:
    opts, args = getopt.getopt(sys.argv[1:], "h",
                               ["help","tablefolder=","cis",
                                "exon","transcript","trans","chr=","cell="])
except getopt.GetoptError:
    usage()
    sys.exit(2)
for opt, arg in opts:
    if opt in ("-h","--help"):
        usage()
        sys.exit()
    elif opt in ("--tablefolder"):
        tablefolder = arg
    elif opt in ("--chr"):
        chromosome = int(arg)
    elif opt in ("--cell"):
        cell_type = arg
    elif opt in ("--cis"):
        cis = True
    elif opt in ("--trans"):
        cis = False
    elif opt in ("--exon"):
        exon = True
    elif opt in ("--transcript"):
        exon = False





table_name = manhattan_helper.determine_table_name(cell_type, exon, cis)
table_loc = os.path.join(tablefolder, table_name)

p_filter = manhattan_helper.determine_filter(cis)
manhattan_helper.multi_filter_helper(table_loc, p_filter, chromosome)
