import os
import sys
import math
import getopt


def split_file(file_loc):
    line1 = True
    base, ext = os.path.splitext(file_loc)
    if not os.path.exists(base):
        os.makedirs(base)
    out_base = os.path.join(base,'chr')
    fileHandles = dict()
    for i in range(1,24):
        fileHandles[str(i)] = open(out_base+str(i)+ext, 'w')
    with open(file_loc, mode="r") as fin:
        for line in fin:
            lsplit = line.strip().split()
            if line1:
                chr_i = lsplit.index('CHR')
                for i in range(1,24):
                    fileHandles[str(i)].write(line)
                line1 = False
            else:
                chromosome = lsplit[chr_i]
                fileHandles[chromosome].write(line)
    for i in range(1,24):
        fileHandles[str(i)].close()

def determine_table_name(cell_type, exon, cis, pca_dict=None):
    if exon:
        exonflag = 'exon'
        t_i=1
    else:
        exonflag = 'transcript'
        t_i=0
    if cis:
        cisflag = 'Cis'
    else:
        cisflag = 'Trans'
    if pca_dict is None:
        pcaflag=''
    else:
        pca = pca_dict[cell_type][t_i]
        pcaflag='_pca'+pca
           
    name = exonflag+cisflag+cell_type+pcaflag+'.txt'
    return name

def determine_table_name_old(cell_type, exon, cis):
    if exon:
        exonflag = 'exon'
        filter_tail = '_filtered_1e-04.txt'
    else:
        exonflag = 'transcript'
        filter_tail = '_filtered_1e-02.txt'
    if cis:
        cisflag = 'Cis'
    else:
        cisflag = 'Trans'
    name = exonflag+cisflag+cell_type+filter_tail
    return name

def determine_filter(cis, exon, chip):
    if chip == 'I':
        if cis:
            p_filter = 0.05/1114
        else:
            p_filter = .05/54821
    if chip == 'E':
        if transcript:
            if cis:
                p_filter= 0.05/710
            else:
                p_filter=0.05/27559
        else:
            if cis:
                p_filter=.05/668
            else:
                p_filter=0.05/24451
    return p_filter

def determine_fudge(pstar):
    #fudge = ((1/200)/(math.floor(-(math.log(pstar,10)))))
    fudge = ((1/10)/(math.floor(-(math.log(pstar,10)))))
    #print("Here's epsilon: {0}".format(fudge))
    return fudge

def filter_table(table_loc, p_filter, fileConn, title=True):
    print table_loc
    snp_dict = dict()
    counter =0
    line1 = True
    #base, ext = os.path.splitext(table_loc)
    #out_loc = base + '_manhattan'+ext
    #print out_loc
    #if title:
    #    out = open(out_loc, mode='w')
    #else:
    #    out=open(out_loc, mode='r')
    read_counter = 0
    with open(table_loc, mode="r") as fin:
        for line in fin:
            lsplit = line.strip().split()
            if line1:
                if title:
                    #out.write(line)
                    fileConn.write(line)
                snp_i = lsplit.index('SNP')
                p_i = lsplit.index('P')
                line1 = False
            else:
                if read_counter%100000 ==0:
                    print("Yes, I'm still reading around line {0}.".format(read_counter))
                snp = lsplit[snp_i]
                pval = lsplit[p_i]
                if snp in snp_dict.keys():
                    snp_dict[snp].append((pval,line))
                else:
                    snp_dict[snp]=[(pval,line)]
            read_counter = read_counter + 1
        print('finished hashing the file')
        print('I have {0} snps'.format(len(snp_dict.keys())))
        for snp in snp_dict.keys():
            snp_dict[snp].sort()
            pvals = [float(item[0]) for item in snp_dict[snp]]
            lines = [str(item[1]) for item in snp_dict[snp]]
            A = [0]
            pstar = pvals[0]
            #print(pstar)
            fudge = determine_fudge(pstar)
            for i in range(1, len(pvals)):
                if pvals[i] < p_filter:
                    A.append(i)
                    pstar = pvals[i]
                    fudge = determine_fudge(pstar)
                elif abs(pvals[i] - pstar) > fudge:
                    A.append(i)
                    pstar = pvals[i]
                    fudge = determine_fudge(pstar)
            #print("For snp {0}, here are the keys:".format(snp))
            #print(A)
            wanted_lines = [lines[index] for index in A]
            #print("I am left with this many lines:"+ str(len(wanted_lines)))
            for line in wanted_lines:
                counter = counter + 1
                #out.write(line)
                fileConn.write(line)
        #out.close()
        print('wrote {0} lines'.format(counter))

def cl_arguments(argv):
    '''Reads arguments from the command line and assigns values to globals

    Keyword arguments:
    argv -- commandline arguments (?)
    '''
    global cis, exon, tablefolder,pca_dict, chip
    tablefolder = None#'/home/jkb4y/ubs/work/data/Achilleas/eQTLs_Feb2013_pcaCorrected'
    exon = False
    cis = True
    pca_dict=None
    chip=None
    try: 
        opts, args = getopt.getopt(argv, "h",
                                   ["help","tablefolder=","cis",
                                    "exon","transcript","trans","pca.list=","chip="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("--tablefolder"):
            tablefolder = arg
        elif opt in ("--chip"):
            chip = str(arg)
        elif opt in ("--cis"):
            cis = True
        elif opt in ("--trans"):
            cis = False
        elif opt in ("--exon"):
            exon = True
        elif opt in ("--transcript"):
            exon = False
        elif opt in ("--pca.list"):
            pca_dict=dict()
            line1=True
            with open(arg,mode="r") as dictfile:
                for line in dictfile:
                    lsplit=line.strip().split()
                    if line1:
                        line1=False
                    else:
                        pca_dict[lsplit[0]]=(lsplit[1],lsplit[2])

def multi_filter_helper(table_loc, p_filter, chromosome):
    print table_loc
    foldername, ext = os.path.splitext(table_loc)
    name_base = os.path.join(foldername, 'chr')
    chr_table_loc = name_base+str(chromosome)+ext
    print chr_table_loc
    title=False
    if int(chromosome) == 1:
        title=True
    filter_table(chr_table_loc, p_filter, title)

def multi_filter(table_loc, p_filter):
    foldername, ext = os.path.splitext(table_loc)
    name_base = os.path.join(foldername, 'chr')
    out_loc = foldername + '_manhattan'+ext
    fileConn = open(out_loc, mode="w")
    for i in range(1,24):
        print("Now filtering chromosome {0} table.".format(str(i)))
        if i==1:
            title = True
        else:
            title = False
        chr_table_loc = name_base+str(i)+ext
        filter_table(chr_table_loc,p_filter, fileConn, title)
    fileConn.close()
            
def non_multifilter(table_loc, p_filter):
    base, ext = os.path.splitext(table_loc)
    out_loc = base + '_manhattan'+ext
    out = open(out_loc, mode="w")
    filter_table(table_loc, p_filter, out, title=True)
    out.close()



def main(argv):
    global cis, exon, tablefolder, pca_dict, chip
    cl_arguments(argv)
    #tablefolder = '/home/jkb4y/ubs/work/data/Achilleas/eQTLs_Feb2013_pcaCorrected'
    #tablefolder = '/home/jkb4y/ubs/work/data/Achilleas/eQTLs_Feb2013'
    #exon = True
    #cis = True
    p_filter = determine_filter(cis, exon, chip)
    for cell_type in ['B','CD4','CD8','MONO','NK']:
        print("About to filter cell_type: "+cell_type)
        table_name = determine_table_name(cell_type, exon, cis, pca_dict)
        table_loc = os.path.join(tablefolder, table_name)
        #multi_filter_helper(table_loc, p_filter, 21)
        #filter_table(table_loc, p_filter)
        split_file(table_loc)
        multi_filter(table_loc, p_filter)
        


if __name__=='__main__':
    main(sys.argv[1:])
