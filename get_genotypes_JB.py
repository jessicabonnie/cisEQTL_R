import sys
# This script takes three arguments
# the first one is the root of the plink file (it uses both the .map and .ped files)
# the second is a plink frequency file produced from the above
# the third one is the name of the output file
#
#
# there no error checking
plink_root = '/m/cphg-expr1/cphg-expr1/temp_data/data/BRI_eqtl_data/intermediate/plinkmap'
fileout = '/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_June2013/data/talk_genotypes.txt'
alleleout = '/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_June2013/data/MmAlleles.txt'
plink_freq = '/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_June2013/data/geno.frq'

#plink_root=sys.argv[1]
plink_map=plink_root+".map"
plink_ped=plink_root+".ped"
#fileout=sys.argv[2]
#snps_table_loc=sys.argv[3]
fout = open(fileout,'w')
index_dict={}

# The list_of_snps below holds all the snps for which we want genotypes
# Note that these are ouput in the file in the same order as in the list
list_of_snps=[]
# The following reads the map file and gets the index for each of the snps 
# The commented line list_of_snps.append creates populates the list
# of snps with all the available snps 
index=1
with open(plink_map) as fin:
    for line in fin:
        snp=line.strip().split()[1]
        index_dict.update({snp:index})
        list_of_snps.append(snp)
        index += 1
#list_of_snps = ['rs2910686']

with open(plink_freq, mode="r") as freq:
    line1=True
    allele_dict = dict()
    for line in freq:
        lsplit = line.strip().split()
        if line1:
            snp_i = lsplit.index('SNP')
            major_i = lsplit.index('A1')
            minor_i = lsplit.index('A2')
            line1=False
        else:
            allele_dict[lsplit[snp_i]] = (lsplit[major_i],lsplit[minor_i])

# major_alleles holds quadruples (quad)
# where the element in the quad are:
# 1 element snp name
# 2 element major allele genotype
# 3 position for the first plink genotype 
# 4 position for the second plink genotype


major_alleles=[]


def get_number(majorAllele,a,b):
    if a == "0" and b == "0":
        return "NA"
    return (2-((a == majorAllele) + (b == majorAllele)))

with open(plink_ped) as fin:
    first_line = fin.next()
    split_line = first_line.strip().split()

        # Note no check in the case that the first line 
        # has a missing allele needs to be changed
        # populate the major_alleles list
    for snp in list_of_snps:
        the_index=index_dict[snp]
        first = 2*(the_index-1)+6
        second = 2*(the_index-1)+6+1
        major_alleles.append((snp,allele_dict[snp][0],first,second))
        # Write out the first line of the output file
    title_line = ['FID','IID', 'SEX','STATUS']
    title_line.extend(list_of_snps)
    fout.write(" ".join([str(item) for item in title_line])+'\n')
    out = [split_line[0]]
    out.extend([split_line[1],split_line[4],split_line[5]])
    out.extend([get_number(quad[1],split_line[quad[2]],split_line[quad[3]]) for quad in major_alleles])
    fout.write(" ".join([str(item) for item in out]))
    fout.write("\n")

        #iterate over the rest of the lines of the ped file
    for line in fin:
        split_line = line.strip().split()
        out = [split_line[0]]
        out.extend([split_line[1],split_line[4],split_line[5]])
        out.extend([get_number(quad[1],split_line[quad[2]],split_line[quad[3]]) for quad in major_alleles])
        fout.write(" ".join([str(item) for item in out]))
        fout.write("\n")

fout.close()

aout = open(alleleout, mode='w')
allele_title_list = ['SNP','A1','A2']
aout.write(" ".join(allele_title_list)+'\n')
for snp in list_of_snps:
    lsplit = [snp, allele_dict[snp][0],allele_dict[snp][1]]
    aout.write(" ".join(lsplit)+'\n')
aout.close()
