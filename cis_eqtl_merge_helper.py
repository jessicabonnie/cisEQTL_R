from collections import namedtuple
import sys
sys.path.insert(0, '/home/jkb4y/h4t1/programs/plink_python/')
import pc_toolbox
import os


def read_r2(r2_loc):
    with open(r2_loc) as si:
        line1 = True
        index_dict = dict()
        storage = dict()
        r2tuple = NamedTuple("r2tuple","SNP_A,SNP_B,R2")
        for line in si:
            line_split = line.strip().split()
            if line1:
                index_dict['SNP_A'] = line_split.index('SNP_A')
                index_dict['SNP_B'] = line_split.index('SNP_B')
                index_dict['R2'] = line_split.index('R2')
                line1 = False
            else:
                SNP_A = line_split[index_dict['SNP_A']]
                SNP_B = line_split[index_dict['SNP_B']]
                R2 = line_split[index_dict['R2']]
                storage[(SNP_A,SNPB)] = r2tuple(SNP_A=SNP_A,SNP_B=SNP_B,R2=R2)
    return storage

def read_si(si_loc):
    si_list = list()
    with open(si_loc, mode="r") as si:
        for line in si:
            si_list.append(line.strip())
    print si_list
    return si_list

def read_literature(literature_loc):
    line1 = True
    index_dict = dict()
    lit_dict = dict()
    with open(literature_loc) as lit:
        for line in lit:
            line_split = line.strip().split('\t')
            if line1:
                index_dict['chr']=line_split.index('Chr_id')
                index_dict['bp']=line_split.index('Chr_pos')
                index_dict['disease'] = line_split.index('Disease/Trait')
                index_dict['link'] = line_split.index('Link')
                print index_dict
                line1 = False
            else:
                if len(line_split) < 2:
                    continue
                gwas_pos = 'chr{0}:{1}'.format(line_split[index_dict['chr']],
                                               line_split[index_dict['bp']])
                if gwas_pos in lit_dict.keys():
                    lit_dict[gwas_pos][0].add(line_split[index_dict['disease']])
                    lit_dict[gwas_pos][1].add(line_split[index_dict['link']])
                    
                else:
                    lit_dict[gwas_pos] = [set([line_split[index_dict['disease']]]),
                                          set([line_split[index_dict['link']]])]
    #print lit_dict
    return lit_dict

def add_lit(merge_loc, lit_dict, out_loc):
    out = open(out_loc, mode="w")
    line1 = True
    with open(merge_loc, mode = "r") as merge:
        for line in merge:
            lsplit = line.strip().split('\t')
            if line1:
                print lsplit
                lsplit.extend(['GWAS_DISEASES','GWAS_LINKS'])
                print lsplit
                #out.write('\t'.join(lsplit)+'\n')
                chro = lsplit.index('CHR')
                bp = lsplit.index('BP')
                line1 = False
            else:
                gwas_pos = 'chr'+lsplit[chro] + ':'+ lsplit[bp]
                try:
                    disease_list = list(lit_dict[gwas_pos][0])
                    diseases = '; '.join(disease_list)
                    link_list = list(lit_dict[gwas_pos][1])
                    links = ';'.join(link_list)
                except KeyError:
                    diseases = 'NA'
                    links = 'NA'
                lsplit.extend([diseases, links])
            out.write('\t'.join(lsplit)+'\n')
    out.close()
        

def main(argv):
    
    basefolder= '/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Oct2012'
    exon = True
    cis = True
    perm = False
    
    if not cis:
        cisfolder = 'trans'
        cis_flag = 'trans_'
        filt = '9.05e-7'
    else:
        cisfolder = 'cis'
        cis_flag = 'cis_'
        filt = '4.45e-5'
    if exon:
        transfolder = "exon"
    else:
        transfolder = "transcript"
    if perm:
        permflag = "_perm"
    else:
        permflag = ""
    outfolder = os.path.join(basefolder, cisfolder, transfolder)
    lit_loc = "/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Oct2012/data/gwascatalog_11262012.txt"
    #r2_loc = '/home/jkb4y/cphgdesk_share/Projects/IMCHIP/Intersect_SNP_list/2012Oct17/eurmeta/eurmeta_LD/all_regions_r2_0.ld'
    #si_loc = '/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Oct2012/data/si_SNP_NoMHC_20121024_Query.txt'
    lit_dict = read_literature(lit_loc)
    #r2_dict = read_r2(r2_loc)
    #si_list = read_si(si_loc)
    mamamerge_loc = os.path.join(outfolder, cis_flag + 'merge_'+ filt +'.tbl')
    thin_merge_loc = os.path.join(outfolder, cis_flag+'thin_merge_'+filt+'.tbl')
    multigene_merge_loc = os.path.join(outfolder, cis_flag+'multigene_region_merge_'+filt+'.tbl')
    for merge_loc in [mamamerge_loc,thin_merge_loc,multigene_merge_loc]:
        base, ext = os.path.splitext(merge_loc)
        out_loc = base + '_lit'+ext
        add_lit(merge_loc, lit_dict, out_loc)
    
if __name__=='__main__':
    main(sys.argv[1:])
