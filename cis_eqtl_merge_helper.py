from collections import namedtuple
import sys
sys.path.insert(0, '/home/jkb4y/h4t1/programs/plink_python/')
import pc_toolbox
import os


def read_literature_additions(lit_add_loc):
    with open(lit_add_loc, mode="r") as lit:
        line1 = True
        index_dict = dict()
        lit_dict = dict()
        for line in lit:
            line_split = line.strip().split("\t")
            if line1:
                index_dict['snp'] = line_split.index('SNP')
                index_dict['pubmed'] = line_split.index('PUBMEDID')
                index_dict['disease'] = line_split.index('Disease')
                index_dict['table'] = line_split.index('Table')
                line1 = False
            else:
                pmed = line_split[index_dict['pubmed']]
                snp = line_split[index_dict['snp']]
                if snp in lit_dict.keys():
                    lit_dict[snp][0].add(line_split[index_dict['disease']])
                    lit_dict[snp][1].add(line_split[index_dict['pubmed']])
                    
                else:
                    lit_dict[snp] = [set([line_split[index_dict['disease']]]),
                                     set([line_split[index_dict['pubmed']]])]
    return lit_dict

def read_r2(r2_loc):
    with open(r2_loc, mode="r") as si:
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
                index_dict['pubmed'] = line_split.index('PUBMEDID')
                print index_dict
                line1 = False
            else:
                if len(line_split) < 2:
                    continue
                gwas_pos = 'chr{0}:{1}'.format(line_split[index_dict['chr']],
                                               line_split[index_dict['bp']])
                if gwas_pos in lit_dict.keys():
                    lit_dict[gwas_pos][0].add(line_split[index_dict['disease']])
                    lit_dict[gwas_pos][1].add(line_split[index_dict['pubmed']])
                    lit_dict[gwas_pos][2].add(line_split[index_dict['link']])
                    
                else:
                    lit_dict[gwas_pos] = [set([line_split[index_dict['disease']]]),
                                          set([line_split[index_dict['pubmed']]]),
                                          set([line_split[index_dict['link']]])]
    #print lit_dict
    return lit_dict

def add_lit(merge_loc, lit_dict, add_lit_dict, out_loc):#, achilleas_key):
    out = open(out_loc, mode="w")
    line1 = True
    with open(merge_loc, mode = "r") as merge:
        for line in merge:
            lsplit = line.strip().split('\t')
            if line1:
                print lsplit
                lsplit.extend(['GWAS_DISEASES','GWAS_PUBMEDIDS','GWAS_LINKS','EXTRA_DISEASES','EXTRA_PUBMEDIDS'])
                print lsplit
                #out.write('\t'.join(lsplit)+'\n')
                chro = lsplit.index('CHR')
                bp = lsplit.index('bp')
                snp_jb_i = lsplit.index('SNP_JB')
                snp_ap_i = lsplit.index('SNP')
                snp_im_i = lsplit.index('SNP_IM')
                line1 = False
            else:
                gwas_pos = 'chr'+lsplit[chro] + ':'+ lsplit[bp]
                snp_ap = lsplit[snp_ap_i]
                snp_jb = lsplit[snp_jb_i]
                snp_im = lsplit[snp_im_i]
                try:
                    disease_list = list(lit_dict[gwas_pos][0])
                    diseases = '; '.join(disease_list)
                    pubmed_list = list(lit_dict[gwas_pos][1])
                    pubmeds = '; '.join(pubmed_list)
                    link_list = list(lit_dict[gwas_pos][2])
                    links = '; '.join(link_list)
                except KeyError:
                    diseases = 'NA'
                    pubmeds = 'NA'
                    links = 'NA'
                try:
                    add_disease_list = list(add_lit_dict[snp_ap][0])
                    add_diseases = '; '.join(add_disease_list)
                    add_pubmed_list = list(add_lit_dict[snp_ap][1])
                    add_pubmeds = '; '.join(add_pubmed_list)
                except KeyError:
                    try:
                        add_disease_list = list(add_lit_dict[snp_jb][0])
                        add_diseases = '; '.join(add_disease_list)
                        add_pubmed_list = list(add_lit_dict[snp_jb][1])
                        add_pubmeds = '; '.join(add_pubmed_list)
                    except KeyError:
                        try:
                            #im = achilleas_key[snp_ap]
                            add_disease_list = list(add_lit_dict[snp_im][0])
                            add_diseases = '; '.join(add_disease_list)
                            add_pubmed_list = list(add_lit_dict[snp_im][1])
                            add_pubmeds = '; '.join(add_pubmed_list)
                        except KeyError:
                            add_diseases = 'NA'
                            add_pubmeds = 'NA'
                lsplit.extend([diseases, pubmeds, links, add_diseases, add_pubmeds])
                
            out.write('\t'.join(lsplit)+'\n')
    out.close()

def build_achilleas_key(key_loc):
    key = dict()
    with open(key_loc, mode="r") as keyfile:
        for line in keyfile:
            lsplit = line.strip().split()
            im = lsplit[0]
            achilleas = lsplit[1]
            key[achilleas] = im
    return key

def main(argv):
    
    basefolder= '/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013_pcaCorrected'
    #basefolder= '/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Oct2012'
    #achilleas_key_loc = '/home/jkb4y/cphgdesk_share/Achilleas/cis-eQTLs_Aug2012/immIDs2currentIDs.tsv'
    #achilleas_key = build_achilleas_key(achilleas_key_loc)
    
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
    outfolder = os.path.join(basefolder, transfolder, cisfolder)
    lit_loc = "/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013_pcaCorrected/data/gwascatalog_20130111.txt"
    add_lit_loc = "/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Feb2013_pcaCorrected/data/additional_journals.txt"
    #r2_loc = '/home/jkb4y/cphgdesk_share/Projects/IMCHIP/Intersect_SNP_list/2012Oct17/eurmeta/eurmeta_LD/all_regions_r2_0.ld'
    #si_loc = '/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Oct2012/data/si_SNP_NoMHC_20121024_Query.txt'
    lit_dict = read_literature(lit_loc)
    add_lit_dict = read_literature_additions(add_lit_loc)
    #r2_dict = read_r2(r2_loc)
    #si_list = read_si(si_loc)
    mamamerge_loc = os.path.join(outfolder, cis_flag + 'merge_'+ filt +'.tbl')
    thin_merge_loc = os.path.join(outfolder, cis_flag+'thin_merge_'+filt+'.tbl')
    multigene_merge_loc = os.path.join(outfolder, cis_flag+'multigene_region_merge_'+filt+'.tbl')
    for merge_loc in [mamamerge_loc,thin_merge_loc,multigene_merge_loc]:
        base, ext = os.path.splitext(merge_loc)
        out_loc = base + '_lit'+ext
        add_lit(merge_loc, lit_dict, add_lit_dict, out_loc)#, achilleas_key)
    
if __name__=='__main__':
    main(sys.argv[1:])
