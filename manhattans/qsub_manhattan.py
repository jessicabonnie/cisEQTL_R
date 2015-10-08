import os
import sys
import subprocess
import manhattan_helper
import getopt

def find_exoncis_filter_flags(exon,cis):
    if cis:
        cisflag = 'cis'
    else:
        cisflag = 'trans'
    if exon:
        exonflag = 'exon'
        filterflag = '_filtered_1e-04'
    else:
        exonflag = 'transcript'
        filterflag = '_filtered_1e-02'
    return (exonflag,cisflag, filterflag)

def write_manhattan_filter_sh(script_loc, cell_type, pca, exon, cis, base_data,folder):
    tablefolder = os.path.join(base_data,folder) 
    exonflag, cisflag, filterflag = find_exoncis_filter_flags(exon,cis)
    if pca:
        pcaflag = '_pca'
    else:
        pcaflag = ''
    tablebase = exonflag + cisflag.title() + cell_type + filterflag
    outlog = os.path.join(tablefolder,tablebase+pcaflag+'_manhattan.log')
    pbs_issue = '{PBS_ARRAY_INDEX}'
    text='''#!/bin/sh
#PBS -l walltime=2:00:00
#PBS -o {0}
#PBS -W group_list=CPHG
#PBS -q cphg
#PBS -j oe
#PBS -m bea
#PBS -M jkb4y@virginia.edu
#PBS -J 1-23


source .kshrc

cd $PBS_O_WORKDIR
module add python26/2.6
echo "MODULE ADDED"

python /home/jkb4y/programs/cisEQTL_R/manhattan_helper_chromosome.py --tablefolder {1} --{2} --{3} --cell {4} --chr ${5}
'''.format(outlog,tablefolder,exonflag,cisflag, cell_type, pbs_issue)
    with open(script_loc, mode="w") as script:
        script.write(text)

def write_concatenate_sh(script_loc,cell_type, pca,exon,cis,base_data, folder,job_id):
    #folder = 'eQTLs_Feb252013'
    tablefolder = os.path.join(base_data,folder) 
    exonflag, cisflag, filterflag = find_exoncis_filter_flags(exon,cis)
    if pca:
        pcaflag = '_pca'
    else:
        pcaflag = ''
    tablebase = exonflag + cisflag.title() + cell_type + filterflag
    subfolder = os.path.join(tablefolder,tablebase)
    outtable = os.path.join(tablefolder, tablebase + '_manhattan.txt')
    outlog = os.path.join(tablefolder, tablebase +pcaflag+'_concatenate.log')
    text='''#!/bin/sh
#PBS -l walltime=8:00:00
#PBS -o {0}

#PBS -W group_list=CPHG
#PBS -q cphg
#PBS -j oe
#PBS -m bea
#PBS -M jkb4y@virginia.edu

source .kshrc

cd $PBS_O_WORKDIR


filebase={2}/chr
filetail=_manhattan.txt
outtable={3}
head -n1 ${4}1${5} > $outtable

for i in ${4}*${5}; do
	awk 'NR > 1' $i >> $outtable;
done

rm -rf {2}
'''.format(outlog, job_id, subfolder,outtable,'{filebase}','{filetail}')
    with open(script_loc, mode="w") as script:
        script.write(text)


def cl_args(argv):
    global cis, exon, pca
    pca = False
    cis = True
    exon = False
    try:
        opts, args = getopt.getopt(argv, "h",
                                   ["help","pca","cis",
                                    "exon","transcript","trans"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("--pca"):
            pca = True
        elif opt in ("--cis"):
            cis = True
        elif opt in ("--trans"):
            cis = False
        elif opt in ("--exon"):
            exon = True
        elif opt in ("--transcript"):
            exon = False


def main(argv):
    global cis, exon, pca
    cl_args(argv)
    
    PBS_folder = '/h4/t1/users/jkb4y/programs/PBS'
    sys.path.append(PBS_folder)
    base_data='/home/jkb4y/ubs/work/data/Achilleas/'
    folder = 'eQTLs_Mar2013'
    if pca:
        folder = folder + '_pcaCorrected'
    exonflag, cisflag, filter_tag = find_exoncis_filter_flags(exon,cis)
    print PBS_folder
    #cell_type = 'B'
    for cell_type in ['B','CD4','CD8','MONO','NK']:
        table_loc = os.path.join(base_data, folder,exonflag+cisflag.title()+cell_type+filter_tag+'.txt')
        manhattan_helper.split_file(table_loc)
        script_base = exonflag + cisflag.title()+cell_type+'manhattan'
        filter_script_loc = os.path.join(PBS_folder,script_base+'_filter.sh')
        concat_script_loc = os.path.join(PBS_folder,script_base+'_concatenate.sh')
        write_manhattan_filter_sh(filter_script_loc, cell_type, pca, exon, cis, base_data, folder)
        #write_concatenate_sh(concat_script_loc,cell_type, pca,exon,cis,base_data,folder)
        #print(filter_script_loc)
        #print(concat_script_loc)
        print(os.path.exists(filter_script_loc))
        print(os.path.exists(concat_script_loc))
        #thing = os.system('qsub {0}'.format(filter_script_loc))
        print('qsub {0}'.format(filter_script_loc))
        sub = subprocess.Popen(['qsub','{0}'.format(filter_script_loc)],stdout=subprocess.PIPE)
        job_id_full = sub.communicate()[0].strip()
        
        #job_id_full = os.system('qsub {0}'.format(filter_script_loc))
        print job_id_full
        job_id_split = job_id_full.split('.')
        print job_id_split
        job_id = job_id_split[0]
        print job_id
        write_concatenate_sh(concat_script_loc,cell_type, pca,exon,cis,base_data,folder,job_id_full)
        #thing = subprocess.Popen(['qsub','{0}'.format(concat_script_loc)])
        cmd ='qsub -W depend=afterok:{0} {1}'.format(job_id_full, concat_script_loc)
        print cmd
        os.system(cmd)

if __name__=='__main__':
    main(sys.argv[1:])
