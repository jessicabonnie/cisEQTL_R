import os
import sys


def find_exoncis_filter_flag(exon,cis):
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

def write_manhattan_filter_sh(script_loc, cell_type, pca, exon, cis, base_data='/home/jkb4y/ubs/work/data/Achilleas/'):
    folder = 'eQTLs_Feb2013'
    tablefolder = os.path.join(base_data,folder) 
    exonflag, cisflag, filterflag = find_exoncis_flags(exon,cis)
    if pca:
        pcaflag = 'pca_'
        tablefolder = tablefolder + 'pcaCorrected'
    else:
        pcaflag = ''
    tablebase = exonflag + cisflag.title() + cell_type + filterflag
    outlog = os.path.join(tablefolder,tablebase+pcaflag+'_manhattan')
    text('''#!/bin/sh
#PBS -l walltime=8:00:00
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

python /home/jkb4y/programs/cisEQTL_R/manhattan_helper_chromosome.py --tablefolder {1} --{2} --{3} --cell {4} --chr ${PBS_ARRAY_INDEX}
'''.format(outlog,tablefolder,exonflag,cisflag, cell_type))
    with open(script_loc, mode="w") as script:
        script.write(text)

def write_concatenate_sh(script_loc,pca,exon,trans):
    folder = 'eQTLs_Feb2013'
    tablefolder = os.path.join(base_data,folder) 
    exonflag, cisflag, filterflag = find_exoncis_flags(exon,cis)
    if pca:
        pcaflag = 'pca_'
        tablefolder = tablefolder + 'pcaCorrected'
    else:
        pcaflag = ''
    tablebase = exonflag + cisflag.title() + cell_type + filterflag
    subfolder = os.path.join(tablefolder,tablebase)
    outtable = os.path.join(tablefolder, tablebase + '_manhattan.txt')
    outlog = os.path.join(tablefolder, tablebase +'_concatenate.log')
    text('''#!/bin/sh
#PBS -l walltime=8:00:00
#PBS -o {0}
#PBS -W group_list=CPHG
#PBS -q cphg
#PBS -j oe
#PBS -m bea
#PBS -M jkb4y@virginia.edu

source .kshrc

cd $PBS_O_WORKDIR


filebase={1}/chr
filetail=_manhattan.txt
outtable={2}
head -n1 ${filebase}1${filetail} > ${outtable}

for i in ${filebase}*${filetail} do
	awk 'NR > 1' $i >> ${outtable}
done
'''.format(outlog, subfolder,outtable))
    with open(script_loc, mode="w") as script:
        script.write(text)



def main(argv):
    exon = False
    cis = True
    
    pca = False
    PBS_folder = '../PBS/{0}'
    cell_type = 'B'
    exonflag, cisflag, thing = find_exoncis_filter_flag(exon,cis)
    script_base = exonflag + cisflag.title()+cell_type+'manhattan'
    filter_script_loc = os.path.join(PBS_folder,script_name+'_filter.sh')
    concat_script_loc = os.path.join(PBS_folder,script_name+'_concatenate.sh')
    array_id = os.popen('qsub {0}'.format(filter_script_loc))
    #command_concatenation(script_loc, array_id):
    os.popen('qsub {0} -W depend=afterokarray:{1}[]'.format(concat_script_loc, array_id))


