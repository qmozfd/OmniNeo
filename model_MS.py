
from tqdm import tqdm
import os
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser


def create_db():
    cmd1='cat ./mut_result/mut_pro_pl.csv.fasta ./mut_result/mut_pro_nl.csv.fasta ./mut_result/output_del.fasta ./mut_result/output_ins.fasta ./mut_result/output_snv.fasta ./mut_result/output_fusion.fasta > ./filter-mass/mut_ref_db.fasta'
    os.system(cmd1)

def xml():
    cmd1='python ./biosoft/gen_mqpar.py ./biosoft/labelfree.xml ./ms -o ./biosoft/mqpar.xml -t 6'
    os.system(cmd1)

def maxquant():
    cmd1= "mono ./biosoft/MaxQuant/bin/MaxQuantCmd.exe ./biosoft/mqpar.xml"
    os.system(cmd1)





def handle04_maxpep():
    peptide = pd.read_table("./mass/combined/txt/peptides.txt", sep='\t', header=0)
    peptide = peptide.dropna(subset=['Proteins'])
    peptide = peptide.reset_index(drop=True)
    output_list = []
    mut_peptide = peptide[~peptide['Proteins'].str.startswith('CON__')]
    mut_peptide['Gene'] = mut_peptide['Proteins'].str.extract(r'([A-Z0-9_]+)[._]')

    output = pd.DataFrame(mut_peptide, columns=['Sequence', 'Gene'])
    output.to_csv("./mass/Maxpep.txt", index=False, sep='\t')


if __name__ == '__main__':
    create_db()
    xml()
    maxquant()
    handle04_maxpep()
    print("Model3 Run Successful!")