import pandas as pd
import os


def handle01_tpm_filter():
    for hla_class in range(1, 3):
        mut_type = mut_types if hla_class == 1 else mut_types[:-1]
        hla = 'I' if hla_class == 1 else 'II'
        input1 = pd.read_table("reference_files/Tran_id_gene.txt", names=['gene', 'EM'], sep=',')
        input1.dropna(axis=0, how='any', inplace=True)
        em_dic = dict(zip(input1['EM'], input1['gene']))
        input2 = pd.read_table("./rna_result/tpm/abundance.tsv", sep='\t')
        input2 = input2[input2['tpm'] > 0]  # 将abundance中tpm>0的留下
        gene_list = []  # 保留dict中存在的gene及其对应tpm
        for i in input2.index:
            target_id = input2.at[i, 'target_id']
            if target_id in em_dic:
                gene_list.append([em_dic[target_id], input2.at[i, 'tpm']])
        gene_df = pd.DataFrame(gene_list, columns=['gene', 'tpm'])
        gene_df.drop_duplicates(inplace=True)
        gene_df.to_csv('./rna_result/kallisto_gene_expression.txt', sep='\t', index=False)  # 保存输出文件1为txt
        tpm_dict = dict(gene_list)
        for mut in mut_type:
            input3 = pd.read_csv('./preneo/filer-tap.csv', sep=',', dtype=str)

            #if mut == 'fusion':  # 融合基因：若有一个tpm>0则保留
            #    input3[['gene1', 'gene2']] = input3['gene'].str.split('--', expand=True)
            #    input3['tpm1'] = input3['gene1'].apply(lambda g: tpm_dict[g] if g in tpm_dict else None)
            #    input3['tpm2'] = input3['gene2'].apply(lambda g: tpm_dict[g] if g in tpm_dict else None)
            #    input3 = input3[~input3[['tpm1', 'tpm2']].isnull().T.all()]  # 保留tpm不全为空的行
            #    input3.drop('gene', axis=1, inplace=True)

            input3['tpm'] = input3['Gene'].apply(lambda g: tpm_dict[g] if g in tpm_dict else None)
            input3.dropna(subset=['Peptide'], inplace=True)
            input3.to_csv('./neofilter/filer-tpm.csv', sep='\t', index=False)  # 保存输出文件2为txt


def handle02_MS_filter():

    input1 = pd.read_table('./neofilter/filer-tpm.csv', sep='\t')
    input2 = pd.read_table("./mass/Maxpep.txt", sep='\t')
    input3 = pd.read_csv("./preneo/noncoding_filer-tap.csv")
    pep_list = input2['Sequence'].values.tolist()
    find_result = []
    for pep in input1['Peptide']:
        find_result.append(find_pep(pep, pep_list))
    input1 = input1[find_result].copy()
    input1.drop_duplicates(inplace=True)
    input1.to_csv('./neofilter/filter-ms.txt', index=False, sep='\t')

    fasta_output = './neofilter/filter-ms.fasta'
    with open(fasta_output, 'w') as fasta_file:
        for index, row in input1.iterrows():
            gene = row['Gene'] if 'Gene' in row else 'UnknownGene'
            sequence = row['Sequence'] if 'Sequence' in row else 'UnknownSequence'
            fasta_file.write(f'>{sequence}\n')
            fasta_file.write(f'{sequence}\n')
    find_result1 = []
    for pep in input3['Peptide']:
        find_result1.append(find_pep(pep, pep_list))
    input3 = input3[find_result1].copy()
    input3.drop_duplicates(inplace=True)
    input3.to_csv('./neofilter/noncondingfilter-ms.txt', index=False, sep='\t')

    fasta_output = './neofilter/noncodingfilter-ms.fasta'
    with open(fasta_output, 'w') as fasta_file:
        for index, row in input1.iterrows():
            gene = row['Gene'] if 'Gene' in row else 'UnknownGene'
            sequence = row['Sequence'] if 'Sequence' in row else 'UnknownSequence'
            fasta_file.write(f'>{sequence}\n')
            fasta_file.write(f'{sequence}\n')
    cmd1 = 'cat ./neofilter/filter-ms.txt ./neofilter/noncodingfilter-ms.txt > ./neofilter/neo_filter-ms.txt'
    os.system(cmd1)

def handle03_filterneo():
    cmd1 = 'cat ./neofilter/filter-ms.fasta ./neofilter/noncodingfilter-ms.fasta > ./neofilter/neo_filter-ms.fasta'
    cmd2 = 'makeblastdb -in ./reference/hcneodb/dbPepNeo.fasta -dbtype prot -title hcneo_db -title "dbPepNeo" -out ./reference/dbPepNeo'
    cmd3 = 'blastp -query ./neofilter/neo_filter-ms.fasta -db ./reference/dbPepNeo -outfmt "6 qacc qseq sacc sseq evalue length pident" -evalue 100000000 -gapopen 11 -gapextend 1  > ./neofilter/filterneo.txt'
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)

if __name__ == '__main__':
    handle01_tpm_filter()
    handle02_MS_filter()
    handle03_filterneo()
    neo_filter_ms_df = pd.read_csv('./neofilter/filter-ms.txt', sep='\t')
    # 读取 BLAST 输出文件
    blast_output_df = pd.read_csv('./neofilter/neo_filter-ms.txt', sep='\t', header=None,
                                  names=['qacc', 'qseq', 'sacc', 'sseq', 'evalue', 'length', 'pident'])
    # 匹配 qacc 列与 Peptide 列，并将 evalue 列数据添加到 neo_filter-ms DataFrame 中
    neo_filter_ms_df['evalue'] = neo_filter_ms_df['Peptide'].map(blast_output_df.set_index('qacc')['evalue'])
    # 将结果保存到新的文件
    neo_filter_ms_df.to_csv('./neofilter/neo_filter-ms_with_blast.csv', sep='\t', index=False)
    # 筛选 evalue 列大于等于 60 的行
    filtered_df = neo_filter_ms_df[neo_filter_ms_df['evalue'] >= 60]
    # 按照 immune_score 列的值从大到小排序
    sorted_df = filtered_df.sort_values(by='immune_score', ascending=False)
    # 将结果保存到新的文件
    sorted_df.to_csv('./neofilter/neo_filter-ms_sorted.csv', sep='\t', index=False)

