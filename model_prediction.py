import pandas as pd
import glob


def hebing():
    cmd1 = 'cat ./wes_mut_result/snv_I.fasta ./wes_mut_result/ins_I.fasta ./wes_mut_result/del_I.fasta > ./wes_mut_result/mut.fasta'
    cmd2 = 'cat ./mut_result/mut_pro_pl.csv.fasta ./mut_result/mut_pro_nl.csv.fasta > ./wes_mut_result/noncoding.fasta'
    # 执行 cmd1
    os.system(cmd1)
    # 检查 ./mut_result/mut_pro_pl.csv.fasta 是否存在
    if os.path.exists('./mut_result/mut_pro_pl.csv.fasta'):
        os.system(cmd2)

def neo_pred_I():
    # 读取 hla_types.txt 文件
    with open('./rna_result/hla_types.txt', 'r') as f:
        hla_types = f.read().splitlines()
      # 将 HLA 分型按照指定形式进行拼接
    formatted_hla_types = ','.join(hla_types)
    cmd1 = "netMHCpan -f ./mut-result/mut.fasta -a " + formatted_hla_types + " -s -BA > ./preneo/netM.csv"
    cmd2 = "netMHCpan -f ./mut-result/noncoding.fasta -a " + formatted_hla_types + " -s -BA > ./preneo/noncoding_netM.csv"
    os.system(cmd1)
    if os.path.exists('./mut_result/noncoding.fasta'):
        os.system(cmd2)

'''
def neo_pred_II():
    # input the hla allele:
    hla_allele2 = input(
        "please input an HLA class II allele like 'DRB1_0101' or multiple alleles like 'DRB1_0102,DRB1_0301,DQB1_0101':")
    hla_allele2 = hla_allele2.replace(' ', '')
    tags = mut_types[:-1]
    ### MHC II:
    command = []
    for tag in tags:
        for i in range(15, 31):
            command.append(
                "netMHCIIpan -f outfile-wes/fasta_files/" + tag + "_MHC_2_" + str(i) + ".fasta -a " + hla_allele2 + " -BA -s -length " + str(i) + " > " + utils.temp_file + "/M2_"HC-II/" + tag + "/" + tag + "_ + str(i) + "_fasta.out")
    os.system(command)
        import pandas as pd
'''
def bindingneo():
    # Step 1: Process netM.csv and create neo1.csv
    filename = "./preneo/netM.csv"
    outFile1 = "./preneo/neo1.csv"

    with open(filename, 'r') as readFile, open(outFile1, 'w') as writeFile:
        for line in readFile:
            ls1 = line.strip().split('\00')
            ls2 = ls1[0].strip().split(' ')
            if "WB" in ls2[-1] or 'SB' in ls2[-1]:
                ls3 = '\t'.join(ls2[:])
                writeFile.write(ls3 + '\n')

    # Step 2: Process neo1.csv and create neo2.csv
    infile = "./preneo/neo1.csv"
    outfile = "./preneo/neo2.csv"

    data = pd.read_table(infile, sep="\s+", header=None)
    temp = data.iloc[:, [1, 2, -8, -6, -3, -1]]
    temp.to_csv(outfile, sep="\t", header=['HLA', 'Peptide', 'Gene', '%Rank', 'Aff(nM)', 'BindLevel'],
                index=False)

    # Step 3: Process neo2.csv and create neo.csv
    filename = "./preneo/neo2.csv"
    outname = "./preneo/neo.csv"

    with open(filename, 'r') as readFile, open(outname, 'w') as writeFile:
        for line in readFile:
            ls = line.strip().split("\t")
            ls1 = list(ls[1])
            if 'X' not in ls1:
                ls2 = ''.join(ls1)
                writeFile.write(
                    ls[0] + '\t' + ls2 + '\t' + ls[2] + '\t' + ls[3] + '\t' + ls[4] + '\t' + ls[5] + '\n')

def bindingneo_nocoding():
    # Step 1: Process netM.csv and create neo1.csv
    filename = "./preneo/noncoding_netM.csv"
    outFile1 = "./preneo/noncoding_neo1.csv"

    with open(filename, 'r') as readFile, open(outFile1, 'w') as writeFile:
        for line in readFile:
            ls1 = line.strip().split('\00')
            ls2 = ls1[0].strip().split(' ')
            if "WB" in ls2[-1] or 'SB' in ls2[-1]:
                ls3 = '\t'.join(ls2[:])
                writeFile.write(ls3 + '\n')

    # Step 2: Process neo1.csv and create neo2.csv
    infile = "./preneo/noncoding_neo1.csv"
    outfile = "./preneo/noncoding_neo2.csv"

    data = pd.read_table(infile, sep="\s+", header=None)
    temp = data.iloc[:, [1, 2, -8, -6, -3, -1]]
    temp.to_csv(outfile, sep="\t", header=['HLA', 'Peptide', 'Gene', '%Rank', 'Aff(nM)', 'BindLevel'],
                index=False)

    # Step 3: Process neo2.csv and create neo.csv
    filename = "./preneo/noncoding_neo2.csv"
    outname = "./preneo/noncoding_neo.csv"

    with open(filename, 'r') as readFile, open(outname, 'w') as writeFile:
        for line in readFile:
            ls = line.strip().split("\t")
            ls1 = list(ls[1])
            if 'X' not in ls1:
                ls2 = ''.join(ls1)
                writeFile.write(
                    ls[0] + '\t' + ls2 + '\t' + ls[2] + '\t' + ls[3] + '\t' + ls[4] + '\t' + ls[5] + '\n')


def getfasta():
    filename = "./preneo/neo.csv"
    outFile = "./preneo/neofa.fasta"
    readFile = open(filename, 'r')
    writeFile = open(outFile, 'w+')

    lines = readFile.readlines()[1:]
    for line in lines:
        ls = line.strip('\n').split('\t')
        writeFile.write(">" + ls[0] + "|" + ls[2] + '\n' + ls[1] + '\n')
    readFile.close()
    writeFile.close()


def neo_TAP():
    #先拼接成长序列接入netctlpan
    input_fasta = './preneo/neofa.fasta'
    with open(input_fasta, 'r') as file:
        content = file.read()
    squence = []
    sequences = content.split('>')[1:]
    for seq in sequences:
        lines = seq.split('\n')
        seq_name = lines[0]
        seq_content = ''.join(lines[1:])
        squence.append(seq_content)
    pre_seq = ''.join(squence)
    with open("./preneo/neo_TAP.fasta", "w") as fasta_file:
        # 写入 FASTA 标识行
        fasta_file.write(">neo_TAP\n")
        # 将整个序列写入单行
        fasta_file.write(pre_seq + "\n")
    # 读取以制表符分隔的文本文件并生成 DataFrame
    cmd1="python ./biosoft/netchop/predict.py --method netctlpan  --length 8 ./preneo/neo_TAP.fasta > ./preneo/neo_8.txt"
    cmd2 = "python ./biosoft/netchop/predict.py --method netctlpan  --length 9 ./preneo/neo_TAP.fasta > ./preneo/neo_9.txt"
    cmd3 = "python ./biosoft/netchop/predict.py --method netctlpan  --length 10 ./preneo/neo_TAP.fasta > ./preneo/neo_10.txt"
    cmd4 = "python ./biosoft/netchop/predict.py --method netctlpan  --length 11 ./preneo/neo_TAP.fasta > ./preneo/neo_11.txt"
    cmd5 = 'cat ./preneo/neo_8.txt ./preneo/neo_9.txt ./preneo/neo_10.txt ./preneo/neo_11.txt > ./preneo/neo_tap.txt'
    command = [cmd1, cmd2, cmd3, cmd4,cmd5]
    for i in command:
        os.system(i)

    data1 = "./preneo/neo_tap.txt"  # 替换为你的文件路径
    data2 = "./preneo/neo.csv"
    df_tap = pd.read_csv(data1, sep='\t')
    df_neo = pd.read_csv(data2, sep='\t')
    df3 = pd.merge(df_tap, df_neo, left_on='peptide', right_on='Peptide', how='inner')
    # 提取 'value1' 列和 'value2' 列
    df3 = df3[['Peptide', 'HLA', 'Gene', '%Rank', 'Aff(nM)', 'BindLevel', 'tap_prediction_score']]
    df3_unique = df3.drop_duplicates(subset='Peptide', keep='first')
    # 筛选列'tap_prediction_score'中大于0的行
    filtered_df = df3_unique[df3_unique['tap_prediction_score'] > 0]
    output = "./preneo/filer-tap.csv"
    filtered_df.to_csv(output, index=False)

def noncoding_getfasta():
    filename = "./preneo/noncoding_neo.csv"
    outFile = "./preneo/noncoding_neofa.fasta"
    readFile = open(filename, 'r')
    writeFile = open(outFile, 'w+')

    lines = readFile.readlines()[1:]
    for line in lines:
        ls = line.strip('\n').split('\t')
        writeFile.write(">" + ls[0] + "|" + ls[2] + '\n' + ls[1] + '\n')
    readFile.close()
    writeFile.close()

def noncoding_neo_TAP():
    #先拼接成长序列接入netctlpan
    input_fasta = './preneo/noncoding_neofa.fasta'
    with open(input_fasta, 'r') as file:
        content = file.read()
    squence = []
    sequences = content.split('>')[1:]
    for seq in sequences:
        lines = seq.split('\n')
        seq_name = lines[0]
        seq_content = ''.join(lines[1:])
        squence.append(seq_content)
    pre_seq = ''.join(squence)
    with open("./preneo/noncoding_neo_TAP.fasta", "w") as fasta_file:
        # 写入 FASTA 标识行
        fasta_file.write(">neo_TAP\n")
        # 将整个序列写入单行
        fasta_file.write(pre_seq + "\n")
    # 读取以制表符分隔的文本文件并生成 DataFrame
    cmd1="python ./biosoft/netchop/predict.py --method netctlpan  --length 8 ./preneo/noncoding_neo_TAP.fasta > ./preneo/noncoding_neo_8.txt"
    cmd2 = "python ./biosoft/netchop/predict.py --method netctlpan  --length 9 ./preneo/noncoding_neo_TAP.fasta > ./preneo/noncoding_neo_9.txt"
    cmd3 = "python ./biosoft/netchop/predict.py --method netctlpan  --length 10 ./preneo/noncoding_neo_TAP.fasta > ./preneo/noncoding_neo_10.txt"
    cmd4 = "python ./biosoft/netchop/predict.py --method netctlpan  --length 11 ./preneo/noncoding_neo_TAP.fasta > ./preneo/noncoding_neo_11.txt"
    cmd5 = 'cat ./preneo/noncoding_neo_8.txt ./preneo/noncoding_neo_9.txt ./preneo/noncoding_neo_10.txt ./preneo/noncoding_neo_11.txt > ./preneo/noncoding_neo_tap.txt'
    command = [cmd1, cmd2, cmd3, cmd4,cmd5]
    for i in command:
        os.system(i)

    data1 = "./preneo/noncoding_neo_tap.txt"  # 替换为你的文件路径
    data2 = "./preneo/noncoding_neo.csv"
    df_tap = pd.read_csv(data1, sep='\t')
    df_neo = pd.read_csv(data2, sep='\t')
    df3 = pd.merge(df_tap, df_neo, left_on='peptide', right_on='Peptide', how='inner')
    # 提取 'value1' 列和 'value2' 列
    df3 = df3[['Peptide', 'HLA', 'Gene', '%Rank', 'Aff(nM)', 'BindLevel', 'tap_prediction_score']]
    df3_unique = df3.drop_duplicates(subset='Peptide', keep='first')
    # 筛选列'tap_prediction_score'中大于0的行
    filtered_df = df3_unique[df3_unique['tap_prediction_score'] > 0]
    output = "./preneo/noncoding_filer-tap.csv"
    filtered_df.to_csv(output, index=False)

if __name__ == '__main__':
    neo_pred_I()
    #neo_pred_II()
    bindingneo()
    getfasta()
    neo_TAP()
    bindingneo_nocoding()
    noncoding_getfasta()
    noncoding_neo_TAP()
    print('Model4 Run Successful!')









