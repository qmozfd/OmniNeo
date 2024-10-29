
import requests
import pandas as pd
from bs4 import BeautifulSoup
import subprocess

# 抗原性的预测
def predict_antigenicity(protein_sequence, threshold=0.5, target_organism="tumour"):
    # 构建 POST 请求的数据
    data = {
        "seq": protein_sequence,
        "threshold": threshold,
        "Target": target_organism,
        # 可能还需要其他参数，具体查看源代码
    }

    # 发送 POST 请求到 VaxiJen
    url = "http://www.ddg-pharmfac.net/vaxijen/scripts/VaxiJen_scripts/VaxiJen3.pl"
    response = requests.post(url, data=data)

    lines = response.text.splitlines()

    # 提取最后几行中的值
    prediction = None
    for line in lines[::-1]:
        if 'Threshold' in line:
            threshold_value = line.split('<FONT COLOR="#0000FF">')[-1].split('</font>')[0]
            parts = line.split('Overall Prediction for the Protective Antigen =<b> ')
            prediction = parts[1].split(' </b> ( Probable ')[0]
            break  # 找到值后即可跳出循环

    return threshold_value


def AllerTOP(sequence):
    payload = {
        'sequence': sequence
    }

    url = 'http://www.ddg-pharmfac.net/AllerTOP/feedback.py'
    response = requests.post(url, data=payload)

    soup = BeautifulSoup(response.text, 'html.parser')
    div_box = soup.find('div', id='box')

    if div_box:
        h4_tags = div_box.find_all( )
        fourth_h4_text = h4_tags[1].text.strip()

        return fourth_h4_text
    else:
        return "Text not found or error occurred"


def toxinpred(input_file, output_file):
    cmd1 = f'python ./toxinpred2-main/toxinpred2-main/toxinpred2.py -i {input_file} -o {output_file}'
    os.system(cmd1)


def run_lineardesign():
    # Define the input and output files
    input_fasta = '../mrna_seq/final_sequence.fasta'  # Adjust the path relative to the cwd
    output_txt = '../mrna_seq/final_mrna.txt'
    lambda_value = 3

    # Construct the command
    cmd = f'cat {input_fasta} | ./lineardesign --lambda {lambda_value} > {output_txt}'

    # Execute the command
    result = subprocess.run(cmd, shell=True, check=True, text=True, cwd='./LinearDesign-main')

    # Check the result
    if result.returncode == 0:
        print(f"Successfully ran lineardesign and saved output to {output_txt}")
    else:
        print(f"Failed to run lineardesign. Return code: {result.returncode}")



if __name__ == '__main__':
    # 读取文件
    neo_filter_ms_sorted_df = pd.read_csv('./neofilter/neo_filter-ms_sorted.csv', sep='\t')
    # 获取 Peptide 列的前十个序列
    peptides = neo_filter_ms_sorted_df['Peptide'].head(10).tolist()
    # 使用 AAY 将前十个序列拼接成一条长序列
    long_sequence = 'AAY'.join(peptides)
    # 指定的前缀序列
    prefix_sequence = 'MIETYNQTSPRSAATGLPISMKIFMYLLTVFLITQMIGSALFAVYLHRRLDKIEDERNLHEDFVFMKTIQRCNTGERSLSLLNCEEIKSQFEGFVKDIMLNKEETKKENSFEMQKGDQNPQIAAHVISEASSKTTSVLQWAEKGYYTMSNNLVTLENGKQLTVKRQGLYYIYAQVTFCSNREASSQAPFIASLCLKSPGRFERILLRAANTHSSAKPCGQQSIHLGGVFELQPGASVFVNVTDPSQVSHGTGFTSFGLLKLGPGPG'
    # 将前缀序列加到拼接后的长序列最前面
    final_sequence = prefix_sequence + long_sequence
    predict_antigenicity(final_sequence)
    AllerTOP(predict_antigenicity)
    # 将final_sequence写入文件
    input_file_path = './mrna_seq/final_sequence.fasta'
    with open(input_file_path, 'w') as file:
        file.write(f'>{final_sequence}\n')
    # 调用toxinpred函数处理序列
    output_file_path = './toxin.csv'
    toxinpred(input_file_path, output_file_path)

    run_lineardesign()









































