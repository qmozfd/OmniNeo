#!/usr/bin/env python
# -*- coding: utf-8 -*-i


"""
Noncoding somatic variant calling and HLA typing:
"""

import os
import sys
import pandas as pd
from pyfasta import Fasta

def wes_processing(input1, input2):
    con_case = input1.split('_')[0]
    cmd1 = 'java -jar ./biosoft/trimmomatic-0.39-1/trimmomatic.jar PE -phred33 ./wes/' + input1 + ' ./wes/' + input2 \
           + ' ./wes_result/' + con_case + '_cut_R1.fastq.gz ./wes_result/' + con_case + '_cut_unpaired_R1.fastq.gz ./wes_result/' + con_case + '_cut_R2.fastq.gz ./wes_result/' + con_case + '_cut_unpaired_R2.fastq.gz ILLUMINACLIP:biosoft/trimmomatic-0.39-1/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 HEADCROP:10 MINLEN:50'
    cmd2 = "bwa mem -t 16 -R '@RG\\tID:foo\\tSM:bar\\tLB:library1' ./reference/hg38/hg38.fa ./wes_result/" + con_case + "_cut_R1.fastq.gz ./wes_result/" + con_case + "_cut_R2.fastq.gz > ./wes_result/" + con_case + "_cut.sam"
    cmd3 = 'samtools fixmate -O bam ./wes_result/' + con_case + '_cut.sam ./wes_result/' + con_case + '_cut_fixmate.bam'

    cmd4 = 'samtools sort -O bam -o ./wes_result/' + con_case + '_cut_sorted.bam -T ./wes_result/' + con_case + '_cut_temp ./wes_result/' + con_case + '_cut_fixmate.bam'

    cmd5 = 'gatk MarkDuplicates -I ./wes_result/' + con_case + '_cut_sorted.bam -O ./wes_result/' + con_case + '_cut_marked_duplicates.bam -M ./wes_result/' + con_case + '_cut_marked_dup_metrics.txt'

    cmd6 = 'gatk AddOrReplaceReadGroups -I ./wes_result/' + con_case + '_cut_marked_duplicates.bam -O ./wes_result/' + con_case + '_cut_marked_duplicates_1.bam -ID 4 -LB lib1 -PL illumina -PU unit1 -SM 20'
    cmd7 = 'gatk BaseRecalibrator -I ./wes_result/' + con_case + '_cut_marked_duplicates_1.bam -R ./reference/hg38/hg38.fa --known-sites ./reference/dbsnp_146.hg38.vcf --known-sites ./reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -O ./wes_result/' + con_case + '_cut_recal_data.table'
    cmd8 = 'gatk ApplyBQSR -R ./reference/hg38/hg38.fa -I ./wes_result/' + con_case + '_cut_marked_duplicates_1.bam --bqsr-recal-file ./wes_result/' + con_case + '_cut_recal_data.table -O ./wes_result/' + con_case + '_cut_recal.bam'

    cmd9 = 'picard AddOrReplaceReadGroups I=./wes_result/' + con_case + '_cut_recal.bam O=./wes_result/' + con_case + '_recal.bam RGID=' + con_case + ' RGLB=library1 RGPL=illumina RGPU=unit1 SORT_ORDER=coordinate RGSM=' + con_case
    cmd10 = 'samtools index ./wes_result/' + con_case + '_recal.bam'
    cmd11 = 'rm ./wes_result/' + con_case + '_cut_unpaired_R1.fastq.gz'
    cmd12 = 'rm ./wes_result/' + con_case + '_cut_unpaired_R2.fastq.gz'
    cmd13 = 'rm ./wes_result/' + con_case + '_cut.sam'
    cmd14 = 'rm ./wes_result/' + con_case + '_cut_fixmate.bam'
    cmd15 = 'rm ./wes_result/' + con_case + '_cut_sorted.bam'
    cmd16 = 'rm ./wes_result/' + con_case + '_cut_marked_duplicates.bam'
    cmd17 = 'rm ./wes_result/' + con_case + '_cut_marked_duplicates_1.bam'
    cmd18 = 'rm ./wes_result/' + con_case + '_cut_recal.bam'
    command = [cmd1, cmd2, cmd3, cmd4, cmd5, cmd6, cmd7, cmd8, cmd9, cmd10, cmd11, cmd12, cmd13, cmd14, cmd15, cmd16,
               cmd17, cmd18]
    for i in command:
        os.system(i)


def call_mutation(input1, input2, input3, input4):
    con = input1.split('_')[0]
    case = input3.split('_')[0]
    combine = con + '_' + case
    cmd1 = 'gatk Mutect2 -R ./reference/hg38/hg38.fa -I ./wes_result/' + con + '_recal.bam -O ./wes_mut_result/normal.vcf'
    cmd2 = 'gatk Mutect2 -R ./reference/hg38/hg38.fa -I ./wes_result/' + case + '_recal.bam -tumor ' + case + ' -I ./wes_result/' + con + '_recal.bam -normal ' + con + ' -pon ./wes_mut_result/normal.vcf -O ./wes_mut_result/' + combine + '_1.vcf'
    cmd3 = 'gatk FilterMutectCalls -R ./reference/hg38/hg38.fa -V ./wes_mut_result/' + combine + '_1.vcf -O ./wes_mut_result/' + combine + '.vcf'
    cmd4 = 'vep -i ./wes_mut_result/' + combine + '.vcf --fork 4 --assembly GRCh38 --cache --cache_version 110 --dir ./software/ensembl-vep-release-110 --offline --fasta ./software/ensembl-vep-release-110/homo_sapiens/110_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --force_overwrite --canonical --symbol -o STDOUT | filter_vep --ontology --filter "Consequence is coding_sequence_variant" -o ./wes_mut_result/result.txt --force_overwrite'

    command = [cmd1, cmd2, cmd3, cmd4]
    for i in command:
        os.system(i)

def nt_reverse(seq):
    reverse_list=[]
    seq_f_list=list(seq)
    for ele in seq_f_list:
        if ele=='A':
            reverse_ele='T'
        elif ele=='T':
            reverse_ele='A'
        elif ele=='C':
            reverse_ele='G'
        else:
            reverse_ele='C'
        reverse_list.append(reverse_ele)
    seq_str=''.join(reverse_list)
    reverse_seq=seq_str[::-1]
    return reverse_seq

def nt2aa(seq):
    aa_list=[]
    seq_len=len(seq)
    aa_num=int(seq_len/3)
    for i in range(aa_num):
        nt_3=seq[i*3:i*3+3]
        nt_3_upper=nt_3.upper()
        aa_single=codon_dic[nt_3_upper]
        if aa_single!='STOP':
            aa_list.append(aa_single)
        else:
            break
    aa_seq=''.join(aa_list)
    return aa_seq



def process_del_fasta(input_file, hum_ref_file, hum_pep_file):
    # 读取VEP数据
    df_vep = pd.read_csv(input_file, delimiter='\t', comment='#', header=None)

    # 处理FASTA文件
    f_fasta = Fasta(human_reference)
    transcript_aa = {}
    with open(human_peptide, 'r') as file:
        for line in file:
            if line.startswith(">"):
                trs_name = line.strip().split(' ')[4][11:26]
                transcript_aa[trs_name] = ''
            else:
                transcript_aa[trs_name] += line.replace('\n', '')

    # 提取相关信息
    location=[]
    allele=[]
    transcript_name=[]
    consequence=[]
    protein_position=[]
    cds_position=[]
    extra=[]
    cdna_position=[]
    animo_acid_change=[]
    cdna_change=[]
    gene_symbol=[]
    strand=[]
    for line in open(input_vep_file):
        if line.startswith('#'):
            continue
        else:
            record = line.strip().split('\t')
            consequence_str = record[6].split(",")[0]
            if (consequence_str=="frameshift_variant" and record[0].split("_")[-1].split("/")[1]=="-") or (consequence_str=="inframe_deletion"):
                loc = record[1]
                alle = record[2]
                tran_n = record[4]
                cons = record[6].split(',')[0]
                pro_pos = record[9]
                cds_pos = record[8]
                ext = record[13]
                cdna_change_p = record[7]
                cdna_c = record[11]
                animo_acid_c = record[10]
                location.append(loc)
                allele.append(alle)
                transcript_name.append(tran_n)
                consequence.append(cons)
                protein_position.append(pro_pos)
                cds_position.append(cds_pos)
                extra.append(ext)
                animo_acid_change.append(animo_acid_c)
                cdna_position.append(cdna_change_p)
                cdna_change.append(cdna_c)
            else:
                continue

            for line in extra:
                element = line.split(';')
                for ele in element:
                    sub_ele = ele.split('=')
                    if sub_ele[0]=="SYMBOL":
                        g_s = sub_ele[1]
                    if sub_ele[0]=="STRAND":
                        st = sub_ele[1]
                gene_symbol.append(g_s)
                strand.append(st)

    transcript_seq=[]
    for name in transcript_name:
        if name not in transcript_aa.keys():
            seq='NULL'
        else:
            seq=transcript_aa[name]
        transcript_seq.append(seq)

    # 生成结果数据框
    df_list = []
    mut_piptide=[]
    wt_piptide=[]
    mt_header=[]
    wt_header=[]


    for i in range(len(consequence)):
        if transcript_seq[i]=="NULL":
            continue
        else:
            if consequence[i]=="inframe_deletion":
                chr_name = location[i].split(':')[0]
                del_start = location[i].split(':')[1].split('-')[0]
                del_end = location[i].split(':')[1].split('-')[-1]
                cds_start = cds_position[i].split('-')[0]
                cds_end = cds_position[i].split('-')[-1]
                protein_change_pos = protein_position[i].split('-')[0]
                strand_n = int(strand[i])
                #    protein_change_pos_end= protein_position[i].split('-')[-1]
                frame_left_num = (int(cds_start) - 1) % 3
                seq = transcript_aa[transcript_name[i]]
                wt_head = '>' + gene_symbol[i] + '_' + cdna_position[i] + '.' + cdna_change[i] + '_' + location[i] + '_' + transcript_name[i]
                mt_head = '>' + gene_symbol[i] + '_' + cdna_position[i] + '.' + cdna_change[i] + '_' + location[i] + '_' + transcript_name[i]
                window_sizes = [8, 9, 10, 11]
                result_sequences = []

                if frame_left_num == 0:
                    del_pro_num = (int(del_end) - int(del_start) + 1) / 3
                    if int(protein_change_pos) <= 10:
                        mt_pt = seq[0:(int(protein_change_pos) - 1)] + seq[int(protein_change_pos) + int(del_pro_num):int(protein_change_pos) + int(del_pro_num) + (22 - int(protein_change_pos))]
                        wt_pt = seq[0:21]
                        for window_size in window_sizes:
                            for a, b in zip(range(0, int(protein_change_pos)),range(window_size, int(protein_change_pos) + window_size)):
                                result_sequences.append(mt_pt[a:b])
                        df = pd.DataFrame({
                            'Result_Sequence': result_sequences,
                            'Additional_Column': [mt_head] * len(result_sequences)
                        })
                    elif (int(protein_change_pos) > 10 and len(seq) - int(protein_change_pos) <= 10):
                        mt_pt = seq[len(seq) - 21 - int(del_pro_num):(int(protein_change_pos) - 1)] + seq[int(protein_change_pos) + int(del_pro_num):len(seq) + 1]
                        wt_pt = seq[len(seq) - 21:len(seq) + 1]
                        protein_change_pos = len(mt_pt) - (len(seq) - int(protein_change_pos))
                        for window_size in window_sizes:
                            for a, b in zip(range(protein_change_pos - window_size, len(mt_pt) - window_size + 1),range(protein_change_pos, len(mt_pt) + 1)):
                                result_sequences.append(mt_pt[a:b])
                        df = pd.DataFrame({
                            'Result_Sequence': result_sequences,
                            'Additional_Column': [mt_head] * len(result_sequences)
                        })
                    else:
                        mt_pt = seq[int(protein_change_pos) - 11:int(protein_change_pos) - 1] + seq[int(protein_change_pos) + int(del_pro_num) - 1:int(protein_change_pos) + int(del_pro_num) + 10]
                        wt_pt = seq[int(protein_change_pos) - 11:int(protein_change_pos) + 10]
                        protein_change_pos = 11
                        for window_size in window_sizes:
                            for i in range(len(mt_pt) - window_size + 1):
                                if i <= int(protein_change_pos) - 1 < i + window_size:
                                    result_sequences.append(mt_pt[i:i + window_size])
                        df = pd.DataFrame({
                            'Result_Sequence': result_sequences,
                            'Additional_Column': [mt_head] * len(result_sequences)
                        })
                    df_list.append(df)

                elif frame_left_num == 1:
                    del_pro_num = (int(del_end) - int(del_start) + 1) / 3
                    if strand_n == 1:
                        change_nt_l = f_fasta[chr_name][int(del_start) - 2]
                        change_nt_r = f_fasta[chr_name][int(del_end):int(del_end) + 2]
                        change_nt = change_nt_l + change_nt_r
                        change_nt_upper = change_nt.upper()
                        change_aa = codon_dic[change_nt_upper]
                    else:
                        change_nt_l = f_fasta[chr_name][int(del_start) - 3:int(del_start) - 1]
                        change_nt_r = f_fasta[chr_name][int(del_end):int(del_end) + 1]
                        change_nt = change_nt_l + change_nt_r
                        change_nt_upper = change_nt.upper()
                        change_nt_upper_reverse = nt_reverse(change_nt_upper)
                        change_aa = codon_dic[change_nt_upper_reverse]
                    if change_aa != 'STOP':
                        if int(protein_change_pos) <= 10:
                            mt_pt = seq[0:int(protein_change_pos) - 1] + change_aa + seq[int(protein_change_pos) + int(del_pro_num):int(protein_change_pos) + int(del_pro_num) + 21 - int(protein_change_pos)]
                            wt_pt = seq[0:21]
                            for window_size in window_sizes:
                                for a, b in zip(range(0, int(protein_change_pos)),range(window_size, int(protein_change_pos) + window_size)):
                                    result_sequences.append(mt_pt[a:b])
                            df = pd.DataFrame({
                                'Result_Sequence': result_sequences,
                                'Additional_Column': [mt_head] * len(result_sequences)
                            })
                        elif (int(protein_change_pos) > 10 and len(seq) - int(protein_change_pos) <= 10):
                            mt_pt = seq[len(seq) - 21 - int(del_pro_num):(int(protein_change_pos) - 1)] + seq[int(protein_change_pos) + int(del_pro_num):len(seq)]
                            wt_pt = seq[len(seq) - 21:len(seq) + 1]
                            protein_change_pos = len(mt_pt) - (len(seq) - int(protein_change_pos))
                            for window_size in window_sizes:
                                for a, b in zip(range(protein_change_pos - window_size, len(mt_pt) - window_size + 1),range(protein_change_pos, len(mt_pt) + 1)):
                                    result_sequences.append(mt_pt[a:b])
                            df = pd.DataFrame({
                                'Result_Sequence': result_sequences,
                                'Additional_Column': [mt_head] * len(result_sequences)
                            })
                        else:
                            mt_pt = seq[int(protein_change_pos) - 11:int(protein_change_pos) - 1] + change_aa + seq[int(protein_change_pos) + int(del_pro_num):int(protein_change_pos) + int(del_pro_num) + 10]
                            wt_pt = seq[int(protein_change_pos) - 11:int(protein_change_pos) + 10]
                            protein_change_pos = 11
                            for window_size in window_sizes:
                                for i in range(len(mt_pt) - window_size + 1):
                                    if i <= int(protein_change_pos) - 1 < i + window_size:
                                        result_sequences.append(mt_pt[i:i + window_size])
                            df = pd.DataFrame({
                                'Result_Sequence': result_sequences,
                                'Additional_Column': [mt_head] * len(result_sequences)
                            })
                        df_list.append(df)

                    else:
                        print('1the deletion result in a STOP codon,so no peptide generate!')
                elif frame_left_num == 2:
                    del_pro_num = (int(del_end) - int(del_start) + 1) / 3
                    if strand_n == 1:
                        change_nt_l = f_fasta[chr_name][int(del_start) - 3:int(del_start) - 1]
                        change_nt_r = f_fasta[chr_name][int(del_end)]
                        change_nt = change_nt_l + change_nt_r
                        change_nt_upper = change_nt.upper()
                        change_aa = codon_dic[change_nt_upper]
                    else:
                        change_nt_l = f_fasta[chr_name][int(del_start) - 2]
                        change_nt_r = f_fasta[chr_name][int(del_end):int(del_end) + 2]
                        change_nt = change_nt_l + change_nt_r
                        change_nt_upper = change_nt.upper()
                        change_nt_upper_reverse = nt_reverse(change_nt_upper)
                        change_aa = codon_dic[change_nt_upper_reverse]
                    if change_aa != 'STOP':
                        if int(protein_change_pos) <= 10:
                            mt_pt = seq[0:int(protein_change_pos) - 1] + change_aa + seq[int(protein_change_pos) + int(del_pro_num):int(protein_change_pos) + int(del_pro_num) + 21 - int(protein_change_pos)]
                            wt_pt = seq[0:21]
                            for window_size in window_sizes:
                                for a, b in zip(range(0, int(protein_change_pos)),range(window_size, int(protein_change_pos) + window_size)):
                                    result_sequences.append(mt_pt[a:b])
                            df = pd.DataFrame({
                                'Result_Sequence': result_sequences,
                                'Additional_Column': [mt_head] * len(result_sequences)
                            })
                        elif (int(protein_change_pos) > 10 and len(seq) - int(protein_change_pos) <= 10):
                            mt_pt = seq[len(seq) - 21 - int(del_pro_num):(int(protein_change_pos) - 1)] + seq[int(protein_change_pos) + int(del_pro_num):len(seq)]
                            wt_pt = seq[len(seq) - 21:len(seq) + 1]
                            protein_change_pos = len(mt_pt) - (len(seq) - int(protein_change_pos))
                            for window_size in window_sizes:
                                for a, b in zip(range(protein_change_pos - window_size, len(mt_pt) - window_size + 1),range(protein_change_pos, len(mt_pt) + 1)):
                                    result_sequences.append(mt_pt[a:b])
                            df = pd.DataFrame({
                                'Result_Sequence': result_sequences,
                                'Additional_Column': [mt_head] * len(result_sequences)
                            })
                        else:
                            mt_pt = seq[int(protein_change_pos) - 11:int(protein_change_pos) - 1] + change_aa + seq[int(protein_change_pos) + int(del_pro_num):int(protein_change_pos) + int(del_pro_num) + 10]
                            wt_pt = seq[int(protein_change_pos) - 11:int(protein_change_pos) + 10]
                            protein_change_pos = 11
                            for window_size in window_sizes:
                                for i in range(len(mt_pt) - window_size + 1):
                                    if i <= int(protein_change_pos) - 1 < i + window_size:
                                        result_sequences.append(mt_pt[i:i + window_size])
                            df = pd.DataFrame({
                                'Result_Sequence': result_sequences,
                                'Additional_Column': [mt_head] * len(result_sequences)
                            })
                        df_list.append(df)

                    else:
                        print('2the deletion result in a STOP codon,so no peptide generate!')
                    # mut_piptide.append(mt_pt)
                    # wt_piptide.append(wt_pt)

            elif consequence[i]=="frameshift_variant":
                chr_name = location[i].split(':')[0]
                del_loc_start = int(location[i].split(':')[1].split('-')[0])
                del_loc_end = int(location[i].split(':')[1].split('-')[-1])
                # print del_loc_start
                strand_n = int(strand[i])
                cds_loc = cds_position[i].split('-')[0]
                # print protein_position[i],consequence[i]
                protein_change_pos_start = int(protein_position[i].split('-')[0])
                frame_left_num = (int(cds_loc) - 1) % 3
                seq = transcript_aa[transcript_name[i]]
                wt_head = '>' + gene_symbol[i] + '_' + cdna_position[i] + '.' + cdna_change[i] + '_' + location[i] + '_' + transcript_name[i]
                mt_head = '>' + gene_symbol[i] + '_' + cdna_position[i] + '.' + cdna_change[i] + '_' + location[i] + '_' + transcript_name[i]
                window_sizes = [8, 9, 10, 11]
                result_sequences = []

                if strand_n == 1:
                    nt_left = f_fasta[chr_name][del_loc_start - 1 - frame_left_num:del_loc_start - 1]
                    if int(protein_change_pos_start) <= 10:
                        nt_right = f_fasta[chr_name][del_loc_end:del_loc_end + 34 - frame_left_num]
                        nt_change = nt_left + nt_right
                        aa_change = nt2aa(nt_change)
                        mt_pt = seq[0:protein_change_pos_start - 1] + aa_change
                        mt_aa_len = len(mt_pt)
                        wt_pt = seq[0:mt_aa_len]

                        for window_size in window_sizes:
                            for a, b in zip(range(0, protein_change_pos_start),range(window_size, protein_change_pos_start + window_size)):
                                result_sequences.append(mt_pt[a:b])
                        df = pd.DataFrame({
                            'Result_Sequence': result_sequences,
                            'Additional_Column': [mt_head] * len(result_sequences)
                        })

                    elif (int(protein_change_pos_start) > 10 and len(seq) - int(protein_change_pos_start) <= 10):
                        right_aa_num = len(seq) - protein_change_pos_start
                        nt_right = f_fasta[chr_name][del_loc_end + 1:del_loc_end + right_aa_num * 3]
                        nt_change = nt_left + nt_right
                        aa_change = nt2aa(nt_change)
                        mt_pt = seq[int(protein_change_pos_start) - 11:protein_change_pos_start - 1] + aa_change
                        wt_pt = seq[int(protein_change_pos_start) - 11:protein_change_pos_start - 1 + len(aa_change)]
                        protein_change_pos_start = len(mt_pt) - (len(seq) - int(protein_change_pos_start))
                        for window_size in window_sizes:
                            for a, b in zip(range(protein_change_pos_start - window_size, len(mt_pt) - window_size + 1),range(protein_change_pos_start, len(mt_pt) + 1)):
                                result_sequences.append(mt_pt[a:b])
                        df = pd.DataFrame({
                            'Result_Sequence': result_sequences,
                            'Additional_Column': [mt_head] * len(result_sequences)
                        })
                    else:
                        nt_right = f_fasta[chr_name][del_loc_end:del_loc_end + 33 - frame_left_num]
                        nt_change = nt_left + nt_right
                        # print nt_change
                        # print len(nt_change)
                        aa_change = nt2aa(nt_change)
                        # print len(aa_change),aa_change
                        mt_pt = seq[int(protein_change_pos_start) - 11:protein_change_pos_start - 1] + aa_change
                        mt_aa_len = len(aa_change)
                        wt_pt = seq[int(protein_change_pos_start) - 11:int(protein_change_pos_start) - 1 + mt_aa_len]
                        protein_change_pos_start = 11
                        for window_size in window_sizes:
                            for i in range(len(mt_pt) - window_size + 1):
                                if i <= int(protein_change_pos_start) - 1 < i + window_size:
                                    result_sequences.append(mt_pt[i:i + window_size])
                        df = pd.DataFrame({
                            'Result_Sequence': result_sequences,
                            'Additional_Column': [mt_head] * len(result_sequences)
                        })
                    df_list.append(df)


                else:
                    change_nt_right_all_reverse = nt_reverse(
                        f_fasta[chr_name][del_loc_start - (3 - frame_left_num) - 31:del_loc_start + frame_left_num])
                    # print change_nt_right_all_reverse
                    if frame_left_num == 0:
                        change_nt = change_nt_right_all_reverse[1:]
                    elif frame_left_num == 1:
                        change_nt = change_nt_right_all_reverse[0] + change_nt_right_all_reverse[2:]
                    else:
                        change_nt = change_nt_right_all_reverse[0:2] + change_nt_right_all_reverse[3:]
                    # print change_nt
                    if int(protein_change_pos_start) <= 10:
                        aa_change = nt2aa(change_nt)
                        mt_pt = seq[0:protein_change_pos_start - 1] + aa_change
                        mt_aa_len = len(mt_pt)
                        wt_pt = seq[0:mt_aa_len]
                        for window_size in window_sizes:
                            for a, b in zip(range(0, protein_change_pos_start),range(window_size, protein_change_pos_start + window_size)):
                                result_sequences.append(mt_pt[a:b])
                        df = pd.DataFrame({
                            'Result_Sequence': result_sequences,
                            'Additional_Column': [mt_head] * len(result_sequences)
                        })
                    elif (int(protein_change_pos_start) > 10 and len(seq) - int(protein_change_pos_start) <= 10):
                        aa_change = nt2aa(change_nt)
                        mt_pt = seq[int(protein_change_pos_start) - 11:protein_change_pos_start] + aa_change
                        wt_pt = seq[int(protein_change_pos_start) - 11:protein_change_pos_start + len(aa_change)]
                        protein_change_pos_start = len(mt_pt) - (len(seq) - int(protein_change_pos_start))
                        for window_size in window_sizes:
                            for a, b in zip(range(protein_change_pos_start - window_size, len(mt_pt) - window_size + 1),range(protein_change_pos_start, len(mt_pt) + 1)):
                                result_sequences.append(mt_pt[a:b])
                        df = pd.DataFrame({
                            'Result_Sequence': result_sequences,
                            'Additional_Column': [mt_head] * len(result_sequences)
                        })
                    else:
                        aa_change = nt2aa(change_nt)
                        # print aa_change
                        mt_pt = seq[int(protein_change_pos_start) - 11:protein_change_pos_start - 1] + aa_change
                        mt_aa_len = len(aa_change)
                        wt_pt = seq[int(protein_change_pos_start) - 11:int(protein_change_pos_start) - 1 + mt_aa_len]
                        protein_change_pos_start = 11
                        for window_size in window_sizes:
                            for i in range(len(mt_pt) - window_size + 1):
                                if i <= int(protein_change_pos_start) - 1 < i + window_size:
                                    result_sequences.append(mt_pt[i:i + window_size])
                        df = pd.DataFrame({
                            'Result_Sequence': result_sequences,
                            'Additional_Column': [mt_head] * len(result_sequences)
                        })
                    mut_piptide.append(mt_pt)
                    wt_piptide.append(wt_pt)
                    mt_header.append(mt_head)
                    wt_header.append(wt_head)
                    df_list.append(df)

            else:
                pass



    result_df = pd.concat(df_list, ignore_index=True)
    result_df = result_df.drop_duplicates(subset=['Result_Sequence'])

    result_df.to_csv("./wes_mut_result/del_I.txt", index=False)
    file_path = "./wes_mut_result/del_I.fasta"
    # 保存成FASTA文件
    with open(file_path, 'w') as file:
        for index, row in result_df.iterrows():
            file.write(row['Additional_Column'] + '\n')
            file.write(row['Result_Sequence'] + '\n')

def process_ins_data(input_file, hum_ref_file, hum_pep_file):
    location = []
    allele = []
    transcript_name = []
    consequence = []
    protein_position = []
    cds_position = []
    ref_animo_acid = []
    alt_animo_acid = []
    extra = []
    cdna_position = []
    animo_acid_change = []
    cdna_change = []

    for line in open(input_file):
        if line.startswith('#'):
            continue
        else:
            record = line.strip().split('\t')
            consequence_str = record[6].split(",")[0]
            if (consequence_str == "frameshift_variant" and record[0].split("_")[-1].split("/")[0] == "-") or (consequence_str == "inframe_insertion"):
                loc = record[1]
                alle = record[2]
                tran_n = record[4]
                cons = record[6].split(',')[0]
                pro_pos = record[9]
                cds_pos = record[8]
                if record[10] == '*':
                    ref_aa = '*'
                    alt_aa = '*'
                else:
                    ref_aa = record[10].split('/')[0]
                    alt_aa = record[10].split('/')[-1]
                ext = record[13]
                cdna_change_p = record[7]
                cdna_c = record[11]
                animo_acid_c = record[10]
                location.append(loc)
                allele.append(alle)
                transcript_name.append(tran_n)
                consequence.append(cons)
                protein_position.append(pro_pos)
                cds_position.append(cds_pos)
                ref_animo_acid.append(ref_aa)
                alt_animo_acid.append(alt_aa)
                extra.append(ext)
                animo_acid_change.append(animo_acid_c)
                cdna_position.append(cdna_change_p)
                cdna_change.append(cdna_c)
            else:
                continue
    gene_symbol = []
    strand = []
    for line in extra:
        element = line.split(';')
        for ele in element:
            sub_ele = ele.split('=')
            if sub_ele[0] == "SYMBOL":
                g_s = sub_ele[1]
            elif sub_ele[0] == "STRAND":
                st = sub_ele[1]
        gene_symbol.append(g_s)
        strand.append(st)
    transcript_seq = []
    for name in transcript_name:
        if name not in transcript_aa.keys():
            seq = 'NULL'
        else:
            seq = transcript_aa[name]
        transcript_seq.append(seq)
    df_list = []
    mut_piptide = []
    wt_piptide = []
    mt_header = []
    wt_header = []
    window_sizes = [8, 9, 10, 11]
    result_sequences = []
    for i in range(len(location)):
        if transcript_seq[i] == "NULL":
            continue
        elif consequence[i] == 'inframe_insertion':
            chr_name = location[i].split(':')[0]
            ins_start = int(location[i].split(':')[1].split('-')[0])
            cds_start = int(cds_position[i].split('-')[0])
            protein_change_pos = int(protein_position[i].split('-')[0])
            ref_aa = ref_animo_acid[i]
            alt_aa = alt_animo_acid[i]
            trans_strand = int(strand[i])
            seq = transcript_aa[transcript_name[i]]
            wt_head = '>' + gene_symbol[i] + '_' + cdna_position[i] + '.' + cdna_change[i] + '_' + location[i] + '_' + transcript_name[i]
            mt_head = '>' + gene_symbol[i] + '_' + cdna_position[i] + '.' + cdna_change[i] + '_' + location[i] + '_' + transcript_name[i]
            alt_aa_len = len(alt_aa)
            if int(protein_change_pos) <= 10:
                wt_pt = seq[0:protein_change_pos + 10]
                mt_pt = seq[0:protein_change_pos - 1] + alt_aa + seq[protein_change_pos:protein_change_pos + 11]
                for window_size in window_sizes:
                    for a, b in zip(range(0, protein_change_pos + alt_aa_len),range(window_size, protein_change_pos + alt_aa_len + window_size)):
                        result_sequences.append(mt_pt[a:b])
                df = pd.DataFrame({
                    'Result_Sequence': result_sequences,
                    'Additional_Column': [mt_head] * len(result_sequences)
                })
            elif protein_change_pos > 10 and len(seq) - protein_change_pos <= 10:
                wt_pt = seq[protein_change_pos - 11:len(seq)]
                mt_pt = seq[protein_change_pos - 11:protein_change_pos - 1] + alt_aa + seq[len(seq) - protein_change_pos - alt_aa_len:len(seq)]
                protein_change_pos = len(mt_pt) - (len(seq) - int(protein_change_pos))
                for window_size in window_sizes:
                    for a, b in zip(range(protein_change_pos - window_size, len(mt_pt) - window_size + 1),range(protein_change_pos, len(mt_pt) + 1)):
                        result_sequences.append(mt_pt[a:b])
                df = pd.DataFrame({
                    'Result_Sequence': result_sequences,
                    'Additional_Column': [mt_head] * len(result_sequences)
                })
            else:
                wt_pt = seq[protein_change_pos - 11:protein_change_pos + 10]
                if ref_aa == '-':
                    mt_pt = seq[protein_change_pos - 11:protein_change_pos - 1] + alt_aa + seq[protein_change_pos - 1:protein_change_pos + 10]
                    protein_change_pos = 11
                    for window_size in window_sizes:
                        for a, b in zip(range(protein_change_pos - window_size, protein_change_pos + alt_aa_len),range(protein_change_pos, len(mt_pt))):
                            result_sequences.append(mt_pt[a:b])
                    df = pd.DataFrame({
                        'Result_Sequence': result_sequences,
                        'Additional_Column': [mt_head] * len(result_sequences)
                    })
                else:
                    mt_pt = seq[protein_change_pos - 11:protein_change_pos - 1] + alt_aa + seq[protein_change_pos:protein_change_pos + 11]
                    protein_change_pos = 11
                    for window_size in window_sizes:
                        for a, b in zip(range(protein_change_pos - window_size, protein_change_pos + alt_aa_len),range(protein_change_pos, len(mt_pt))):
                            result_sequences.append(mt_pt[a:b])
                    df = pd.DataFrame({
                        'Result_Sequence': result_sequences,
                        'Additional_Column': [mt_head] * len(result_sequences)
                    })
                df_list.append(df)

        elif consequence[i] == 'frameshift_variant':
            chr_name = location[i].split(':')[0]
            ins_start = int(location[i].split(':')[1].split('-')[0])
            cds_start = int(cds_position[i].split('-')[0])
            protein_change_pos = int(protein_position[i].split('-')[0])
            allele_ins = allele[i]
            alt_aa = alt_animo_acid[i]
            alt_aa_len = len(alt_aa)
            trans_strand = int(strand[i])
            frame_left_num = int(cds_start) % 3
            seq = transcript_aa[transcript_name[i]]
            if protein_change_pos < 11:
                aa_ten_before = seq[0:protein_change_pos - 1]
            else:
                aa_ten_before = seq[protein_change_pos - 11:protein_change_pos - 1]
            wt_head = '>' + gene_symbol[i] + '_' + cdna_position[i] + '.' + cdna_change[i] + '_' + location[i] + '_' + transcript_name[i]
            mt_head = '>' + gene_symbol[i] + '_' + cdna_position[i] + '.' + cdna_change[i] + '_' + location[i] + '_' + transcript_name[i]
            # print wt_head
            if trans_strand == 1:
                nt_all = f_fasta[chr_name][ins_start - frame_left_num:ins_start] + allele_ins + f_fasta[chr_name][ins_start:ins_start + 33 - frame_left_num]
                # print len(nt_all)
                # print aa_ten_before,nt2aa(nt_all)
                mt_pt = aa_ten_before + nt2aa(nt_all)
                wt_pt = seq[protein_change_pos - 11:protein_change_pos - 1 + len(nt2aa(nt_all))]
                protein_change_pos = 11
                for window_size in window_sizes:
                    for a, b in zip(range(protein_change_pos - window_size, protein_change_pos + alt_aa_len),range(protein_change_pos, len(mt_pt))):
                        result_sequences.append(mt_pt[a:b])
                df = pd.DataFrame({
                    'Result_Sequence': result_sequences,
                    'Additional_Column': [mt_head] * len(result_sequences)
                })
            else:
                reverse_ins = nt_reverse(allele_ins)
                nt_all = nt_reverse(f_fasta[chr_name][ins_start:ins_start + frame_left_num]) + reverse_ins + nt_reverse(f_fasta[chr_name][ins_start + frame_left_num - 34:ins_start - 1])
                # print len(nt_all),nt_all,nt2aa(nt_all)
                mt_ptl = aa_ten_before + nt2aa(nt_all)
                protein_change_pos = 11
                for window_size in window_sizes:
                    for a, b in zip(range(protein_change_pos - window_size, protein_change_pos + alt_aa_len),range(protein_change_pos, len(mt_pt))):
                        result_sequences.append(mt_pt[a:b])
                df = pd.DataFrame({
                    'Result_Sequence': result_sequences,
                    'Additional_Column': [mt_head] * len(result_sequences)
                })
            if mt_pt != '':
                df_list.append(df)
                mut_piptide.append(mt_pt)
                wt_piptide.append(wt_pt)
                mt_header.append(mt_head)
                wt_header.append(wt_head)
            else:
                continue
        else:
            continue

    result_df = pd.concat(df_list, ignore_index=True)
    result_df = result_df.drop_duplicates(subset=['Result_Sequence'])
    result_df.to_csv("./wes_mut_result/ins_I.txt", index=False)


    file_path = "./wes_mut_result/ins_I.fasta"
    # 保存成FASTA文件
    with open(file_path, 'w') as file:
        for index, row in result_df.iterrows():
            file.write(row['Additional_Column'] + '\n')
            file.write(row['Result_Sequence'] + '\n')


def process_snv_data(input_file, hum_ref_file, hum_pep_file):
    # 读取FASTA文件
    f_fasta = Fasta(hum_ref_file)

    # 读取人类肽序列
    transcript_aa = {}
    with open(hum_pep_file, 'r') as file:
        transcript_name = None
        for line in file:
            if line.startswith(">"):
                transcript_name = line.strip().split(' ')[4][11:26]
                transcript_aa[transcript_name] = ''
            else:
                transcript_aa[transcript_name] += line.replace('\n', '')

    # 解析输入文件
    protein_position = []
    cdna_position = []
    extra = []
    trans_name = []
    ref_animo_acid = []
    alt_animo_acid = []
    ref_nucleotide = []
    alt_nucleotide = []
    chrom_pos = []

    with open(input_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            elif line.strip().split('\t')[6] != "missense_variant":
                continue
            else:
                record = line.strip().split('\t')
                chr_p = record[1]
                tran_n = record[4]
                pro_pos = record[9].split('-')[-1]
                ext = record[13]
                alt_aa = record[10].split('/')[1]
                ref_aa = record[10].split('/')[0]
                cdna_p = record[7]
                ref_n = re.findall('[A-Z]', record[11].split('/')[0])[0]
                alt_n = re.findall('[A-Z]', record[11].split('/')[1])[0]
                alt_animo_acid.append(alt_aa)
                ref_animo_acid.append(ref_aa)
                trans_name.append(tran_n)
                protein_position.append(pro_pos)
                extra.append(ext)
                cdna_position.append(cdna_p)
                ref_nucleotide.append(ref_n)
                alt_nucleotide.append(alt_n)
                chrom_pos.append(chr_p)

    # 提取基因符号
    gene_symbol = []
    for line in extra:
        element = line.split(';')
        for ele in element:
            sub_ele = ele.split('=')
            if sub_ele[0] == "SYMBOL":
                g_s = sub_ele[1]
        gene_symbol.append(g_s)

    # 生成转录本序列
    transcript_seq = []
    for name in trans_name:
        if name not in transcript_aa.keys():
            seq = 'NULL'
        else:
            seq = transcript_aa[name]
        transcript_seq.append(seq)

    # 生成突变肽序列
    window_sizes = [8, 9, 10, 11]
    result_sequences = []
    df_list = []
    mut_peptide = []
    wt_peptide = []
    mt_header = []
    wt_header = []

    for i in range(len(trans_name)):
        if transcript_seq[i] == "NULL":
            continue
        else:
            pro_change_pos = int(protein_position[i])
            wt_head = f'>{gene_symbol[i]}_{alt_animo_acid[i]}_{ref_nucleotide[i]}{cdna_position[i]}{alt_nucleotide[i]}_{chrom_pos[i]}_{trans_name[i]}'
            mt_head = f'>{gene_symbol[i]}_{alt_animo_acid[i]}_{ref_nucleotide[i]}{cdna_position[i]}{alt_nucleotide[i]}_{chrom_pos[i]}_{trans_name[i]}'
            ref_animo_acid_seq = transcript_seq[i]

            if pro_change_pos <= 10:
                wt_pt = ref_animo_acid_seq[0:21]
                mt_pt = ref_animo_acid_seq[0:pro_change_pos - 1] + alt_animo_acid[i] + ref_animo_acid_seq[pro_change_pos:21]
                for window_size in window_sizes:
                    for a, b in zip(range(0, int(pro_change_pos)),
                                    range(window_size, int(pro_change_pos) + window_size)):
                        result_sequences.append(mt_pt[a:b])
                df = pd.DataFrame({
                    'Result_Sequence': result_sequences,
                    'Additional_Column': [mt_head] * len(result_sequences)
                })
            elif pro_change_pos > 10 and len(ref_animo_acid_seq) - pro_change_pos <= 10:
                wt_pt = ref_animo_acid_seq[len(ref_animo_acid_seq) - 21:len(ref_animo_acid_seq)]
                mt_pt = ref_animo_acid_seq[len(ref_animo_acid_seq) - 21:pro_change_pos - 1] + alt_animo_acid[
                    i] + ref_animo_acid_seq[pro_change_pos:len(ref_animo_acid_seq)]
                protein_change_pos = len(mt_pt) - (len(seq) - int(protein_change_pos))
                for window_size in window_sizes:
                    for a, b in zip(range(protein_change_pos - window_size, len(mt_pt) - window_size + 1),
                                    range(protein_change_pos, len(mt_pt) + 1)):
                        result_sequences.append(mt_pt[a:b])
                df = pd.DataFrame({
                    'Result_Sequence': result_sequences,
                    'Additional_Column': [mt_head] * len(result_sequences)
                })
            else:
                wt_pt = ref_animo_acid_seq[pro_change_pos - 11:pro_change_pos + 10]
                mt_pt = ref_animo_acid_seq[pro_change_pos - 11:pro_change_pos - 1] + alt_animo_acid[
                    i] + ref_animo_acid_seq[pro_change_pos:pro_change_pos + 10]
                protein_change_pos = 11
                for window_size in window_sizes:
                    for j in range(len(mt_pt) - window_size + 1):
                        if j <= int(protein_change_pos) - 1 < j + window_size:
                            result_sequences.append(mt_pt[j:j + window_size])
                df = pd.DataFrame({
                    'Result_Sequence': result_sequences,
                    'Additional_Column': [mt_head] * len(result_sequences)
                })
            df_list.append(df)
            mt_header.append(mt_head)
            wt_header.append(wt_head)
            mut_peptide.append(mt_pt)
            wt_peptide.append(wt_pt)

    # 保存结果
    result_df = pd.concat(df_list, ignore_index=True)
    result_df = result_df.drop_duplicates(subset=['Result_Sequence'])
    result_df.to_csv("./wes_mut_result/snv_I.txt", index=False)

    with open("./wes_mut_result/snv_I.fasta", 'w') as file:
        for index, row in result_df.iterrows():
            file.write(row['Additional_Column'] + '\n')
            file.write(row['Result_Sequence'] + '\n')

def process_delII_fasta(input_file, hum_ref_file, hum_pep_file):
    # 读取VEP数据
    df_vep = pd.read_csv(input_file, delimiter='\t', comment='#', header=None)

    # 处理FASTA文件
    f_fasta = Fasta(human_reference)
    transcript_aa = {}
    with open(human_peptide, 'r') as file:
        for line in file:
            if line.startswith(">"):
                trs_name = line.strip().split(' ')[4][11:26]
                transcript_aa[trs_name] = ''
            else:
                transcript_aa[trs_name] += line.replace('\n', '')

    # 提取相关信息
    location=[]
    allele=[]
    transcript_name=[]
    consequence=[]
    protein_position=[]
    cds_position=[]
    extra=[]
    cdna_position=[]
    animo_acid_change=[]
    cdna_change=[]
    gene_symbol=[]
    strand=[]
    for line in open(input_vep_file):
        if line.startswith('#'):
            continue
        else:
            record = line.strip().split('\t')
            consequence_str = record[6].split(",")[0]
            if (consequence_str=="frameshift_variant" and record[0].split("_")[-1].split("/")[1]=="-") or (consequence_str=="inframe_deletion"):
                loc = record[1]
                alle = record[2]
                tran_n = record[4]
                cons = record[6].split(',')[0]
                pro_pos = record[9]
                cds_pos = record[8]
                ext = record[13]
                cdna_change_p = record[7]
                cdna_c = record[11]
                animo_acid_c = record[10]
                location.append(loc)
                allele.append(alle)
                transcript_name.append(tran_n)
                consequence.append(cons)
                protein_position.append(pro_pos)
                cds_position.append(cds_pos)
                extra.append(ext)
                animo_acid_change.append(animo_acid_c)
                cdna_position.append(cdna_change_p)
                cdna_change.append(cdna_c)
            else:
                continue

            for line in extra:
                element = line.split(';')
                for ele in element:
                    sub_ele = ele.split('=')
                    if sub_ele[0]=="SYMBOL":
                        g_s = sub_ele[1]
                    if sub_ele[0]=="STRAND":
                        st = sub_ele[1]
                gene_symbol.append(g_s)
                strand.append(st)

    transcript_seq=[]
    for name in transcript_name:
        if name not in transcript_aa.keys():
            seq='NULL'
        else:
            seq=transcript_aa[name]
        transcript_seq.append(seq)

    # 生成结果数据框
    df_list = []
    mut_piptide=[]
    wt_piptide=[]
    mt_header=[]
    wt_header=[]


    for i in range(len(consequence)):
        if transcript_seq[i]=="NULL":
            continue
        else:
            if consequence[i]=="inframe_deletion":
                chr_name = location[i].split(':')[0]
                del_start = location[i].split(':')[1].split('-')[0]
                del_end = location[i].split(':')[1].split('-')[-1]
                cds_start = cds_position[i].split('-')[0]
                cds_end = cds_position[i].split('-')[-1]
                protein_change_pos = protein_position[i].split('-')[0]
                strand_n = int(strand[i])
                #    protein_change_pos_end= protein_position[i].split('-')[-1]
                frame_left_num = (int(cds_start) - 1) % 3
                seq = transcript_aa[transcript_name[i]]
                wt_head = '>' + gene_symbol[i] + '_' + cdna_position[i] + '.' + cdna_change[i] + '_' + location[i] + '_' + transcript_name[i]
                mt_head = '>' + gene_symbol[i] + '_' + cdna_position[i] + '.' + cdna_change[i] + '_' + location[i] + '_' + transcript_name[i]
                window_sizes = [15, 16, 17, 18]
                result_sequences = []

                if frame_left_num == 0:
                    del_pro_num = (int(del_end) - int(del_start) + 1) / 3
                    if int(protein_change_pos) <= 10:
                        mt_pt = seq[0:(int(protein_change_pos) - 1)] + seq[int(protein_change_pos) + int(del_pro_num):int(protein_change_pos) + int(del_pro_num) + (36 - int(protein_change_pos))]
                        wt_pt = seq[0:56]
                        for window_size in window_sizes:
                            for a, b in zip(range(0, int(protein_change_pos)),range(window_size, int(protein_change_pos) + window_size)):
                                result_sequences.append(mt_pt[a:b])
                        df = pd.DataFrame({
                            'Result_Sequence': result_sequences,
                            'Additional_Column': [mt_head] * len(result_sequences)
                        })
                    elif (int(protein_change_pos) > 10 and len(seq) - int(protein_change_pos) <= 10):
                        mt_pt = seq[len(seq) - 56 - int(del_pro_num):(int(protein_change_pos) - 1)] + seq[int(protein_change_pos) + int(del_pro_num):len(seq) + 1]
                        wt_pt = seq[len(seq) - 56:len(seq) + 1]
                        protein_change_pos = len(mt_pt) - (len(seq) - int(protein_change_pos))
                        for window_size in window_sizes:
                            for a, b in zip(range(protein_change_pos - window_size, len(mt_pt) - window_size + 1),range(protein_change_pos, len(mt_pt) + 1)):
                                result_sequences.append(mt_pt[a:b])
                        df = pd.DataFrame({
                            'Result_Sequence': result_sequences,
                            'Additional_Column': [mt_head] * len(result_sequences)
                        })
                    else:
                        mt_pt = seq[int(protein_change_pos) - 11:int(protein_change_pos) - 1] + seq[int(protein_change_pos) + int(del_pro_num) - 1:int(protein_change_pos) + int(del_pro_num) + 10]
                        wt_pt = seq[int(protein_change_pos) - 11:int(protein_change_pos) + 10]
                        protein_change_pos = 11
                        for window_size in window_sizes:
                            for i in range(len(mt_pt) - window_size + 1):
                                if i <= int(protein_change_pos) - 1 < i + window_size:
                                    result_sequences.append(mt_pt[i:i + window_size])
                        df = pd.DataFrame({
                            'Result_Sequence': result_sequences,
                            'Additional_Column': [mt_head] * len(result_sequences)
                        })
                    df_list.append(df)

                elif frame_left_num == 1:
                    del_pro_num = (int(del_end) - int(del_start) + 1) / 3
                    if strand_n == 1:
                        change_nt_l = f_fasta[chr_name][int(del_start) - 2]
                        change_nt_r = f_fasta[chr_name][int(del_end):int(del_end) + 2]
                        change_nt = change_nt_l + change_nt_r
                        change_nt_upper = change_nt.upper()
                        change_aa = codon_dic[change_nt_upper]
                    else:
                        change_nt_l = f_fasta[chr_name][int(del_start) - 3:int(del_start) - 1]
                        change_nt_r = f_fasta[chr_name][int(del_end):int(del_end) + 1]
                        change_nt = change_nt_l + change_nt_r
                        change_nt_upper = change_nt.upper()
                        change_nt_upper_reverse = nt_reverse(change_nt_upper)
                        change_aa = codon_dic[change_nt_upper_reverse]
                    if change_aa != 'STOP':
                        if int(protein_change_pos) <= 10:
                            mt_pt = seq[0:int(protein_change_pos) - 1] + change_aa + seq[int(protein_change_pos) + int(del_pro_num):int(protein_change_pos) + int(del_pro_num) + 56 - int(protein_change_pos)]
                            wt_pt = seq[0:56]
                            for window_size in window_sizes:
                                for a, b in zip(range(0, int(protein_change_pos)),range(window_size, int(protein_change_pos) + window_size)):
                                    result_sequences.append(mt_pt[a:b])
                            df = pd.DataFrame({
                                'Result_Sequence': result_sequences,
                                'Additional_Column': [mt_head] * len(result_sequences)
                            })
                        elif (int(protein_change_pos) > 10 and len(seq) - int(protein_change_pos) <= 10):
                            mt_pt = seq[len(seq) - 56 - int(del_pro_num):(int(protein_change_pos) - 1)] + seq[int(protein_change_pos) + int(del_pro_num):len(seq)]
                            wt_pt = seq[len(seq) - 56:len(seq) + 1]
                            protein_change_pos = len(mt_pt) - (len(seq) - int(protein_change_pos))
                            for window_size in window_sizes:
                                for a, b in zip(range(protein_change_pos - window_size, len(mt_pt) - window_size + 1),range(protein_change_pos, len(mt_pt) + 1)):
                                    result_sequences.append(mt_pt[a:b])
                            df = pd.DataFrame({
                                'Result_Sequence': result_sequences,
                                'Additional_Column': [mt_head] * len(result_sequences)
                            })
                        else:
                            mt_pt = seq[int(protein_change_pos) - 11:int(protein_change_pos) - 1] + change_aa + seq[int(protein_change_pos) + int(del_pro_num):int(protein_change_pos) + int(del_pro_num) + 10]
                            wt_pt = seq[int(protein_change_pos) - 11:int(protein_change_pos) + 10]
                            protein_change_pos = 11
                            for window_size in window_sizes:
                                for i in range(len(mt_pt) - window_size + 1):
                                    if i <= int(protein_change_pos) - 1 < i + window_size:
                                        result_sequences.append(mt_pt[i:i + window_size])
                            df = pd.DataFrame({
                                'Result_Sequence': result_sequences,
                                'Additional_Column': [mt_head] * len(result_sequences)
                            })
                        df_list.append(df)

                    else:
                        print('1the deletion result in a STOP codon,so no peptide generate!')
                elif frame_left_num == 2:
                    del_pro_num = (int(del_end) - int(del_start) + 1) / 3
                    if strand_n == 1:
                        change_nt_l = f_fasta[chr_name][int(del_start) - 3:int(del_start) - 1]
                        change_nt_r = f_fasta[chr_name][int(del_end)]
                        change_nt = change_nt_l + change_nt_r
                        change_nt_upper = change_nt.upper()
                        change_aa = codon_dic[change_nt_upper]
                    else:
                        change_nt_l = f_fasta[chr_name][int(del_start) - 2]
                        change_nt_r = f_fasta[chr_name][int(del_end):int(del_end) + 2]
                        change_nt = change_nt_l + change_nt_r
                        change_nt_upper = change_nt.upper()
                        change_nt_upper_reverse = nt_reverse(change_nt_upper)
                        change_aa = codon_dic[change_nt_upper_reverse]
                    if change_aa != 'STOP':
                        if int(protein_change_pos) <= 10:
                            mt_pt = seq[0:int(protein_change_pos) - 1] + change_aa + seq[int(protein_change_pos) + int(del_pro_num):int(protein_change_pos) + int(del_pro_num) + 56 - int(protein_change_pos)]
                            wt_pt = seq[0:56]
                            for window_size in window_sizes:
                                for a, b in zip(range(0, int(protein_change_pos)),range(window_size, int(protein_change_pos) + window_size)):
                                    result_sequences.append(mt_pt[a:b])
                            df = pd.DataFrame({
                                'Result_Sequence': result_sequences,
                                'Additional_Column': [mt_head] * len(result_sequences)
                            })
                        elif (int(protein_change_pos) > 10 and len(seq) - int(protein_change_pos) <= 10):
                            mt_pt = seq[len(seq) - 56 - int(del_pro_num):(int(protein_change_pos) - 1)] + seq[int(protein_change_pos) + int(del_pro_num):len(seq)]
                            wt_pt = seq[len(seq) - 56:len(seq) + 1]
                            protein_change_pos = len(mt_pt) - (len(seq) - int(protein_change_pos))
                            for window_size in window_sizes:
                                for a, b in zip(range(protein_change_pos - window_size, len(mt_pt) - window_size + 1),range(protein_change_pos, len(mt_pt) + 1)):
                                    result_sequences.append(mt_pt[a:b])
                            df = pd.DataFrame({
                                'Result_Sequence': result_sequences,
                                'Additional_Column': [mt_head] * len(result_sequences)
                            })
                        else:
                            mt_pt = seq[int(protein_change_pos) - 11:int(protein_change_pos) - 1] + change_aa + seq[int(protein_change_pos) + int(del_pro_num):int(protein_change_pos) + int(del_pro_num) + 10]
                            wt_pt = seq[int(protein_change_pos) - 11:int(protein_change_pos) + 10]
                            protein_change_pos = 11
                            for window_size in window_sizes:
                                for i in range(len(mt_pt) - window_size + 1):
                                    if i <= int(protein_change_pos) - 1 < i + window_size:
                                        result_sequences.append(mt_pt[i:i + window_size])
                            df = pd.DataFrame({
                                'Result_Sequence': result_sequences,
                                'Additional_Column': [mt_head] * len(result_sequences)
                            })
                        df_list.append(df)

                    else:
                        print('2the deletion result in a STOP codon,so no peptide generate!')
                    # mut_piptide.append(mt_pt)
                    # wt_piptide.append(wt_pt)

            elif consequence[i]=="frameshift_variant":
                chr_name = location[i].split(':')[0]
                del_loc_start = int(location[i].split(':')[1].split('-')[0])
                del_loc_end = int(location[i].split(':')[1].split('-')[-1])
                # print del_loc_start
                strand_n = int(strand[i])
                cds_loc = cds_position[i].split('-')[0]
                # print protein_position[i],consequence[i]
                protein_change_pos_start = int(protein_position[i].split('-')[0])
                frame_left_num = (int(cds_loc) - 1) % 3
                seq = transcript_aa[transcript_name[i]]
                wt_head = '>' + gene_symbol[i] + '_' + cdna_position[i] + '.' + cdna_change[i] + '_' + location[i] + '_' + transcript_name[i]
                mt_head = '>' + gene_symbol[i] + '_' + cdna_position[i] + '.' + cdna_change[i] + '_' + location[i] + '_' + transcript_name[i]
                window_sizes = [15, 16, 17, 18]
                result_sequences = []

                if strand_n == 1:
                    nt_left = f_fasta[chr_name][del_loc_start - 1 - frame_left_num:del_loc_start - 1]
                    if int(protein_change_pos_start) <= 10:
                        nt_right = f_fasta[chr_name][del_loc_end:del_loc_end + 55 - frame_left_num]
                        nt_change = nt_left + nt_right
                        aa_change = nt2aa(nt_change)
                        mt_pt = seq[0:protein_change_pos_start - 1] + aa_change
                        mt_aa_len = len(mt_pt)
                        wt_pt = seq[0:mt_aa_len]

                        for window_size in window_sizes:
                            for a, b in zip(range(0, protein_change_pos_start),range(window_size, protein_change_pos_start + window_size)):
                                result_sequences.append(mt_pt[a:b])
                        df = pd.DataFrame({
                            'Result_Sequence': result_sequences,
                            'Additional_Column': [mt_head] * len(result_sequences)
                        })

                    elif (int(protein_change_pos_start) > 10 and len(seq) - int(protein_change_pos_start) <= 10):
                        right_aa_num = len(seq) - protein_change_pos_start
                        nt_right = f_fasta[chr_name][del_loc_end + 1:del_loc_end + right_aa_num * 3]
                        nt_change = nt_left + nt_right
                        aa_change = nt2aa(nt_change)
                        mt_pt = seq[int(protein_change_pos_start) - 11:protein_change_pos_start - 1] + aa_change
                        wt_pt = seq[int(protein_change_pos_start) - 11:protein_change_pos_start - 1 + len(aa_change)]
                        protein_change_pos_start = len(mt_pt) - (len(seq) - int(protein_change_pos_start))
                        for window_size in window_sizes:
                            for a, b in zip(range(protein_change_pos_start - window_size, len(mt_pt) - window_size + 1),range(protein_change_pos_start, len(mt_pt) + 1)):
                                result_sequences.append(mt_pt[a:b])
                        df = pd.DataFrame({
                            'Result_Sequence': result_sequences,
                            'Additional_Column': [mt_head] * len(result_sequences)
                        })
                    else:
                        nt_right = f_fasta[chr_name][del_loc_end:del_loc_end + 54 - frame_left_num]
                        nt_change = nt_left + nt_right
                        # print nt_change
                        # print len(nt_change)
                        aa_change = nt2aa(nt_change)
                        # print len(aa_change),aa_change
                        mt_pt = seq[int(protein_change_pos_start) - 11:protein_change_pos_start - 1] + aa_change
                        mt_aa_len = len(aa_change)
                        wt_pt = seq[int(protein_change_pos_start) - 11:int(protein_change_pos_start) - 1 + mt_aa_len]
                        protein_change_pos_start = 11
                        for window_size in window_sizes:
                            for i in range(len(mt_pt) - window_size + 1):
                                if i <= int(protein_change_pos_start) - 1 < i + window_size:
                                    result_sequences.append(mt_pt[i:i + window_size])
                        df = pd.DataFrame({
                            'Result_Sequence': result_sequences,
                            'Additional_Column': [mt_head] * len(result_sequences)
                        })
                    df_list.append(df)


                else:
                    change_nt_right_all_reverse = nt_reverse(
                        f_fasta[chr_name][del_loc_start - (3 - frame_left_num) - 52:del_loc_start + frame_left_num])
                    # print change_nt_right_all_reverse
                    if frame_left_num == 0:
                        change_nt = change_nt_right_all_reverse[1:]
                    elif frame_left_num == 1:
                        change_nt = change_nt_right_all_reverse[0] + change_nt_right_all_reverse[2:]
                    else:
                        change_nt = change_nt_right_all_reverse[0:2] + change_nt_right_all_reverse[3:]
                    # print change_nt
                    if int(protein_change_pos_start) <= 10:
                        aa_change = nt2aa(change_nt)
                        mt_pt = seq[0:protein_change_pos_start - 1] + aa_change
                        mt_aa_len = len(mt_pt)
                        wt_pt = seq[0:mt_aa_len]
                        for window_size in window_sizes:
                            for a, b in zip(range(0, protein_change_pos_start),range(window_size, protein_change_pos_start + window_size)):
                                result_sequences.append(mt_pt[a:b])
                        df = pd.DataFrame({
                            'Result_Sequence': result_sequences,
                            'Additional_Column': [mt_head] * len(result_sequences)
                        })
                    elif (int(protein_change_pos_start) > 10 and len(seq) - int(protein_change_pos_start) <= 10):
                        aa_change = nt2aa(change_nt)
                        mt_pt = seq[int(protein_change_pos_start) - 11:protein_change_pos_start] + aa_change
                        wt_pt = seq[int(protein_change_pos_start) - 11:protein_change_pos_start + len(aa_change)]
                        protein_change_pos_start = len(mt_pt) - (len(seq) - int(protein_change_pos_start))
                        for window_size in window_sizes:
                            for a, b in zip(range(protein_change_pos_start - window_size, len(mt_pt) - window_size + 1),range(protein_change_pos_start, len(mt_pt) + 1)):
                                result_sequences.append(mt_pt[a:b])
                        df = pd.DataFrame({
                            'Result_Sequence': result_sequences,
                            'Additional_Column': [mt_head] * len(result_sequences)
                        })
                    else:
                        aa_change = nt2aa(change_nt)
                        # print aa_change
                        mt_pt = seq[int(protein_change_pos_start) - 11:protein_change_pos_start - 1] + aa_change
                        mt_aa_len = len(aa_change)
                        wt_pt = seq[int(protein_change_pos_start) - 11:int(protein_change_pos_start) - 1 + mt_aa_len]
                        protein_change_pos_start = 11
                        for window_size in window_sizes:
                            for i in range(len(mt_pt) - window_size + 1):
                                if i <= int(protein_change_pos_start) - 1 < i + window_size:
                                    result_sequences.append(mt_pt[i:i + window_size])
                        df = pd.DataFrame({
                            'Result_Sequence': result_sequences,
                            'Additional_Column': [mt_head] * len(result_sequences)
                        })
                    mut_piptide.append(mt_pt)
                    wt_piptide.append(wt_pt)
                    mt_header.append(mt_head)
                    wt_header.append(wt_head)
                    df_list.append(df)

            else:
                pass



    result_df = pd.concat(df_list, ignore_index=True)
    result_df = result_df.drop_duplicates(subset=['Result_Sequence'])

    result_df.to_csv("./wes_mut_result/del_I.txt", index=False)
    file_path = "./wes_mut_result/del_I.fasta"
    # 保存成FASTA文件
    with open(file_path, 'w') as file:
        for index, row in result_df.iterrows():
            file.write(row['Additional_Column'] + '\n')
            file.write(row['Result_Sequence'] + '\n')

def process_insII_data(input_file, hum_ref_file, hum_pep_file):
    location = []
    allele = []
    transcript_name = []
    consequence = []
    protein_position = []
    cds_position = []
    ref_animo_acid = []
    alt_animo_acid = []
    extra = []
    cdna_position = []
    animo_acid_change = []
    cdna_change = []

    for line in open(input_file):
        if line.startswith('#'):
            continue
        else:
            record = line.strip().split('\t')
            consequence_str = record[6].split(",")[0]
            if (consequence_str == "frameshift_variant" and record[0].split("_")[-1].split("/")[0] == "-") or (consequence_str == "inframe_insertion"):
                loc = record[1]
                alle = record[2]
                tran_n = record[4]
                cons = record[6].split(',')[0]
                pro_pos = record[9]
                cds_pos = record[8]
                if record[10] == '*':
                    ref_aa = '*'
                    alt_aa = '*'
                else:
                    ref_aa = record[10].split('/')[0]
                    alt_aa = record[10].split('/')[-1]
                ext = record[13]
                cdna_change_p = record[7]
                cdna_c = record[11]
                animo_acid_c = record[10]
                location.append(loc)
                allele.append(alle)
                transcript_name.append(tran_n)
                consequence.append(cons)
                protein_position.append(pro_pos)
                cds_position.append(cds_pos)
                ref_animo_acid.append(ref_aa)
                alt_animo_acid.append(alt_aa)
                extra.append(ext)
                animo_acid_change.append(animo_acid_c)
                cdna_position.append(cdna_change_p)
                cdna_change.append(cdna_c)
            else:
                continue
    gene_symbol = []
    strand = []
    for line in extra:
        element = line.split(';')
        for ele in element:
            sub_ele = ele.split('=')
            if sub_ele[0] == "SYMBOL":
                g_s = sub_ele[1]
            elif sub_ele[0] == "STRAND":
                st = sub_ele[1]
        gene_symbol.append(g_s)
        strand.append(st)
    transcript_seq = []
    for name in transcript_name:
        if name not in transcript_aa.keys():
            seq = 'NULL'
        else:
            seq = transcript_aa[name]
        transcript_seq.append(seq)
    df_list = []
    mut_piptide = []
    wt_piptide = []
    mt_header = []
    wt_header = []
    window_sizes = [15, 16, 17, 18]
    result_sequences = []
    for i in range(len(location)):
        if transcript_seq[i] == "NULL":
            continue
        elif consequence[i] == 'inframe_insertion':
            chr_name = location[i].split(':')[0]
            ins_start = int(location[i].split(':')[1].split('-')[0])
            cds_start = int(cds_position[i].split('-')[0])
            protein_change_pos = int(protein_position[i].split('-')[0])
            ref_aa = ref_animo_acid[i]
            alt_aa = alt_animo_acid[i]
            trans_strand = int(strand[i])
            seq = transcript_aa[transcript_name[i]]
            wt_head = '>' + gene_symbol[i] + '_' + cdna_position[i] + '.' + cdna_change[i] + '_' + location[i] + '_' + transcript_name[i]
            mt_head = '>' + gene_symbol[i] + '_' + cdna_position[i] + '.' + cdna_change[i] + '_' + location[i] + '_' + transcript_name[i]
            alt_aa_len = len(alt_aa)
            if int(protein_change_pos) <= 10:
                wt_pt = seq[0:protein_change_pos + 10]
                mt_pt = seq[0:protein_change_pos - 1] + alt_aa + seq[protein_change_pos:protein_change_pos + 11]
                for window_size in window_sizes:
                    for a, b in zip(range(0, protein_change_pos + alt_aa_len),range(window_size, protein_change_pos + alt_aa_len + window_size)):
                        result_sequences.append(mt_pt[a:b])
                df = pd.DataFrame({
                    'Result_Sequence': result_sequences,
                    'Additional_Column': [mt_head] * len(result_sequences)
                })
            elif protein_change_pos > 10 and len(seq) - protein_change_pos <= 10:
                wt_pt = seq[protein_change_pos - 11:len(seq)]
                mt_pt = seq[protein_change_pos - 11:protein_change_pos - 1] + alt_aa + seq[len(seq) - protein_change_pos - alt_aa_len:len(seq)]
                protein_change_pos = len(mt_pt) - (len(seq) - int(protein_change_pos))
                for window_size in window_sizes:
                    for a, b in zip(range(protein_change_pos - window_size, len(mt_pt) - window_size + 1),range(protein_change_pos, len(mt_pt) + 1)):
                        result_sequences.append(mt_pt[a:b])
                df = pd.DataFrame({
                    'Result_Sequence': result_sequences,
                    'Additional_Column': [mt_head] * len(result_sequences)
                })
            else:
                wt_pt = seq[protein_change_pos - 11:protein_change_pos + 10]
                if ref_aa == '-':
                    mt_pt = seq[protein_change_pos - 11:protein_change_pos - 1] + alt_aa + seq[protein_change_pos - 1:protein_change_pos + 10]
                    protein_change_pos = 11
                    for window_size in window_sizes:
                        for a, b in zip(range(protein_change_pos - window_size, protein_change_pos + alt_aa_len),range(protein_change_pos, len(mt_pt))):
                            result_sequences.append(mt_pt[a:b])
                    df = pd.DataFrame({
                        'Result_Sequence': result_sequences,
                        'Additional_Column': [mt_head] * len(result_sequences)
                    })
                else:
                    mt_pt = seq[protein_change_pos - 11:protein_change_pos - 1] + alt_aa + seq[protein_change_pos:protein_change_pos + 11]
                    protein_change_pos = 11
                    for window_size in window_sizes:
                        for a, b in zip(range(protein_change_pos - window_size, protein_change_pos + alt_aa_len),range(protein_change_pos, len(mt_pt))):
                            result_sequences.append(mt_pt[a:b])
                    df = pd.DataFrame({
                        'Result_Sequence': result_sequences,
                        'Additional_Column': [mt_head] * len(result_sequences)
                    })
                df_list.append(df)

        elif consequence[i] == 'frameshift_variant':
            chr_name = location[i].split(':')[0]
            ins_start = int(location[i].split(':')[1].split('-')[0])
            cds_start = int(cds_position[i].split('-')[0])
            protein_change_pos = int(protein_position[i].split('-')[0])
            allele_ins = allele[i]
            alt_aa = alt_animo_acid[i]
            alt_aa_len = len(alt_aa)
            trans_strand = int(strand[i])
            frame_left_num = int(cds_start) % 3
            seq = transcript_aa[transcript_name[i]]
            if protein_change_pos < 11:
                aa_ten_before = seq[0:protein_change_pos - 1]
            else:
                aa_ten_before = seq[protein_change_pos - 11:protein_change_pos - 1]
            wt_head = '>' + gene_symbol[i] + '_' + cdna_position[i] + '.' + cdna_change[i] + '_' + location[i] + '_' + transcript_name[i]
            mt_head = '>' + gene_symbol[i] + '_' + cdna_position[i] + '.' + cdna_change[i] + '_' + location[i] + '_' + transcript_name[i]
            # print wt_head
            if trans_strand == 1:
                nt_all = f_fasta[chr_name][ins_start - frame_left_num:ins_start] + allele_ins + f_fasta[chr_name][ins_start:ins_start + 54 - frame_left_num]
                # print len(nt_all)
                # print aa_ten_before,nt2aa(nt_all)
                mt_pt = aa_ten_before + nt2aa(nt_all)
                wt_pt = seq[protein_change_pos - 11:protein_change_pos - 1 + len(nt2aa(nt_all))]
                protein_change_pos = 11
                for window_size in window_sizes:
                    for a, b in zip(range(protein_change_pos - window_size, protein_change_pos + alt_aa_len),range(protein_change_pos, len(mt_pt))):
                        result_sequences.append(mt_pt[a:b])
                df = pd.DataFrame({
                    'Result_Sequence': result_sequences,
                    'Additional_Column': [mt_head] * len(result_sequences)
                })
            else:
                reverse_ins = nt_reverse(allele_ins)
                nt_all = nt_reverse(f_fasta[chr_name][ins_start:ins_start + frame_left_num]) + reverse_ins + nt_reverse(f_fasta[chr_name][ins_start + frame_left_num - 55:ins_start - 1])
                # print len(nt_all),nt_all,nt2aa(nt_all)
                mt_ptl = aa_ten_before + nt2aa(nt_all)
                protein_change_pos = 11
                for window_size in window_sizes:
                    for a, b in zip(range(protein_change_pos - window_size, protein_change_pos + alt_aa_len),range(protein_change_pos, len(mt_pt))):
                        result_sequences.append(mt_pt[a:b])
                df = pd.DataFrame({
                    'Result_Sequence': result_sequences,
                    'Additional_Column': [mt_head] * len(result_sequences)
                })
            if mt_pt != '':
                df_list.append(df)
                mut_piptide.append(mt_pt)
                wt_piptide.append(wt_pt)
                mt_header.append(mt_head)
                wt_header.append(wt_head)
            else:
                continue
        else:
            continue

    result_df = pd.concat(df_list, ignore_index=True)
    result_df = result_df.drop_duplicates(subset=['Result_Sequence'])
    result_df.to_csv("./wes_mut_result/ins_II.txt", index=False)


    file_path = "./wes_mut_result/ins_II.fasta"
    # 保存成FASTA文件
    with open(file_path, 'w') as file:
        for index, row in result_df.iterrows():
            file.write(row['Additional_Column'] + '\n')
            file.write(row['Result_Sequence'] + '\n')


def process_snvII_data(input_file, hum_ref_file, hum_pep_file):
    # 读取FASTA文件
    f_fasta = Fasta(hum_ref_file)

    # 读取人类肽序列
    transcript_aa = {}
    with open(hum_pep_file, 'r') as file:
        transcript_name = None
        for line in file:
            if line.startswith(">"):
                transcript_name = line.strip().split(' ')[4][11:26]
                transcript_aa[transcript_name] = ''
            else:
                transcript_aa[transcript_name] += line.replace('\n', '')

    # 解析输入文件
    protein_position = []
    cdna_position = []
    extra = []
    trans_name = []
    ref_animo_acid = []
    alt_animo_acid = []
    ref_nucleotide = []
    alt_nucleotide = []
    chrom_pos = []

    with open(input_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            elif line.strip().split('\t')[6] != "missense_variant":
                continue
            else:
                record = line.strip().split('\t')
                chr_p = record[1]
                tran_n = record[4]
                pro_pos = record[9].split('-')[-1]
                ext = record[13]
                alt_aa = record[10].split('/')[1]
                ref_aa = record[10].split('/')[0]
                cdna_p = record[7]
                ref_n = re.findall('[A-Z]', record[11].split('/')[0])[0]
                alt_n = re.findall('[A-Z]', record[11].split('/')[1])[0]
                alt_animo_acid.append(alt_aa)
                ref_animo_acid.append(ref_aa)
                trans_name.append(tran_n)
                protein_position.append(pro_pos)
                extra.append(ext)
                cdna_position.append(cdna_p)
                ref_nucleotide.append(ref_n)
                alt_nucleotide.append(alt_n)
                chrom_pos.append(chr_p)

    # 提取基因符号
    gene_symbol = []
    for line in extra:
        element = line.split(';')
        for ele in element:
            sub_ele = ele.split('=')
            if sub_ele[0] == "SYMBOL":
                g_s = sub_ele[1]
        gene_symbol.append(g_s)

    # 生成转录本序列
    transcript_seq = []
    for name in trans_name:
        if name not in transcript_aa.keys():
            seq = 'NULL'
        else:
            seq = transcript_aa[name]
        transcript_seq.append(seq)

    # 生成突变肽序列
    window_sizes = [15, 16, 17, 18]
    result_sequences = []
    df_list = []
    mut_peptide = []
    wt_peptide = []
    mt_header = []
    wt_header = []

    for i in range(len(trans_name)):
        if transcript_seq[i] == "NULL":
            continue
        else:
            pro_change_pos = int(protein_position[i])
            wt_head = f'>{gene_symbol[i]}_{alt_animo_acid[i]}_{ref_nucleotide[i]}{cdna_position[i]}{alt_nucleotide[i]}_{chrom_pos[i]}_{trans_name[i]}'
            mt_head = f'>{gene_symbol[i]}_{alt_animo_acid[i]}_{ref_nucleotide[i]}{cdna_position[i]}{alt_nucleotide[i]}_{chrom_pos[i]}_{trans_name[i]}'
            ref_animo_acid_seq = transcript_seq[i]

            if pro_change_pos <= 10:
                wt_pt = ref_animo_acid_seq[0:56]
                mt_pt = ref_animo_acid_seq[0:pro_change_pos - 1] + alt_animo_acid[i] + ref_animo_acid_seq[pro_change_pos:56]
                for window_size in window_sizes:
                    for a, b in zip(range(0, int(pro_change_pos)),
                                    range(window_size, int(pro_change_pos) + window_size)):
                        result_sequences.append(mt_pt[a:b])
                df = pd.DataFrame({
                    'Result_Sequence': result_sequences,
                    'Additional_Column': [mt_head] * len(result_sequences)
                })
            elif pro_change_pos > 10 and len(ref_animo_acid_seq) - pro_change_pos <= 10:
                wt_pt = ref_animo_acid_seq[len(ref_animo_acid_seq) - 56:len(ref_animo_acid_seq)]
                mt_pt = ref_animo_acid_seq[len(ref_animo_acid_seq) - 56:pro_change_pos - 1] + alt_animo_acid[
                    i] + ref_animo_acid_seq[pro_change_pos:len(ref_animo_acid_seq)]
                protein_change_pos = len(mt_pt) - (len(seq) - int(protein_change_pos))
                for window_size in window_sizes:
                    for a, b in zip(range(protein_change_pos - window_size, len(mt_pt) - window_size + 1),
                                    range(protein_change_pos, len(mt_pt) + 1)):
                        result_sequences.append(mt_pt[a:b])
                df = pd.DataFrame({
                    'Result_Sequence': result_sequences,
                    'Additional_Column': [mt_head] * len(result_sequences)
                })
            else:
                wt_pt = ref_animo_acid_seq[pro_change_pos - 11:pro_change_pos + 10]
                mt_pt = ref_animo_acid_seq[pro_change_pos - 11:pro_change_pos - 1] + alt_animo_acid[
                    i] + ref_animo_acid_seq[pro_change_pos:pro_change_pos + 10]
                protein_change_pos = 11
                for window_size in window_sizes:
                    for j in range(len(mt_pt) - window_size + 1):
                        if j <= int(protein_change_pos) - 1 < j + window_size:
                            result_sequences.append(mt_pt[j:j + window_size])
                df = pd.DataFrame({
                    'Result_Sequence': result_sequences,
                    'Additional_Column': [mt_head] * len(result_sequences)
                })
            df_list.append(df)
            mt_header.append(mt_head)
            wt_header.append(wt_head)
            mut_peptide.append(mt_pt)
            wt_peptide.append(wt_pt)

    # 保存结果
    result_df = pd.concat(df_list, ignore_index=True)
    result_df = result_df.drop_duplicates(subset=['Result_Sequence'])
    result_df.to_csv("./wes_mut_result/snv_II.txt", index=False)

    with open("./wes_mut_result/snv_II.fasta", 'w') as file:
        for index, row in result_df.iterrows():
            file.write(row['Additional_Column'] + '\n')
            file.write(row['Result_Sequence'] + '\n')



if __name__ == '__main__':
    con_input1 = sys.argv[1]
    con_input2 = sys.argv[2]
    case_input1 = sys.argv[3]
    case_input2 = sys.argv[4]
    input_file = "./wes_mut_result/result.txt"
    hum_ref_file = "./reference/hg38/hg38.fa"
    hum_pep_file = "./reference/hg38/Homo_sapiens.GRCh38.pep.all.fa"

    wes_processing(con_input1, con_input2)
    wes_processing(case_input1, case_input2)
    call_mutation(con_input1, con_input2, case_input1, case_input2)
    process_snv_data(input_file, hum_ref_file, hum_pep_file)
    process_ins_data(input_file, hum_ref_file, hum_pep_file)
    process_del_data(input_file, hum_ref_file, hum_pep_file)
    process_snvII_data(input_file, hum_ref_file, hum_pep_file)
    process_insII_data(input_file, hum_ref_file, hum_pep_file)
    process_delII_data(input_file, hum_ref_file, hum_pep_file)
    print('Model1 Run Successful!')
    ### python model1_rnaseq_mutation_hla.py con_R1.fastq.gz con_R2.fastq.gz case_R1.fastq.gz case_R2.fastq.gz



