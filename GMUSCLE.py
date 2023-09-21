# python3.8

__author__ = "Peng Zhang"
__copyright__ = "Copyright 2023, HGID, The Rockefeller University"
__license__ = "CC BY-NC-ND 4.0"
__version__ = "08-24-2023"

import os
import time
import gzip
import argparse
import pandas
import plotnine
import random
from plotnine import *
from Bio import SeqIO
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

print('*************************************************')
print('  ###   #   #  #   #   ####   ####  #      ##### ')
print(' #      ## ##  #   #  #      #      #      #     ')
print(' #  ##  # # #  #   #   ###   #      #      ####  ')
print(' #   #  #   #  #   #      #  #      #      #     ')
print('  ###   #   #   ###   ####    ####  #####  ##### ')
print('*************************************************\n')


###
# input parameters
###
TIME_START = time.time()

parser = argparse.ArgumentParser(description="GMUSCLE")
parser.add_argument("--dir", type=str, help="directory of the fastq files")
parser.add_argument("--sample", type=str, help="file of sample list (include path, if not in the GMUSCLE folder)")
parser.add_argument("--readtype", type=str, default='se', choices=['se','pe'], help="single-end or pair-end [default: se]")
parser.add_argument("--genome", type=str, help="file of reference genome (include path, if not in the folder of .py file)")
parser.add_argument("--chr", type=str, help="chromosome (consistent format with the reference genome)")
parser.add_argument("--start", type=int, help="start position of the sequencing region")
parser.add_argument("--end", type=int, help="end position of the sequencing region")
parser.add_argument("--count_cutoff", type=int, default=30, help="read count cuoff [default: 30]")
parser.add_argument("--ratio_cutoff", type=float, default=0.01, help="genotype ratio cuoff [default: 0.01]")
parser.add_argument("--prefix", type=str, default='', help="prefix for the output files [default: none]")

args = parser.parse_args()
DIR = args.dir
SAMPLES = args.sample
READ_TYPE = args.readtype
GENOME = args.genome
CHR = args.chr
START = args.start
END = args.end
COUNT_CUTOFF = args.count_cutoff
RATIO_CUTOFF = args.ratio_cutoff
PREFIX = args.prefix

file_summary = open(PREFIX + 'GMUSCLE_all_summary.txt', 'w')
file_genotype = open(PREFIX + 'GMUSCLE_genotype_matrix.csv', 'w')
file_alignment = open(PREFIX + 'GMUSCLE_genotype_alignment.csv', 'w')

file_summary.write('* * * * * * * * * * * * * * * * * * * * \n*\n')
file_summary.write('* GMUSCLE parameters:\n*\n')
file_summary.write('* FASTQ directory: ' + DIR + '\n')
file_summary.write('* Sample file: ' + SAMPLES + '\n')
file_summary.write('* Read type: ' + READ_TYPE + '\n')
file_summary.write('* Genome file: ' + GENOME + '\n')
file_summary.write('* CHR: ' + CHR + '\n')
file_summary.write('* START: ' + str(START) + '\n')
file_summary.write('* END: ' + str(END) + '\n')
file_summary.write('* Read-count cutoff: ' + str(COUNT_CUTOFF) + '\n')
file_summary.write('* Genotype-ratio cutoff: ' + str(RATIO_CUTOFF) + '\n')
file_summary.write('* Output prefix: ' + PREFIX + '\n')
file_summary.write('*\n* * * * * * * * * * * * * * * * * * * * \n\n\n')


###
# reverse sequence function
###
base_pairing = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
def rev_seq(fwd_seq):
    output_seq = ''
    for fwd_i in range(0, len(fwd_seq)):
        output_seq += base_pairing[fwd_seq[len(fwd_seq) - fwd_i - 1]]
    return output_seq


error_message = ''
try:
    try:
        ###
        # load genome sequence, extract reference sequence, and build BLAST database for reference sequence
        ###
        file_genome = SeqIO.parse(GENOME, 'fasta')
        ref_seq = ''
        for fasta in file_genome:
            chrom, sequence = fasta.id, str(fasta.seq)
            if chrom == CHR:
                ref_seq = sequence[START-1: END]

        file_ref_seq_fasta = open('ref_seq.fasta', 'w')
        file_ref_seq_fasta.write('>ref_seq\n' + ref_seq + '\n')
        file_ref_seq_fasta.close()

        command_blast_makedb_ref_seq = 'makeblastdb -in ref_seq.fasta -dbtype nucl'
        os.system(command_blast_makedb_ref_seq + ' > /dev/null')

        if ref_seq == '':
            error_message = 'error in genomic region.'
    except:
        error_message = 'error in genomic region.'


    ###
    # load fastq files for each sample, and identify major genotypes
    ###
    global_genotype_muttype_dict = defaultdict(str)
    global_genotype_sample_ratio_dict = defaultdict(lambda: defaultdict(str))
    global_genotype_pos_list = list()
    global_genotype_ratio_list = list()
    file_sample_list = open(SAMPLES, 'r')
    sample_list = list()
    for each_sample in file_sample_list:
        sample = each_sample.strip()
        sample_list.append(sample)

        try:
            ###
            # single-end read
            ###
            if READ_TYPE == 'se':
                # fastq format
                filename_fastq = DIR + sample + '.fastq'
                if os.path.exists(filename_fastq):
                    file_fastq = open(filename_fastq, 'r')

                # fastq.gz format
                else:
                    filename_fastq_gz = DIR + sample + '.fastq.gz'
                    if os.path.exists(filename_fastq_gz):
                        file_fastq = gzip.open(filename_fastq_gz, 'rt')
                    else:
                        print('files do not exist.\n')
                        break

            ###
            # pair-end read
            ###
            elif READ_TYPE == 'pe':
                # fastq format
                filename_fastq_r1 = DIR + sample + '_R1.fastq'
                filename_fastq_r2 = DIR + sample + '_R2.fastq'
                filename_fastq_join = DIR + sample + '.fastq'
                if os.path.exists(filename_fastq_r1) and os.path.exists(filename_fastq_r2):
                    command_fastq_join = 'fastq-join ' + filename_fastq_r1 + ' ' + filename_fastq_r2 + ' -o ' + filename_fastq_join
                    os.system(command_fastq_join)
                    file_fastq = open(filename_fastq_join + 'join', 'r')

                # fastq.gz format
                else:
                    filename_fastq_r1_gz = DIR + sample + '_R1.fastq.gz'
                    filename_fastq_r2_gz = DIR + sample + '_R2.fastq.gz'
                    filename_fastq_join = DIR + sample + '.fastq'
                    if os.path.exists(filename_fastq_r1_gz) and os.path.exists(filename_fastq_r2_gz):
                        command_fastq_join = 'fastq-join ' + filename_fastq_r1_gz + ' ' + filename_fastq_r2_gz + ' -o ' + filename_fastq_join
                        os.system(command_fastq_join)
                        file_fastq = open(filename_fastq_join + 'join', 'r')
                    else:
                        print('files do not exist.\n')
                        break
        except:
            error_message = 'error in fastq file. Please check your fastq, and confirm the filenames match with the fastq.'


        ###
        # load raw reads
        ###
        file_fasta = open('temp_' + sample + '_reads.fasta', 'w')
        print('> Sample:', sample)
        total_read_count = 0
        read_count_dict = defaultdict(int)
        for eachline in file_fastq:
            if eachline[0] == '@':
                line_header = eachline.strip()
                line_read = file_fastq.readline().strip()
                line_plus = file_fastq.readline().strip()
                line_qc = file_fastq.readline().strip()
                total_read_count += 1

                try:
                    read_direction = line_header.split(' ')[1][0]
                    # store read sequences in forward strand
                    if read_direction == '1':
                        read_seq = line_read
                    else:
                        read_seq = rev_seq(line_read)
                    read_count_dict[read_seq] += 1

                except:
                    total_read_count += 1
                    read_seq1 = line_read
                    read_seq2 = rev_seq(line_read)
                    read_count_dict[read_seq1] += 1
                    read_count_dict[read_seq2] += 1

        ###
        # sort and filter raw reads by the count, and write remaining reads into fasta
        ###
        read_count_dict_sorted = dict(sorted(read_count_dict.items(), key=lambda item: item[1], reverse=True))
        read_id_seq_dict = defaultdict(str)
        read_index = 0
        passed_read_count = 0
        for each_read in read_count_dict_sorted.keys():
            each_read_count = read_count_dict_sorted[each_read]
            # read filtering by the count
            if each_read_count > COUNT_CUTOFF:
                read_index += 1
                fasta_header = sample + '|read_' + str(read_index) + '|count:' + str(read_count_dict_sorted[each_read])
                read_id_seq_dict[fasta_header] = each_read
                passed_read_count += each_read_count
                file_fasta.write('>' + fasta_header + '\n' + each_read + '\n')

        total_read_count_unique = len(read_count_dict.keys())
        passed_read_count_unique = read_index - 1
        file_summary.write('>> Sample:\t' + sample + '\n\n')
        file_summary.write('# total reads (raw):\t' + str(total_read_count) + '\n')
        file_summary.write('# total reads (pass):\t' + str(passed_read_count) + '\n')
        file_fastq.close()
        file_fasta.close()


        ###
        # BLAST: remaining read sequences to reference sequence, with gap-finding parameters
        ###
        command_blast_reads_to_ref_seq = 'blastn -query temp_' + sample + '_reads.fasta -db ref_seq.fasta '\
                                         '-reward 5 -penalty -4 -gapopen 10 -gapextend 6 ' \
                                         '-outfmt \"6 qseqid gaps mismatch length qstart qend qseq sstart send sseq\" ' \
                                         '-out temp_' + sample + '_reads_blast.txt'
        os.system(command_blast_reads_to_ref_seq)


        ###
        # load BLAST result, and identify major genotypes
        ###
        file_blast = open('temp_' + sample + '_reads_blast.txt', 'r')
        wt_read_count = 0
        genotype_count_dict = defaultdict(int)
        for eachline in file_blast:
            column = eachline.strip().split('\t')
            read_id = column[0]
            gap_size = int(column[1])
            mismatch_size = int(column[2])
            length = int(column[3])
            read_start = int(column[4])
            read_end = int(column[5])
            read_alignment = column[6]
            ref_start = int(column[7])
            ref_end = int(column[8])
            ref_alignment = column[9]
            read_count = int(read_id.split('|')[2].split(':')[1])

            if gap_size == 0 and mismatch_size == 0:
                wt_read_count += 1

            elif gap_size != 0 or mismatch_size != 0:
                # generate genotype
                try:
                    REF = ALT = ''
                    flag = -1
                    for i in range(0, length):
                        if (read_alignment[i] != ref_alignment[i]) and (flag == -1):
                            flag = i
                            REF += ref_alignment[i]
                            ALT += read_alignment[i]
                        else:
                            if flag != -1:
                                if (read_alignment[i] == ref_alignment[i]) and \
                                   (read_alignment[i+1] == ref_alignment[i+1]) and \
                                   (read_alignment[i+2] == ref_alignment[i+2]) and \
                                   (read_alignment[i+3] == ref_alignment[i+3]) and \
                                   (i >= flag + gap_size + mismatch_size):
                                    break
                                else:
                                    REF += ref_alignment[i]
                                    ALT += read_alignment[i]
                except IndexError:
                    pass

                ###
                # generate VCF format for each genotype
                ###
                REF = REF.replace('-', '')
                if not REF:
                    flag = flag-1
                    REF = ref_alignment[flag]
                    ALT = read_alignment[flag] + ALT

                ALT = ALT.replace('-', '')
                if not ALT:
                    flag = flag-1
                    REF = ref_alignment[flag] + REF
                    ALT = read_alignment[flag]

                POS = START + (ref_start-1) + flag
                genotype = CHR + '-' + str(POS) + '-' + REF + '-' + ALT
                genotype_count_dict[genotype] += read_count

                ###
                # mut type
                ###
                muttype = '.'
                if len(REF) == 1 and len(ALT) == 1:
                    muttype = 'snv'
                elif len(REF) == 1 and len(ALT) > 1:
                    muttype = 'ins-' + str(len(ALT) - 1) + 'nt'
                elif len(REF) > 1 and len(ALT) == 1:
                    muttype = 'del-' + str(len(REF) - 1) + 'nt'
                    else:
                    if len(REF) >= len(ALT):
                        muttype = 'indel-' + str(len(REF)) + 'nt'
                    else:
                        muttype = 'indel-' + str(len(ALT)) + 'nt'

                if genotype not in global_genotype_muttype_dict.keys():
                    global_genotype_muttype_dict[genotype] = muttype


        ###
        # sort genotype by count, and output
        ###
        file_plot = open('temp_' + sample + '_plot.txt', 'w')
        file_plot.write('POS\tTYPE\tSIZE\tRATIO\n')
        genotype_count_dict_sorted = dict(sorted(genotype_count_dict.items(), key=lambda item: item[1], reverse=True))
        genotype_index = 0
        genotype_output = ''
        genotype_pos_set = set()
        for each_genotype in genotype_count_dict_sorted.keys():
            genotype_count = genotype_count_dict_sorted[each_genotype]
            genotype_ratio = genotype_count / total_read_count
            muttype = global_genotype_muttype_dict[each_genotype]

            if genotype_ratio > RATIO_CUTOFF:
                genotype_index += 1
                genotype_ratio_str = '%.4f' % genotype_ratio
                global_genotype_sample_ratio_dict[each_genotype][sample] = genotype_ratio_str
                genotype_output += 'genotype_' + str(genotype_index) + '\t' + each_genotype + '\t' + muttype + '\t' + genotype_ratio_str + '\n'

                genotype_item = each_genotype.split('-')
                POS = int(genotype_item[1])

                if muttype == 'snv':
                    indel_type = 'snv'
                    indel_size = '1'
                else:
                    muttype_item = muttype.split('-')
                    indel_type = muttype_item[0]
                    indel_size = muttype_item[1].replace('nt', '')

                global_genotype_pos_list.append(POS)
                global_genotype_ratio_list.append(float(genotype_ratio_str))

                if genotype_ratio >= 0.01:
                    if POS not in genotype_pos_set:
                        genotype_pos_set.add(POS)
                    else:
                        while(1):
                            new_POS = POS + random.random()
                            if new_POS not in genotype_pos_set:
                                POS = new_POS
                                genotype_pos_set.add(POS)
                                break
                    file_plot.write(str(POS) + '\t' + indel_type + '\t' + indel_size + '\t' + str(int(100*genotype_ratio)) + '\n')

        print('> # genotypes:', str(genotype_index), '\n')
        file_summary.write('# wt reads:\t' + str(wt_read_count) + '\n')
        file_summary.write('# genotypes:\t' + str(genotype_index) + '\n')
        file_summary.write(genotype_output + '\n\n')
        file_plot.close()
        os.system('rm temp_' + sample + '_reads.fasta')
        os.system('rm temp_' + sample + '_reads_blast.txt')


    ###
    # output genotypes and their ratios in each sample
    ###
    file_genotype.write('Genotype,Indel,' + ','.join(sample_list) + '\n')
    global_genotype_list = list(global_genotype_sample_ratio_dict.keys())
    global_genotype_list.sort()
    for each_genotype in global_genotype_list:
        global_genotype_output = each_genotype+','+global_genotype_muttype_dict[each_genotype]+','
        for each_sample in sample_list:
            genotype_ratio = global_genotype_sample_ratio_dict[each_genotype][each_sample]
            if not genotype_ratio:
                genotype_ratio = '.'
            global_genotype_output += genotype_ratio + ','
        file_genotype.write(global_genotype_output[0:-1] + '\n')


    ###
    # output genotypes alignment for each sample
    ###
    ref_seq_list = list(ref_seq)

    file_alignment.write('Chrom,' + CHR + ',,Position,')
    for i in range(START, END+1):
        file_alignment.write(str(i)+',')
    file_alignment.write('\n')

    file_alignment.write(',,,Sequence,')
    for j in range(0, len(ref_seq_list)):
        file_alignment.write(str(ref_seq_list[j])+',')
    file_alignment.write('\n\n')

    file_alignment.write('Sample,Genotype,Indel,Proportion\n')
    for each_sample in sample_list:
        for each_genotype in global_genotype_list:
            if global_genotype_sample_ratio_dict[each_genotype][each_sample]:
                genotype_muttype = global_genotype_muttype_dict[each_genotype]
                genotype_ratio = global_genotype_sample_ratio_dict[each_genotype][each_sample]
                genotype_item = each_genotype.split('-')
                POS = int(genotype_item[1])
                REF = genotype_item[2]
                ALT = genotype_item[3]
                mut_seq_list = list(ref_seq)

                try:
                    mut_index = POS - START
                    for m in range(0, len(REF)):
                        mut_seq_list[mut_index + m] = '.'

                    if len(ALT) <= len(REF):
                        for n in range(0, len(ALT)):
                            mut_seq_list[mut_index + n] = ALT[n]
                    else:
                        for n in range(0, len(REF)-1):
                            mut_seq_list[mut_index + n] = ALT[n]
                        mut_seq_list[mut_index + len(REF) - 1] = ALT[n+1:]

                    file_alignment.write(each_sample + ',' + each_genotype + ',' + genotype_muttype + ',' + genotype_ratio + ',')
                    for k in range(0, len(mut_seq_list)):
                        file_alignment.write(str(mut_seq_list[k]) + ',')
                    file_alignment.write('\n')

                except IndexError:
                    pass


    if global_genotype_pos_list:
        ###
        # plot genotypes for each sample (plotnine)
        ###
        POS_min = min(global_genotype_pos_list)-2
        POS_max = max(global_genotype_pos_list)+2
        ratio_max = 100*max(global_genotype_ratio_list)+5
        for each_sample in sample_list:
            filename_plot_in = 'temp_' + each_sample + '_plot.txt'
            filename_plot_out = PREFIX + 'GMUSCLE_genotype_plot_' + each_sample + '.png'

            plot_data = pandas.read_csv(filename_plot_in, sep='\t')
            plot_figure = (plotnine.ggplot(plot_data)
            + geom_col(aes(x='POS', y='RATIO', fill='TYPE'), width=0.3, stat='identity', color='darkgrey', size=0.3, position='dodge')
            + geom_text(aes(x='POS', y='RATIO', label='SIZE'), size=9, angle=22, nudge_y=2)
            + scale_fill_manual({'snv':'#E74C3C', 'del':'#4D6FE6', 'ins':'#E8C943', 'indel':'#16A085'})
            + scale_x_continuous(limits=(POS_min, POS_max))
            + scale_y_continuous(limits=(0, ratio_max))
            + labs(title='Sample: '+each_sample+'', x='Genomic Position ('+CHR+')', y='Proportion (%)', fill='Type (nt)')
            + theme_538()
            + theme(plot_background=element_rect(fill='white'),
                    panel_background=element_rect(fill='white'),
                    title=element_text(size=12, color='black'),
                    axis_title_x=element_text(size=12, color='black'),
                    axis_text_x=element_text(size=9, color='black', angle=22, vjust=1),
                    axis_title_y=element_text(size=12, color='black'),
                    axis_text_y=element_text(size=10, color='black'),
                    legend_title=element_text(size=9, color='black'),
                    legend_text=element_text(size=9, color='black'))
            )
            plot_figure.save(filename_plot_out, width=6, height=3, dpi=300)

            os.system('rm ' + filename_plot_in)


        ###
        # print out summary
        ###
        TIME_COST = int(time.time() - TIME_START)
        print('# Total samples:', len(sample_list))
        print('# Total genotypes:', len(global_genotype_list))
        print('Time cost:', TIME_COST, 'seconds\n')

        file_summary.write('>> TOTAL SUMMARY\n\n')
        file_summary.write('# total samples:\t' + str(len(sample_list)) + '\n')
        file_summary.write('# total genotypes:\t' + str(len(global_genotype_list)) + '\n')
        file_summary.write('time cost:\t' + str(TIME_COST) + ' seconds\n')

    else:
        print('Zero genotype identified, please try to: check the genome assembly, genomic region and fastq; '
              'lower the read-count cutoff and the genotype proportion cutoff, even to zero.\n')

except:
    if error_message != '':
        print('Error encountered: ' + error_message + '\n')
    else:
        print('Error encountered. Please check your input parameters and files.\n')
    print('If error continues, please contact us for feedback and troubleshooting.\n')
