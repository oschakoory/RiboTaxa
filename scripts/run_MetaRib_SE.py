#!/usr/bin/python

# -*- coding: utf-8 -*-
__author__ = "Yaxin Xue"
__license__ = "GPL"
__version__ = "3.0"
__affiliation__ = "CBU, University of Bergen"
__email__ = "xue.ethan@gmail.com, yaxin.xue@uib.no"

# import modules

import ConfigParser, argparse
import os, sys, re, shutil
import random
import pandas as pd

class MyParser(argparse.ArgumentParser):
   def error(self, message):
      sys.stderr.write('error: Check your Parameters!\n')
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)

def parse_arg():
    example_text = '''Example:
    python2 run_MetaRib.py -cfg MetaRib.cfg -p path_to_data -b path_to/bwt_indexes -l path_to/reference database
    '''
    parser = argparse.ArgumentParser(description='Constructing ribosomal genes from large scale total RNA meta-transcriptomic data\n',epilog=example_text, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser._optionals.title = 'Mandatory Arguments'
    parser.add_argument('-cfg', required=True, help='MetaRib configure file')
    #parser.add_argument('-n', required=True, help="Number of subsampling to be used", dest="subsampling")
    parser.add_argument('-p', required=True, help="path to data", dest="data")
    #parser.add_argument('-1', required=True, help="Forward file in fastq format", dest="forward")
    #parser.add_argument('-2', required=True, help="Reverse file in fastq format", dest="reverse")
    parser.add_argument('-b', required=True, help="Bowtie indexed database files", dest="bwt")
    parser.add_argument('-l', required=True, help="Reference database file in fasta format", dest="ref")
    parser.add_argument('-o', required=True, help="output of MetaRib", dest="output")
    global PROJECT_DIR, data_dir, EM_REF, EM_BT
    args = parser.parse_args()
    #SAMPLING_NUM = args.subsampling
    #forward_file = args.forward
    #print(forward_file)
    #reverse_file = args.reverse
    data_dir = args.data
    EM_BT = args.bwt
    EM_REF = args.ref
    PROJECT_DIR = args.output
    #print('args.forward = ', args.forward)
    return(parser)

def parse_cfg(config):
    # BASE
    global SAMPLING_NUM, THREAD
    #DATA_DIR = config.get('BASE', 'DATA_DIR')
    #PROJECT_DIR = config.get('BASE', 'OUTPUT')
    SAMPLING_NUM = config.get('METARIB', 'SAMPLING_NUM')
    THREAD = config.getint('BASE','THREAD')
    # EMIRGE
    global MAX_LENGTH, IDENTITY, MEAN_INSERT_SIZE, STD_DEV, EMIRGE_DB, NUM_ITERATION, RAM, MIN_COV
    #EM_PATH = config.get('EMIRGE', 'EM_PATH')
    #EM_PARA = config.get('EMIRGE', 'EM_PARA')
    #EM_REF = config.get('METARIB', 'EM_REF')
    #EM_BT = config.get('METARIB', 'EM_BT')
    MAX_LENGTH = config.get('EMIRGE', 'MAX_LENGTH')
    IDENTITY = config.get('EMIRGE', 'IDENTITY')
    NUM_ITERATION = config.get('EMIRGE', 'NUM_ITERATION')
    MEAN_INSERT_SIZE = config.get('EMIRGE', 'MEAN_INSERT_SIZE')
    STD_DEV = config.get('EMIRGE', 'STD_DEV')
    EMIRGE_DB = config.get('EMIRGE', 'EMIRGE_DB')
    RAM = config.get('BBMAP', 'RAM')
    MIN_COV = config.get('EMIRGE', 'MIN_COV')
   
    # BBTOOL
    #global BM_PATH, MAP_PARA, CLS_PARA
    #BM_PATH = config.get('BBTOOL', 'BB_PATH')
    #MAP_PARA = config.get('BBTOOL', 'MAP_PARA')
    #CLS_PARA = config.get('BBTOOL', 'CLS_PARA')
    # FILTER
    #global MIN_COV, MIN_PER
    #MIN_COV = config.get('METARIB', 'MIN_COV')
    #MIN_PER = config.get('METARIB', 'MIN_PER')

    return(1)


def init(config):
    #data_dir = str(config.get('BASE', 'DATA_DIR'))
    samples_list_path = PROJECT_DIR+'/samples.list.txt'
    samples_list = []
    samples_fq1_path = {}

    for i in open(samples_list_path):
    	sample_id = i.strip()
        samples_list.append(sample_id)
    	all_fq1 = data_dir+'/'+sample_id+'_trimmed.fastq'
    
        fq1_path = data_dir+'/'+sample_id+'_trimmed.fastq'
        samples_fq1_path[sample_id] = fq1_path

    return(samples_list, samples_fq1_path, all_fq1)

def cal_fastq_num(fastq):
    fastq = fastq
    total_number_of_file = 0
    with open(fastq) as f:
        total_number_of_file = sum(1 for _ in f)
    number_of_fq = int(total_number_of_file)/4.0
    return (number_of_fq)

def cal_fa_num(in_fa):
    fa_file = open(in_fa, 'r')
    num_fa = 0
    for inp in fa_file:
        inp = inp.strip()
        if re.match('>', inp):
            num_fa +=1
    return(num_fa)

def parse_fa_ids(fa_file):
    fa_h = open(fa_file)
    fa_ids = []
    for inp in fa_h:
        inp = inp.strip()
        if inp.startswith('>'):
            header = inp.split('>')[1]
            fa_ids.append(str(header))
        else:
            continue
    return(fa_ids)

def subsampling_reads(unmap_fq1):
    curr_dir = os.getcwd()
    sub_fq1 = curr_dir+'/sub.1.fq'
    sampling_num = int(SAMPLING_NUM)
    max_reads = 1000 * sampling_num
    seeds = random.randint(1, 100)
    cmd = ' '.join(['reformat.sh', '-Xmx'+RAM+'g', 'in='+unmap_fq1, '-Xmx'+RAM+'g','out='+sub_fq1, 'sample='+str(sampling_num), 'sampleseed='+str(seeds),'ow=t', 'reads='+str(max_reads), '2> subsample.log'])
    os.system(cmd)
    return(sub_fq1)

def dedup_contig(old_fa, new_fa):
    work_dir = os.getcwd()
    #step1: join fasta
    cur_mg_fa = work_dir+'/current.merged.fasta'
    cmd = 'cat '+old_fa+ ' '+new_fa+' >'+cur_mg_fa
    os.system(cmd)
    #step2: sort by length
    cur_sort_fa = work_dir+'/current.sorted.fasta'
    cmd = ' '.join(['sortbyname.sh', '-Xmx'+RAM+'g', 'in='+cur_mg_fa, 'out='+cur_sort_fa, 'length descending', '2> sort.log'])
    os.system(cmd)
    #step3: keep uniq id
    cur_uniq_fa = work_dir+'/current.uniqname.fasta'
    cmd = ' '.join(['reformat.sh', '-Xmx'+RAM+'g', 'in='+cur_sort_fa, 'out='+cur_uniq_fa, 'uniquenames', '2> rename.log'])
    os.system(cmd)
    #step4: dedup fasta
    all_dedup_fa = work_dir+'/all.dedup.fasta'
    all_dup_fa = work_dir+'/all.dup.fasta'
    cmd = ' '.join(['dedupe.sh', '-Xmx'+RAM+'g', 'in='+cur_uniq_fa, 'out='+all_dedup_fa, 'outd='+all_dup_fa, 'fo=t', 'ow=t', 'c=t', 'mcs=1', 'e=5', 'mid=99', '2> dedupe.log'])
    os.system(cmd)
    return(all_dedup_fa)

def run_align_bbmap(current_iter_fa, unmap_fq1):
    ref = current_iter_fa
    # run bbmap alignment, since we may have duplicates, cannot calculate stats
    cmd = ' '.join(['bbmap.sh', '-Xmx'+RAM+'g', 'in='+unmap_fq1, 'ref='+ref,
    'threads='+str(THREAD), 'minid=0.96', 'maxindel=1', 'minhits=2', 'idfilter=0.98', 'outu=bbmap.unmap.fq', '32bit=t', 'ow=t', 'statsfile=bbmap.statsfile.txt','sortscafs=t', 'scafstats=bbmap.scafstats.txt', 'covstats=bbmap.covstats.txt','2> bbmap.log'])
    os.system(cmd)
    # reformat to two fastq files
    new_unmap_fq1 = os.getcwd()+'/bbmap.unmaped.1.fq'
    cmd = ' '.join(['reformat.sh', '-Xmx'+RAM+'g', 'in=bbmap.unmap.fq', 'out='+new_unmap_fq1,
    '2> deinterleave.log'])
    os.system(cmd)
    # remove unmapped fq
    cmd = 'rm bbmap.unmap.fq'
    os.system(cmd)
    return(new_unmap_fq1)

def run_emirge_and_dedup(sub_fq1, dedup_fa, iter_time):
    cmd = ' '.join(['emirge_amplicon.py', 'emirge_amp/', '-1', sub_fq1,
    '--phred33', '-l', MAX_LENGTH, '-i', MEAN_INSERT_SIZE,'-j', IDENTITY, '-s', STD_DEV, '-c', MIN_COV, '-a', str(THREAD), '-n', NUM_ITERATION, '-f', EM_REF, '-b', EM_BT, '>> iter_'+str(iter_time)+'_emirge.log','2>> iter_'+str(iter_time)+'_emirge.log'])
   # print(cmd)
    os.system(cmd)
    # change to last iteration folder in EMIRGE
    os.chdir('emirge_amp')
    dirs = [d for d in os.listdir('.') if os.path.isdir(d)]
    last_cycle = sorted(dirs, key=lambda x: os.path.getctime(x), reverse=True)[0]
    os.chdir(last_cycle)
    iter_fa = os.getcwd()+'/'+last_cycle+'.cons.fasta'
    # check if it has novel contigs
    all_dedup_fa = dedup_contig(dedup_fa, iter_fa)
    return(all_dedup_fa, iter_fa)

def run_iteration(unmap_fq1, dedup_fa, iter_time, keep_running):
    iter_dir = '/'.join([PROJECT_DIR, 'Iteration', 'iter_'+str(iter_time)])
    if not os.path.isdir(iter_dir):
        os.mkdir(iter_dir)
    os.chdir(iter_dir)
    sub_fq1= ''
    (sub_fq1) = subsampling_reads(unmap_fq1)
    # run emirge and dedup
    all_dedup_fa, iter_fa = run_emirge_and_dedup(sub_fq1, dedup_fa, iter_time)
    # EMIRGE stop: no more new contigs
    if os.stat(iter_fa).st_size == 0:
        keep_running = 0
        new_iter_time = iter_time + 1
        return(unmap_fq1, all_dedup_fa, new_iter_time, keep_running)
    # print fasta stats
    prev_fa_num = cal_fa_num(dedup_fa)
    cur_fa_num = cal_fa_num(all_dedup_fa)
    new_fa_num = cur_fa_num - prev_fa_num
    fa_stat = 'Iteration: '+str(iter_time)+'\tTotal contigs: '+str(cur_fa_num)+'\tNew contigs: '+str(new_fa_num)
    (new_unmap_fq1) = run_align_bbmap(all_dedup_fa, unmap_fq1)
    # calculate unmapped fq size (MB)
    new_unmap_fq_size = os.stat(new_unmap_fq1).st_size/(1024.0*1024)
    old_unmap_fq_size = os.stat(unmap_fq1).st_size/(1024.0*1024)
    iter_stat = fa_stat+'\tunmapped fastq (MB): '+str(round(new_unmap_fq_size,2))
    print(iter_stat)
    curr_unmap_fq_size = old_unmap_fq_size - new_unmap_fq_size
#    if curr_unmap_fq_size <= (0.01*new_unmap_fq_size):
#        keep_running = 0
    # case2: less than 1% novel contigs
    new_unmap_fq_num = cal_fastq_num(new_unmap_fq1)
    old_unmap_fq_num = cal_fastq_num(unmap_fq1)
    curr_unmap_fq_num = old_unmap_fq_num - new_unmap_fq_num
    if curr_unmap_fq_num <= (0.01*new_unmap_fq_num):
        keep_running = 0
    # case3: maximum iteration
    if iter_time == 10:
        keep_running = 0
    new_iter_time = iter_time + 1
    iter_dir = '/'.join([PROJECT_DIR, '/Iteration'])
    os.chdir(iter_dir)
    return (new_unmap_fq1, all_dedup_fa, new_iter_time, keep_running)

def run_last_iteration(unmap_fq1, dedup_fa, iter_time, keep_running):
    iter_dir = '/'.join([PROJECT_DIR,'Iteration','iter_'+str(iter_time)+'_L'])
    if not os.path.isdir(iter_dir):
        os.mkdir(iter_dir)
    os.chdir(iter_dir)
    # case1: rest unmaped reads <= subsamping reads
    if (keep_running == 1):
        # use all unmaped reads
        sub_fq1= unmap_fq1
        # run emrige_amp and dedup
        all_dedup_fa, iter_fa = run_emirge_and_dedup(sub_fq1,dedup_fa, iter_time)

    # case2 and 3: only a few new contigs or reach the last iteration
    if (keep_running == 0):
        num_unmap_fq = cal_fastq_num(unmap_fq1)
        if (num_unmap_fq <= 2.0*float(SAMPLING_NUM)):
            # use all unmaped reads
            sub_fq1 = unmap_fq1
            # run emrige_amp and dedup
            all_dedup_fa, iter_fa = run_emirge_and_dedup(sub_fq1, dedup_fa, iter_time)
        else:
            # increase subsampling reads * 2
            work_dir = os.getcwd()
            sub_fq1 = work_dir+'/sub.1.fq'

            new_sampling_num = 2.0 * float(SAMPLING_NUM)
            max_reads = 100 * new_sampling_num
            cmd = ' '.join(['reformat.sh', '-Xmx'+RAM+'g', 'in='+unmap_fq1, 'out='+sub_fq1, 'sample='+str(new_sampling_num), 'ow=t', 'reads='+str(max_reads), '2> subsample.log'])
            os.system(cmd)
            # run emrige_amp and dedup
            all_dedup_fa, iter_fa = run_emirge_and_dedup(sub_fq1, dedup_fa, iter_time)
    # print iter fasta stat
    prev_fa_num = cal_fa_num(dedup_fa)
    cur_fa_num = cal_fa_num(all_dedup_fa)
    new_fa_num = cur_fa_num - prev_fa_num
    fa_stat = 'Iteration: '+str(iter_time)+'\tTotal contigs: '+str(cur_fa_num)+'\tNew contigs: '+str(new_fa_num)
    print(fa_stat)
    return(all_dedup_fa)


def cal_mapping_stats(samples_list, samples_fq1_path, all_dedup_fa):
    ab_dir = PROJECT_DIR+'/Abundance/'
    if not os.path.isdir(ab_dir):
        os.mkdir(ab_dir)
    else:
        shutil.rmtree(ab_dir)
        os.mkdir(ab_dir)
    os.chdir(ab_dir)
    dedup_ref = ab_dir+('all.dedup.fasta')
    shutil.copy(all_dedup_fa, dedup_ref)
    all_scafstats_path = {}
    all_covstats_path = {}
    # calculate coverage for each samples
    for idx, val in enumerate(samples_list):
        sample_idx = idx
        sample_name = str(val)
        reads1 = samples_fq1_path[sample_name]
        statsfile = sample_name+'.statsfile.txt'
        scafstats = sample_name+'.scafstats.txt'
        covstats = sample_name+'.covstats.txt'
        # run bbmap alignment, but only we only need statistics file, set ozo=f to print all cov info
        cmd = ' '.join(['bbmap.sh', '-Xmx'+RAM+'g', 'in='+reads1, 'ref='+dedup_ref,
        'threads='+str(THREAD), 'minid=0.96', 'maxindel=1', 'minhits=2', 'idfilter=0.98', 'ow=t', '32bit=t', 'statsfile='+statsfile, 'nzo=f',
        'sortscafs=t', 'scafstats='+scafstats, 'covstats='+covstats, '2> run.'+sample_name+'.log'])
        scafstats = os.getcwd()+'/'+scafstats
        covstats = os.getcwd()+'/'+covstats
        os.system(cmd)
        # sava scafsats file path
        all_scafstats_path[sample_name] = scafstats
        all_covstats_path[sample_name] = covstats
    return(all_scafstats_path, all_covstats_path, dedup_ref)

def generate_and_filter_abundance_table(samples_list, all_scafstats_path, all_covstats_path, dedup_ref):
    ab_dir = PROJECT_DIR+'/Abundance/'
    os.chdir(ab_dir)
    fa_ids = parse_fa_ids(dedup_ref)
    all_ab_df = pd.DataFrame()
    all_ab_df['Contig_ID'] = fa_ids
    all_ab_df.set_index('Contig_ID')
    all_keeped_ctgs = []
    for sample_name in samples_list:
        sample_df = pd.read_csv(all_scafstats_path[sample_name], sep = '\t')
        saved_df = pd.DataFrame()
        saved_df['Contig_ID'] = sample_df['#name']
        saved_df['Contig_ID'] = saved_df.Contig_ID.astype(str)
        saved_df.set_index('Contig_ID')
        # extract the % of amb and unamb reads, and sum as the abundance
        per_unamb = sample_df['%unambiguousReads']
        per_amb = sample_df['%ambiguousReads']
        per_all = per_unamb + per_amb
        saved_df[sample_name+'_estab'] = per_all
        # merge two df based on the all_ab column keys
        all_ab_df = all_ab_df.merge(saved_df, 'left')
        # filter fasta by coverage information
        # parse coverage info
        cov_df = pd.read_csv(all_covstats_path[sample_name], sep= '\t')
        min_cov = float(2)
        #min_per = float(80)
        # filter low avg fold and low covered percent
        cov_df_filter = cov_df.loc[(cov_df['Avg_fold'] >=min_cov)]
        keeped_ids = cov_df_filter['#ID'].tolist()
        for ids in keeped_ids:
            if ids not in all_keeped_ctgs:
                all_keeped_ctgs.append(ids)
    # save filtered fasta
    filter_dedup_fa = os.getcwd()+'/all.dedup.filtered.fasta'
    keep_id_f = open('all.keeped.ids.txt','w')
    for inp in all_keeped_ctgs:
        inp = inp.strip()
        keep_id_f.write(inp+'\n')
    keep_id_f.close()
    cmd = ' '.join(['filterbyname.sh', '-Xmx'+RAM+'g', 'in='+dedup_ref, 'names=all.keeped.ids.txt','out='+filter_dedup_fa,'include=t', '2>ft.log'])
    os.system(cmd)
    # remove id text file
    os.remove('all.keeped.ids.txt')
    # save abundance file
    # save filtered abundance file
    raw_fa_num = cal_fa_num(os.getcwd()+'/all.dedup.fasta')
    ft_fa_num = cal_fa_num(os.getcwd()+'/all.dedup.filtered.fasta')
    fa_stat = 'Number of reconstructed contig:'+str(raw_fa_num)
    print(fa_stat)
    filter_ab_df = all_ab_df.loc[all_ab_df['Contig_ID'].isin(all_keeped_ctgs)]
    filter_ab_file = os.getcwd()+'/all.dedup.filtered.est.ab.txt'
    filter_ab_df.to_csv(filter_ab_file, sep='\t', header=True, index = False, float_format='%.5f', na_rep='NaN')
    os.chdir(PROJECT_DIR)
    return(filter_ab_file)

def main():
    # parse config file
    parser = parse_arg()
    args = parser.parse_args()
    config_file = args.cfg
    config = ConfigParser.ConfigParser()
    config.read(config_file)
    run_cfg = parse_cfg(config)
    # INIT
    samples_list, samples_fq1_path, all_fq1 = init(config)
    # build work folder
    work_dir = PROJECT_DIR +'/'
    if not os.path.isdir(work_dir):
        os.mkdir(work_dir)
    else:
        shutil.rmtree(work_dir)
        os.mkdir(work_dir)
    os.chdir(work_dir)
    dedup_fa = os.getcwd()+'/dedup_contigs.fasta'
    open(dedup_fa,'w').close()
    orig_dedup_fa = dedup_fa
    # iteration parameters
    unmap_fq1, iter_time = all_fq1, 1
    keep_running = 1
    max_iter = 10 # set maximum 10 iterations
    keep_running = 1
    iteration_dir = work_dir+'/Iteration/'
    if not os.path.isdir(iteration_dir):
        os.mkdir(iteration_dir)
    # main iterative steps by while loop
    while (keep_running == 1 and iter_time <= max_iter):
        curr_iter_time = iter_time
        print('====START ITERATION '+str(curr_iter_time)+'====')
        (unmap_fq1, dedup_fa, iter_time, keep_running) = run_iteration(unmap_fq1,dedup_fa, iter_time, keep_running)
        print('====FINISH ITERATION '+str(curr_iter_time)+'====')
    # it reaches maximum iteration
    if (curr_iter_time == max_iter):
        keeep_running = 0
    #  run last iteration
    curr_iter_time = iter_time
    print('====START LAST ITERATION '+str(curr_iter_time)+'====')
    all_dedup_fa = run_last_iteration(unmap_fq1, dedup_fa, iter_time, keep_running)
    print('====FINISH ITERATION '+str(curr_iter_time)+'====')
    #print('====START POSTPROCESSING====')
    # calculate mapping stats for each sample
    all_scafstats_path, all_covstats_path, dedup_ref = cal_mapping_stats(samples_list, samples_fq1_path, all_dedup_fa)
    # generate abundance table based on scafstats, and filter by coverage info
    #(filter_ab_file)  = generate_and_filter_abundance_table(samples_list, all_scafstats_path, all_covstats_path, dedup_ref)
    #print('====FINISH POSTPROCESSING====')
    # print final fasta stat
    os.remove(orig_dedup_fa)
    print('====PROGRAM FINISHED!====')
    return()

if __name__ == "__main__":
   main()
