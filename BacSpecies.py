# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 14:04:53 2023

@author: Genglin Guo
"""
import argparse
import sys
import pathlib
import multiprocessing
import subprocess
import time

__version__ = '1.0'

def get_argument():
    # Parsers
    parser = argparse.ArgumentParser(description = 'BacSpecies', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser_group_1 = parser.add_argument_group('Input and Output')
    parser_group_2 = parser.add_argument_group('Parameters')

    # Input and output
    parser_group_1.add_argument('-i', '--input', required = True, nargs = '+', type = str, 
                                help = 'Input FASTA file')
    parser_group_1.add_argument('-o', '--outdir', required = False, type = str, default = 'BS_results.txt',
                              help = 'Output directory')
    parser_group_1.add_argument('-r', '--refs', required = False, type = str, default = 'reference_database',
                              help = 'fasta file contain all 16S rRNA sequence from NCBI')

    # Parameters
    parser_group_2.add_argument('-t', '--threads', required = False, type = int, default = min(multiprocessing.cpu_count(), 4), 
                        help = 'Threads to use for BLAST searches')
    parser_group_2.add_argument('-v', '--version', action = 'version', version = 'BacSpecies v' + __version__,
                        help = "Show version number and exit")
    return parser

def check_dependencies():
    # Checks dependencies are available
    dependencies = ['makeblastdb', 'blastn', 'tblastn']
    for i in dependencies:
        try:
            subprocess.check_call(['which', i], stdout = subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            print('Error: could not find %s tool' % i, file=sys.stderr)
            sys.exit(1)

def process_reference(refs):
    # Generate a dict, genbank_number : bacteria_name
    database = {}
    with open(refs, 'rt') as file:
    	for line in file:
    	    if line.startswith('>'):
    	        line = line[1:]
    	        info = line.split(' ')
    	        id = info[0]
    	        name = info[1] + ' ' + info[2]
    	        database[id] = name
    	    else:
    	        continue
    return database

def get_best_result(inputfile, refs, threads, database):
    # Perform blast, find the best match in reference_database
    inpa = pathlib.Path(inputfile).resolve()
    repa = pathlib.Path(refs).resolve()
    command = ['blastn', '-query', repa, '-subject', inpa, '-num_threads', str(threads), '-evalue', '0.00001', '-perc_identity', '99.5',
    		'-outfmt', '6 qseqid sseqid qstart qend sstart send evalue bitscore length pident qlen qseq']
    process = subprocess.run(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out = process.stdout.decode()
    if out == 0:
        print('The sequence file {} has no blast hit, please check the fasta file'.format(inputfile))
    ''' for test
    with open('123.txt', 'wt') as file:
    	file.write(out)
    	'''
    blast_hits = []
    for line in line_iterator(out):
        blast_hits.append(BlastResult(line))
    best_match = []
    best_cov = 0.0
    best_pident = 0.0
    for hit in blast_hits:
        if hit.length < 1000:
            continue
        elif hit.pident >= best_pident and hit.query_cov >= best_cov:
            best_pident = hit.pident
            best_cov = hit.query_cov
            if database[hit.qseqid] in best_match:
                continue
            else:
                best_match.append(database[hit.qseqid])
    return best_match

def line_iterator(line_breaks):
    # Handle the BLAST output and remove the line breaks 
    line = -1
    while True:
        nextline = line_breaks.find('\n', line + 1)
        if nextline < 0:
            break
        yield line_breaks[line + 1:nextline]
        line = nextline

class BlastResult(object):
    # Handle the BLAST output
    def __init__(self, hit_string):
        parts = hit_string.split('\t')
        self.qseqid = parts[0]
        self.length = int(parts[8])
        self.pident = float(parts[9])
        self.query_cov = 100.0 * len(parts[11]) / float(parts[10])
            
def generate_output(outdir):
    # Generate a blank output table file
    if pathlib.Path(outdir).is_file():
        return
    headers = ['Sequence', 'Best_match', 'Problems', 'Extra informations']
    with open(outdir, 'wt') as file:
        file.write('\t'.join(headers))
        file.write('\n')

def output(outdir, best_match, inputfile):
    # Generate output
    if len(best_match) == 1:
        simple_output = str(inputfile) + ': ' + best_match[0]
        line = [str(inputfile), best_match[0]]
    elif len(best_match) == 0:
        simple_output = str(inputfile) + ': ' + 'Unavailable'
        line = [str(inputfile), 'NA', 'Unavailable']
    else:
        simple_output = str(inputfile) + ': ' + 'Multiple results, please check the output file'
        line = [str(inputfile), best_match[0], 'Multiple results']
        for i in best_match[1:]:
            line.append(i)
    print(simple_output)
    table = open(outdir, 'at')
    table.write('\t'.join(line))
    table.write('\n')
    table.close()
        
def main():
    # Initialize
    starttime = time.perf_counter()
    args = get_argument().parse_args()
    check_dependencies()
    # Processing database
    database = process_reference(args.refs)
    generate_output(args.outdir)
    #Start BLAST
    for inputfile in args.input:
        best_match = get_best_result(inputfile, args.refs, args.threads, database)
    #Generate output
        output(args.outdir, best_match, inputfile)
    endtime = time.perf_counter() - starttime
    print('Total time consumed : {:.1f}h{:.1f}m{:.1f}s'.format(endtime // 3600, endtime % 3600 // 60, endtime % 60))
   
main()
