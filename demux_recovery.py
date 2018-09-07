#!/usr/bin/env python3
from Bio import pairwise2
from Bio import SeqIO
from Bio import Seq
import argparse
import gzip
import logging
import os
import sys
import yaml

# MATCH_SCORE = 2
# NON_MATCH_SCORE = 1
# GAP_SCORE = -3
# GAP_EXTENSION_SCORE = -1
# MATCHING_SCORE_FACTOR = 0.85

# I7_ADAPTER = 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

# REQUIRED_MIN_BARCODE_DIST = 2

# REQUIRED_MIN_BARCODE_LENGTH = 3


if __name__ == "__main__":

    ############################ Parameter ############################
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', '--input_fastq', required=True,
                    default=None, help='Input fastq file')
    ap.add_argument('-b', '--samples_file', required=True,
                    default=None, help='Yaml file with barcodes, sample names and adapter')
    ap.add_argument('-o', '--output_dir', required=False,
                    default='.', help='Output directory for the demultiplexed files .')
    ap.add_argument('-c', '--config_file', required=False,
                    default=os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),
                                         'demux_recovery.yaml'))
    ap.add_argument('-l', '--log_level', required=False,
                    help ='log level, DEBUG or INFO default is INFO',
                    default='INFO')
    args = vars(ap.parse_args())

    if args['log_level'] == 'INFO':
        log_level = logging.INFO
    else:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level)
    logger = logging.getLogger('log')

    if not os.path.isdir(args['output_dir']):
        try:
            os.makedirs(args['output_dir'])
        except OSError:
            logger.critical('Could not create directory {}!'.format(args['output_dir']))
            sys.exit(1)
    if not os.path.isfile(args['input_fastq']):
        logger.critical('Input file {} does not exist!'.format(args['input_fastq']))
        sys.exit(1)
    if not os.path.isfile(args['samples_file']):
        logger.critical('Demux setup file {} does not exist!'.format(args['samples_file']))

    # Read barcodes and adapter
    with open(args['samples_file'], 'r') as stream:
        f = yaml.load(stream)
        samples = f['samples']
        I7_ADAPTER = f['i7_adapter']

    # Read config file
    with open(args['config_file'], 'r') as stream:
        cfg = yaml.load(stream)
    # Determine minimal required length of barcodes that is required for demultiplexing
    BARCODE_LENGTH = None
    for barcode_length in range(cfg['REQUIRED_MIN_BARCODE_LENGTH'], len(samples[0]['i7_barcode'])):
        distances = []
        for i in range(0, len(samples)):
            for j in range(i + 1, len(samples)):
                distances.append(
                    sum([samples[i]['i7_barcode'][k] != samples[j]['i7_barcode'][k]
                         for k in range(0, min(len(samples[i]['i7_barcode']),
                                                   barcode_length
                                               )
                                        )
                         ])
                )
        min_dist = min(distances)
        if min_dist >= cfg['REQUIRED_MIN_BARCODE_DIST']:
            BARCODE_LENGTH = barcode_length
            break
    logger.info('Determined minimal barcode length of {}, with a barcode distance of {}'.format(BARCODE_LENGTH, cfg['REQUIRED_MIN_BARCODE_DIST']))



    # Input setup
    if args['input_fastq'][-2:] == 'gz':
        handle = gzip.open(args['input_fastq'], 'rt')
        seqp = SeqIO.parse(handle, 'fastq')
    else:
        seqp = SeqIO.parse(args['input_fastq'], 'fastq')

    # Output setup
    for s in samples:
        filename = os.path.join(args['output_dir'], '{}_{}.fastq'.format(s['name'], s['i7_barcode']))
        s['fastq'] = open(filename, 'wt')
    unassigned_samples = {}
    unassigned_samples['fastq'] = open(os.path.join(args['output_dir'],
                                                    'unassigned.fastq'
                                                    ), 'wt')


    ############################ Processing #############################
    record_counter = 0
    adapter_found_counter = 0
    adapter_position_usable_counter = 0
    barcode_assigned_counter = 0
    for record in seqp:
        record_counter += 1
        al = pairwise2.align.localms(str(record.seq), I7_ADAPTER, cfg['MATCH_SCORE'], cfg['NON_MATCH_SCORE'], cfg['GAP_SCORE'], cfg['GAP_EXTENSION_SCORE'], one_alignment_only=True)
        alignment_score = al[0][2]
        start_position =  al[0][3]
        end_position = al[0][4]
        alignment_visualization_sequence = al[0][0]
        alignment_visualization_i7_adapter = al[0][1]
        if alignment_score >= int(cfg['MATCHING_SCORE_FACTOR'] * cfg['MATCH_SCORE'] * len(I7_ADAPTER)):
            adapter_found_counter += 1
            logger.debug('Alignment score: {}, start: {}'.format(alignment_score, start_position))
            logger.debug(alignment_visualization_sequence)
            logger.debug(alignment_visualization_i7_adapter)
            # Is there enough sequence left once the adapter is cut off?
            if end_position <= len(record.seq) - BARCODE_LENGTH:
                adapter_position_usable_counter += 1
                # Get the barcode
                barcode = str(record.seq[end_position:end_position + BARCODE_LENGTH])
                # Put the barcode into the header
                record.description = record.description[0:-1] + barcode
                # Trim sequence
                #_new_record = SeqRecord(Seq(_sequence, generic_dna), id=_record.id, name=_record.name, description=_record.description, annotations=_record.annotations, letter_annotations=_record.letter_annotations)

                letter_annotations = record.letter_annotations.copy()
                letter_annotations['phred_quality'] = letter_annotations['phred_quality'][0:start_position]
                record.letter_annotations = {}
                record.seq = Seq.Seq(str(record.seq[0:start_position]))
                record.letter_annotations = letter_annotations
                # Matches?
                matches = [s['i7_barcode'][0:BARCODE_LENGTH] == barcode for s in samples]
                try:
                    barcode_index = matches.index(True)
                    assert(matches.count(True) <= 1)
                    fastq_file = samples[barcode_index]['fastq']
                    barcode_assigned_counter += 1
                    logger.debug('Found valid barcode {} in position {}'.format(barcode, end_position))
                except ValueError:
                    # to unassigned
                    fastq_file = unassigned_samples['fastq']
                SeqIO.write(record, fastq_file, 'fastq')
        if record_counter % 10000 == 0:
            logger.info('Sequences processed: {}, with adapters: {} ({}%), demultiplexed: {} ({}%), barcode assigned: {} ({}%)'.format(
        record_counter, adapter_found_counter, adapter_found_counter/record_counter * 100,
        adapter_position_usable_counter, adapter_position_usable_counter/record_counter * 100,
        barcode_assigned_counter, barcode_assigned_counter/record_counter * 100
        ))

    logger.info(
        'Sequences processed: {}, with adapters: {} ({}%), demultiplexed: {} ({}%), barcode assigned: {} ({}%)'.format(
            record_counter, adapter_found_counter, adapter_found_counter / record_counter * 100,
            adapter_position_usable_counter, adapter_position_usable_counter / record_counter * 100,
            barcode_assigned_counter, barcode_assigned_counter / record_counter * 100
        ))


