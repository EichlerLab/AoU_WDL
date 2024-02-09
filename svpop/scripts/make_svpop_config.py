#!/bin/env python

import json
import argparse


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter) 

parser.add_argument('--input', '-i', type=str, required=True, help='Sample list')

args = parser.parse_args()

with open(args.input, "w") as infile:
    sample_list = [line.rstrip() for line in infile]

config = {}
config['reference'] = "link_data/ref/ref.fa"
config['samplelist'] = {}
config['samplelist']['aou'] = sample_list
config['sampleset'] = {}
config['sampleset']['aou'] = {}
config['sampleset']['aou']['sourcetype'] = "caller" 
config['sampleset']['aou']['sourcename'] = "svpop"
config['sampleset']['aou']['merge'] = {}  
config['sampleset']['aou']['merge']["svindel:ins,del"] = "nr::ro:szro(0.8):exact:match"
config['sampleset']['aou']['name'] = "AoU LR"
config['sampleset']['aou']['description'] = "AoU SVPOP merge"


with open("config.json", "w") as outfile:
    json.dump(config, outfile)