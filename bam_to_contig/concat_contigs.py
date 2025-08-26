# Replace contig names with sample_hap and concatenate
import sys
import os

out_fname = sys.argv[1]
contig_fpaths = sys.argv[2:]
    
# read in fasta
with open(out_fname, 'w') as out_f:
    for contig_fpath in contig_fpaths:
        contig_name = os.path.basename(contig_fpath).replace('.fa', '')
        sequences = {}
        current_seq_name = None
        current_seq = []

        # read in fasta
        with open(contig_fpath) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_seq_name and current_seq:
                        sequences[current_seq_name] = ''.join(current_seq).upper()
                    current_seq_name = line.strip().strip('>')
                    current_seq = []
                else:
                    current_seq.append(line)
                if current_seq_name and current_seq: # last sequence
                    sequences[current_seq_name] = ''.join(current_seq).upper()

        for seqname in sequences:
            out_f.write(f">{contig_name}:{seqname}\n{sequences[seqname]}\n")

