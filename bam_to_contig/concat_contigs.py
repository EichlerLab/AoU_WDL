# Replace contig names with sample_hap and concatenate
import sys
import os

out_fname = sys.argv[1]
contig_fpaths = sys.argv[2:]
        
with open(out_fname, 'w') as out_f:
  for contig_fpath in contig_fpaths:
    contig_name = os.path.basename(contig_fpath).replace('.fa', '')
    print(contig_fpath)
    with open(contig_fpath) as f:
      for line in f:
        if line.startswith('>'):
          if '/rc' in line:  # Reverse strand
              out_f.write(f">{contig_name}/rc:{line.strip('>')}")
          else:
              out_f.write(f">{contig_name}:{line.strip('>')}")
        else:
          out_f.write(line)
