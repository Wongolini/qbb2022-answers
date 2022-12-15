import os
import sys 

# usage: python bwa_mem.py /path/to/assembly.fasta /dir/to/fastq/
# assumes you have r1 and r2 for all samples 

if __name__ == "__main__":
    assembly = sys.argv[1]
    print(assembly)
    fastqs = os.listdir(sys.argv[2])
    fastq_path = os.path.abspath(sys.argv[2])
    cmd = 'bwa mem -t 4 {} {} {} > {}'

    fastqs.sort()
   
    for i,f in enumerate(fastqs[::2]):
        f1 = fastq_path+'/'+f
        f2 = fastq_path+'/'+fastqs[i+1]
        sample = f.replace('_1.fastq','.paired.sam')
        cmd_sys = cmd.format(assembly, f1, f2, sample)
        print(cmd_sys)
        os.system(cmd_sys)
