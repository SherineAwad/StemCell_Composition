


rule all:
      input: 
         expand("{path}/CT_conversion/genome_mfa.CT_conversion.fa",  path = config['BISULFITE_PATH']),
         expand("{path}/GA_conversion/genome_mfa.GA_conversion.fa", path = config['BISULFITE_PATH']), 
         expand("{sample}_trimmed_bismark_bt2.bam", sample = config['SAMPLE']), 
         expand("CpG_context_{sample}_trimmed_bismark_bt2.txt", sample = config['SAMPLE']), 
         expand("CHG_context_{sample}_trimmed_bismark_bt2.txt", sample = config['SAMPLE']), 
         expand("CHH_context_{sample}_trimmed_bismark_bt2.txt", sample = config['SAMPLE']), 

rule genome_preparation:
     output: 
         expand("{path}/CT_conversion/genome_mfa.CT_conversion.fa",  path = config['BISULFITE_PATH']), 
         expand("{path}/GA_conversion/genome_mfa.GA_conversion.fa", path = config['BISULFITE_PATH'])    
     params: 
        genome_path = config['GENOME_PATH'] ,
        bisulfite_path = config['BISULFITE_PATH'], 
        threads = config['THREADS'] 
     log: "logs/genome_prepare.log"
     benchmark: "logs/prepare.benchmark"
     shell: 
         """
         bismark_genome_preparation {params.genome_path} --parallel {params.threads} path_to_genome_folder {params.bisulfite_path} 
         """

if config['SINGLE'] == "TRUE": 
       rule trim:
           input:
              "{sample}.fastq",
           params:
              threads= config['THREADS']
           output:
              "galore/{sample}_trimmed.fq.gz",
           conda: 'env/env-trim.yaml'
           shell:
              """
              mkdir -p galore
              mkdir -p fastqc
              trim_galore --gzip --fastqc --fastqc_args "--outdir fastqc" -o galore {input} --cores {params.threads} --rrbs 
              """
       rule run_bismark:
            input:
                 "galore/{sample}_trimmed.fq.gz"
            output:
                 "{sample}_trimmed_bismark_bt2.bam"
            shell:
               """
               bismark . {input}
               """ 
       rule methylation_extractor:
          input:
               "{sample}_trimmed_bismark_bt2.bam"
          output:
              "CpG_context_{sample}_trimmed_bismark_bt2.txt",
              "CHG_context_{sample}_trimmed_bismark_bt2.txt",
              "CHH_context_{sample}_trimmed_bismark_bt2.txt"
          shell:
              """
              bismark_methylation_extractor --comprehensive {input}
              """

 
 
else: 
        rule run_bismark:
            input:
                 "{sample}_r1.fastq",
                 "{sample}_r2.fastq"
            output:
                 "{sample}_bismark_bt2.bam"
            shell:
               """
               bismark . {input}
               """

