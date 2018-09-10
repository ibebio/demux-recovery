configfile: "config.yaml"

rule all:
    input: directory("Demux")


rule join:
    input: dynamic(directory("Demux.{splitid}"))
    output: directory("Demux")
    shell: """
        echo {input}
        
        """


rule demux_recovery:
    input: undemultiplexed=config['SPLIT_OUT'] + "/undemultiplexed.fastq.{splitid}"
    params: demux_recovery=config['DEMUX_RECOVERY'], samples_file=config['SAMPLES_FILE']
    output: directory("Demux.{splitid}")
    shell: """
    {params.demux_recovery} -i {input} --samples_file {params.samples_file} -o {output}
    """
        


rule split:
    input: config['INPUT_FILE']
    output:
        dynamic(config['SPLIT_OUT'] + "/" + "undemultiplexed.fastq.{splitid}")
    params:
        lines_per_file=str(config['LINES_PER_FILE']), out_dir=config['SPLIT_OUT']
    shell: """
    mkdir -p params.out_dir
    zcat {input} | split -d -l {params.lines_per_file} - {params.out_dir}/undemultiplexed.fastq.
    """    
