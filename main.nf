import groovy.json.JsonSlurper
import groovy.json.JsonOutput


// read config
def jsonSlurper = new JsonSlurper()
config_file = jsonSlurper.parseText(file(params.config).text)

process demux {
    echo true
    cpus 30
    container "nkrumm/nextflow-demux:latest"
    input:
        val run_id from Channel.from(params.run_id)
        file(samplesheet) from Channel.fromPath(params.samplesheet)
    output:
        file("output/**.fastq.gz") into demux_fastq_out_ch
        file("output/Reports")
        file("output/Stats")

        publishDir params.output_path, pattern: 'output/Reports', mode: 'copy', saveAs: {f -> f.replaceFirst("output/", "")}
        publishDir params.output_path, pattern: 'output/Stats', mode: 'copy', saveAs: {f -> f.replaceFirst("output/", "")}
        publishDir params.output_path, pattern: 'output/*.fastq.gz', mode: 'copy' // these are the "Undetermined" fastq.gz files


    script:
        rundir = "inputs/$run_id"
        """
        mkdir -p ${rundir}
        aws s3 sync --only-show-errors s3://uwlm-personal/nkrumm/${run_id} ${rundir}

        if [ -f ${rundir}/Data.tar ]; then
         tar xf ${rundir}/Data.tar -C ${rundir}/
         rm ${rundir}/Data.tar
        fi
        
        bcl2fastq \
            --runfolder-dir ${rundir} \
            --output-dir output/ \
            --sample-sheet ${samplesheet} \
            --min-log-level WARNING \
            --ignore-missing-positions \
            --barcode-mismatches 0,0 \
            --create-fastq-for-index-reads \
            --use-bases-mask ${config_file.basemask} \
            --mask-short-adapter-reads 0 \
            #--tiles s_1_1101
        """
}


demux_fastq_out_ch.flatMap()
    // process output directories from bcl2fastq and re-associate them with config entries
    // we need a key of (lane, project_name, sample_id) to map back to the nested config file.
    .filter { path -> !("${path.getName()}" =~ /^Undetermined_S0_L/) } // first filter ignore Undetermined read files
    .filter { path -> !("${path.getName()}" =~ /I\d_001.fastq.gz$/) } // filter out indexing reads
    .map { path -> 
          def filename = path.getName() // gets filename
          def tokens = path.toString().tokenize('/') // tokenize path
          println (tokens)
          def lane_re = filename =~ /(?!L00)(\d)(?=_R\d_001.fastq.gz)/ // matches L00? in the filename for the lane.
          def lane = lane_re[0][0] // weird groovy way to get first match
          def project = tokens.get(tokens.size() - 2) // Project_Name
          //def sample_id = tokens.get(tokens.size() - 2) // Sample_ID
          def sample_id_re = filename =~ /(.*)(?=_S[\d]*_L00\d_.*.fastq.gz)/
          def sample_id = sample_id_re[0][0]
          println tuple(lane, project, sample_id)
          def key = tuple(lane, project, sample_id)
          return tuple(key, path)
         }
    .groupTuple() // group FASTQ files by key
    .map { key, files -> 
          // attach config information
          def config = config_file.lanes[key[0]][key[1]][key[2]]
          [key, files, config] 
         } 
    .view{ JsonOutput.prettyPrint(JsonOutput.toJson(it[2])) } // diagnostic print
    .branch { 
        // send files to postprocess_ch.umi_true if is_umi is set.
        umi_true: it[2].is_umi
        umi_false: it[2].is_umi==false
    }
    .set{postprocess_ch}


process postprocess_umi {
    echo true
    input:
        set key, file(fastqs), config from postprocess_ch.umi_true
    output:
        set key, file(fastqs), config into umi_out_ch
    script:
        """
        echo UMI $fastqs
        """
}

// route samples to trim process if needed
umi_out_ch.mix(postprocess_ch.umi_false)
    .branch {
        trim_true: it[2].fwd_adapter || it[2].rev_adapter
        trim_false: true // default, if no fwd or rev adapters are specified
    }
    .set { trim_in_ch }


process trim {
    cpus 8
    container 'nkrumm/atropos-paper:latest'
    input:
        set key, file(fastqs), config from trim_in_ch.trim_true
        file(params.adapters)
    output:
        set key, file("trimmed/*.fastq.gz"), config into trim_out_ch
    
    script:
        // TODO: consider using --stats and --report-file here to get report of read trimming
        def trim_options = config.additional_trim_options ? config.additional_trim_options : ""
        """
        # # -F file://${params.adapters} --no-default-adapters --no-cache-adapters
        mkdir -p trimmed/
        atropos \
          -T ${task.cpus} \
          -a ${config.fwd_adapter} \
          -A ${config.rev_adapter} \
          ${trim_options} \
          -pe1 ${fastqs[0]} \
          -pe2 ${fastqs[1]} \
          -o trimmed/1.fastq.gz \
          -p trimmed/2.fastq.gz \
          > report.txt
        """
}

process fastqc {
    cpus 2
    container 'quay.io/biocontainers/fastqc:0.11.8--1'
    input:
        set key, file(fastqs), config from trim_out_ch.mix(trim_in_ch.trim_false)
    output:
        path "fastqc/*", type:"dir" into fastqc_report_ch
        path("libraries/**.fastq.gz")

    publishDir params.output_path, pattern: "*.html", mode: "copy"
    script:
        lane = key[0]
        readgroup = "${params.fcid}.${lane}.${config.index}-${config.index2}"
        library_path = "libraries/${config.Sample_Name}/${config.library_type}/${readgroup}/"
        fastqc_path = "fastqc/${config.Sample_Name}/${config.library_type}/${readgroup}/"
        """
        mkdir -p ${fastqc_path}
        fastqc --quiet -o ${fastqc_path} ${fastqs[0]} ${fastqs[1]}

        # save files using appropriate paths
        mkdir -p ${library_path}
        mv ${fastqs[0]} ${library_path}/1.fastq.gz
        mv ${fastqs[1]} ${library_path}/2.fastq.gz
        """
}

process multiqc {
    echo true
    container 'quay.io/biocontainers/multiqc:1.7--py_4'
    memory '4 GB'
    cpus 4
    input:
        path('fastqc.*') from fastqc_report_ch.flatMap().collect()
        //file('*') from trim_reads_report_ch.collect()
    output:
       file "multiqc_report.${params.fcid}.html"
    publishDir params.output_path, saveAs: {f -> "multiqc/${f}"}, mode: "copy"
    script:
        """
        multiqc -v -d --filename "multiqc_report.${params.fcid}.html" .
        """
}