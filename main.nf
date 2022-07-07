#!/usr/bin/env nextflow

// If the user uses the --help flag, print the help text below
params.help = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Blast sequences against a database
    
    Required Arguments:
      --seedfile    list of s3 paths to *.tar.gz files from IMG (not required with filepath)
      --genome      s3 path to an individual *.tar.gz file from IMG (not required with seedfile)
      --prefix      Output prefix (default: ${params.prefix})
      --project     Folder to place analysis outputs (default: ${params.project})
    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Show help message if the user misses a required argument
if ((params.seedfile  == null) && (params.genome  == null)){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 1, "Input genome must be defined using either '--genome' or '--seedfile' parameter. Please choose at least one (but not both)"
}


//Creates working dir
// fixes #1
workingpath = params.outdir + "/" + params.project
workingdir = file(workingpath)

if( !workingdir.exists() ) {
    if( !workingdir.mkdirs() )     {
        exit 1, "Cannot create working directory: $workingpath"
    } 
}    

if(params.prefix){
    workingpath = workingpath + "/" + params.prefix
}


//
if (params.genome && params.seedfile){
   // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 1, "Input genome must be defined using either '--genome' or '--seedfile' parameter. Please choose at least one (but not both)"
}

if(params.seedfile){
    Channel
        .fromPath(params.seedfile)
        .ifEmpty { exit 1, "Cannot find any seed file matching: ${params.seedfile}." }
        .splitCsv(header:true)
        .map{ row -> tuple(row.genome_id, file(row.genome_path))}
        .set { img_genome_tarball_ch }
} else {
    Channel
        .fromPath(params.genome)
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}." }
        .map {it -> tuple(params.prefix, it)}
        .set { img_genome_tarball_ch }
}

/* 
 * 
 */
process img_parser {
    tag "${genome_id}"

    container params.job_container

    publishDir "${workingpath}", mode: 'copy', pattern: "${genome_id}*.csv.gz"

    input:
    tuple val(genome_id), path(genome_tarball) from img_genome_tarball_ch

    output:
    path("${genome_id}*.csv.gz")

    script:
    """
    tar xzf ${genome_tarball}
    rm -rf ${genome_tarball}
    
    # some gffs contain a fasta portion as well
    # split the gff into gff + fasta
    # sort the gff by col1 (contig name) and col 4,5 (start, stop)
    setup_inputs.sh ${genome_id}
    
    img_to_neptune_via_gremlin.py \\
      --gff ${genome_id}/${genome_id}.gff.sorted \\
      --pfam ${genome_id}/${genome_id}.pfam.tab.txt \\
      --prefix ${genome_id}/${genome_id}
    """
}
