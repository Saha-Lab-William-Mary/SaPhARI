params.data_dir        = ''
params.code_dir        = ''
params.out_dir         = ''
params.pdbpath         = ''
params.ndbpath         = ''
params.dimevalue       = '0.0001'
params.dimmatchamount  = '10'
params.dimcov          = '20'
params.dimpi           = '20'
params.nbevalue        = '0.0001'
params.nbmatchamount   = '10'
params.nbcov           = '20'
params.nbpi            = '20'


process protorf {
    publishDir "${params.out_dir}/amino_acid_orf", mode: 'copy'

    input:
      path contigfile

    output:
      path '*_prot_orfed.fna'

    script:
    """
    prodigal -i ${contigfile} -a ${contigfile.baseName}_prot_orfed.fna
    """
}

process protorfmeta {
    publishDir "${params.out_dir}/amino_acid_orf", mode: 'copy'

    input:
      path contigfile

    output:
      path '*_prot_orfed.fna'

    script:
    """
    prodigal -i ${contigfile} -a ${contigfile.baseName}_prot_orfed.fna -p meta
    """
}

process nucorf {
    publishDir "${params.out_dir}/nucleotide_orf", mode: 'copy'

    input:
      path contigfile

    output:
      path '*_nuc_orfed.fna'

    script:
    """
    prodigal -i ${contigfile} -d ${contigfile.baseName}_nuc_orfed.fna
    """
}

process nucorfmeta {
    publishDir "${params.out_dir}/nucleotide_orf", mode: 'copy'

    input:
      path contigfile

    output:
      path '*_nuc_orfed.fna'

    script:
    """
    prodigal -i ${contigfile} -d ${contigfile.baseName}_nuc_orfed.fna -p meta
    """
}

process formatHeaders {
    publishDir "${params.out_dir}/cleaned_orfs", mode: 'copy'

    input:
      path orf_file

    output:
      path '*_cleaned.fna'

    script:
    """
    awk '/^>/{gsub(/#/, "|"); gsub(/ /, ""); print} !/^>/' ${orf_file} > ${orf_file.baseName}_cleaned.fna
    """
}

process diamond {
    publishDir "${params.out_dir}/protein_DIAMOND_BLAST", mode: 'copy'

    input:
      path orfedfile

    output:
      path '*_b_diamond.tsv'

    script:
    """
    diamond blastp \
      -d ${params.pdbpath} \
      -q ${orfedfile} \
      -o ${orfedfile.baseName}_b_diamond.tsv \
      -f 6 qseqid stitle pident qcovhsp evalue \
      -k ${params.dimmatchamount} \
      -e ${params.dimevalue} \
      --id ${params.dimpi} \
      --query-cover ${params.dimcov}
    """
}

process blastn {
    publishDir "${params.out_dir}/BLASTn", mode: 'copy'

    input:
      path orf_nuc_file

    output:
      path '*_nuc_blasted'

    script:
    """
    blastn \
      -db ${params.ndbpath} \
      -query ${orf_nuc_file} \
      -out ${orf_nuc_file.baseName}_nuc_blasted \
      -outfmt '6 qseqid stitle pident qcovhsp evalue' \
      -max_target_seqs ${params.nbmatchamount} \
      -evalue ${params.nbevalue} \
      -perc_identity ${params.nbpi} \
      -qcov_hsp_perc ${params.nbcov}
    """
}

process proteincleaner {
    publishDir "${params.out_dir}/annotated", mode: 'copy'

    input:
      path blasted_file

    output:
      path '*.tsv'

    script:
    """
    python3 ${params.code_dir}/utils/extract.py --input ${blasted_file}
    """
}

process formatGenbank {
    publishDir "${params.out_dir}/genbanktofasta", mode: 'copy'

    input:
      path genbank_file

    output:
      path '*_genbank.fna'

    script:
    """
    python3 ${params.code_dir}/utils/genbanktofasta.py --input ${genbank_file}
    """
}

workflow n_t_a {
    /* Nucleotide FASTA → ORFs → DIAMOND → cleanup */
    sequence_ch = Channel.fromPath("${params.data_dir}/*")

    porfs = protorf(sequence_ch)
    formatted_orfs = formatHeaders(porfs)
    blasts = diamond(formatted_orfs)
    proteincleaner(blasts)
}

workflow m_n_t_a {
    /* Nucleotide FASTA (5kb–20kb) → ORFs (meta) → DIAMOND → cleanup */
    sequence_ch = Channel.fromPath("${params.data_dir}/*")

    porfs = protorfmeta(sequence_ch)
    formatted_orfs = formatHeaders(porfs)
    blasts = diamond(formatted_orfs)
    proteincleaner(blasts)
}

workflow a_t_a {
    /* GenBank CDS → convert to FASTA → DIAMOND → cleanup */
    sequence_ch = Channel.fromPath("${params.data_dir}/*")

    formatted_genbank = formatGenbank(sequence_ch)
    blasts = diamond(formatted_genbank)
    proteincleaner(blasts)
}

workflow n_b_a {
    /* Nucleotide FASTA → ORFs → BLASTn → cleanup */
    sequence_ch = Channel.fromPath("${params.data_dir}/*")

    norfs = nucorf(sequence_ch)
    formatted_orfs = formatHeaders(norfs)
    blasts = blastn(formatted_orfs)
    proteincleaner(blasts)
}

workflow m_n_b_a {
    /* Nucleotide FASTA (5kb–20kb) → ORFs (meta) → BLASTn → cleanup */
    sequence_ch = Channel.fromPath("${params.data_dir}/*")

    norfs = nucorfmeta(sequence_ch)
    formatted_orfs = formatHeaders(norfs)
    blasts = blastn(formatted_orfs)
    proteincleaner(blasts)
}
