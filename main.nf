#!/usr/bin/env nextflow

params {
    reads: Path
    reference: Path
    primers: Path?
    freyja_barcodes: Path?
    freyja_lineages: Path?
    min_depth: Integer = 20
    min_read_length: Integer = 0
    max_read_length: Integer = 0
    medaka_model: String = 'auto'
}

process PREPARE {
    label 'medaka'

    input:
    path reference, stageAs: 'input/*'
    path primers, stageAs: 'input/*'

    output:
    tuple path('reference.fasta'), path('reference.fasta.fai'), emit: reference
    tuple path('primers.bed'), path('amplicons.bed'), emit: scheme

    script:
    def scheme = primers
        ? "primers_to_bed.py ${primers} --fai reference.fasta.fai --primers-out primers.bed --amplicons-out amplicons.bed"
        : "touch primers.bed amplicons.bed"
    """
    case "${reference}" in
        *.gz) gzip -dc ${reference} > reference.fasta ;;
        *) cp -L ${reference} reference.fasta ;;
    esac
    sed -i -E 's/^>([^[:space:]]+).*/>\\1/' reference.fasta
    samtools faidx reference.fasta
    ${scheme}
    """
}

process READS {
    tag "${id}"
    label 'seqkit'
    cpus 2

    input:
    tuple val(id), path(reads, stageAs: 'input/*')
    val min_length
    val max_length

    output:
    tuple val(id), path("${id}.fastq.gz"), emit: reads
    tuple val(id), path("${id}.raw_stats.tsv"), path("${id}.filtered_stats.tsv"), emit: stats

    script:
    def min_arg = min_length > 0 ? "--min-len ${min_length}" : ''
    def max_arg = max_length > 0 ? "--max-len ${max_length}" : ''
    """
    seqkit seq -j ${task.cpus} input/* | seqkit stats -a -T -j ${task.cpus} > ${id}.raw_stats.tsv
    seqkit seq -j ${task.cpus} ${min_arg} ${max_arg} input/* -o ${id}.fastq.gz
    seqkit stats -a -T ${id}.fastq.gz > ${id}.filtered_stats.tsv
    """
}

process ALIGN {
    tag "${id}"
    label 'medaka'
    cpus 4

    input:
    tuple val(id), path(reads)
    tuple path(reference), path(fai)
    tuple path(primers_bed), path(amplicons_bed)

    output:
    tuple val(id), path("${id}.bam"), path("${id}.bam.bai"), emit: bam
    tuple val(id), path("${id}.counts.tsv"), emit: counts

    script:
    """
    minimap2 -a -x map-ont --secondary=no -t ${task.cpus} ${reference} ${reads} | samtools sort -@ ${task.cpus} -o aligned.bam -
    samtools view -b -F 0x904 -o mapped.bam aligned.bam
    printf 'mapped_reads\\t%s\\n' "\$(samtools view -c mapped.bam)" > ${id}.counts.tsv
    if [ -s ${primers_bed} ]; then
        samtools ampliconclip -@ ${task.cpus} --both-ends --soft-clip --no-excluded --tolerance 10 --filter-len 30 \\
            -b ${primers_bed} -o clipped.bam mapped.bam
        samtools sort -@ ${task.cpus} -o ${id}.bam clipped.bam
        printf 'trimmed_reads\\t%s\\n' "\$(samtools view -c ${id}.bam)" >> ${id}.counts.tsv
    else
        mv mapped.bam ${id}.bam
    fi
    samtools index ${id}.bam
    rm -f aligned.bam mapped.bam clipped.bam
    """
}

process MEDAKA {
    tag "${id}"
    label 'medaka'
    cpus 2
    memory 8.GB

    input:
    tuple val(id), path(bam), path(bai), path(probe, stageAs: 'probe/*')
    tuple path(reference), path(fai)
    val model
    val min_depth

    output:
    tuple val(id), path("${id}.consensus.fasta"), emit: consensus
    tuple val(id), path("${id}.depth.tsv.gz"), path("${id}.variants.tsv"), path("${id}.model.txt"), emit: stats

    script:
    """
    if [ ! -w "\${HOME:-/nonexistent}" ]; then export HOME="\$PWD"; fi
    if [ "${model}" = "auto" ]; then
        model_path=\$(medaka tools resolve_model --auto_model variant ${probe} 2>/dev/null | tail -n1 || true)
        if [ -z "\$model_path" ]; then
            fallback=\$(medaka tools list_models | awk -F':' '/^Default variant/ {gsub(/ /, "", \$2); print \$2}')
            model_path=\$(medaka tools resolve_model --model "\$fallback" | tail -n1)
        fi
    else
        model_path=\$(medaka tools resolve_model --model "${model}" | tail -n1)
    fi
    basename "\$model_path" | sed -E 's/(_model_pt)?\\.tar\\.gz\$//' > ${id}.model.txt

    if [ "\$(samtools view -c ${bam})" -gt 0 ]; then
        normalise_depth.py ${bam} medaka.bam --target 200
        medaka inference medaka.bam probs.hdf --model "\$model_path" --threads 2 --min_mapq 0
        medaka vcf probs.hdf ${reference} raw.vcf
        bcftools sort raw.vcf -o sorted.vcf
        medaka tools annotate sorted.vcf ${reference} ${bam} annotated.vcf
    else
        {
            echo '##fileformat=VCFv4.1'
            awk '{printf "##contig=<ID=%s,length=%s>\\n", \$1, \$2}' ${fai}
            printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tSAMPLE\\n'
        } > annotated.vcf
    fi

    samtools depth -aa -J ${bam} | gzip -c > ${id}.depth.tsv.gz
    depth_mask.py ${id}.depth.tsv.gz --fai ${fai} --min-depth ${min_depth} --output low_depth.bed
    filter_variants.py annotated.vcf --mask low_depth.bed --min-depth ${min_depth} --min-qual 20 \\
        --out-vcf filtered.vcf --out-tsv ${id}.variants.tsv --out-mask low_qual.bed
    cat low_depth.bed low_qual.bed | sort -k1,1 -k2,2n > mask.bed
    bcftools view -Oz -o filtered.vcf.gz filtered.vcf
    bcftools index -t filtered.vcf.gz
    bcftools consensus -f ${reference} -m mask.bed filtered.vcf.gz \\
        | awk -v s="${id}" -v n="\$(grep -c '^>' ${reference})" '/^>/ { print (n == 1 ? ">" s : ">" s "|" substr(\$1, 2)); next } { print }' \\
        > ${id}.consensus.fasta
    """
}

process FREYJA {
    tag "${id}"
    label 'freyja'

    input:
    tuple val(id), path(bam), path(bai)
    tuple path(reference), path(fai)
    tuple path(barcodes, stageAs: 'db/*'), path(lineages, stageAs: 'db/*')

    output:
    tuple val(id), path("${id}.freyja.tsv"), emit: demix

    script:
    def lineage_arg = lineages ? "--lineageyml ${lineages} --depthcutoff 1" : ''
    """
    out=${id}.freyja.tsv
    status=ok
    if refname=\$(check_barcode_reference.py ${barcodes} ${reference} 2> refcheck.log); then
        freyja variants ${bam} --ref ${reference} --refname "\$refname" --minq 0 --variants variants.tsv --depths depths.tsv
        if ! awk '\$4 >= 1 { found = 1; exit } END { exit !found }' depths.tsv; then
            status='insufficient coverage'
        elif ! freyja demix variants.tsv depths.tsv --barcodes ${barcodes} ${lineage_arg} --adapt 0.035 --output \$out > demix.log 2>&1 || [ ! -s \$out ]; then
            if grep -q 'Not enough coverage' demix.log; then
                status='insufficient coverage'
            else
                status="failed: \$(grep -iE 'error' demix.log | tail -n1 | cut -c1-200)"
            fi
        fi
    else
        status="skipped: \$(tail -n1 refcheck.log)"
    fi
    if [ "\$status" != ok ]; then printf 'lineages\\t\\nabundances\\t\\n' > \$out; fi
    printf 'status\\t%s\\n' "\$status" >> \$out
    """
}

process STATS {
    tag "${id}"
    label 'medaka'

    input:
    tuple val(id), path(raw_stats), path(filtered_stats), path(counts), path(depth), path(variants), path(model), path(consensus), path(freyja, stageAs: 'freyja/*')
    tuple path(reference), path(fai)
    tuple path(primers_bed), path(amplicons_bed)
    val min_depth

    output:
    path "${id}.summary.csv", emit: summary
    path "${id}.amplicons.csv", optional: true, emit: amplicons
    path "${id}.lineages.csv", optional: true, emit: lineages

    script:
    def freyja_arg = freyja ? "--freyja ${freyja}" : ''
    """
    sample_stats.py --sample ${id} --prefix ${id} --min-depth ${min_depth} \\
        --raw-stats ${raw_stats} --filtered-stats ${filtered_stats} --counts ${counts} \\
        --fai ${fai} --depth ${depth} --amplicons ${amplicons_bed} \\
        --variants ${variants} --model ${model} --consensus ${consensus} ${freyja_arg}
    """
}

process MERGE {
    label 'medaka'

    input:
    path summaries, stageAs: 'summary/*'
    path amplicons, stageAs: 'amplicons/*'
    path lineages, stageAs: 'lineages/*'

    output:
    path '*.csv'

    script:
    """
    merge_tables.py summary.csv summary/*.csv
    if compgen -G 'amplicons/*.csv' > /dev/null; then
        merge_tables.py amplicon_coverage.csv amplicons/*.csv
    fi
    if compgen -G 'lineages/*.csv' > /dev/null; then
        merge_tables.py freyja_lineages.csv lineages/*.csv
    fi
    """
}

workflow {
    main:
    ensureExecutableScripts()
    def samples = findSamples(params.reads)
    if( !samples ) {
        error("No FASTQ files found in ${params.reads} (expected one sub-directory of FASTQ files per sample)")
    }
    log.info("mpxv-freyja: ${samples.size()} sample(s), depth threshold ${params.min_depth}x, primer trimming ${params.primers ? 'on' : 'off'}, Freyja ${params.freyja_barcodes ? 'on' : 'off'}")

    prep = PREPARE(params.reference, params.primers ?: [])
    reads = READS(channel.fromList(samples), params.min_read_length, params.max_read_length)
    aln = ALIGN(reads.reads, prep.reference, prep.scheme)

    ch_probe = channel.fromList(samples.collect { s -> [s[0], s[1][0]] })
    med = MEDAKA(aln.bam.join(ch_probe), prep.reference, params.medaka_model, params.min_depth)

    ch_barcodes = params.freyja_barcodes != null
        ? channel.value([params.freyja_barcodes, params.freyja_lineages ?: []])
        : channel.empty()
    fr = FREYJA(params.freyja_barcodes != null ? aln.bam : channel.empty(), prep.reference, ch_barcodes)

    ch_stats = reads.stats
        .join(aln.counts)
        .join(med.stats)
        .join(med.consensus)
        .join(fr.demix, remainder: true)
        .filter { _id, raw, _filt, _counts, _depth, _variants, _model, consensus, _freyja -> raw != null && consensus != null }
        .map { id, raw, filt, counts, depth, variants, model, consensus, freyja -> [id, raw, filt, counts, depth, variants, model, consensus, freyja ?: []] }
    st = STATS(ch_stats, prep.reference, prep.scheme, params.min_depth)

    merged = MERGE(st.summary.collect(), st.amplicons.collect().ifEmpty([]), st.lineages.collect().ifEmpty([]))

    publish:
    consensus = med.consensus.map { _id, fasta -> fasta }
    tables = merged

    onComplete:
    if( workflow.success ) {
        log.info("Done: ${workflow.outputDir}/summary.csv")
    }
}

output {
    consensus {
        path 'consensus'
    }
    tables {
        path '.'
    }
}

def findSamples(dir) {
    def fastq = /(?i).+\.(fastq|fq)(\.gz)?/
    def samples = []
    files("${dir}/*", type: 'dir').sort { d -> d.name }.each { d ->
        def reads = files("${d}/*", type: 'file').findAll { f -> f.name ==~ fastq }.sort { f -> f.name }
        if( reads && d.name != 'unclassified' ) {
            samples.add([d.name, reads])
        }
    }
    files("${dir}/*", type: 'file').findAll { f -> f.name ==~ fastq }.sort { f -> f.name }.each { f ->
        samples.add([f.name.replaceAll(/(?i)\.(fastq|fq)(\.gz)?$/, ''), [f]])
    }
    return samples
}

def ensureExecutableScripts() {
    files("${projectDir}/bin/*.py").findAll { s -> !s.toFile().canExecute() }.each { s ->
        s.toFile().setExecutable(true, false)
        if( !s.toFile().canExecute() ) {
            error("${s} is not executable; run: chmod +x ${projectDir}/bin/*.py")
        }
    }
}
