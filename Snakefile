import pandas

def get_second_in_pair(library_id, ss):
    # Get the second-in-pair fastq file for this library_id. Return None if
    # library is single-end.
    if 'fastq_r2' not in ss.columns:
        return None
    r2 = ss[ss.library_id == library_id].fastq_r2.iloc[0]
    if pandas.isna(r2):
        return None
    if not os.path.exists(r2):
        sys.exit('Second-in-pair file %s does not exist' % r2)
    return r2

genomes= pandas.read_csv(config['genomes'], sep= '\t')
ss = pandas.read_csv(config['ss'], sep= '\t', comment='#')

for genome in ss.genome:
    if genome not in list(genomes.genome):
        sys.exit('Genome "%s" not found in genome table %s' % (genome, config['genomes']))

assert len(list(ss.library_id)) == len(set(list(ss.library_id)))

wildcard_constraints:
    library_id= '|'.join([re.escape(x) for x in ss.library_id]),
    strand= '|'.join([re.escape(str(x)) for x in [0, 1, 2]]),
    genome= '|'.join([re.escape(str(x)) for x in ss.genome]),
    feature_type= '|'.join([re.escape(str(x)) for x in ['gene', 'macs', 'notexon']]),


rule all:
    input:
        os.path.join(workflow.basedir, 'results/rRNA_read_summary.tsv'),
        os.path.join(workflow.basedir, 'results/exon_read_summary.tsv'),
        os.path.join(workflow.basedir, 'results/gene.pct_reads_on_strand.pdf'),
        os.path.join(workflow.basedir, 'results/macs.pct_reads_on_strand.pdf'),
        os.path.join(workflow.basedir, 'results/notexon.pct_reads_on_strand.pdf'),
        os.path.join(workflow.basedir, 'results/dedup.pdf'),
        os.path.join(workflow.basedir, 'results/peak_ranks_pileup_rpm.pdf'),
        expand('macs/{library_id}.bed', library_id= ss.library_id),
        expand('bigwig/{library_id}.bigwig', library_id= ss.library_id),
        expand('edger/{genome}.rpkm.tsv.gz', genome= 'TbruceiLister427_2018'),


rule fastqc:
    # TODO
    input:
        # ...

rule download_genome:
    output:
        fa= 'ref/{genome}.fasta',
        fai= 'ref/{genome}.fasta.fai',
    params:
        ftp= lambda wc: genomes[genomes.genome == wc.genome].fasta.iloc[0],
    run:
        if params.ftp.endswith('.gz'):
            gunzip = '| gunzip'
        else:
            gunzip = ''
        shell(r"""
              curl -s -L {params.ftp} %s > {output.fa}
              samtools faidx {output.fa}
              """ % gunzip)

rule download_gff:
    output:
        gff= 'ref/{genome}.gff',
    params:
        ftp= lambda wc: genomes[genomes.genome == wc.genome].gff.iloc[0],
    run:
        if params.ftp.endswith('.gz'):
            gunzip = '| gunzip'
        else:
            gunzip = ''
        shell(r"""
              curl -s -L {params.ftp} %s \
              | {workflow.basedir}/scripts/addGeneIdToGff.py --geneid gene_id -tss '' > {output.gff}
              """ % gunzip)


def get_cutadapt_umi_barcode(sample_sheet, library_id):
    row = sample_sheet[sample_sheet.library_id == library_id]
    umi_position = row.umi_position
    assert len(umi_position) == 1
    umi_position = umi_position.iloc[0].split(':')
    assert len(umi_position) == 2
    umi_position = [int(x) for x in umi_position]
    umi_len = umi_position[1] - umi_position[0]
    n_bar = '^' + 'N' * umi_len
    barcode = row.barcode.iloc[0]
    return n_bar + barcode


rule detect_barcode:
    input:
        fq= lambda wc: ss[ss.library_id == wc.library_id].fastq_r1,
    output:
        info= temp('cutadapt/{library_id}.R1.info.gz'),
        xlog= 'cutadapt/{library_id}.log',
    params:
        umi_barcode= lambda wc: get_cutadapt_umi_barcode(ss, wc.library_id),
    shell:
        r"""
        cutadapt -g '{params.umi_barcode}' --error-rate 0.125 --cores 16 --info-file {output.info} {input.fq} > /dev/null 2> {output.xlog}
        """


rule cutadaptInfoToFastq:
    input:
        info= 'cutadapt/{library_id}.R1.info.gz',
    output:
        fq= temp('fastq/{library_id}.R1.fastq.gz'),
    params:
        umi_pos= lambda wc: ss[ss.library_id == wc.library_id].umi_position.iloc[0].split(':'),
    shell:
        r"""
        pigz -cd {input.info} \
        | {workflow.basedir}/scripts/cutadaptInfoToFastq.py -u {params.umi_pos} -s '#' -i - \
        | pigz > {output.fq}
        """


rule repair_read2:
    input:
        R1= 'fastq/{library_id}.R1.fastq.gz',
    output:
        R2= temp('fastq/{library_id}.R2.fastq.gz'),
    params:
        R2= lambda wc: get_second_in_pair(wc.library_id, ss),
    run:
        if params.R2 is None:
            shell('touch {output.R2}')
        else:
            shell('{workflow.basedir}/scripts/repairFastq.py -r <(pigz -cd {input.R1}) -i <(pigz -cd {params.R2}) -s "#" | pigz > {output.R2}')


rule bwa_index:
    input:
        fa= 'ref/{genome}.fasta',
    output:
        bwt= 'ref/{genome}.fasta.bwt',
    shell:
        r"""
        bwa index {input.fa}
        """


rule bwa:
    input:
        bwt= lambda wc: 'ref/{genome}.fasta.bwt'.format(genome= ss[ss.library_id == wc.library_id].genome.iloc[0]),
        R1= 'fastq/{library_id}.R1.fastq.gz',
        R2= 'fastq/{library_id}.R2.fastq.gz',
    output:
        bam= 'bwa/{library_id}.bam',
    run:
        if os.path.getsize(input.R2) > 0:
            R2 = input.R2
        else:
            R2 = ''
        shell(
            r"""
            genome=`echo {input.bwt} | sed 's/.bwt$//'` 

            bwa mem -t 16 $genome {input.R1} %s \
            | samtools sort -@ 8 > {output.bam}
            samtools index -@ 4 {output.bam}
            """ % R2)


rule dedup:
    input:
        bam= 'bwa/{library_id}.bam',
    output:
        bam= 'dedup/{library_id}.bam',
        xlog= 'dedup/{library_id}.log',
    params:
        stats_prefix= 'dedup/{library_id}',
        paired= lambda wc: '' if get_second_in_pair(wc.library_id, ss) is None else '--paired',
    shell:
        r"""
        umi_tools dedup {params.paired} -I {input.bam} -S {output.bam} --temp-dir . \
            --umi-separator '#' -L {output.xlog} --output-stats {params.stats_prefix}
        samtools index -@ 2 {output.bam}
        """

rule dedup_parse_log:
    input:
        xlog= 'dedup/{library_id}.log',
    output:
        xlog=temp('dedup/{library_id}.tmp'),
    shell:
        r"""
        in=`grep 'Reads: Input Reads: ' {input.xlog} | sed 's/.*: //'`
        out=`grep 'Number of reads out: ' {input.xlog} | sed 's/.*: //'`
        echo "{wildcards.library_id} $in $out" | tr ' ' '\t' > {output.xlog}
        """


rule dedup_summary:
    input:
        xlog= expand('dedup/{library_id}.tmp', library_id=set(ss.library_id)),
        ss=config['ss'],
    output:
        smry=os.path.join(workflow.basedir, 'results/dedup.pdf'),
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(ggplot2)

ss <- fread('{input.ss}')
log <- fread(cmd='cat {input.xlog}', header=FALSE, col.names=c('library_id', 'reads_in', 'reads_out'))
log[, pct_out := 100 * (reads_out/reads_in)]

log <- merge(log, ss[, list(library_id, dataset)], by='library_id')
xdat <- log[, list(pct_out=mean(pct_out)), dataset][order(pct_out)]$dataset
log[, dataset := factor(dataset, xdat)]
xord <- log[order(dataset, pct_out)]$library_id
log[, library_id := factor(library_id, xord)]

lab <- log[, .SD[which.max(pct_out)], dataset]

gg <- ggplot(data=log, aes(x=library_id, y=pct_out, fill=dataset)) +
    geom_col() +
    geom_text(data=lab, aes(label=dataset, colour=dataset), hjust=-0.05) +
    coord_flip() +
    xlab('') +
    ylab('% reads left after deduplication') +
    ylim(0, 100) +
    theme_light() +
    theme(legend.position='none')
ggsave('{output.smry}', width=18, height=16, units='cm')

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """

rule feature_read_count:
    input:
        bam= 'dedup/{library_id}.bam',
        gff= lambda wc: 'ref/%s.gff' % ss[ss.library_id == wc.library_id].genome.iloc[0],
    output:
        cnt= temp('{feature}/{library_id}.tsv'),
        smry= temp('{feature}/{library_id}.tsv.summary'),
    params:
        paired= lambda wc: '' if get_second_in_pair(wc.library_id, ss) is None else '-p',
    shell:
        r"""
        featureCounts -a {input.gff} -o {output.cnt} {params.paired} \
            -t {wildcards.feature} -g Parent --largestOverlap -Q 0 -T 8 {input.bam}
        """

rule feature_size:
    input:
        gff= lambda wc: 'ref/{genome}.gff',
    output:
        txt=temp('{feature}/{genome}.txt'),
    shell:
        r"""
        awk '$3 == "{wildcards.feature}"' {input.gff} \
        | sortBed \
        | mergeBed \
        | awk '{{tot += $3 - $2}}END{{print tot}}' > {output.txt}
        """

rule feature_count_summary:
    input:
        smry= expand('{{feature}}/{library_id}.tsv.summary', library_id= ss.library_id),
        rrna_size=expand('{{feature}}/{genome}.txt', genome=set(ss.genome)),
        ss=config['ss'],
        fai=expand('ref/{genome}.fasta.fai', genome=set(ss.genome)),
    output:
        smry= os.path.join(workflow.basedir, 'results/{feature}_read_summary.tsv'),
    run:
        smry = []
        for fin in input.smry:
            tsv = pandas.read_csv(fin, sep= '\t')
            assert len(tsv.columns)
            bam_id = tsv.columns[1]
            tot = sum(tsv[bam_id]) 
            tsv['total_reads'] = tot
            tsv['pct'] = 100 * tsv[bam_id] / tot
            tsv['library_id'] = re.sub('.bam$', '', os.path.basename(bam_id))
            tsv.rename(columns= {bam_id: 'count'}, inplace= True)
            tsv = tsv[tsv.Status == 'Assigned']
            smry.append(tsv[['library_id', 'total_reads', 'count', 'pct']])
        smry = pandas.concat(smry)
        ss = pandas.read_csv(input.ss, sep='\t', comment='#')
        feature_size = []
        genome_size = []
        for idx, row in smry.iterrows():
            genome = ss[ss.library_id == row.library_id].genome.iloc[0]
            rs = open('%s/%s.txt' % (wildcards.feature, genome)).readlines()
            assert len(rs) == 1
            rs = int(rs[0].strip())
            feature_size.append(rs)
            
            gs = pandas.read_csv('ref/%s.fasta.fai' % genome, header=None, sep='\t')
            gs = sum(gs[1])
            genome_size.append(gs)

        smry['size'] = feature_size
        smry['genome_size'] = genome_size
        smry['pct_genome'] = 100 * smry['size'] / smry['genome_size']
        smry['enrichment'] = smry['pct'] / smry['pct_genome'] 

        smry.sort_values('pct', inplace= True, ascending= False)
        smry.to_csv(output.smry, sep= '\t', index= False)


rule macs_bed_to_gff:
    input:
        bed='macs/{library_id}.bed',
    output:
        gff=temp('macs/{library_id}.gff'),
    shell:
        r"""
        awk -v OFS='\t' -v FS='\t' '$1 !~ "^#" && $6 > 10 && $8 > 6 {{
            print $1, "macs", "peak", $2+1, $3, ".", "+", ".", "ID="$10 ";pileup="$6 ";fold_enrichment=" $8
        }}' {input.bed} > {output.gff}
        """

rule macs_not_exon_to_gff:
    input:
        bed='macs/{library_id}.bed',
        gff=lambda wc: 'ref/%s.gff' % ss[ss.library_id == wc.library_id].genome.iloc[0],
    output:
        gff=temp('notexon/{library_id}.gff'),
    shell:
        r"""
        awk '$3 == "exon"' {input.gff} \
        | intersectBed -a {input.bed} -b - -v \
        | awk -v OFS='\t' -v FS='\t' '$1 !~ "^#" && $6 > 10 && $8 > 6 {{
            print $1, "macs", "peak", $2+1, $3, ".", "+", ".", "ID="$10 ";pileup="$6 ";fold_enrichment=" $8
        }}' > {output.gff}
        """

rule strand_mapping:
    input:
        bam= 'dedup/{library_id}.bam',
        gff= lambda wc: 'ref/%s.gff' % ss[ss.library_id == wc.library_id].genome.iloc[0],
    output:
        cnt= temp('gene/{library_id}.{strand}.tsv'),
        smry= temp('gene/{library_id}.{strand}.tsv.summary'),
    params:
        paired= lambda wc: '' if get_second_in_pair(wc.library_id, ss) is None else '-p',
    shell:
        r"""
        featureCounts -a {input.gff} -o {output.cnt} -s {wildcards.strand} {params.paired} \
            -t exon -g gene_id -Q 0 -T 8 {input.bam}
        """


rule macs_strand_mapping:
    input:
        bam= 'dedup/{library_id}.bam',
        gff= 'macs/{library_id}.gff',
    output:
        cnt= temp('macs/{library_id}.{strand}.tsv'),
        smry= temp('macs/{library_id}.{strand}.tsv.summary'),
    params:
        paired= lambda wc: '' if get_second_in_pair(wc.library_id, ss) is None else '-p',
    shell:
        r"""
        featureCounts -a {input.gff} -o {output.cnt} -s {wildcards.strand} {params.paired} \
            -t peak -g ID -Q 0 -T 8 {input.bam}
        """

rule macs_notexon_strand_mapping:
    input:
        bam= 'dedup/{library_id}.bam',
        gff= 'notexon/{library_id}.gff',
    output:
        cnt= temp('notexon/{library_id}.{strand}.tsv'),
        smry= temp('notexon/{library_id}.{strand}.tsv.summary'),
    params:
        paired= lambda wc: '' if get_second_in_pair(wc.library_id, ss) is None else '-p',
    shell:
        r"""
        featureCounts -a {input.gff} -o {output.cnt} -s {wildcards.strand} {params.paired} \
            -t peak -g ID -Q 0 -T 8 {input.bam}
        """

rule merge_strand_mapping:
    input:
        tsv= expand('{{feature_type}}/{library_id}.{strand}.tsv', library_id= ss.library_id, strand= [1, 2]),
    output:
        tsv= temp('{feature_type}.strand.cnt.tsv'),
    run:
        tsv = []
        for x in input.tsv:
            cnt = pandas.read_csv(x, sep= '\t', comment= '#')
            bam_id = cnt.columns[len(cnt.columns) - 1]
            cnt['library_id'] = re.sub('.bam$', '', os.path.basename(bam_id))
            strand = int(os.path.basename(x).split('.')[-2])
            if strand == 1:
                strand = 'same_strand'
            elif strand == 2:
                strand = 'opposite_strand'
            elif strand == 0:
                strand = 'unstranded'
            else:
                raise Exception()
            cnt['strand_count'] = strand
            cnt['strand'] = [x.split(';')[0] for x in cnt.Strand]
            cnt['chrom'] = [x.split(';')[0] for x in cnt.Chr]
            cnt.rename(columns= {bam_id: 'count'}, inplace= True)
            tsv.append(cnt[['library_id', 'Geneid', 'chrom', 'count', 'strand', 'strand_count']])

        tsv = pandas.concat(tsv)
        tsv = pandas.pivot_table(tsv, index= ['library_id', 'Geneid', 'chrom', 'strand'], columns= ['strand_count'], values= 'count').reset_index()

        tsv['tot'] = tsv.same_strand + tsv.opposite_strand
        tsv = tsv[tsv.tot > 0]
        tsv['same_strand_pct'] = 100 * tsv.same_strand / tsv.tot
        tsv.to_csv(output.tsv, index= False, sep= '\t') 


rule strand_mapping_summary:
    input:
        tsv= '{feature_type}.strand.cnt.tsv',
        ss= config['ss'],
    output:
        hist= os.path.join(workflow.basedir, 'results/{feature_type}.pct_reads_on_strand.pdf'),
    params:
        min_depth= config['min_depth'],
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(ggplot2)

ss <- fread('{input.ss}')
tsv <- fread('{input.tsv}')

tsv <- merge(tsv, unique(ss[, list(library_id, dataset)]), by= 'library_id')
xord <- unique(tsv[order(dataset, library_id)]$library_id)
tsv[, library_id := factor(library_id, xord)]

stsv <- tsv[tot > {params.min_depth}]

lab <- stsv[, .SD[which.min(library_id)], dataset]

if('{wildcards.feature_type}' == 'gene') {{
    xlab <- '% reads on gene on the same strand as the gene' 
}} else if('{wildcards.feature_type}' == 'macs') {{
    xlab <- '% reads on forward strand in macs peaks'
}} else if('{wildcards.feature_type}' == 'notexon') {{
    xlab <- '% reads on forward strand in macs peaks not in exons'
}} else {{
    stop(sprintf('Unexpected feature type: %s', '{wildcards.feature_type}'))
}}

gg <- ggplot(data= stsv, aes(x= same_strand_pct, fill= dataset)) +
    geom_histogram(bins= 20) +
    geom_text(data=lab, aes(label=dataset, x=0, y=Inf, colour=dataset), hjust=0.05, vjust=1.5, size=4) +
    facet_wrap(~library_id, scales= 'free_y') +
    xlab(xlab) +
    ylab('N. genes') +
    theme_light() +
    theme(strip.text= element_text(colour= 'black'), legend.position='none')
ggsave('{output.hist}', width= 30, height= 18, units= 'cm')

# stsv[, qq := qqnorm(same_strand_pct, plot.it=FALSE)$x, by=library_id]
# 
# expected_quantiles <- function(x, shape=10){{
#     # ord <- order(x)
#     n <- length(x)
#     qq <- qbeta(seq(1/n, 1-1/n, length.out=n), shape1=shape, shape2=shape)
#     # qq <- qq[ord]
#     return(qq)
# }}
# 
# dummy <- data.table(
#     library_id='dummy',
#     Geneid=1:10000,
#     chrom=NA,
#     strand=NA,
#     tot=rpois(n=10000, lambda=20),
#     genome='dummy'
# )
# dummy[, same_strand := rbinom(n=length(tot), prob=0.5, size=tot)]
# dummy[, opposite_strand := tot - same_strand]
# dummy[, same_strand_pct := 100 * (same_strand/tot)]
# 
# stsv <- rbind(stsv, dummy)
# 
# stsv <- stsv[order(library_id, same_strand_pct)]
# stsv[, qb := expected_quantiles(x=same_strand_pct, shape=median(.SD$tot)), library_id]
# stsv[, shape1 := 10]
# stsv[, shape2 := 10]
# gg <- ggplot(data= stsv, aes(sample=same_strand_pct/100, color=genome, shape1=shape1)) +
#     geom_qq2(distribution=qbeta) + 
#     # stat_qq_line(distribution=qbeta, dparams=list(shape1=10, shape2=10), colour='grey50') +
#     facet_wrap(~library_id) +
#     theme_light() +
#     theme(strip.text= element_text(colour= 'black'))
# ggsave('tmp.png', width= 26, height= 18, units= 'cm')
# 
# gg <- ggplot(data= stsv, aes(x= same_strand_pct/100, y=qb, color=genome)) +
#     geom_abline(slope=1, intercept=0, colour='grey50') +
#     geom_point(pch='.') +
#     facet_wrap(~library_id, scales= 'free_y') +
#     theme_light() +
#     theme(strip.text= element_text(colour= 'black'))
# ggsave('tmp.pdf', width= 26, height= 18, units= 'cm')

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


rule macs_bam:
    input:
        bam= 'dedup/{library_id}.bam',
    output:
        bam= temp('macs/{library_id}.bam'),
    shell:
        r"""
        samtools view -@ 4 -b -q 10 -F 2048 {input.bam} > {output.bam}
        """


rule macs_callpeak:
    input:
        bam= 'macs/{library_id}.bam',
    output:
        bed= 'macs/{library_id}.bed',
    params:
        gs= lambda wc: genomes[genomes.genome == ss[ss.library_id == wc.library_id].genome.iloc[0]].macs2_genome_size.iloc[0],
        fmt= lambda wc: 'AUTO' if get_second_in_pair(wc.library_id, ss) is None else 'BAMPE',
    shell:
        r"""
        macs2 callpeak -f {params.fmt} -t {input.bam} -g {params.gs} --keep-dup all --nomodel --extsize 80 -n macs/{wildcards.library_id}
        {workflow.basedir}/scripts/reformatMacsXls.py macs/{wildcards.library_id}_peaks.xls > {output.bed}

        rm macs/{wildcards.library_id}_peaks.xls 
        rm macs/{wildcards.library_id}_summits.bed
        rm macs/{wildcards.library_id}_peaks.narrowPeak
        """


rule bigwig:
    input:
        bam= 'dedup/{library_id}.bam',
    output:
        bw= 'bigwig/{library_id}.bigwig',
    shell:
        r"""
        bamCoverage --bam {input.bam} -o {output.bw} \
            --binSize 25 \
            --numberOfProcessors 8 \
            --normalizeUsing CPM
        """


rule plot_peak_ranks:
    input:
        W= 'macs/W.bed',
        Mu= 'macs/Mu.bed',
    output:
        rank= os.path.join(workflow.basedir, 'results/peak_ranks.pdf'),
        pvalue= os.path.join(workflow.basedir, 'results/peak_ranks_pvalue.pdf'),
        pileup_rpm= os.path.join(workflow.basedir, 'results/peak_ranks_pileup_rpm.pdf'),
    script:
        os.path.join(workflow.basedir, 'scripts/plot_peak_ranks.R')

rule makeWindowsGff:
    input:
        fai= 'ref/{genome}.fasta.fai',
    output:
        gff= 'ref/{genome}.windows.gff',
    shell:
        r"""
        bedtools makewindows -g {input.fai} -w 500 -i srcwinnum \
        | awk -v OFS='\t' '{{print $1, ".", "window", $2+1, $3, ".", ".", ".", "ID=" $4}}' > {output.gff} 
        """


rule featureCountsWindows:
    input:
        bam= 'dedup/{library_id}.bam',
        gff= lambda wc: 'ref/%s.windows.gff' % ss[ss.library_id == wc.library_id].genome.iloc[0],
    output:
        cnt= temp('featureCounts/{library_id}.windows.tsv'),
    shell:
        r"""
        featureCounts -a {input.gff} -o {output.cnt} \
            -t window -g ID --largestOverlap -Q 0 -T 8 {input.bam}
        """


rule normaliseWindows:
    input:
        cnt= lambda wc:
            expand('featureCounts/{library_id}.windows.tsv', library_id= ss[ss.genome == wc.genome].library_id),
    output:
        rpkm= 'edger/{genome}.rpkm.tsv.gz',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(edgeR)

flist <- strsplit('{input.cnt}', ' ')[[1]]

tsv <- list()
for(x in flist) {{
    cnt <- fread(x)
    col_id <- names(cnt)[ncol(cnt)]
    id <- sub('\\.bam$', '', basename(col_id))
    setnames(cnt, col_id, 'count')
    stopifnot(!grepl('#', cnt$Chr))
    tsv[[length(tsv) + 1]] <- cnt[, list(locus_id= paste(Chr, Start, End, sep= '#'), Length, library_id= id, count)]
}}
tsv <- rbindlist(tsv)
tsv <- dcast(data= tsv, locus_id + Length ~ library_id, value.var= 'count')
Length <- tsv[, list(locus_id, Length)]
tsv[, Length := NULL]
mat <- as.matrix(tsv, rownames= 'locus_id')

y <- DGEList(mat)
y <- calcNormFactors(y)
rpkm <- rpkm(y, gene.length= Length$Length)
coords <- as.data.table(do.call(rbind, strsplit(rownames(rpkm), '#')))
names(coords) <- c('chrom', 'start', 'end')
coords[, start := as.integer(start)]
coords[, end := as.integer(end)]
rpkm <- cbind(coords, rpkm)

rpkm <- rpkm[order(chrom, start, end)]
gz <- gzfile('{output.rpkm}', 'w')
write.table(rpkm, gz, sep= '\t', row.names= FALSE, quote= FALSE)
close(gz)

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


rule get_exons:
    input:
        gff='ref/{genome}.gff',
    output:
        bed=temp('ref/{genome}.exon.bed'),
    shell:
        r"""
        awk '$3 == "exon"' {input.gff} \
        | sortBed \
        | mergeBed > {output.bed}
        """
