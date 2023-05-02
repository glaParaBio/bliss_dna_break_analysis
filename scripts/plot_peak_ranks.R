library(data.table)
library(ggplot2)

get_total_tags <- function(macs_file) {
    tag_line <- system(sprintf("grep '## total tags in treatment:' %s", macs_file), intern=TRUE)
    n_tags <- as.integer(sub('.*: ', '', tag_line))
    return(n_tags)
}

w <- fread(snakemake@input[['W']])
mu <- fread(snakemake@input[['Mu']])

w_n_tags <- get_total_tags(snakemake@input[['W']])
mu_n_tags <- get_total_tags(snakemake@input[['Mu']])

macs <- rbind(mu, w)
macs[, library_id := sub('_peak_\\d+$', '', basename(name))]
macs[, n_tags := ifelse(library_id == 'W', w_n_tags, 
    ifelse(library_id == 'Mu', mu_n_tags, NA))]
stopifnot(!is.na(macs$n_tags))

hl <- list(chrom='BES1_Tb427v10', start=76490, end=76930)
macs[, highlight := ifelse(`#chrom` == hl$chrom & start >= hl$start & end <= hl$end, TRUE, FALSE)]

# Only one peak selected to be highlighted
stopifnot(nrow(macs[highlight == TRUE, .N, by=library_id][N == 1]) == 2)

macs <- macs[order(library_id, -fold_enrichment)]
macs[, order := 1:nrow(.SD), by=library_id]
macs[, n_peaks := .N, library_id]
macs[, percentile := 100 * order/n_peaks]

gg <- ggplot(data=macs, aes(x=percentile, y=fold_enrichment)) +
    geom_point(pch='.') +
    geom_hline(data=macs[highlight == TRUE], aes(yintercept=fold_enrichment), colour='blue', linetype='dashed', size=0.1) +
    geom_point(data=macs[highlight == TRUE], colour='red') +
    geom_label(data=macs[highlight == TRUE], aes(label=sprintf('#%s of %s peaks\nTop %.2f%%', order, n_peaks, percentile)), hjust=-0.05, size=2.5) +
    facet_wrap(~library_id, scale='free_x') +
    ggtitle(sprintf('Ranking of bless peaks with top peak on %s highlighted', hl$chrom)) +
    theme_light() +
    theme(strip.text=element_text(colour='black'), plot.title=element_text(size=8))
ggsave(snakemake@output[['rank']], width=16, height=8, units='cm')

macs <- macs[order(library_id, -`-log10(pvalue)`)]
macs[, order := 1:nrow(.SD), by=library_id]
macs[, n_peaks := .N, library_id]
macs[, percentile := 100 * order/n_peaks]

gg <- ggplot(data=macs, aes(x=percentile, y=`-log10(pvalue)`)) +
    geom_point(pch='.') +
    geom_hline(data=macs[highlight == TRUE], aes(yintercept=`-log10(pvalue)`), colour='blue', linetype='dashed', size=0.1) +
    geom_point(data=macs[highlight == TRUE], colour='red') +
    geom_label(data=macs[highlight == TRUE], aes(label=sprintf('#%s of %s peaks\nTop %.2f%%', order, n_peaks, percentile)), hjust=-0.05, vjust=-0.1, size=2.5) +
    ylab('-log10(pvalue)') +
    facet_wrap(~library_id, scale='free_x') +
    ggtitle(sprintf('Ranking of bless peaks with top peak on %s highlighted', hl$chrom)) +
    theme_light() +
    theme(strip.text=element_text(colour='black'), plot.title=element_text(size=8))
ggsave(snakemake@output[['pvalue']], width=16, height=8, units='cm')

macs <- macs[, pileup_rpm := (pileup / n_tags) * 1e6, by=library_id]
macs <- macs[order(library_id, -pileup_rpm)]
macs[, order := 1:nrow(.SD), by=library_id]
macs[, n_peaks := .N, library_id]
macs[, percentile := 100 * order/n_peaks]

gg <- ggplot(data=macs, aes(x=percentile, y=pileup_rpm)) +
    geom_point(pch='.') +
    geom_hline(data=macs[highlight == TRUE], aes(yintercept=pileup_rpm), colour='blue', linetype='dashed', size=0.1) +
    geom_point(data=macs[highlight == TRUE], colour='red') +
    geom_label(data=macs[highlight == TRUE], aes(label=sprintf('#%s of %s peaks\nTop %.2f%%', order, n_peaks, percentile)), hjust=-0.05, vjust=-0.1, size=2.5) +
    ylab('N reads at pileup') +
    facet_wrap(~library_id, scale='free_x') +
    ggtitle(sprintf('Ranking of bless peaks with top peak on %s highlighted', hl$chrom)) +
    theme_light() +
    theme(strip.text=element_text(colour='black'), plot.title=element_text(size=8))
ggsave(snakemake@output[['pileup_rpm']], width=16, height=8, units='cm')
