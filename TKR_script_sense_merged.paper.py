import pandas as pd
import pysam
from dataclasses import dataclass,field
import os, argparse
from typing import List
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import matplotlib.colors as clr
import scipy as scp
import tqdm
import pathlib
import warnings
import colorsys
import numpy as np

warnings.simplefilter(action='ignore', category=FutureWarning)

@dataclass
class Exon:
    chrom: str
    start: int
    stop: int
    name: int
    strand: str
    is_trk: str
    coverage: list=field(default_factory=lambda: [0,0])

def calc_coverage(samfile: pysam.AlignmentFile, exon: Exon, paired=False,cov_coef=1, stream="up") -> float:
    """
    Calculate exon coverage
    :param samfile: Alignment data object
    :param exon: Exon object
    :return: coverage of exon
    """
    overlap = 0
    
    for read in samfile.fetch(exon.chrom, exon.start, exon.stop):
        condition_qual = read.mapping_quality > 100 and not read.is_qcfail and read.get_overlap(exon.start,
                                                                                                exon.stop) > 10

        if paired:
            if exon.strand == '+':
                condition_strand = read.is_reverse and read.is_read1 or read.is_forward and read.is_read2
            elif exon.strand == '-':
                condition_strand = read.is_reverse and read.is_read2 or read.is_forward and read.is_read1
        else:
            if exon.strand == '+':
                condition_strand = read.is_reverse
            elif exon.strand == '-':
                condition_strand = read.is_forward

        if condition_qual:
            if  (stream=="up" and condition_strand) or (stream=="down" and not condition_strand):
                overlap += read.get_overlap(exon.start, exon.stop)
    
    coverage = overlap / (exon.stop - exon.start)
    coverage /= cov_coef
    
    return coverage


def calc_gene_coverage(samfile: pysam.AlignmentFile, exon_list: List[Exon],
                       paired=False, cov_coef=1,strandness="both"):
    """
    Calculates coverage for whole gene
    :param samfile: AlignmentFile object for patient
    :param exon_list: list with the Exon objects for gene analyzed
    :return: dict with coverage for each exon in gene
    
    Added normalisation coeficient with total coverage
    """
    res_dict = {}
    for exon in exon_list:
        cov_up,cov_down = 0, 0
        if strandness == "up" or strandness == "both":
            cov_up = calc_coverage(samfile, exon, paired, cov_coef, stream="up")
        if strandness == "down" or strandness == "both":
            cov_down = calc_coverage(samfile, exon, paired, cov_coef, stream="down")
        exon.coverage = [cov_up, cov_down]
        res_dict[exon.name] = [cov_up, cov_down]
    return res_dict


def exons_from_table(table_path, gene_list, no_intron=False):
    """
    Creates dicts with lists of Exons classes
    :param table_path: path to table with exons data
    :return: dict with gene as key and list of Exons as values
    """

    exons_data = pd.read_csv(table_path, sep='\t', index_col=0)
    genes_list = gene_list if gene_list else list(exons_data.gene.unique())
    genes_dict = {}

    for gene in genes_list:

        exons_data_gene = exons_data[exons_data.gene == gene]
        exons_list = []

        i = 1
        for exon in exons_data_gene.iterrows():
            exon = exon[1]
            
            if no_intron and exon.is_TRK not in ["True", "False"]:
                continue

            exons_list.append(Exon(str(exon.chr),
                                   exon.start,
                                   exon.finish,
                                   i,
                                   exon.strand,
                                   exon.is_TRK))
            i += 1
        genes_dict[gene] = exons_list

    return genes_dict


def plot_coverage(gene, name, trk_border):
    """
    Plots coverage of particular gene
    :param gene:
    """

    fig, ax = plt.subplots(figsize=(14/2.54, 10/2.54))

    color_dict = {k:v for k,v in zip(["True", "False"],
                                    ["#D9153F", "#485DF6"])}
    color_dict_sense = {k:v for k,v in zip(["True", "False"],
                                    ["#F7C1C1", "#7BC2FD"])}

    # draw exon borders
    ax.axvspan(*trk_border[:2], alpha=0.2, color='gray',zorder=0,linewidth=0)
    ax.axvspan(*trk_border[2:], alpha=0.2, color='gray',zorder=0,linewidth=0) 

    # for upstream
    cov_up=[exon.coverage[0] for exon in genes_dict[gene]]
    ax.bar(x=[exon.name for exon in genes_dict[gene]],
           height=cov_up,
           color=[color_dict[exon.is_trk] for exon in genes_dict[gene]],zorder=3)

    # for downstream
    cov_down = [exon.coverage[1] for exon in genes_dict[gene]]
    ax.bar(x=[exon.name for exon in genes_dict[gene]],
           height=[-i for i in cov_down],
           color=[color_dict_sense[exon.is_trk] for exon in genes_dict[gene]],zorder=3)
    
    title = f'{gene} coverage plot'
    title += f"    "

    if ttest_res is not None:
        title += f'P_val = {round(ttest_res, 3)} (U test)\n'

    title += f'non-TK/TK coverage = {round(cov_non_trk_raw,3)}/{round(cov_trk_raw,3)}  (antisense reads)'
    title = name + "    " + title if name else title

    plt.title(title)
    ax.xaxis.set_major_locator(tck.MultipleLocator())
    plt.ylabel('Coverage')
    plt.xlabel('Exon number')
    
    #y_min = -0.05*1.1 if max(cov_down) < 0.05 else -max(cov_down)*1.1
    #y_max = 0.5*1.1 if max(cov_up) < 0.5 else max(cov_up)*1.1

    y_min, y_max = -0.2*1.1, 0.5*1.1
    plt.xlim(0.5, len(genes_dict[gene]) + 0.5)

    plt.axhline(0, color='black', linewidth = 1,zorder=0)

    red_patch = plt.Rectangle((0, 0), 0.5, 0.5, fc="#D9153F")
    blue_patch = plt.Rectangle((0, 0), 0.5, 0.5, fc="#485DF6")

    x_labs = [exon.name for exon in genes_dict[gene]]
    leg = ax.legend([red_patch, blue_patch],
        ['TK', 'Non-TK'], handlelength=0.75, handleheight=0.75,
        ncol=1, loc='upper left', labelspacing=0.5, handletextpad=0.5, borderpad=0.25)
    leg.get_frame().set_linewidth(0.0)

    font_items = ([ax.xaxis.label, ax.yaxis.label] +
                  ax.get_xticklabels() + ax.get_yticklabels() +
                  ax.get_legend().get_texts())

    plt.xticks([exon.name for exon in genes_dict[gene]], x_labs, rotation='vertical', fontsize=12)
    plt.yticks(np.arange(round(y_min,1), round(y_max,1), 0.1), fontsize=12)
    ax.title.set_fontsize(16)

    for item in font_items:
        item.set_fontsize(12)

    plt.savefig(f'{args.output}/plots/{gene}/{name}_{gene}.jpg', dpi=300, bbox_inches='tight')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Input and output for transposition calculator')
    parser.add_argument('-i',
                        '--input',
                        type=str,
                        default='.',
                        help='input .bam file'
                        )

    parser.add_argument('-o',
                        '--output',
                        type=str,
                        default='./',
                        help='output folder for TKR detection results'
                        )
     
    parser.add_argument('-p',
                        '--plot_all',
                        action='store_true',
                        help='make plot even if no significant results'
                        )
    parser.add_argument('-s',
                        '--stats',
                        type=str,
                        nargs='?',
                        default=None,
                        help='create file with intron coverate stats. Append, if file exists'
                        )
    parser.add_argument('-n',
                        '--names',
                        type=str,
                        nargs='?',
                        default=None,
                        help='file with sample names'
                        )
    parser.add_argument('-st',
                        '--strandness',
                        type=str,
                        nargs='?',
                        choices=['both', 'sense', 'antisense'],
                        default="both",
                        help='strandness of calculated transcripts'
                        )

    parser.add_argument('-t',
                        '--table',
                        type=str,
                        nargs='?',
                        default=None,
                        help='create file with info on exons coverage. Append, if file exists'
                        )

    parser.add_argument('-g',
                        '--genes',
                        #type=str,
                        nargs='+',
                        default=None,
                        help='calculate asymmetry only for these genes'
                        )

    parser.add_argument('--original',
                       type=str,
                       nargs='?',
                       default=None,
                       help='original BAM file before filtration'
                        )

    args = parser.parse_args()
    script_path = os.path.dirname(os.path.abspath(__file__))
    gene_list = args.genes
    genes_dict = exons_from_table(f'{script_path}/all_tkr_exons_with_flag.merged.tsv', gene_list, no_intron=True)
    samfile = pysam.AlignmentFile(args.input, "rb")
    samfile_original = pysam.AlignmentFile(args.original, "rb") if args.original else samfile
    
    # Get short sample name
    if args.names:
        names = pd.read_csv(args.names, index_col=0, header=None,sep=";").to_dict()[1]
        name = names[args.input.split("/")[-1].replace("Aligned.sortedByCoord.out.bam", "")]
    else:
        name == ""

    # calculate number of reads
    total_reads = samfile_original.mapped
    cov_coef = total_reads / 1000000
    
    samfile_original = None
    res_dict = {}
    paired = '_R2' in samfile.header.get('PG')[0]['CL']
    print(f'Sample is paired: {paired}')

    for gene in gene_list:
        os.makedirs(f'{args.output}/plots/{gene}', exist_ok=True)

        calc_gene_coverage(samfile, genes_dict[gene], paired, cov_coef, args.strandness)
        gene_df = pd.DataFrame(genes_dict[gene])
        gene_df = gene_df.loc[gene_df["is_trk"].isin(["True", "False"])]
        
        gene_df["coverage_sep"] = gene_df["coverage"].copy()
        gene_df["coverage"] = gene_df["coverage"].apply(sum)
        gene_df["coverage_asense"] = gene_df["coverage_sep"].str[0]

        non_trk_exons = []
        trk_exons = []
    
        met_trk = False
        for exon in genes_dict[gene]:
            # if we reached non-TK exons after TK
            if exon.is_trk == "False" and met_trk:
                break
            if exon.is_trk == "False":
                non_trk_exons.append(exon.name -1)
            if exon.is_trk == "True":
                trk_exons.append(exon.name -1)
                met_trk = True

        # if found TK domain
        if len(trk_exons) > 0:
            # selected top-5 of each TK and non-TK
            non_trk_exons = non_trk_exons[1:6]
            trk_exons = trk_exons[:5]
            cov_non_trk = gene_df.loc[non_trk_exons].coverage_asense
            cov_trk = gene_df.loc[trk_exons].coverage_asense

            #ttest_res = scp.stats.ttest_ind(cov_non_trk, cov_trk, alternative='less').pvalue
            ttest_res = scp.stats.mannwhitneyu(cov_non_trk, cov_trk, alternative='less').pvalue
        else:
            ttest_res = None

        res_dict[gene] = False

        # calc mean coverage of TK and non-TK exons
        cov_trk_raw, cov_non_trk_raw = 0, 0
        trk_len, non_trk_len = 0, 0
        for exon in genes_dict[gene]:
            exon_i = exon.name - 1
            if exon_i in trk_exons:
                cov_trk_raw += gene_df.loc[exon_i].coverage_asense * (exon.stop - exon.start) * cov_coef
                trk_len += exon.stop - exon.start
            elif exon_i in non_trk_exons:
                cov_non_trk_raw += gene_df.loc[exon_i].coverage_asense * (exon.stop - exon.start) * cov_coef
                non_trk_len += exon.stop - exon.start
        cov_trk_raw /= trk_len
        cov_non_trk_raw /= non_trk_len

        mean_cov = gene_df.coverage.mean() * cov_coef 
        mean_cov_non_trk, mean_cov_trk = sum(cov_non_trk)/len(cov_non_trk), sum(cov_trk)/len(cov_trk)

        # Потом сделать порог по cov_trk_raw
        if ttest_res is not None and ttest_res < 0.05 and cov_trk_raw > 0:
            if ttest_res < 0.05:
                print(f'Transposition detected in gene {gene}')
                res_dict[gene] = True

        if (ttest_res is not None and ttest_res < 0.05 and gene_df.coverage.mean() > 2/cov_coef) or args.plot_all:
            if ttest_res is not None:
                trk_border = [non_trk_exons[0] + 0.5, non_trk_exons[-1] + 1.5,
                           trk_exons[0] + 0.5, trk_exons[-1] + 1.5]
            else:
                trk_border = [0,0,0,0]

            plot_coverage(gene, name, trk_border)
    
        if args.table:
            gene_stats = pd.DataFrame(genes_dict[gene])
            gene_stats["coverage_sep"] = gene_stats["coverage"].copy()
            gene_stats["coverage"] = gene_stats["coverage"].apply(sum)

            gene_stats.index = [str(i) for i in range(1,len(gene_stats)+1)]
            gene_stats[['cov_antisense','cov_sense']] = pd.DataFrame(gene_stats["coverage_sep"].tolist(),
                                                                     index= gene_stats.index)
            gene_stats = gene_stats[['is_trk', "coverage", 'cov_antisense','cov_sense']]
            gene_stats["sample"] = name or args.input.split("/")[-1].split(".")[0]
            gene_stats["gene"] = gene
            gene_stats.to_csv(args.table, mode='a', header=not pathlib.Path(args.table).exists())

        # table with coverage stats
        if args.stats:
            sample_name = name or args.input.split("/")[-1].split(".")[0]
            data = pd.DataFrame(dict(sample=[sample_name], gene=[gene], p_val=[ttest_res],
                                    mean_cov=[mean_cov],
                                    cov_trk=[cov_trk_raw],
                                    cov_non_trk=[cov_non_trk_raw],
                                    total_reads=[total_reads]))
            data.to_csv(args.stats, mode='a',index=False, header=not pathlib.Path(args.stats).exists())

    res = pd.Series(res_dict)
    res.name = 'Transposition_status'
    res.index.name = 'Gene'
    res.to_csv(f'{args.output}/trk_transposition_res.csv')
