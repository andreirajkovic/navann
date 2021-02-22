"""
vcf annotator: version 0.0.1
Author: Andrei Rajkovic
Contact: rajkovic.1@osu.edu
"""
import os
import pandas as pd
import subprocess
import sys
import argparse
import requests
from dpcontracts import require, ensure
import numpy as np

# load the csq data
csq_order = ["transcript_ablation","splice_acceptor_variant",
             "splice_donor_variant","stop_gained","frameshift_variant",
             "stop_lost","start_lost","initiator_codon_variant",
             "transcript_amplification", "inframe_insertion", 
             "inframe_deletion","missense_variant",
             "protein_altering_variant","splice_region_variant",
             "incomplete_terminal_codon_variant","stop_retained_variant",
             "synonymous_variant","coding_sequence_variant",
             "mature_miRNA_variant","5_prime_UTR_variant",
             "3_prime_UTR_variant","non_coding_transcript_exon_variant",
             "non_coding_exon_variant","intron_variant",
             "NMD_transcript_variant","non_coding_transcript_variant",
             "nc_transcript_variant","upstream_gene_variant",
             "downstream_gene_variant","TFBS_ablation",
             "TFBS_amplification","TF_binding_site_variant",
             "regulatory_region_ablation","regulatory_region_amplification",
             "feature_elongation","regulatory_region_variant",
             "feature_truncation","intergenic_variant"]

csq_pd = pd.DataFrame({'csq_order': csq_order})
csq_pd.reset_index(inplace=True)
csq_pd.set_index('csq_order', inplace=True)


def runSnpEff(jar_path: str, input_vcf: str,
              output_vcf: str, version: str) -> None:
    """
    Run a local version of snpEff if it exists

    jar_path - path of the snpEff jar file
    input_vcf - location of the vcf, this will be after the file has been annotated
    output_vcf - same as the input but gets processed to have a .snpEff. ext added 
    version - genome assembly
    """

    # check to make sure jar_path is vaild
    if os.path.isfile(jar_path) is False:
        raise RuntimeError(f"{jar_path} could not be located")
        sys.exit()
    output_vcf, ext = os.path.splitext(output_vcf)
    output_vcf = output_vcf+'snpEff.vcf'
    if os.path.isfile(output_vcf):
        raise RuntimeError(f"{output_vcf} is already a file")
        sys.exit()
    with open(output_vcf, 'w') as vcf_out:
        return_code = subprocess.call(['java', '-Xmx4g', '-jar', jar_path,
                                       version, input_vcf], stdout=vcf_out)
    if return_code == 0:
        print("vcf annontated...")
    else:
        print("vcf annontatation failed exiting...")
        sys.exit()


@require("`sample_name` must be in vcf_pd.columns",
         lambda args: args.sample_name in args.vcf_pd.columns)
def parseSampleField(vcf_pd: pd.core.frame.DataFrame,
                     sample_name: str) -> pd.core.frame.DataFrame:
    """
    Extract the data from the sample columns

    vcf_pd - vcf converted into a pandas dataframe
    sample_name - sample name typically found after field #9 in a vcf
    """
    sample_pd = vcf_pd[sample_name].str.split(':', expand=True)
    sample_pd.columns = vcf_pd['FORMAT'][0].split(':')[:len(sample_pd.columns)]
    return sample_pd


@require("`sample_pd` is not empty",
         lambda args: args.sample_pd.empty is False)
@require("`sample_pd` has DPR in sample_pd.columns",
         lambda args: 'DPR' in args.sample_pd.columns)
@require("`sample_pd` DPR column is obj",
         lambda args: isinstance(args.sample_pd.DPR,
                                 np.object))
@ensure("`results` ref column is np.int64",
        lambda args, results: isinstance(results[0].values[0],
                                   np.int64))
def parseDepthData(sample_pd: pd.core.frame.DataFrame):
    """
    Use the data found in the sample field to get the 
    read counts. The DPR column will be used to do this.

    sample_pd - pandas dataframe with ploidy status, counts ect
    """
    # parse depth from the sample field
    read_count_pd = sample_pd['DPR'].str.split(',', expand=True)
    read_count_pd.fillna(0, inplace=True)
    read_count_pd = read_count_pd.astype(int)
    return read_count_pd


def computeTotalDepth(read_count_pd: pd.core.frame.DataFrame,
                      sample: str, vcf_pd, add_format: bool = True) -> None:
    """
    Compute total depth by summing all the REF and ALT
    read count values.

    read_count_pd - dataframe of the read counts
    sample - sample name
    add_format -whether to modify the FORMAT field
    vcf_pd - vcf converted into a pandas dataframe
    """
    # Depth of sequence coverage at the site of variation.
    total_dp = read_count_pd.sum(axis=1)
    total_dp = ':'+total_dp.astype(str)
    if add_format is True:
        vcf_pd['FORMAT'] = vcf_pd['FORMAT'] + ':TD'
    vcf_pd[sample] = vcf_pd[sample] + total_dp


def parseAltSeqDepth(read_count_pd: pd.core.frame.DataFrame,
                     sample: str, vcf_pd, add_format: bool = True) -> None:

    """
    Parse the ALT sequence by simply pulling out the
    values that are not REF counts.

    read_count_pd - dataframe of the read counts
    sample - sample name
    add_format -whether to modify the FORMAT field
    vcf_pd - vcf converted into a pandas dataframe
    """
    # Number of reads supporting the variant.
    alt_counts = read_count_pd.iloc[:, 1:].astype(str).T
    alt_counts = alt_counts.apply(lambda x: ",".join(x))
    # remove zeros if variant doesn't exist; possible edge case 12,0,0,1
    alt_counts = alt_counts.str.replace('[^1-9]*,0$', '')
    if add_format is True:
        vcf_pd['FORMAT'] = vcf_pd['FORMAT'] + ':AS'
    alt_counts = ':' + alt_counts
    vcf_pd[sample] = vcf_pd[sample] + alt_counts


def computeVAF(read_count_pd: pd.core.frame.DataFrame,
               sample: str, vcf_pd, add_format: bool = True) -> None:
    """
    Perform a simple computation for the VAF by
    dividing the ALT and REF reads by the total.
    The output is converted into a string and
    regex is used to parse out potentially
    complicated situations arising due to mutiallelic
    variants.

    read_count_pd - dataframe of the read counts
    sample - sample name
    add_format -whether to modify the FORMAT field
    vcf_pd - vcf converted into a pandas dataframe
    """
    # Percentage of reads supporting the variant
    # versus those supporting reference reads.
    total_reads = read_count_pd.sum(axis=1)
    vaf_pd = read_count_pd.div(total_reads, axis=0).round(4).astype(str)
    vaf_start = vaf_pd.iloc[:, :2].T.apply(lambda x: ",".join(x))
    vaf_end = ',' + vaf_pd.iloc[:, 2:].T.apply(lambda x: ",".join(x))
    vaf_end = vaf_end.str.replace(r'(?:,0\.0){1,}?$', '')
    vaf_pd = vaf_start + vaf_end
    vaf_pd = ':'+vaf_pd
    if add_format is True:
        vcf_pd['FORMAT'] = vcf_pd['FORMAT'] + ':VAF'
    vcf_pd[sample] = vcf_pd[sample] + vaf_pd


def restExAC(vcf_pd, ext: str = "/rest/bulk/variant/variant") -> tuple:
    """
    POST the vcf as a bulk request and retrieve the data.
    We perform this on biallelic and multiallelic data.
    Multiallelic data is handled by iterating through
    each of the alleles and adding that to the bulk
    request.
    ext - extension to use in case we need to use a different
          API
    vcf_pd - vcf converted into a pandas dataframe
    """
    server = "http://exac.hms.harvard.edu/"
    # handle multi-allelic and biallelic separately
    track_idx = []
    key = []
    split_allelic = vcf_pd['ALT'].str.split(',', expand=True)
    multi_allelic = split_allelic[~split_allelic[1].isna()]
    # parse data for the biallelic data
    biallelic_idx = split_allelic[split_allelic[1].isna()].index
    track_idx += list(biallelic_idx.values)
    biallelic_pd = vcf_pd.loc[biallelic_idx, ['#CHROM', 'POS', 'REF', 'ALT']]
    biallelic_pd = biallelic_pd.astype(str).T
    bi_vars = list(biallelic_pd.apply(lambda x: "-".join(x)).values)
    key += bi_vars
    str_vars = '\"' + '\",\"'.join(bi_vars) + '\"'
    for c in multi_allelic.columns:
        alts = multi_allelic[~multi_allelic[c].isna()][c]  # no missing
        track_idx += list(alts.index)  # capture original index
        var_info = vcf_pd.loc[alts.index][['#CHROM', 'POS', 'REF']]
        var_info['ALT'] = alts
        var_info = var_info.astype(str).T
        mvars = list(var_info.apply(lambda x: "-".join(x)).values)
        key += mvars
        str_mvar = ',\"' + '\",\"'.join(mvars) + '\"'
        str_vars = str_vars + str_mvar  # combine the strings
    str_vars = '[' + str_vars + ']'
    # pull out the data
    r = requests.post(server+ext, data=str_vars,
                      headers={"Content-Type": "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    return (decoded, track_idx, key)


def addExACFreq(decoded: dict, track_idx: list, key: list, vcf_pd):
    """
    Parse the ExAC json for the allele frequences found
    in the variant json. If no allele frequency is present,
    report it as zero.

    decoded - return json from the ExAC server
    track_idx - index of the data so as to keep track when biallelic
                and multiallelic sites get mixed up during the annotation
                part
    key - the CHROM-POS-REF-ALT input when posting to the ExAC server
    vcf_pd - vcf converted into a pandas dataframe
    """
    freq = []
    exac_key = []
    if len(decoded) != len(track_idx):
        raise RuntimeError("The output from ExAC does not match the expected output")
    for x in decoded:
        exac_key.append(x)
        try:
            freq.append(decoded[x]['allele_freq'])
        except KeyError:
            freq.append(0)
    # coerce the data back into a dataframe
    lookup_pd = pd.DataFrame({"o_idx": track_idx, 'key': key})
    exac_pd = pd.DataFrame({"freq": freq, 'key': exac_key})
    freq_pd = pd.merge(lookup_pd, exac_pd, on=['key'])
    # annontate the vcf wit the the ExAC freq data
    for g, r in freq_pd.groupby('o_idx'):
        info = ";ExAC_AF="+','.join([f"{x:.2E}" for x in r['freq']])
        vcf_pd.loc[g, 'INFO'] = vcf_pd.loc[g, 'INFO'] + info
    # format the data so that the info field holds the proper data


def order_csqs(vep_annotations):
    """
    compute the order and select the most impactful consequence of the variant
    vep_annotations - the dictionary of transcripts for a variant and their potential impact
    """
    consequences = []
    for ann in vep_annotations:
        mjc = ann['major_consequence']
        consequences.append(mjc)
    try:
        most_severe = csq_pd.loc[consequences].idxmin().values[0]
    except Exception:
        print(f"failed to annotate {vep_annotations}")
        most_severe = None
    return most_severe


def addExACcscq(decoded: dict, track_idx: list, key: list, vcf_pd):
    """
    Parse the ExAC json for the consequences found in the
    VEP annotations. Report only the highest impact. The
    order of impact can be found here:
    https://github.com/konradjk/exac_browser/blob/a212465c5b75752abe8990cf6aa581295835ab58/utils.py#L211

    decoded - return json from the ExAC server
    track_idx - index of the data so as to keep track when biallelic
                and multiallelic sites get mixed up during the annotation
                part
    key - the CHROM-POS-REF-ALT input when posting to the ExAC server
    vcf_pd - vcf converted into a pandas dataframe
    """
    csqs = []
    exac_key = []
    if len(decoded) != len(track_idx):
        raise RuntimeError("The output from ExAC does not match the expected output")
    for x in decoded:
        exac_key.append(x)
        try:
            vep_annotations = decoded[x]['vep_annotations']
            if len(vep_annotations) < 1:
                csqs.append(None)
                print(f"failed to annotate {x} due to an empty vep list")
            else:
                ms = order_csqs(vep_annotations)
                csqs.append(ms)
        except KeyError:
            csqs.append(None)
    # coerce the data back into a dataframe
    lookup_pd = pd.DataFrame({"o_idx": track_idx, 'key': key})
    exac_pd = pd.DataFrame({"csqs": csqs, 'key': exac_key})
    csqs_pd = pd.merge(lookup_pd, exac_pd, on=['key'])
    # annontate the vcf wit the the ExAC csqs data
    for g, r in csqs_pd.groupby('o_idx'):
        info = ";ExAC_CSQS="+','.join([f"{x}" for x in r['csqs']])
        vcf_pd.loc[g, 'INFO'] = vcf_pd.loc[g, 'INFO'] + info


def extract_vcf_header(vcf_path: str, outpath: str):
    """
    Reads the original vcf in and then extracts the header.
    The header is written out to a new file that will serve
    as the input for writing out the file after the annotations
    have been added.

    vcf_path - path to the vcf in question being annotated
    outpath - path to were the user would like to save the
              annotated vcf.
    """
    # parse the header from the vcf
    fout = open(outpath, 'a')
    with open(vcf_path, 'r') as fin:
        for line in fin:
            line = line.strip()
            if line.startswith('##'):  # info and filter lines
                fout.write(line+'\n')
            elif line.startswith('#CHROM'):  # header line
                header = line.split('\t')
                return header, fout


def main():
    """
    Handles the execution of the code
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--snpEff_jar", default=None, required=False,
                        help="path of local snpEff jar")
    parser.add_argument("--vcf", type=str, required=True,
                        help="path of vcf to annotate")
    parser.add_argument("--outpath", type=str,
                        help="path to write out the vcf")
    parser.add_argument('--ver', default='GRCh37.75', required=False,
                        help='genome assembly')
    parser.add_argument("--all", required=False, action='store_true',
                        help='run all the annotations')
    parser.add_argument("--var_eff", required=False, action='store_true',
                        help='parse variant effect')
    parser.add_argument("--total_depth", required=False, action='store_true',
                        help='compute total depth for variants')
    parser.add_argument("--alt_depth", required=False, action='store_true',
                        help='compute alt depth for variants')
    parser.add_argument("--vaf", required=False, action='store_true',
                        help='compute vaf')
    parser.add_argument("--allele_freq", required=False, action='store_true',
                        help='ExaAC allele_frequency')
    parser.add_argument("--var_type", required=False, action='store_true',
                        help='annotate variant type')
    parser.add_argument("--verbose", required=False, action='store_true',
                        help='enable verbose output')
    parser.add_argument("--simple_output", required=False, default=None,
                        help='outpath for simplified report of data')
    args = parser.parse_args()    
    if args.verbose:
        def verboseprint(message):
            print(message)
    else:
        verboseprint = lambda *a: None
    # parse the header
    header_str, header_file = extract_vcf_header(args.vcf, args.outpath)
    verboseprint(f"{','.join(header_str)} extracted from {args.vcf}")
    # load the annotated vcf
    vcf_pd = pd.read_csv(args.vcf, sep='\t', comment='#',
                         names=header_str)
    verboseprint("Loaded vcf...")
    # extract the var effect
    if args.allele_freq or args.var_eff or args.all:
        verboseprint("Running ExAC annontation...")
        decoded, track_idx, key = restExAC(vcf_pd)
        if args.var_eff or args.all:
            addExACcscq(decoded, track_idx, key, vcf_pd)
            verboseprint("ExAC variant consequences added...")
        if args.allele_freq or args.all:
            addExACFreq(decoded, track_idx, key, vcf_pd)
            verboseprint("ExAC variant frequencies added...")
    if args.total_depth or args.all or args.alt_depth or args.vaf:
        verboseprint("Computing depth metrics...")
        add_format = True
        for sample in vcf_pd.columns[9:]:
            sample_pd = parseSampleField(vcf_pd, sample)
            read_count_pd = parseDepthData(sample_pd)
            if args.total_depth or args.all:
                computeTotalDepth(read_count_pd, sample, vcf_pd, add_format)
                verboseprint(f"Total depth added to {sample}...")
            if args.alt_depth or args.all:
                parseAltSeqDepth(read_count_pd, sample, vcf_pd, add_format)
                verboseprint(f"Alt depth added to {sample}...")
            if args.vaf or args.all:
                computeVAF(read_count_pd, sample, vcf_pd, add_format)
                verboseprint(f"VAF added to {sample}...")
            add_format = False
    if args.outpath:
        vcf_pd.to_csv(header_file, sep='\t', index=False)
        verboseprint(f"annotated vcf was saved to {args.outpath}")
    else:
        verboseprint("File was annotated, but not written out as \
              --outpath was missing")
    header_file.close()
    # annotate the vcf with snpEff
    if args.snpEff_jar:
        verboseprint("Running snpEff...")
        runSnpEff(jar_path=args.snpEff_jar, input_vcf=args.outpath,
                  output_vcf=args.outpath, version=args.ver)
    if type(args.simple_output) == str:
        verboseprint(f"Generating simplified report {args.simple_output}")
        simple_df = vcf_pd[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']]
        if args.var_eff or args.all:
            simple_df['ExAC_CSQS'] = vcf_pd.INFO.str.extract('ExAC_CSQS=(.*?);')
            verboseprint(f"added variant effect to {args.simple_output}")
        if args.allele_freq or args.all:
            simple_df['ExAC_AF'] = vcf_pd.INFO.str.extract('(?:ExAC_AF=)(.*)')
            verboseprint(f"added variant frequency to {args.simple_output}")
        if args.total_depth or args.all:
            loc = vcf_pd.FORMAT[0].split(':').index('TD')
            for sample in vcf_pd.columns[9:]:
                simple_df[f'TD_{sample}'] = vcf_pd[sample].str.split(':', expand=True)[loc]
                verboseprint(f"added variant depth to {args.simple_output}")
        if args.alt_depth or args.all:
            loc = vcf_pd.FORMAT[0].split(':').index('AS')
            for sample in vcf_pd.columns[9:]:
                simple_df[f'AS_{sample}'] = vcf_pd[sample].str.split(':', expand=True)[loc]
                verboseprint(f"added alternative depth to {args.simple_output}")
        if args.vaf or args.all:
            loc = vcf_pd.FORMAT[0].split(':').index('VAF')
            for sample in vcf_pd.columns[9:]:
                simple_df[f'VAF{sample}'] = vcf_pd[sample].str.split(':', expand=True)[loc]
                verboseprint(f"added vaf to {args.simple_output}")
        if args.var_type or args.all:
            simple_df['TYPE'] = vcf_pd.INFO.str.extract('TYPE=(.*?);')
            verboseprint(f"added variant type to {args.simple_output}")
        simple_df.to_csv(args.simple_output, index=False)


if __name__ == "__main__":
    main()
