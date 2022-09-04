# modules
import sys
import warnings
import argparse

def make_wide(formatter, w=120, h=36):
    """Return a wider HelpFormatter, if possible."""
    try:
        # https://stackoverflow.com/a/5464440
        # beware: "Only the name of this class is considered a public API."
        kwargs = {'width': w, 'max_help_position': h}
        formatter(None, **kwargs)
        return lambda prog: formatter(prog, **kwargs)
    except TypeError:
        warnings.warn("argparse help formatter failed, falling back.")
        return formatter

# argparse
def parse_args(program_version, arguments=sys.argv[1:]):
    # main_arguments
    parser = argparse.ArgumentParser(
        add_help=True,
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter),
        description="hapmash calls crossovers and gene conversions"
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=program_version),
    )
    # subcommands: init
    subparsers = parser.add_subparsers(dest="sub", metavar="")

    # subcommands: call
    parser_call = subparsers.add_parser(
        "call",
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter, w=180, h=60),
        help="detect crossovers and gene conversions using PacBio circular consensus seuqence reads"
    )
    parser_call.add_argument(
        "-i",
        "--bam",
        type=str,
        required=True,
        help="minimap2 (parameters: -ax map-hifi --cs=short) aligned SAM/BAM files"
    )
    parser_call.add_argument(
        "--region",
        type=str,
        required=False,
        help="target chromosome",
    )
    parser_call.add_argument(
        "--region_list",
        type=str,
        required=False,
        help="list of target chromosomes separated by new line"
    )
    parser_call.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="phased deepvariant VCF file" 
    )
    parser_call.add_argument(
        "--min_bq",
        type=int,
        default=93,
        required=False,
        help="minimum base quality score threshold"
    )
    parser_call.add_argument(
        "--min_mapq",
        type=int,
        default=60,
        required=False,
        help="minimum mapping quality score"
    )
    parser_call.add_argument(
        "--min_trim",
        type=float,
        default=0.01,
        required=False,
        help="minimum proportion of bases to be trimmed from the start and end of the read"
    )
    parser_call.add_argument(
        "--min_alignment_proportion",
        type=float,
        default=0.80,
        required=False,
        help="minimum proportion of aligned bases"
    )
    parser_call.add_argument(
        "--min_sequence_identity",
        type=float,
        default=0.95,
        required=False,
        help="minimum blast sequence identity threshold"
    )
    parser_call.add_argument(
        "--mismatch_window",
        type=int,
        default=20,
        required=False,
        help="mismatch window size"
    )
    parser_call.add_argument(
        "--max_mismatch_count",
        type=int,
        default=0,
        required=False,
        help="maximum number of mismatches within the mismatch window"
    )
    parser_call.add_argument(
        "--threads",
        type=int,
        default=1,
        required=False,
        help="maximum number of threads to be used"
    )
    parser_call.add_argument(
        "-o",
        "--out",
        type=str,
        required=True,
        help="TSV file to write the gene conversions"
    )
    # subcommands: phase
    parser_phase = subparsers.add_parser(
        "phase",
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter, w=180, h=60),
        help="returns phased hetsnps",
    )
    parser_phase.add_argument(
        "-i",
        "--bam",
        type=str,
        required=True,
        help="minimap2 (parameters: -ax map-hifi --cs=short) aligned SAM/BAM files"
    )
    parser_phase.add_argument(
        "--region",
        type=str,
        required=False,
        help="target chromosome",
    )
    parser_phase.add_argument(
        "--region_list",
        type=str,
        required=False,
        help="list of target chromosomes separated by new line"
    )
    parser_phase.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="deepvariant VCF file with germline mutations",
    )
    parser_phase.add_argument(
        "--min_bq",
        type=int,
        default=20,
        required=False,
        help="minimum base quality score threshold"
    )
    parser_phase.add_argument(
        "--min_mapq",
        type=int,
        default=60,
        required=False,
        help="minimum mapping quality score"
    )
    parser_phase.add_argument(
        "--min_qual",
        type=int,
        default=30,
        required=False,
        help="minimum quality (QUAL) score for hetsnps",
    )
    parser_phase.add_argument(
        "--min_phase_proportion",
        type=float,
        default=0.2,
        required=False,
        help="minimum proportion of phase consistent edges"
    )
    parser_phase.add_argument(
        "--threads",
        type=int,
        default=1,
        required=False,
        help="maximum number of threads to be used"
    )
    parser_phase.add_argument(
        "-o",
        "--out",
        type=str,
        required=True,
        help="VCF file to write phased hetsnps"
    )
    if len(arguments) == 0: 
        parser.print_help()
        parser.exit()
    else:
        return parser, parser.parse_args(arguments)
