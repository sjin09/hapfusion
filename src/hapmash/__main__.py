#!/usr/bin/env python3
__version__ = "0.0.1"
__author__ = "Sangjin Lee"

## modules
import hapmash.caller
import hapmash.phaselib
from hapmash.parse_args import parse_args

def main():
    parser, options = parse_args(program_version=__version__)
    if options.sub == "call": # return gene conversions
        hapmash.util.check_num_threads(options.threads)
        hapmash.caller.call_recombinantion(
            options.bam, # input # bam_file
            options.vcf, # deepvariant phased germline mutations
            options.region, # region 
            options.region_list, #
            options.min_mapq, # int: 0-60
            options.min_sequence_identity, # float: 0.0 - 1.0
            options.min_alignment_proportion, # float: 0.0 - 1.0
            options.min_bq, # minimum base quality score: int
            options.min_trim, # float: 0.0 1.0
            options.mismatch_window, # int: 20
            options.max_mismatch_count, # int: 0 
            options.threads, # maxminum number of threads
            __version__, # str
            options.out, # output # hapmash vcf file
        )
    elif options.sub == "phase": # returns phased hetsnps
        hapmash.phaselib.get_chrom_hblock(
            options.bam, 
            options.vcf, 
            options.region, 
            options.region_list,
            options.min_qual, 
            options.min_bq, 
            options.min_mapq,
            options.min_phase_proportion, 
            options.threads, 
            __version__, 
            options.out, 
        )
    else:
        print("The subcommand does not exist!\n")
        parser.print_help()
        parser.exit()


if __name__ == "__main__":
    main(prog_name="hapmash")  # pragma: no cover
    hapmash.util.exit()
