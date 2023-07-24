import os
import math
import time
import pysam
import bisect
import natsort
import svgwrite
import himut.bamlib
import hapfusion.bamlib
import hapfusion.caller
import hapfusion.haplib
import hapfusion.vcflib
import multiprocessing as mp
from svglib.svglib import svg2rlg
from collections import defaultdict
from typing import Dict, List, Tuple
from reportlab.graphics import renderPDF
hap_inconsistent_fill = {"0": "black", "1": "white", "2": "blue", "-": "grey"}
hap0_recombination_fill = {"0": "black", "1": "red", "2": "blue", "-": "grey"}
hap1_recombination_fill = {"0": "red", "1": "white", "2": "blue", "-": "grey"}


class REGION:

    def __init__(self, arr):
        self.tname = arr[0] # chrom
        self.start = arr[1] # start # x-coordinate: 50
        self.end = arr[2] # end # x-coordinate: 950
        self.len = arr[3] # len # width margiin: 900
        self.event = arr[4] # CO, NCO, CNCO
        self.pos = arr[5] # PDF is centered around recombination
        self.center = self.len/2 # PDF is centered around recombination
        self.ccs2hbit = arr[6]  
        self.ccs2coord = arr[7] # k: qname, v: [start, end]
        self.ccs2state = arr[8] 
        self.hap2ccs_lst = arr[9]
        self.ccs2hpos_lst = arr[10] # hpos := hetsnp position
        self.prefix = "{}-{}-{}_{}".format(self.tname, self.start, self.end, self.event)


class TITLE():
    def __init__(self):
        self.x = 25 # title x-coordinate
        self.y = 50 # title y-coordinate
        self.stroke = "black"
        self.font_size = "15" # title font size
        self.font_family="Helvetica" # title font family

    def get_title(self, sample, region):

        if region.event == "CO":
            state = "crossover"
        elif region.event == "NCO":
            state = "gene conversion"
        elif region.event == "CNCO":
            state = "complex gene conversion"
        elif region.event == "NCO_candidate":
            state = "gene conversion candidate"
        elif region.event == "CO_NCO_candidate":
            state = "crossover and gene conversion candidate"
        else: 
            state = "."

        self.title = "sample: {},\t coordinate: {}:{}-{},\t event: {}".format( # title
            sample, 
            region.tname, 
            region.start, 
            region.end,
            state            
        ) 

    def write_title(self, dwg, sample, bundle, canvas):
        self.get_title(sample, bundle)
        canvas.write_text(dwg, self.x, self.y, self.title, self.font_size, self.font_family) 

        
class RULER():
    def __init__(self, canvas):
        self.dist = 20
        self.font_size = "12"         
        self.font_family="Helvetica"
        self.x1 = canvas.width_x1 # default: 50 # ruler x-coordinate
        self.x2 = canvas.width_x2 # default: 950
        self.y1 = self.y2 = canvas.height - canvas.margin # ruler y-coordinate
        self.stroke = "black" # ruler colour
        self.stroke_width = 1 # ruler line stroke

    def draw_ruler(self, dwg, region, canvas):
        canvas.draw_line(dwg, self.x1, self.y1, self.x2, self.y2)
        self.draw_ruler_ticks(dwg, region, canvas)

    def get_ruler_ticks(self, region, canvas):
        
        self.ruler_tick_xcoord_lst = []
        breakpoints = [100, 1000, 10000] 
        d = bisect.bisect(breakpoints, region.len) # d: 0, 1, 2, 3 # number of digits
        region_tick_step = int(round(region.len/10, -d)) # round to nearest 10, 100, 1000, 10000, 100000, 1000000, 10000000  
        udist_to_center = region.pos - region.start
        ddist_to_center = region.end - region.pos
        utick_count = int(math.floor(udist_to_center/region_tick_step)) 
        dtick_count = int(math.ceil(ddist_to_center/region_tick_step))
        region_tick_start = region.pos - (utick_count * region_tick_step)
        region_tick_end = region.pos + (dtick_count * region_tick_step)
        for tpos in range(region_tick_start, region_tick_end, region_tick_step): # furthest from center # to center # to furtherst from center
            tick_dist_to_center = tpos - region.pos 
            if tick_dist_to_center < 0:
                x1 = canvas.get_xcoord(tpos, region) 
                tick_text = "{} bp".format(tick_dist_to_center) 
                self.ruler_tick_xcoord_lst.append((x1, tick_text))
            elif tick_dist_to_center == 0:
                x1 = canvas.get_xcoord(tpos, region) 
                self.ruler_tick_xcoord_lst.append((x1, "0"))
            elif tick_dist_to_center > 0:
                x1 = canvas.get_xcoord(tpos, region) 
                tick_text = "+{} bp".format(tick_dist_to_center)
                self.ruler_tick_xcoord_lst.append((x1, tick_text))
        
    def draw_ruler_ticks(self, dwg, region, canvas): 

        y1 = self.y1
        y2 = self.y1 - 5
        y3 = self.y1 - 10
        self.get_ruler_ticks(region, canvas)
        for (tick_xcoord, tick_text) in self.ruler_tick_xcoord_lst:
            canvas.draw_line(dwg, tick_xcoord, y1, tick_xcoord, y2)
            canvas.write_text(dwg, tick_xcoord, y3, tick_text, self.font_size, self.font_family)


class CANVAS():
    def __init__(self):
        self.width = 1200 # width
        self.height = 700 # height
        self.center = 500 # center
        self.margin = 50
        self.width_x1 = 50
        self.width_x2 = self.width - self.width_x1
        self.width_xlen = self.width_x2 - self.width_x1
        self.height_y1 = 2 * self.margin # height start
        self.height_y2 = self.height - self.margin # height stop
        self.height_ylen = self.height_y2 - self.height_y1 # height length
        self.stroke = "black" # line stroke
        self.stroke_width = 1 # line stroke width
        self.text_xdist = 120 
        self.text_size = "6"
        self.radius = 2.5  # circle radius

    def add_cross(
        dwg,
        svg,
        x,
        y,
    ):
        dwg.add(
            dwg.line(
                start=(x - 3, y - 3), 
                end=(x + 3, y + 3), 
                stroke=svg.lcolour, 
                stroke_width=svg.lstroke_width
                )
            )
        dwg.add(
            dwg.line(
                start=(x + 3, y - 3), 
                end=(x - 3, y + 3), 
                stroke=svg.lcolour, 
                stroke_width=svg.lstroke_width
                )
            )

    def draw_line(
        self, 
        dwg,
        x1,  
        y1,
        x2,
        y2,
    ):
        dwg.add(
            dwg.line(
                start=(x1, y1), 
                end=(x2, y2), 
                stroke=self.stroke, 
                stroke_width=self.stroke_width
                )
            )

    def draw_circle(
        self,
        dwg,
        x1: int,
        y1: int,
        fill: str,
    ):
        dwg.add(
            dwg.circle(
                fill=fill,
                r=self.radius,
                center=(x1,  y1),
                stroke=self.stroke,
            )
        )

    def get_ccs_yd(self, region):
        self.get_ccs_count(region)
        self.ccs_yd = math.floor(self.height_ylen/self.ccs_count) 
        
    def get_ccs_count(self, region):
        self.ccs_count = 3
        for ccs_lst in region.hap2ccs_lst.values():
            self.ccs_count += len(ccs_lst) 
        return self.ccs_count

    def get_xcoord(self, pos, region):
        x1 = self.width_x1 + ((pos  - region.start)/region.len) * self.width_xlen
        return x1

    def write_text(
        self, 
        dwg, 
        x1, 
        y1, 
        text,
        font_size,
        font_family
    ):
        dwg.add(
            dwg.text
            (
                text,
                insert=(x1, y1),
                stroke=self.stroke,
                font_size=font_size,
                font_family=font_family
            )
        )


def load_recombination(recomb_file: str):

    chrom2recomb_lst = defaultdict(list)
    chrom2recomb_candidate_lst = defaultdict(list)
    for line in open(recomb_file).readlines():
        if line.startswith("#"):
            continue
        fusion = hapfusion.util.FUSION(line)
        if fusion.is_pass:
            chrom2recomb_lst[fusion.tname].append((fusion.qname, fusion.tname, int(fusion.tstart), int(fusion.tend), fusion.phase_set, fusion.event))
        else:
            chrom2recomb_candidate_lst[fusion.tname].append((fusion.qname, fusion.tname, int(fusion.tstart), int(fusion.tend), fusion.phase_set, fusion.event))
    return chrom2recomb_lst, chrom2recomb_candidate_lst


def load_recombination_region_ccs(
    bam_file: str,
    recomb_lst: List[Tuple[str, str, int, int]],
    ps2hbit_lst: Dict[str, List[str]],
    ps2hpos_lst: Dict[str, List[int]],
    ps2hetsnp_lst: Dict[str, List[Tuple[int, str, str]]],
    chrom2recomb_region_ccs_lst
):

    recomb_region_ccs_lst = [] 
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for (qname, tname, tstart, tend, phase_set, event) in recomb_lst:
        ccs_lst = []
        tend_lst = []
        tstart_lst = [] 
        hetsnp_lst = []
        ccs2hbit = {}
        ccs2coord = {} 
        ccs2state = {}
        ccs2hpos_lst = {}
        hap2ccs_lst = defaultdict(list)
        ps_hbit_lst = ps2hbit_lst[phase_set] 
        ps_hpos_lst = ps2hpos_lst[phase_set] 
        ps_hetsnp_lst = ps2hetsnp_lst[phase_set] 
        for i in alignments.fetch(tname, tstart, tend): ## collect candidates
            ccs = hapfusion.bamlib.BAM(i)
            if ccs.mapq == 0: 
                continue
            if not ccs.is_primary: 
                continue
            ccs.get_cs2tpos2qbase() # o: ccs.rpos2qpos, ccs.tpos2qbase, ccs.mismatch_lst 
            hapfusion.haplib.get_ccs_hap(ccs, ps_hbit_lst, ps_hpos_lst, ps_hetsnp_lst) # o: ccs.hap 
            if ccs.hap == ".":
                continue
            ccs_lst.append(ccs.qname)
            ccs2hbit[ccs.qname] = ccs.hbit
            hetsnp_lst.extend(ccs.hetsnp_lst)
            ccs2coord[ccs.qname] = (ccs.tstart, ccs.tend)
            ccs2hpos_lst[ccs.qname] = [hetsnp[0] for hetsnp in ccs.hetsnp_lst] 
            if hapfusion.caller.is_ccs_phased(ccs.hap):
                ccs2state[ccs.qname] = "hap_consistent"
                hap2ccs_lst[ccs.hap].append(ccs.qname)
            else:
                if ccs.hap == "2":
                    ccs2state[ccs.qname] = "hap_ambiguous"
                    hap2ccs_lst[ccs.hap].append(ccs.qname)
                else:
                    if qname != ccs.qname:
                        ccs2state[ccs.qname] = "hap_inconsistent"
                    else:
                        ccs2state[ccs.qname] = "recombination"
                        q, r = divmod(len(ccs.hetsnp_lst), 2)
                        if r == 0:
                            dpos = ccs.hetsnp_lst[q][0]
                            upos = ccs.hetsnp_lst[q - 1][0]
                            grid_center = int((upos + dpos)/2)
                        else:
                            grid_center = ccs.hetsnp_lst[q][0]
                    hap2ccs_lst[ccs.hap[0]].append(ccs.qname)
            tend_lst.append(ccs.tend)
            tstart_lst.append(ccs.tstart)
            
        hetsnp_lst = natsort.natsorted(list(set(hetsnp_lst)))
        grid_end = max(tend_lst)
        grid_start = min(tstart_lst)
        grid_len = grid_end - grid_start
        recomb_region_ccs_lst.append(
            [
                tname,
                grid_start,
                grid_end,
                grid_len,
                event,
                grid_center,
                ccs2hbit,
                ccs2coord,
                ccs2state,
                hap2ccs_lst,
                ccs2hpos_lst,
            ] 
        )
    chrom2recomb_region_ccs_lst[tname] = recomb_region_ccs_lst


def parallel_load_recombination_coordinates(
    bam_file: str,
    vcf_file: str,
    recomb_file: str,
    chrom_lst: List[str],
    tname2tsize: Dict[str, int],
    threads: int
):

    (
        chrom2ps2hbit_lst,
        chrom2ps2hpos_lst,
        chrom2ps2hetsnp_lst,
        chrom2chunkloci_lst,
    ) = himut.vcflib.load_phased_hetsnps(vcf_file, chrom_lst, tname2tsize) # load
    chrom2recomb_lst, chrom2recomb_candidate_lst = load_recombination(recomb_file)
    del chrom2chunkloci_lst

    # parallel
    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2recomb_region_ccs_lst = manager.dict()
    chrom2recomb_candidate_region_ccs_lst = manager.dict()
    load_recombination_region_ccs_arg_lst = [
        (
            bam_file,
            chrom2recomb_lst[chrom],
            chrom2ps2hbit_lst[chrom],
            chrom2ps2hpos_lst[chrom],
            chrom2ps2hetsnp_lst[chrom], 
            chrom2recomb_region_ccs_lst
        )
        for chrom in chrom_lst
    ]
    # load_recombination_candidate_region_ccs_arg_lst = [ ## EVALUATION
    #     (
    #         bam_file,
    #         chrom2recomb_lst[chrom],
    #         chrom2ps2hbit_lst[chrom],
    #         chrom2ps2hpos_lst[chrom],
    #         chrom2ps2hetsnp_lst[chrom], 
    #         chrom2recomb_candidate_region_ccs_lst
    #     )
    #     for chrom in chrom_lst
    # ]
    p.starmap(
        load_recombination_region_ccs, load_recombination_region_ccs_arg_lst,
    )
    # p.starmap(
    #     load_recombination_region_ccs, load_recombination_candidate_region_ccs_arg_lst,
    # )
    p.close()
    p.join()
    return dict(chrom2recomb_region_ccs_lst), dict(chrom2recomb_candidate_region_ccs_lst)


def get_dwg(width, height, svg_out):
    dwg = svgwrite.Drawing(svg_out, size=(width, height), profile="full")
    return dwg

def svg2pdf(svgin, pdfout): 
    drawing = svg2rlg(svgin)
    renderPDF.drawToFile(drawing, pdfout)


def get_ccs_hpos_xcoord(
    hap: str, 
    state: str, 
    ccs_hbit: str, 
    ccs_hpos_lst: List[int],
    canvas: CANVAS, 
    region: REGION, 
):

    if hap == "0":
        if state == "recombination":
            ccs_hpos_x1_lst = [(canvas.get_xcoord(hpos, region), hap0_recombination_fill[hbit]) for (hbit, hpos) in zip(ccs_hbit, ccs_hpos_lst)]
        elif state == "hap_consistent": 
            ccs_hpos_x1_lst = [(canvas.get_xcoord(hpos, region), "black") for hpos in ccs_hpos_lst]
        elif state == "hap_inconsistent":
            ccs_hpos_x1_lst = [(canvas.get_xcoord(hpos, region), hap_inconsistent_fill[hbit]) for (hbit, hpos) in zip(ccs_hbit, ccs_hpos_lst)]
    else: # hap == 1
        if state == "recombination":
            ccs_hpos_x1_lst = [(canvas.get_xcoord(hpos, region), hap1_recombination_fill[hbit]) for (hbit, hpos) in zip(ccs_hbit, ccs_hpos_lst)]
        elif state == "hap_consistent": 
            ccs_hpos_x1_lst = [(canvas.get_xcoord(hpos, region), "white") for hpos in ccs_hpos_lst]
        elif state == "hap_inconsistent":
            ccs_hpos_x1_lst = [(canvas.get_xcoord(hpos, region), hap_inconsistent_fill[hbit]) for (hbit, hpos) in zip(ccs_hbit, ccs_hpos_lst)]
    return ccs_hpos_x1_lst


def dump_recombination_plot(
    sample: str,
    pdf_dir: str, 
    recomb_region_ccs_lst
):

    for i in recomb_region_ccs_lst: ## iterate through each recombination event
        title = TITLE()
        canvas = CANVAS()
        region = REGION(i)
        ruler = RULER(canvas)
        canvas.get_ccs_yd(region)
        svg_out = "{}/{}.svg".format(pdf_dir, region.prefix)
        pdf_out = "{}/{}.pdf".format(pdf_dir, region.prefix)
        dwg = get_dwg(canvas.width, canvas.height, svg_out) # get a blankpage 
        ruler.draw_ruler(dwg, region, canvas) 
        title.write_title(dwg, sample, region, canvas)

        counter = 0 
        for i, hap in enumerate(hapfusion.haplib.hap_lst): # plotting CCS reads
            counter += i
            if hap == "2":
                continue 
            for ccs in region.hap2ccs_lst[hap]:  
                state = region.ccs2state[ccs]  # hap_consistent, hap_inconsistent, hap_ambiguous, recombination
                ccs_start, ccs_end  = region.ccs2coord[ccs] 
                ccs_x1 = canvas.get_xcoord(ccs_start, region) 
                ccs_x2 = canvas.get_xcoord(ccs_end, region) 
                ccs_y1 = ccs_y2 = hpos_y1 = (counter * canvas.ccs_yd) + canvas.height_y1
                canvas.draw_line(dwg, ccs_x1, ccs_y1, ccs_x2, ccs_y2)
                ccs_hpos_x1_lst = get_ccs_hpos_xcoord(
                    hap, 
                    state, 
                    region.ccs2hbit[ccs],
                    region.ccs2hpos_lst[ccs], 
                    canvas, 
                    region
                ) 
                for (hpos_x1, hpos_fill) in ccs_hpos_x1_lst:
                    canvas.draw_circle(dwg, hpos_x1, hpos_y1, hpos_fill)
        dwg.save()
        svg2pdf(svg_out, pdf_out)     



def parallel_dump_recombination_plot(
    bam_file: str,
    bed_file: str,
    vcf_file: str, 
    recomb_file: str, 
    region: str, 
    region_list: str,
    threads: int,
    debug: bool, ## TODO
    pdf_dir: str
):

    # init
    cwd = os.getcwd()
    cpu_start = time.time() / 60
    os.makedirs("{}/{}".format(cwd, pdf_dir), exist_ok=True)
        
    # load
    sample = hapfusion.bamlib.get_sample(bam_file)
    tname2tsize = himut.bamlib.get_tname2tsize(bam_file)[1]
    chrom_lst = himut.util.load_loci(region, region_list, tname2tsize)[0]
    hapfusion.util.check_plot_input_exists(
        bam_file,
        vcf_file,
        recomb_file,
        chrom_lst,
        tname2tsize,
    )
    chrom2recomb_coordinate_lst, chrom2recomb_candidate_coordinate_lst = parallel_load_recombination_coordinates(
        bam_file,
        vcf_file,
        recomb_file,
        chrom_lst,
        tname2tsize,
        threads,
    )
    p = mp.Pool(threads)
    dump_recombination_plot_arg_lst = [
        (
            sample,
            pdf_dir,
            chrom2recomb_coordinate_lst[chrom]
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        dump_recombination_plot, dump_recombination_plot_arg_lst,
    )

    # if debug:
    #     dump_hapfusion_candidate_plot_arg_lst = [
    #         (
    #             bam_file,
    #             chrom2recomb_candidate_lst[chrom],
    #             chrom2ps2hbit_lst[chrom],
    #             chrom2ps2hpos_lst[chrom],
    #             chrom2ps2hetsnp_lst[chrom], 
    #             chrom2chunkloci_lst[chrom],
    #             outdir
    #         )
    #         for chrom in chrom_lst
    #     ]
    #     p.starmap(
    #         dump_recombination_plot, dump_recombination_candidate_plot_arg_lst,
    #     )
    p.close()
    p.join()
    # end

    print("hapfusion finished returning crossover and gene conversions plots")
    cpu_end = time.time() / 60
    duration = cpu_end - cpu_start
    print(
        "hapfusion plot generation took {} minutes".format(
            duration
        )
    )
    hapfusion.util.exit()
