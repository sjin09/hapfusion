import os
import math
import svgwrite
import hapsmash.bamlib
import hapsmash.haplib
import hapsmash.vcflib
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF
from typing import List, Tuple


class SMASH:
    def __init__(self, arr):
        self.center = arr[0] 
        self.ccs2mapq = arr[1]
        self.ccs2state = arr[2]
        self.ccs2coord = arr[3]
        self.coord2ccs = {v: k for k,v in arr[3].items()} 
        self.ccs2hetpos_lst = arr[5]
        self.ccs2denovo_sbs_lst = arr[6]
        self.ccs2denovo_indel_lst = arr[7]
        self.ccs2smash_set = arr[4]
        self.hap2ccs_lst = arr[8]
        self.homsnp_lst = arr[9]
        self.hetsnp_lst = arr[10]
        self.hetindel_lst = arr[11]
        self.tstart, self.tend = arr[12]
        self.tlen = self.tend - self.tstart 
        self.tdist = int(round((self.tlen/10), -3))
        umut = int(math.floor((self.center - self.tstart)/self.tdist))
        dmut = int(math.ceil((self.tend - self.center)/self.tdist))
        self.tick_lst = list(range(self.center - umut * self.tdist, self.center + dmut * self.tdist, self.tdist))


class SVG():
     def __init__(self):
        self.width = 1000
        self.height = 500
        self.center = 500
        self.title_x = 25
        self.title_y = 50
        self.title_font_size = "15"
        self.lx_dist = 450
        self.ly_init = 100
        self.ruler_x1 = self.center - self.lx_dist
        self.ruler_x2 = self.center + self.lx_dist
        self.ruler_xlen = self.ruler_x2 - self.ruler_x1
        self.ruler_y = 450
        self.ruler_dist = 20
        self.circle_radius = 2
        self.text_xdist = 120 
        self.text_size = "6"
        self.smash_fill = "red"
        self.line_colour = "black"



HAP_LST = ["0", "1", "."] 
def get_ccs_count(smash):
    ccs_count = sum([len(smash.hap2ccs_lst[hap]) for hap in HAP_LST])
    return ccs_count 
       

def get_ly_dist(
    svg, 
    smash,
):
    ccs_count = sum([len(smash.hap2ccs_lst[hap]) for hap in HAP_LST])
    svg.ly_dist = math.floor((svg.ruler_y - (svg.ruler_dist +  svg.ly_init))/ccs_count)
       
        
def get_dwg(width, height, svgout):
    dwg = svgwrite.Drawing(svgout, size=(width, height), profile="full")
    return dwg


def add_text(
    dwg,
    x: int, 
    y: int,
    text: str,
    font_size: str
):            
    dwg.add(
        dwg.text
        (
            text,
            insert=(x, y),
            stroke='black',
            font_size=font_size,
            font_family="Helvetica"
        )
    )


def get_initial_ccs_coordinates(
    svg: SVG,
    smash: SMASH,
) -> Tuple[int, int, int, int]:

    counter = 0
    ccs2lxy = {}  
    for hap in HAP_LST:
        for ccs in smash.hap2ccs_lst[hap]:                        
            ccs_start, ccs_end = smash.ccs2coord[ccs]
            lx1 = (svg.ruler_x1 + (svg.ruler_xlen * ((ccs_start - smash.tstart)/smash.tlen)))
            lx2 = (svg.ruler_x1 + (svg.ruler_xlen * ((ccs_end - smash.tstart)/smash.tlen)))    
            ly = svg.ly_init + (counter * svg.ly_dist)
            ccs2lxy[ccs] = (lx1, ly, lx2, ly) 
            counter += 1
    return ccs2lxy


def get_ccs_coordinates(svg, smash):

    counter = 0
    ccs2lxy = get_initial_ccs_coordinates(svg, smash)
    for hap in HAP_LST:
        seen = set()
        counter += 1
        for i in smash.hap2ccs_lst[hap]:
            if i in seen:
                continue
            qlx1, _, qlx2, _ = ccs2lxy[i]
            ly = svg.ly_init + (counter * svg.ly_dist)
            for j in smash.hap2ccs_lst[hap]: 
                if i == j: 
                    continue
                if j in seen:
                    continue
                tlx1, _, tlx2, _ = ccs2lxy[j]
                if (tlx1 < qlx1 and qlx1 < tlx2) and (tlx1 < qlx2 and qlx2 < tlx2):
                    continue
                elif (qlx1 < tlx1 and tlx1 < qlx2) and (qlx1 < tlx2 and tlx2 < qlx2):                     
                    continue 
                elif (tlx1 < qlx2 and qlx2 < tlx2):
                    continue
                elif (tlx1 < qlx1 and qlx1 < tlx2):
                    continue
                else:
                    seen.add(i)
                    seen.add(j)
                    ccs2lxy[j] = (tlx1, ly, tlx2, ly)
                    break
            counter += 1
            ccs2lxy[i] = (qlx1, ly, qlx2, ly)
    return ccs2lxy


def add_line(
    dwg, 
    x1: int, 
    y1: int, 
    x2: int, 
    y2: int, 
    lcolour: str, 
    lwidth: str
):
    dwg.add(
        dwg.line(
            start=(x1, y1), 
            end=(x2, y2), 
            stroke=lcolour, 
            stroke_width=lwidth
            )
        )


def get_hetsnp_coordinate(
    pos,
    svg,
    smash, 
) -> int:

    cx =  svg.ruler_x1 + (svg.ruler_xlen * (pos - smash.tstart)/smash.tlen)
    return cx


def add_circle(
    dwg,
    cx: int,
    cy: int,
    radius: int,
    cfill: str 
):
    dwg.add(
        dwg.circle(center=(cx, cy),
        r=radius,
        stroke=svgwrite.rgb(15, 15, 15, '%'),
        fill=cfill)
    )


def add_ruler(
    dwg,
    svg, 
    smash
):
    add_line(dwg, svg.ruler_x1, svg.ruler_y, svg.ruler_x2, svg.ruler_y, svg.line_colour, 1) # draw ruler
    for tick in smash.tick_lst:
        tick_x = svg.ruler_x1 + (svg.ruler_xlen * (tick - smash.tstart)/smash.tlen)
        if tick < smash.center:
            tick_txt = "-{} bp".format(smash.center - tick)
        elif tick == smash.center:
            tick_txt = 0
        else:
            tick_txt = "+{} bp".format(tick - smash.center)
        add_line(dwg, tick_x, svg.ruler_y, tick_x, svg.ruler_y - 5, svg.line_colour, 1) # draw ruler
        add_text(
            dwg,
            tick_x,
            svg.ruler_y - 10,
            tick_txt,
            svg.text_size
    ) 


def svg2pdf(svgin, pdfout): 
    drawing = svg2rlg(svgin)
    renderPDF.drawToFile(drawing, pdfout)


def dump_hapsmash_pdf(
    sample: str,
    chrom_lst: List[str],
    chrom2recomb_lst,
    pdf_dir: str,
):
  
    if not os.path.exists(pdf_dir):
        os.mkdir(pdf_dir)
        
    for chrom in chrom_lst:
        for recomb in chrom2recomb_lst[chrom]:
            svg = SVG()
            smash = SMASH(recomb)
            get_ly_dist(svg, smash)
            svgout = "{}/{}_{}_{}.svg".format(pdf_dir, chrom, smash.tstart, smash.tend)
            pdfout = "{}/{}_{}_{}.pdf".format(pdf_dir, chrom, smash.tstart, smash.tend)
            dwg = get_dwg(svg.width, svg.height, svgout) 
            add_ruler(dwg, svg, smash)
            add_text(
                dwg,
                svg.title_x,
                svg.title_y,
                "{}:{}-{}\t({})".format(chrom, smash.tstart, smash.tend, sample),
                svg.title_font_size,
            )

            ccs2lxy = get_ccs_coordinates(svg, smash)
            for hap in HAP_LST:
                for ccs in smash.hap2ccs_lst[hap]:                        
                    state = smash.ccs2state[ccs]  
                    lx1, ly1, lx2, ly2 = ccs2lxy[ccs]
                    circle_fill = "black" if hap == "0" else "white" 
                    add_line(dwg, lx1, ly1, lx2, ly2, svg.line_colour, 1) # draw reads
                    if state: 
                        for i, hetpos in enumerate(smash.ccs2hetpos_lst[ccs]):
                            cx = get_hetsnp_coordinate(hetpos, svg, smash)
                            if i in smash.ccs2smash_set[ccs]:
                                add_circle(dwg, cx, ly1, svg.circle_radius, svg.smash_fill)
                            else:
                                add_circle(dwg, cx, ly1, svg.circle_radius, circle_fill)
                        add_text(
                            dwg,
                            lx1 - svg.text_xdist,
                            ly1 + (svg.ly_dist * 1/3),
                            ccs,
                            svg.text_size
                        )
                    else:
                        for hetpos in smash.ccs2hetpos_lst[ccs]:
                            cx = get_hetsnp_coordinate(hetpos, svg, smash)
                            add_circle(dwg, cx, ly1, svg.circle_radius, circle_fill)
            dwg.save()
            svg2pdf(svgout, pdfout) 

