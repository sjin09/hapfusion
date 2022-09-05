

# class PDF:
#     def __init__(self):
       
       ## >>> from PIL import Image
## >>> from PIL import Image, ImageDraw
## >>> im = Image.new('RGBA', (100, 100), (0, 0, 0, 0))
## >>> draw = ImageDraw.Draw(im)
## >>> draw.rectangle((10, 10, 90, 90), fill="yellow", outline="red")
## >>> draw.save("out.png")

## 1. crossover and gene conversion image generation:
## 2. select all CCS reads that intersect with CCS read with the crossover or gene conversion
## 3. determine the haplotype of CCS read without the crossover or gene conversion
## 4. call SNP and indel from all the CCS reads in the target region
## 5. red circle should indicate reference allele
## 6. blue circle should indicate alternative allele
## 7. grey rectangle indicates CCS reads without the crossover or gene conversion
## 8. white rectangle should indicate CCS read with crossover or gene conversion
## 9. rectangle should indicate undeleted section of CCS read
## 10. there should be sufficient space between each CCS read.
## 11. alpha * MAPQ for each CCS read to control rectangle colour transparency
