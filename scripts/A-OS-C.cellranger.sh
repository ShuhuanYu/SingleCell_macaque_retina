#!/bin/sh

cd /opt/yushuhuan/macaque_retina/result/0.cellranger_output
cellranger-arc count --id=A-OS-C \
                       --reference=/opt/yushuhuan/macaque_retina/result/0.cellranger_output/ARC/macaque_genome \
                       --libraries=/opt/yushuhuan/macaque_retina/result/0.cellranger_output/scripts/A-OS-C.libraries.csv \
                       --localcores=16 \
                       --localmem=64
