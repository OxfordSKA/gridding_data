#!/bin/bash

pointings=(0 1)
# snapshots=(0 1)
snapshots=(0 1 2)
freqs=(0 1)

for p in "${pointings[@]}"; do
    for s in "${snapshots[@]}"; do
        for f in "${freqs[@]}"; do
            ini=$(printf "setup_p%02i_s%02i_f%02i.ini" "$p" "$s" "$f")
            oskar_sim_interferometer "$ini"
            # ms=$(printf "test_p%02i_s%02i_f%02i.ms" $p $s $f)
            # casapy --nogui --nologger --log2term -c casa_make_image.py $ms
        done
    done
done
#oskar_sim_interferometer setup_p00_s00_f00.ini
