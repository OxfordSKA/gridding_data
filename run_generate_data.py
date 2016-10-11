#!/usr/bin/python -u
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from subprocess import call, PIPE
import math

cmd = 'sbatch'
cmd_exists = call("type " + cmd, shell=True, stdout=PIPE, stderr=PIPE) == 0
script = 'slurm_generate_data.wilkes'
num_channels_total = 1500
num_channels_per_sim = 1

# Loop over simulations.
num_sims = int(math.ceil(num_channels_total / num_channels_per_sim))
for i in range(num_sims):
    start_channel = i * num_channels_per_sim
    end_channel = (i + 1) * num_channels_per_sim - 1
    if end_channel >= num_channels_total:
        end_channel = num_channels_total - 1
    num_channels = 1 + end_channel - start_channel
    print('Running channels %d-%d' % (start_channel, end_channel))
    if not cmd_exists:
        continue
    call([cmd, script, '%i' % start_channel, '%i' % num_channels])