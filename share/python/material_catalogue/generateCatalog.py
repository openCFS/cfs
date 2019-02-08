#!/usr/bin/python
from __future__ import print_function

import numpy as np
import argparse
import itertools
import shutil


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--res', type=int, required=True,
                        help='resolution of catalog')
    parser.add_argument('-x', type=int, default=20,
                        help='resolution of microcell in x-direction')
    parser.add_argument('-y', type=int, default=20,
                        help='resolution of microcell in y-direction')
    parser.add_argument('-z', type=int, default=20,
                        help='resolution of microcell in y-direction')
    parser.add_argument('-c', '--catalogfile', default='detailed_stats',
                        help='catalogfile to write')
    parser.add_argument('-f', '--generatorfun', dest='function_name',
                        help='function to generate micro cell')

    args = parser.parse_args()

    np = 2

    p = range(args.res)
    params = []
    for _ in range(np):
        params.append(p)
    
    with open(args.catalogfile, 'wb') as f:
        for _ in range(np):
            f.write('{:5d} '.format(args.res))
        for _ in range(6):
            f.write('{:e} '.format(0))
        f.write('\n') 

    with open('jobs', 'wb') as f:
        pass

    for combination in itertools.product(*params):
        str = 'tensor' + '_{:d}' * np + '.sh'
        qsubfile = str.format(*combination)
        shutil.copyfile('qsub_template.sh', qsubfile)
        with open(qsubfile, 'ab') as f:
            command = 'python ./calculateTensor.py ' + '{:d} ' * np
            command = command.format(*combination)
            command += '--res {:d} -c {:s} -f {:s}\n'.format(args.res, args.catalogfile, args.function_name)
            f.write(command)
        with open('jobs', 'ab') as f:
            f.write('qsub ' + qsubfile + '\n')
