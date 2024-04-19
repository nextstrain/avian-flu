"""
This script runs the entire suite of avian flu builds via AWS batch
One batch run is created per subtype x segment combination
"""

import subprocess
import argparse
import os

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run flu builds on aws')
    parser.add_argument('--subtypes', nargs='+', type=str, help ="avian flu subtypes to include")
    parser.add_argument('--segments', nargs='+', type=str, help ="avian flu segments to include")
    params = parser.parse_args()

    subtypes = ["h5n1", "h7n9"]
    segments = ["pb2", "pb1", "pa", "ha", "np", "na", "mp", "ns"]

    if not os.path.exists("logs"):
        os.makedirs("logs")

    if params.subtypes is None:
        params.subtypes = subtypes

    if params.segments is None:
        params.segments = segments

    for subtype in params.subtypes:
        for segment in params.segments:
            call = ['nextstrain', 'build', '--aws-batch', '.', '-j 1']
            targets = []
            targets.append('auspice/avian-flu_%s_%s_tree.json'%(subtype, segment))
            targets.append('auspice/avian-flu_%s_%s_meta.json'%(subtype, segment))
            call.extend(targets)
            print(' '.join(call))
            log = open('logs/%s_%s.txt'%(subtype, segment), 'w')
            # run in background with subprocess.Popen
            pro = subprocess.Popen(call, stdout=log, stderr=log)
