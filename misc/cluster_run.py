#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 15:36:09 2022

@author: insauer
"""
import argparse
import os
import sys
#import imp
import glob
import numpy
from shutil import copyfile

"""The scripts starts a run for each combination of climate forcing and GHM
   and permits damage generation for each combination by calling schedule_sim.py
"""

parser = argparse.ArgumentParser(
    description='schedules climada runs for different parameter combinations')
# parser.add_argument(
#     '--parameters', type=str, default="parameters.py",
#     help='parameters file')
parser.add_argument(
    '--dry', action="store_true",
    help='dry run (do not run Climada)')
parser.add_argument(
    '--shared', action="store_true",
    help='share nodes on cluster')
parser.add_argument(
    '--notify', action="store_true",
    help='notify per mail when done')
parser.add_argument(
    '--minutes', type=int, default=0,
    help='maximal minutes to run on cluster (< 60)')
parser.add_argument(
    '--hours', type=int, default=24,
    help='maximal hours to run on cluster (168=week, 720=month)')
parser.add_argument(
    '--threads', type=int, default=8,
    help='maximal number of threads on cluster (<= 16)')
parser.add_argument(
    '--mem_per_cpu', type=int, default=8000,
    help='number of memory per CPU (3584 is MaxMemPerCPU on cluster)')
parser.add_argument(
    '--largemem', action="store_true",
    help='use ram_gpu partition')
parser.add_argument(
    '--verbose', action="store_true",
    help='be verbose')

args = parser.parse_args()

indices = []
sys.dont_write_bytecode = True

# shock_files

# run_names = ['shocks_0',
#  'shocks_1',
#  'shocks_2',
#  'shocks_3',
#  'shocks_4',
#   'shocks_5',
#   'shocks_6',
#   'shocks_7',
#   'shocks_8',
#   'shocks_9',
#   'shocks_10',
#   'shocks_11',
#   'shocks_12',
#   'shocks_13',
#   'shocks_14',
#   'shocks_15',
#   'shocks_16',
#   'shocks_17',
#   'shocks_18',
#   'shocks_19',
#   'shocks_20',
#   'shocks_21',
#   'shocks_22',
#   'shocks_23',
#   'shocks_24',
#   'shocks_25',
#   'shocks_26',
#   'shocks_27',
#   'shocks_28',
#   'shocks_29',
#   'shocks_30',
#   'shocks_31',
#   'shocks_32',
#   'shocks_33',
#   'shocks_34',
#   'shocks_35',
#   'shocks_36',
#   'shocks_37',
#   'shocks_38',
#   'shocks_39',
#   'shocks_40',
#   'shocks_41',
#   'shocks_42',
#   'shocks_43',
#   'shocks_44',
#   'shocks_45',
#   'shocks_46',
#   'shocks_47',
#   'shocks_48',
#   'shocks_49']

# run_names=['shocks_syn_0',
#   'shocks_syn_1',
#   'shocks_syn_2',
#   'shocks_syn_3',
#   'shocks_syn_4',
#   'shocks_syn_5',
#   'shocks_syn_6',
#   'shocks_syn_7',
#   'shocks_syn_8',
#   'shocks_syn_9',
#   'shocks_syn_10',
#   'shocks_syn_11',
#   'shocks_syn_12',
#   'shocks_syn_13',
#   'shocks_syn_14',
#   'shocks_syn_15',
#   'shocks_syn_16',
#   'shocks_syn_17',
#   'shocks_syn_18',
#   'shocks_syn_19',
#   'shocks_syn_20',
#   'shocks_syn_21',
#   'shocks_syn_22',
#   'shocks_syn_23',
#   'shocks_syn_24',
#   'shocks_syn_25',
#   'shocks_syn_26',
#   'shocks_syn_27',
#   'shocks_syn_28',
#   'shocks_syn_29',
#   'shocks_syn_30',
#   'shocks_syn_31',
#   'shocks_syn_32',
#   'shocks_syn_33',
#   'shocks_syn_34',
#   'shocks_syn_35',
#   'shocks_syn_36',
#   'shocks_syn_37',
#   'shocks_syn_38',
#   'shocks_syn_39',
#   'shocks_syn_40',
#   'shocks_syn_41',
#   'shocks_syn_42',
#   'shocks_syn_43',
#   'shocks_syn_44',
#   'shocks_syn_45',
#   'shocks_syn_46',
#   'shocks_syn_47',
#   'shocks_syn_48',
#   'shocks_syn_49']

seeds = [ 0,  1,  2,  3,  4,
        5,  6,  7,  8,  9, 10]
       #  , 11, 12, 13, 14, 15, 16,
       # 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
       # 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49]

run_names= ['shocks_7',
    'shocks_19',
    'shocks_20',
    'shocks_26',
    'shocks_28',
    'shocks_32',
    'shocks_34',
    'shocks_38',
    'shocks_47',
    'shocks_48',
    'shocks_53',
    'shocks_59',
    'shocks_64',
    'shocks_66',
    'shocks_76',
    'shocks_80',
    'shocks_83',
    'shocks_86',
    'shocks_91',
    'shocks_92',
    'shocks_99',
    'shocks_101',
    'shocks_110',
    'shocks_117',
    'shocks_118',
    'shocks_123',
    'shocks_127',
    'shocks_130',
    'shocks_132',
    'shocks_133',
    'shocks_142',
    'shocks_144',
    'shocks_162',
    'shocks_166',
    'shocks_169',
    'shocks_180',
    'shocks_182',
    'shocks_201']
  

run_times=[
    162,
    162
    ]


def schedule_run(run_nb,flag, run_name, run_time, seed):
    if not flag:
        run_label = "run_%s_%s" %(run_name, str(seed))
        if os.path.exists(run_label):
        #    run_id += 1
            return
        os.mkdir(run_label)
 #       desc = run_description()
 #       f = open("%s/parameters.txt" % run_label, 'w')
 #       f.write(desc)
 #       f.close()
 #       run_index.write(run_description_csv(run_label))
 #       run_index.write("\n")
        # with open("%s/settings.yml" % run_label, 'w') as f:
        #     f.write(pyaml.dump(settings_yml))
        #     for nc in glob.glob('*.nc'):
        #         copyfile(nc,"%s/%s" % (run_label,nc))
        #run_id += 1
    else:
        run_label = "."
    if args.dry:
        return
    else:
        if (int(args.hours) <= 24):
            _class = "short"
        elif (int(args.hours) <= 24 * 7):
            _class = "medium"
        else:
            _class = "long"

        run_params = {
            "job_name": "%s/%s" % (os.path.basename(os.getcwd()), run_label),
            "minutes": args.minutes,
            "hours": args.hours,
            "class": _class,
            "initialdir": run_label,
            "node_usage": "share" if args.shared else "exclusive",
            "notification": "END,FAIL,TIME_LIMIT" if args.notify else "FAIL,TIME_LIMIT",
            "comment": "%s/%s" % (os.getcwd(), run_label),
            "environment": "ALL",
            "executable": 'cluster_sim.py',
            "options": "--run_name %s --run_time %i --seed %i"%(run_name, run_time, seed),
            "num_threads": args.threads,
            "mem_per_cpu": args.mem_per_cpu if not args.largemem else 15360,   # if mem_per_cpu is larger than MaxMemPerCPU then num_threads is reduced
            "other": "" if args.largemem else ""
            #"other": "#SBATCH --partition=ram_gpu" if args.largemem else ""
        }

        cmd = """echo "#!/bin/sh
#SBATCH --job-name=\\\"%(job_name)s\\\"
#SBATCH --comment=\\\"%(comment)s\\\"
#SBATCH --time=%(hours)02d:%(minutes)02d:00
#SBATCH --qos=%(class)s
#SBATCH --output=output.txt
#SBATCH --error=errors.txt
#SBATCH --export=%(environment)s
#SBATCH --mail-type=%(notification)s
#SBATCH --%(node_usage)s
#SBATCH --account=ebm        
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=%(num_threads)1d
#SBATCH --mem-per-cpu=%(mem_per_cpu)1d
#SBATCH --workdir=%(initialdir)s
%(other)s
export OMP_PROC_BIND=true
export OMP_NUM_THREADS=%(num_threads)1d
source activate inga_base        
ulimit -c unlimited
%(executable)s %(options)s
" | sbatch -Q""" % run_params

        if args.verbose:
            print(cmd)

        os.system(cmd)
        #run_cnt += 1

num = 1
num *= len(run_names)*len(seeds)

single = True if num == 1 else False
if num > 1:
    print("Number of runs to be scheduled: %s" % num)
    sys.stdout.write('Run? y/N : ')
    if sys.version_info >= (3, 0):
        if input() != "y":
            exit("Aborted")
    else:
        if raw_input() != "y":
            exit("Aborted")

enum = 0

for r, run_name in enumerate(run_names):
    for s, seed in enumerate(seeds):
        schedule_run(run_nb=enum,flag=single, run_name=run_name, run_time=run_times[1], seed = seed)
        enum += 1
if num > 1:
    print("Scheduled %s runs" % enum)

# def set_in_yml(paths, value):
#     global yml_nodes
#     for p in paths:
#         node = yml_nodes
#         nodes = p.split(".")
#         for n in nodes[:-1]:
#             try:
#                 n = int(n)
#             except ValueError:
#                 if not n in node:
#                     exit("Path '%s' not found!" % p)
#             node = node[n]
#         n = nodes[-1]
#         if not n in node:
#             exit("Path '%s' not found!" % p)
#         node[n] = value


# def next_step():
#     for i, ind in enumerate(indices):
#         indices[i] += 1
#         if indices[i] < len(parameters[i]["values"]):
#             set_in_yml(
#                 parameters[i]["paths"], parameters[i]["values"][indices[i]])
#             return True
#         else:
#             indices[i] = 0
#             set_in_yml(
#                 parameters[i]["paths"], parameters[i]["values"][indices[i]])



# if os.path.exists(args.parameters):
#     parameters = imp.load_source("parameters", args.parameters).parameters
#     single = False
# else:
#     single = True
# run_cnt = 0
# def run_description():
#     res = ""
#     for i, ind in enumerate(indices):
#         res += "%s = %s\n" % (parameters[i]
#                               ["name"], parameters[i]["values"][ind])
#     return res


# def run_description_csv(run_label):
#     res = "\"%s\"" % run_label
#     for i, ind in enumerate(indices):
#         res += ",\"%s\"" % parameters[i]["values"][ind]
#     return res
