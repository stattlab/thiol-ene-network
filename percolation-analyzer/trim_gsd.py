#!/usr/bin/env python3
import numpy as np
import sys, re, os, gsd
import gsd.hoomd
import argparse
import math

#python3 trim_gsd.py -i /home/bj21/simulations/nvt_vs_npt/out/nvt_deform_finalLam_5_lamella_611_x2.gsd -o nvt_deform_finalLam_5_lamella_611_x2_trimmed.gsd -p 5

def convert_size(size_bytes):
    ''' small helper function to convert bytes to human readable format '''
    if size_bytes == 0:
        return "0 B"
    size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
    i = int(math.floor(math.log(size_bytes, 1024)))
    p = math.pow(1024, i)
    s = round(size_bytes / p, 2)
    return "%s %s" % (s, size_name[i])

############################### PARSE ##################################
parser = argparse.ArgumentParser(description='Selectively delete frames from a gsd-file')

non_opt = parser.add_argument_group('mandatory arguments')

non_opt.add_argument('-o','--out',\
                metavar='<gsd-file>', dest='output_file',type=str,\
                required=False,help='Name for the thinned out gsd file (not needed if --info is used)')

non_opt.add_argument('-i','--in',\
                metavar='<gsd-file>', dest='input_file',type=str,\
                required=True,help='input trajectory *.gsd file')

parser.add_argument('--info', dest='info', action='store_true',\
                    help=" Only print file size information and potential period values and exit without writing a new file")

parser.add_argument('-p','--p', metavar='<int>', type=int, default=1, required=False,dest='period',
    help='Period of frames to take, -p 10 means take every 10th frame, DEFAULT=1/every')

parser.add_argument('-s','--skip',metavar='<int>', dest='skip',type=int,required=False,default=0,
    help='Number of frames to skip at the beginning (positive s) or take from the end (negative s) of the trajectory. DEFAULT=0/no skip')

parser.add_argument('-l','--list',\
                metavar='<ints>', dest='to_del', nargs='+', type=int,\
                required=False,help='frame numbers to delete, DEFAULT=none')

args = parser.parse_args()
input_file = args.input_file
to_del = args.to_del

if not os.path.isfile(input_file):
    print("ERROR: file does not exist! %s"%(input_file))
    exit(1)
if args.info==True:
    # print some info about the file and exit without writing a new file
    print("Only printing information about gsd file: ", input_file)
    trajectory = gsd.hoomd.open(name=input_file, mode='r')
    print("Total number of frames: ", len(trajectory))
    if len(trajectory)==0:
        print("ERROR: trajectory does not seem to contain any snapshots, check file")
        exit(1)
    # more info about file content
    bytes = os.path.getsize(input_file)
    print("Total file size:",convert_size(bytes))
    s_per_snap = (bytes/len(trajectory))
    print("Size per snapshot", convert_size(s_per_snap))
    print("potential values for p:")
    print("period  frames   size")
    print("--------------------------")
    for p in [2,3,5,10]:
        slice = trajectory[::p]
        new_s_per_snap = s_per_snap*len(slice)
        print(p,"    \t",len(slice),"\t",convert_size(new_s_per_snap))
    exit(0)

if args.info==False and (args.input_file is None or args.output_file is None):
    parser.error("Output file is required, specifiy -o/--out <file> or use --info")

output_file = args.output_file
if input_file==output_file:
    parser.error("Input and output file cannot be the SAME file! %s"%(output_file))
period = args.period
if period<1:
    parser.error("Period needs to be >= 1! Current value %s"%period)
skip = args.skip

trajectory = gsd.hoomd.open(name=input_file, mode='r')
N_frames = len(trajectory)
# figure out which frames to analyze
if np.abs(skip)>=N_frames:
    print("ERROR: trajectory only has %d frames, skip=%d is too large!"%(N_frames,skip))
    exit(1)
if skip>=0:
    frames_to_write=np.arange(skip,N_frames,1).astype(int)
else:
    frames_to_write=np.arange(N_frames+skip,N_frames,1).astype(int)

frames_to_write=frames_to_write[::period]

# Return the unique values in frames_to_write that are not in to_del.
frames_to_write = np.setdiff1d(frames_to_write,to_del)

if len(frames_to_write)<=0:
    print("ERROR: frames to write seems to be <=0! %s "%(frames_to_write))
    exit(1)
print("total number of frames ",N_frames)
print("writing frames ",frames_to_write[0],"...",frames_to_write[-1])
print("writing every nth frame: ",period)

#gsd.fl.create(name=output_file, application='gsd.hoomd ' + gsd.__version__, schema='hoomd', schema_version=[1,1])

try:
    t =  gsd.hoomd.open(name=f'{os.path.dirname(input_file)}/{output_file}', mode='w')
except:
    t =  gsd.hoomd.open(name=output_file, mode='w')

# Iterate over all frames of the sequence.
for i in frames_to_write:
    sys.stdout.write("\rCurrent frame %3d/%3d"%(i,N_frames))
    sys.stdout.flush()
    t.append(trajectory[int(i)])


sys.stdout.write("\nDONE\n")
sys.stdout.flush()
