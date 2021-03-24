from casatasks import listobs
import casatools
import sys
import os
import argparse

def get_observation_name(data):
    listobs(vis=data, listfile="getTelescope.out")
    # check for the Telescope used here
    with open ("getTelescope.out", "r") as fout:
        for line in fout:
            if "Observation:" in line:
                print("--------")
                print(line)
                print("Dataset: ", data)
                print("--------")

    os.remove("getTelescope.out")

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help='name of input table to read', required=True)

    args, unknownArgs = parser.parse_known_args()
    print(args)
    
    input_table = ''
    if args.input is not None:
        input_table = args.input
    else:
        print('Please provide the name of a dataset')
        sys.exit()

    try:
        get_observation_name(input_table)
    except:
        print('Cannot get the name of Observatory')
        traceback.print_exc()
