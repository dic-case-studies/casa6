
import argparse
import os
import traceback
from casatools import table
############################################################
##
##
## TAQL Documentation:
## https://casacore.github.io/casacore-notes/199.html
## https://casa.nrao.edu/docs/CasaRef/table.taql.html
## https://casa.nrao.edu/aips2_docs/notes/199/node3.html
##
##
############################################################

def set_taql_query(itable, otable = None, lrows='2', overwrite = False):
    '''Select the limit of rows to keep from an input table and copy them to an output table'''

    if otable is None:
        otable = os.path.basename(itable) + '_OUT'

    tb=table()
    
    # TaQL query command 
    selq="SELECT from "+itable+" LIMIT "+ str(lrows)
    #selq="SELECT from "+msfile+" WHERE FLAG_ROW==TRUE && NOT ALL(FLAG)==TRUE LIMIT 10"
    
    # Open the table
    try:
        tb.open(itable, nomodify=True)
        print('Running TaQL selection on {}'.format(itable))
        print('Original table has {} rows'.format(tb.nrows()))
    
        # Query the table. It returns an instance of the tb tool
        seltbl = tb.taql(selq)
        nrow = seltbl.nrows()
        print('Selection has {} rows'.format(nrow))
        
        if (nrow > 0):
            if os.path.exists(otable) and overwrite == False:
                print('ERROR: Output table {} already exist. Will not overwrite it'.format(otable))
                return
                #os.system('rm -rf '+outtable)
                
            print('Saving selection to {}'.format(otable))
            out = seltbl.copy(otable,deep=True,valuecopy=True)
            out.close() 
        else:
            print('ERROR: {} has {} rows. Will not write to {}'.format(itable,nrow, otable))
            return
    #        seltbl.flush()
               
    finally:
        tb.close()
        seltbl.close()
        del seltbl

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help='name of input table to read select from', required=True)
    parser.add_argument("-o", "--output",help='name of output table to create', required=False)
    parser.add_argument("-n", "--nrows",help="number of rows to save in output", required=False)
     
    args, unknownArgs = parser.parse_known_args()
    print(args)
    
    input_table = ''
    if args.input is not None:
        input_table = args.input
    else:
        print('Please provide the name of a CASA table')
        sys.exit()
    
    output_table = os.path.basename(input_table) + '_OUT'
    if args.output is not None:
        output_table = args.output
         
    limit_rows = '2'
    if args.nrows is not None:
        limit_rows = str(args.nrows)
             
    try:
        set_taql_query(input_table, output_table, limit_rows)
    except:
        print('Cannot create output TaQL table')
        traceback.print_exc()

