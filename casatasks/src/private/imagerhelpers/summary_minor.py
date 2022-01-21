import copy
import numpy as np
import os

# from casatasks import casalog

class SummaryMinor:
    """Gathers the information together from the tclean return value in a way that makes it easier to query for exactly what you want.

    The structure for this nested dictionary is:
        {
            multi-field id: {
                "chans": [chan0, chan1, ...],
                "stoks": [stok0, stok1, ...],
                "ncycs": number of major/minor cycles
                channel id: {
                    stokes id: {
                        summary key: {
                            cycle: value
                        }
                    }
                }
            },
            "matrix": original numpy ndarray matrix returned from the cpp code
        }

    Examples:
    
        1. To get the number of iterations done on the first field, fifth channel, first stokes plane, during the middle minor cycle:
            field0 = summ_min[0]
            chan   = field0['chans'][5] # channel index doesn't necessarily start at 0
            stok   = field0['stoks'][0] # stokes index doesn't necessarily start at 0
            cycle  = field0['ncycs']/2
            chan5iters = field0[field][chan][stok]['iterDone'][cycle]
    
        2. To get the number of available channels, and the ids of those channels:
            nchans = len(summ_min[0]["chans"])
            avail_chans = summ_min[0]["chans"]

        3. To get the available minor cycle summary statistics:
            field0 = summ_min[0]
            summaryKeys = field0[field0['chans'][0]][field0['stoks'][0]].keys()
    """
    #                           0           1          2            3              4          5       6      7                  8                9               10                11 "No Mask"      12           13         14           15         16         17              18
    rowDescriptionsOldOrder = ["iterDone", "peakRes", "modelFlux", "cycleThresh", "deconId", "chan", "stok", "cycleStartIters", "startIterDone", "startPeakRes", "startModelFlux", "startPeakResNM", "peakResNM", "masksum", "mpiServer", "peakMem", "runtime", "multifieldId", "stopCode"]
    rowDescriptions13683    = ["iterDone", "peakRes", "modelFlux", "cycleThresh", "deconId", "chan"]
    # rowDescriptions does not include {"multifieldId", "chan", "stok", "deconId"}, and so the returned dictionary will not include those values in the summary keys
    rowDescriptions         = ["startIterDone", "iterDone", "startPeakRes", "peakRes", "startModelFlux", "modelFlux", "startPeakResNM", "peakResNM", "cycleThresh", "cycleStartIters", "masksum", "mpiServer", "peakMem", "runtime", "stopCode"]
    rowStartDescs           = ["startIterDone",             "startPeakRes",            "startModelFlux",              "startPeakResNM"]

    def convertMatrix(summaryminor_matrix, calc_iterdone_deltas=None, keep_startvals=None):
        ret = {}

        # edge case: no iterations
        if summaryminor_matrix.shape[1] == 0:
            return { 0: {
                'chans': [],
                'stoks': [],
                'ncycs': 0
            }}

        # get individual dictionaries for each field id
        field_ids = SummaryMinor._getFieldIds(summaryminor_matrix)
        if len(field_ids) > 1:
            for fieldId in field_ids:
                singleFieldMatrix = SummaryMinor._getSingleFieldMatrix(summaryminor_matrix, field_ids[0])
                ret[fieldId] = SummaryMinor._convertSingleFieldMatrix(singleFieldMatrix, calc_iterdone_deltas, keep_startvals)
        elif len(field_ids) == 1:
            ret[field_ids[0]] = SummaryMinor._convertSingleFieldMatrix(summaryminor_matrix, calc_iterdone_deltas, keep_startvals)
        else:
            raise RuntimeError("No multifield ids were found. Failed to parse summary minor matrix after tclean finished running.")

        # insert convenience information
        for field_id in field_ids:
            chan_ids = sorted(list( ret[field_id].keys() ))
            stok_ids  = []
            nCycles = 0

            if len(chan_ids) > 0:
                chan0  = chan_ids[0]
                stok_ids = sorted(list( ret[field_id][chan0].keys() ))
                if len(stok_ids) > 0:
                    stok0    = stok_ids[0]
                    nCycles = len( ret[field_id][chan0][stok0]['iterDone'] )

            ret[field_id]['chans'] = chan_ids
            ret[field_id]['stoks'] = stok_ids
            ret[field_id]['ncycs'] = nCycles

        return ret

    def _convertSingleFieldMatrix(single_field_matrix, calc_iterdone_deltas=None, keep_startvals=None):
        summaryminor_dict = SummaryMinor.indexMinorCycleSummaryBySubimage(single_field_matrix)
        percycleiters_dict = SummaryMinor._getPerCycleDict(copy.deepcopy(summaryminor_dict), calc_iterdone_deltas, keep_startvals)
        return percycleiters_dict

    def _getFieldIds(matrix):
        """ Get a sorted list of available outlier field ids in the given matrix """
        availRows = SummaryMinor.getRowDescriptionsOldOrder()
        if not "multifieldId" in availRows:
            return [0]
        multifieldIdx = availRows.index("multifieldId")
        nrows = matrix.shape[0]
        ncols = matrix.shape[1]
        fieldIds = sorted(np.unique(matrix[multifieldIdx,:]).tolist())
        fieldIds = list(map(lambda x: int(x), fieldIds))
        return fieldIds

    def _getSingleFieldMatrix(matrixIn, fieldId):
        """ Create a new matrix to hold all the values of the given matrix, but only for the given outlier field id """
        availRows = SummaryMinor.getRowDescriptionsOldOrder()
        if not "multifieldId" in availRows:
            return matrixIn
        multifieldIdx = availRows.index("multifieldId")
        nrowsIn = matrixIn.shape[0]
        ncolsIn = matrixIn.shape[1]
        nrowsOut = nrowsIn
        ncolsOut = matrixIn[multifieldIdx,:].tolist().count(fieldId)

        matrixOut = np.zeros((nrowsOut, ncolsOut))
        colOut = 0
        maxColOut = 0
        maxRowOut = 0
        for colIn in range(ncolsIn):
            if matrixIn[multifieldIdx,colIn] != fieldId:
                continue
            for rowIn in range(nrowsIn):
                rowOut = rowIn
                matrixOut[rowOut,colOut] = matrixIn[rowIn,colIn]
                maxRowOut = max(rowOut, maxRowOut)
            maxColOut = colOut
            colOut += 1

        return matrixOut

    def useSmallSummaryminor(ignored_parameter=None):
        """Temporary CAS-13683 workaround"""
        if ('USE_SMALL_SUMMARYMINOR' in os.environ):
            uss = os.environ['USE_SMALL_SUMMARYMINOR'].lower()
            if (uss == "true"):
                return True
        return False

    def _getRowDescriptionsOldOrder(useSmallSummaryminor):
        """Temporary CAS-13683 workaround"""
        if (useSmallSummaryminor):
            return SummaryMinor.rowDescriptions13683
        return SummaryMinor.rowDescriptionsOldOrder

    def getRowDescriptionsOldOrder():
        """ Retrieves brief descriptions of the available minor cycle summary rows, in the old (matrix) order. """
        return SummaryMinor._getRowDescriptionsOldOrder(SummaryMinor.useSmallSummaryminor())

    def _getRowDescriptions(useSmallSummaryminor):
        """Temporary CAS-13683 workaround"""
        ret = SummaryMinor.rowDescriptions
        availRows = SummaryMinor._getRowDescriptionsOldOrder(useSmallSummaryminor)
        ret = list(filter(lambda x: x in availRows, ret))
        return ret

    def getRowDescriptions():
        """ Retrieves brief descriptions of the available minor cycle summary rows """
        return SummaryMinor._getRowDescriptions(SummaryMinor.useSmallSummaryminor())

    def _getRowStartDescs(useSmallSummaryminor):
        """Temporary CAS-13683 workaround"""
        ret = SummaryMinor.rowStartDescs
        availRows = SummaryMinor._getRowDescriptionsOldOrder(useSmallSummaryminor)
        ret = list(filter(lambda x: x in availRows, ret))
        return ret

    def getRowStartDescs():
        """ Retrieves abreviated names of the available minor cycle summary "start" rows.

        These are the rows that catalog the values at the beggining of a minor cycle (pre-deconvolution). """
        return SummaryMinor._getRowStartDescs(SummaryMinor.useSmallSummaryminor())

    def indexMinorCycleSummaryBySubimage(matrix):
        """ Re-indexes matrix from [row,column] to [channel,stokes,row,cycle]. 

        Param matrix: the original matrix to convert.
        """
        # get some properties of the summary_minor matrix
        nrows = matrix.shape[0]
        ncols = matrix.shape[1]
        uss = SummaryMinor.useSmallSummaryminor() # Temporary CAS-13683 workaround
        import sys
        oldChanIdx = SummaryMinor.getRowDescriptionsOldOrder().index("chan")
        if not uss:
            oldStokIdx  = SummaryMinor.getRowDescriptionsOldOrder().index("stok")
        chans = list(np.sort(np.unique(matrix[oldChanIdx])))
        chans = [int(x) for x in chans]
        if uss:
            stoks = [0]
        else:
            stoks = list(np.sort(np.unique(matrix[oldStokIdx])))
            stoks = [int(x) for x in stoks]
        ncycles = 0
        if len(chans) > 0 and len(stoks) > 0:
            ncycles = int( ncols / (len(chans)*len(stoks)) )
            if uss:
                try:
                    from casampi.MPIEnvironment import MPIEnvironment
                    if MPIEnvironment.is_mpi_enabled:
                        # This is necessary because we may have an odd number of "channels" due to each process getting only a subchunk.
                        # Example:
                        #     Process 1 gets stokes 0-1, process 2 gets stokes 2
                        #     Each of them assigns channel id = chan + stok * nsubstoks
                        #     Process 1 assigns channel ids [0,2], Process 2 assigns channel id 0.
                        # This hack is not needed when not using a small summary minor because we have the extra knowledge of the stokes, instead of mapping stokes + channels onto chunks.
                        chanslist = matrix[oldChanIdx].tolist()
                        for chan in chans:
                            singlechan_occurances = list(filter(lambda x: x == chan, chanslist))
                            ncycles = max(ncycles, len(singlechan_occurances))
                except ModuleNotFoundError as e:
                    raise

        # ret is the return dictionary[chans][stoks][rows][cycles]
        # cummulativeCnt counts how many cols we've read for each channel/stokes/row
        ret = {desc:[0]*ncycles for desc in SummaryMinor.getRowDescriptions()}
        ret = {stok:copy.deepcopy(ret) for stok in stoks}
        ret = {chan:copy.deepcopy(ret) for chan in chans}
        cummulativeCnt = copy.deepcopy(ret) # copy ret's structure

        # reindex based on subimage index (aka chan/stok index)
        for desc in SummaryMinor.getRowDescriptions():
            oldRowIdx = SummaryMinor.getRowDescriptionsOldOrder().index(desc)
            for colIdx in range(ncols):
                chan = int(matrix[oldChanIdx][colIdx])
                if uss:
                    stok = 0
                else:
                    stok = int(matrix[oldStokIdx][colIdx])
                val = matrix[oldRowIdx][colIdx]
                cummulativeCol = int(cummulativeCnt[chan][stok][desc][0]) # const 0: cummulativeCnt doesn't make use of 'cycle' index from copied ret structure
                ret[chan][stok][desc][cummulativeCol] = val
                cummulativeCnt[chan][stok][desc][0] += 1

        return ret

    def _getPerCycleDict(summaryminor_dict, calc_iterdone_deltas=None, keep_startvals=None):
        calc_iterdone_deltas = True if (calc_iterdone_deltas == None) else calc_iterdone_deltas
        keep_startvals       = True if (keep_startvals == None)       else keep_startvals
        ret = summaryminor_dict
        availRows = SummaryMinor.getRowDescriptionsOldOrder()

        if (calc_iterdone_deltas) and ("startIterDone" in availRows):
            for chan in ret:
                for stok in ret[chan]:
                    for cyc in range(len(ret[chan][stok]["startIterDone"])):
                        ret[chan][stok]["iterDone"][cyc] -= ret[chan][stok]["startIterDone"][cyc]
        if not keep_startvals:
            for chan in ret:
                for stok in ret[chan]:
                    for desc in SummaryMinor.getRowStartDescs():
                        del ret[chan][stok][desc]

        return ret