import copy
import numpy as np
import os

# from casatasks import casalog

class SummaryMinor:
    """Gathers the information together from the tclean return value in a way that makes it easier to query for exactly what you want.

    The structure for this nested dictionary is:
        {
            multi-field id: {
                channel id: {
                    stokes id: {
                        summary key: {
                            cycle: value
                        }
                    }
                }
            }
        }

    Examples:
    
        1. To get the number of available channels, and the ids of those channels:
            nchans = len(summ_min[0].keys())
            avail_chans = summ_min[0].keys()
    
        2. To get the number of iterations done on the main field, fifth channel, first stokes plane, during the middle minor cycle:
            field0 = summ_min[0] # field 0 is the main field
            chan = field0.keys()[5] # channel index doesn't necessarily start at 0
            stoke = field0[chan].keys()[0] # stokes index doesn't necessarily start at 0
            ncycles = len(field0[chan][stoke]['iterDone'])
            itersDone = field0[chan][stoke]['iterDone'][ncycles/2]

        3. To get the available minor cycle summary statistics:
            field0 = summ_min[0]
            chan0 = field0.keys()[0]
            stoke0 = field0[chan0].keys()[0]
            availSummStats = field0[field0][stoke0].keys()
    """
    #                           0           1          2            3              4          5       6      7                  8                9               10                11 "No Mask"      12           13         14           15         16         17              18
    rowDescriptionsOldOrder = ["iterDone", "peakRes", "modelFlux", "cycleThresh", "deconId", "chan", "stoke", "cycleStartIters", "startIterDone", "startPeakRes", "startModelFlux", "startPeakResNM", "peakResNM", "masksum", "mpiServer", "peakMem", "runtime", "multifieldId", "stopCode"]
    rowDescriptions13683    = ["iterDone", "peakRes", "modelFlux", "cycleThresh", "deconId", "chan"]
    # rowDescriptions does not include {"multifieldId", "chan", "stoke", "deconId"}, and so the returned dictionary will not include those values in the summary keys
    rowDescriptions         = ["startIterDone", "iterDone", "startPeakRes", "peakRes", "startModelFlux", "modelFlux", "startPeakResNM", "peakResNM", "cycleThresh", "cycleStartIters", "masksum", "mpiServer", "peakMem", "runtime", "stopCode"]
    rowStartDescs           = ["startIterDone",             "startPeakRes",            "startModelFlux",              "startPeakResNM"]

    def convertMatrix(summaryminor_matrix, calc_iterdone_deltas=None, keep_startvals=None):
        # casalog.post(summaryminor_matrix, "SEVERE")
        ret = {}

        # edge case: no iterations were done (eg threshold < model flux)
        if summaryminor_matrix.shape[1] == 0:
            return { 0: {} }

        # get individual dictionaries for each field id
        field_ids = SummaryMinor._getFieldIds(summaryminor_matrix)
        if len(field_ids) > 1:
            for fieldId in field_ids:
                singleFieldMatrix = SummaryMinor._getSingleFieldMatrix(summaryminor_matrix, field_ids[fieldId])
                ret[fieldId] = SummaryMinor._convertSingleFieldMatrix(singleFieldMatrix, calc_iterdone_deltas, keep_startvals)
        elif len(field_ids) == 1:
            ret[field_ids[0]] = SummaryMinor._convertSingleFieldMatrix(summaryminor_matrix, calc_iterdone_deltas, keep_startvals)
        else:
            raise RuntimeError("No multifield ids were found. Failed to parse summary minor matrix after tclean finished running.")

        return ret

    def _convertSingleFieldMatrix(single_field_matrix, calc_iterdone_deltas=None, keep_startvals=None):
        # edge case: no iterations were done (eg threshold < model flux)
        if single_field_matrix.shape[1] == 0:
            return {}

        summaryminor_dict = SummaryMinor.indexMinorCycleSummaryBySubimage(single_field_matrix)
        percycleiters_dict = SummaryMinor._getPerCycleDict(copy.deepcopy(summaryminor_dict), calc_iterdone_deltas, keep_startvals)
        return percycleiters_dict

    def _getFieldIds(matrix):
        """ Get a sorted list of available outlier field ids in the given matrix """

        # edge case: running with MPI and CAS-13683 hasn't been fixed yet
        availRows = SummaryMinor.getRowDescriptionsOldOrder()
        if not "multifieldId" in availRows:
            return [0] # can't differentiate multiple fields from available data, assume one field

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
            oldStokeIdx  = SummaryMinor.getRowDescriptionsOldOrder().index("stoke")
        chans = list(np.sort(np.unique(matrix[oldChanIdx])))
        chans = [int(x) for x in chans]
        if uss:
            stokes = [0]
        else:
            stokes = list(np.sort(np.unique(matrix[oldStokeIdx])))
            stokes = [int(x) for x in stokes]
        ncycles = 0
        if len(chans) > 0 and len(stokes) > 0:
            ncycles = int( ncols / (len(chans)*len(stokes)) )
            if uss:
                try:
                    from casampi.MPIEnvironment import MPIEnvironment
                    if MPIEnvironment.is_mpi_enabled:
                        # This is necessary because we may have an odd number of "channels" due to each process getting only a subchunk.
                        # Example:
                        #     Process 1 gets stokes 0-1, process 2 gets stokes 2
                        #     Each of them assigns channel id = chan + stoke * nsubstokes
                        #     Process 1 assigns channel ids [0,2], Process 2 assigns channel id 0.
                        # This hack is not needed when not using a small summary minor because we have the extra knowledge of the stokes, instead of mapping stokes + channels onto chunks.
                        chanslist = matrix[oldChanIdx].tolist()
                        for chan in chans:
                            singlechan_occurances = list(filter(lambda x: x == chan, chanslist))
                            ncycles = max(ncycles, len(singlechan_occurances))
                except ModuleNotFoundError as e:
                    raise

        # ret is the return dictionary[chans][stokes][rows][cycles]
        # cummulativeCnt counts how many cols we've read for each channel/stokes/row
        ret = {desc:[0]*ncycles for desc in SummaryMinor.getRowDescriptions()}
        ret = {stoke:copy.deepcopy(ret) for stoke in stokes}
        ret = {chan:copy.deepcopy(ret) for chan in chans}
        cummulativeCnt = copy.deepcopy(ret) # copy ret's structure

        # reindex based on subimage index (aka chan/stoke index)
        for desc in SummaryMinor.getRowDescriptions():
            oldRowIdx = SummaryMinor.getRowDescriptionsOldOrder().index(desc)
            for colIdx in range(ncols):
                chan = int(matrix[oldChanIdx][colIdx])
                if uss:
                    stoke = 0
                else:
                    stoke = int(matrix[oldStokeIdx][colIdx])
                val = matrix[oldRowIdx][colIdx]
                cummulativeCol = int(cummulativeCnt[chan][stoke][desc][0]) # const 0: cummulativeCnt doesn't make use of 'cycle' index from copied ret structure
                ret[chan][stoke][desc][cummulativeCol] = val
                cummulativeCnt[chan][stoke][desc][0] += 1

        return ret

    def _getPerCycleDict(summaryminor_dict, calc_iterdone_deltas=None, keep_startvals=None):
        calc_iterdone_deltas = True if (calc_iterdone_deltas == None) else calc_iterdone_deltas
        keep_startvals       = True if (keep_startvals == None)       else keep_startvals
        ret = summaryminor_dict
        availRows = SummaryMinor.getRowDescriptionsOldOrder()

        if (calc_iterdone_deltas) and ("startIterDone" in availRows):
            for chan in ret:
                for stoke in ret[chan]:
                    for cyc in range(len(ret[chan][stoke]["startIterDone"])):
                        ret[chan][stoke]["iterDone"][cyc] -= ret[chan][stoke]["startIterDone"][cyc]
        if not keep_startvals:
            for chan in ret:
                for stoke in ret[chan]:
                    for desc in SummaryMinor.getRowStartDescs():
                        del ret[chan][stoke][desc]

        return ret

