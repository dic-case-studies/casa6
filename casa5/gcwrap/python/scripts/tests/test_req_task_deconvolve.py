##########################################################################
##########################################################################
#
# Test programs for the minor cycle deconvolve task :  test_req_task_deconvolve
#
# Each of the following categories (classes) has a set of tests within it.
#
#  test_onefield                 # basic tests, deconvolution algorithms
#  test_iterbot                  # iteration control options for mfs and cube
#  test_multifield               # multiple fields of same type and with different shapes/deconvolvers/gridders
#  test_stokes                   # multiple stokes planes, imaging with flagged correlations..
#  test_cube                     # all things cube. Spectral frame setup, handling empty channels, etc
#  test_mask                     # input mask options : regridding, mask file, automasking, etc
#  test_multirun                 # running deconvolve multiple times in a row
#  test_imgval                   # validation of input images and copying/converting startmodel images for the hogbom deconvolver
#  test_mtmfsimgval              # validation of input images and copying/converting startmodel images for the mtmfs deconvolver
#  test_residual_update          # .residual images should get updated between consecutive runs of deconvolve
#
# To run from within casapy :  
#
#  runUnitTest.main(['test_req_task_deconvolve'])                                              # Run all tests
#  runUnitTest.main(['test_req_task_deconvolve[test_onefield]'])                               # Run tests from test_onefield
#  runUnitTest.main(['test_req_task_deconvolve[test_onefield_mtmfs]'])                         # Run one specific test
#  runUnitTest.main(['test_req_task_deconvolve[test_onefield_mtmfs,test_onefield_hogbom]'])    # Multiple specific tests
#
# To see the full list of tests :   grep "\"\"\" \[" test_req_task_deconvolve.py
#
#  These tests need data stored in casa-data/regression/unittest/clean/refimager
#  The test_mask_autobox_multithresh_standard_cube_eph test needs data stored in asa-data-req/stakeholders/alma
#
#  For a developer build, to get the datasets locally
#  - https://open-bitbucket.nrao.edu/projects/CASA/repos/casa-data/browse
#  - https://open-bitbucket.nrao.edu/scm/casa/casa-data-req.git
#
##########################################################################
#
#  Datasets
#
#  refim_twochan.ms : 2 channels, one 1Jy point source with spectral index of -1.0
#  refim_point.ms : 1-2 GHz, 20 channels, 1 spw, one 1Jy point source with spectral index -1.0.
#  refim_point_linRL.ms : I=1, Q=2, U=3, V=4  in circular pol basis.
#  2017.1.00750.T_tclean_exe*.ms: testing data from casa-data-req/stakeholders/alma, for the pipeline test test_mask_autobox_multithresh_standard_cube_eph
#
#  task_tclean is used to generate the input tables to task_deconvolve from these data, then
#  the test is run on task_deconvolve. It is probable, therefore, that if task_tclean tests
#  start failing, then task_deconvolve tests will also start failing.
#
##########################################################################
#
#  Tests
#
#Onefield tests (mimicing tclean)
#1. Onefield Hogbom: Should produce the same results as the test by the same name for tclean.
#imsize:100, cell:'8.0arcsec', deconvolver:'hogbom', niter:10
#testname: test_onefield_hogbom
#
#2. Onefield Clark: Should produce the same results as the test by the same name for tclean.
#imsize:100, cell:'8.0arcsec', deconvolver:'clark', niter:10
#testname: test_onefield_clark
#
#3. Onefield Multiscale: Should produce the same results as the test by the same name for tclean.
#imsize:200, cell:'8.0arcsec', deconvolver:'multiscale', scales:[0,20,40,100], niter:10
#testname: test_onefield_multiscale
#
#4. Onefield mfmfs: Should produce the same results as the test by the same name for tclean.
#imsize:100, cell:'8.0arcsec', deconvolver:'mtmfs', niter:10
#testname: test_onefield_mtmfs
#
#5. onefield_rectangular_pixels(self): Should produce the same results as the test by the same name for tclean.
#imsize:100, cell:['10.0arcsec','30.0arcsec'], niter:10
#testname: test_onefield_rectangular_pixels
#
#
#
#Iterbot tests (mimicing tclean)
#6. Iterbot Clark mfs: Should produce the same results as the test by the same name for tclean.
#imsize:100,cell:'8.0arcsec',deconvolver:'clark',gain:0.15, niter=20
#testname: test_iterbot_mfs_4
#
#7. Iterbot Clark mfs: Should produce the same results as the test by the same name for tclean.
#imsize:100,cell:'8.0arcsec',deconvolver:'clark',threshold:'0.5Jy',gain:0.15, niter=10, threshold='0.5Jy'
#testname: test_iterbot_mfs_5
#
#
#
#Stokes tests (mimicing tclean)
#8. Stokes I mfs: Should produce the same results as the test by the same name for tclean.
#imsize:100,cell:'8.0arcsec',stokes:'I', niter=10
#testname: test_stokes_mfs_I
#
#9. Stokes IQUV mtmfs: Should produce the same results as the test by the same name for tclean.
#imsize:100,cell:'8.0arcsec',stokes:'IQUV',deconvolver:'mtmfs',nterms:2, niter=10
#testname: test_stokes_mtmfs_IQUV
#
#
#
#Cube tests (mimicing tclean)
#10. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': 0, 'width': 1, 'veltype': 'radio', 'outframe': 'LSRK', 'interpolation': 'linear'
#testname: test_cube_0
#
#11. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': 0, 'width': 1, 'veltype': 'radio', 'outframe': 'TOPO', 'interpolation': 'linear'
#testname: test_cube_1
#
#12. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': 0, 'width': 2, 'veltype': 'radio', 'outframe': 'LSRK', 'interpolation': 'linear'
#testname: test_cube_2
#
#13. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': 5, 'width': 1, 'veltype': 'radio', 'outframe': 'LSRK', 'interpolation': 'linear'
#testname: test_cube_3
#
#14. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0:5~19', 'start': 0, 'width': 1, 'veltype': 'radio', 'outframe': 'LSRK', 'interpolation': 'linear'
#testname: test_cube_4
#
#15. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': '', 'width': '100MHz', 'veltype': 'radio', 'outframe': 'TOPO', 'interpolation': 'linear'
#testname: test_cube_5
#
#16. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': '1.1GHz', 'width': '', 'veltype': 'radio', 'outframe': 'TOPO', 'interpolation': 'linear'
#testname: test_cube_6
#
#17. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0:4~19', 'start': '1.1GHz', 'width': '', 'veltype': 'radio', 'outframe': 'TOPO', 'interpolation': 'linear'
#testname: test_cube_7
#
#18. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': '1.5GHz', 'width': '-50MHz', 'veltype': 'radio', 'outframe': 'TOPO', 'interpolation': 'linear'
#testname: test_cube_8
#
#19. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': '', 'width': '23983.4km/s', 'veltype': 'radio', 'outframe': 'TOPO', 'interpolation': 'linear'
#testname: test_cube_9
#
#20. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': '', 'width': '-23983.4km/s', 'veltype': 'radio', 'outframe': 'TOPO', 'interpolation': 'linear'
#testname: test_cube_10
#
#21. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': '11991.7km/s', 'width': '', 'veltype': 'radio', 'outframe': 'TOPO', 'interpolation': 'linear'
#testname: test_cube_11
#
#22. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': '11977.6km/s', 'width': '', 'veltype': 'radio', 'outframe': 'BARY', 'interpolation': 'linear'
#testname: test_cube_12
#
#23. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', specmode': 'cube', 'nchan': 8, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': '-41347.8km/s', 'width': '20000km/s', 'veltype': 'optical', 'outframe': 'LSRK'
#testname: test_cube_13
#
#24. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': {'unit': 'GHz', 'value': 1.2}, 'width': '', 'veltype': 'radio', 'outframe': '', 'interpolation': 'linear'
#testname: test_cube_14
#
#25. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': {'m0': {'unit': 'Hz', 'value': 1199989000.0}, 'refer': 'LSRK', 'type': 'frequency'}, 'width': '', 'veltype': 'radio', 'outframe': '', 'interpolation': 'linear'
#testname: test_cube_15
#
#26. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': {'unit': 'km/s', 'value': 11991.7}, 'width': '', 'veltype': 'radio', 'outframe': 'TOPO', 'interpolation': 'linear'
#testname: test_cube_16
#
#27. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': {'m0': {'unit': 'm/s', 'value': 11977600.0}, 'refer': 'BARY', 'type': 'radialvelocity'}, 'width': '', 'veltype': 'radio', 'outframe': 'TOPO', 'interpolation': 'linear'
#testname: test_cube_17
#
#28. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': '', 'width': {'unit': 'km/s', 'value': 11991.7}, 'veltype': 'radio', 'outframe': 'TOPO', 'interpolation': 'linear'
#testname: test_cube_18
#
#29. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': '', 'width': {'m0': {'unit': 'm/s', 'value': 11991700.0}, 'refer': 'TOPO', 'type': 'radialvelocity'}, 'veltype': 'radio', 'outframe': 'TOPO', 'interpolation': 'linear'
#testname: test_cube_19
#
#30. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0', 'start': {'m0': {'unit': 'm/s', 'value': 11994336.49363042}, 'refer': 'RADIO', 'type': 'doppler'}, 'width': '', 'veltype': 'radio', 'outframe': 'LSRK', 'interpolation': 'linear'
#testname: test_cube_20
#
#31. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0:4~9;12~14', 'start': 4, 'width': '', 'veltype': 'radio', 'outframe': 'LSRK', 'interpolation': 'nearest'
#testname: test_cube_21
#
#32. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0:0~10^2', 'start': 0, 'width': '', 'veltype': 'radio', 'outframe': 'LSRK', 'interpolation': 'nearest'
#testname: test_cube_22
#
#33. Cube: Should produce the same results as the test by the same name for tclean.
#'field': '0', 'imsize': 100, 'cell': '8.0arcsec', 'specmode': 'cube', 'nchan': 10, 'restfreq': ['1.25GHz'], 'phasecenter': 'J2000 19:59:28.500 +40.44.01.50', 'deconvolver': 'hogbom', 'spw': '0:4~13', 'start': '', 'width': '', 'veltype': 'radio', 'outframe': 'TOPO', 'interpolation': 'linear'
#testname: test_cube_23
#
#34. Cube, 'chanchunks' Auto: Should produce the same results as the test by the same name for tclean.
#specmode:'cube',imsize:100,cell:'10.0arcsec',deconvolver:'hogbom',niter=10, deconvolver='hogbom'
#testname: test_cube_chanchunks_auto
#
#
#
#Masking tests (mimicing tclean)
#35. User Mask: Should produce the same results as the test by the same name for tclean.
#imsize:100,cell:'8.0arcsec',deconvolver:'hogbom',usemask:'user',mask:self.img+'.mask.txt', niter=10
#deconvolve runs: mask='tst.mask.txt', mask=mstr
#testname: test_mask_1
#
#36. User Mask: Should produce the same results as the test by the same name for tclean.
#imsize:100,cell:'8.0arcsec',deconvolver:'hogbom',specmode:'cube',interactive:0,usemask:'user', niter=10
#two mask runs: mask:'tst.mask.txt', mask:mstr
#testname: test_mask_2
#
#37. Missing Mask File: tst.mask is sometimes required
#mask='tst.model.txt'
#rm 'tst.model.txt', mask='tst.model.txt'
#testname: test_mask_missingfile
#
#38. Basic PB Mask: create .mask from .pb with default pbmask value
#rm 'tst.mask', usemask='pb', don't set pbmask
#verify the mask file exists and all pixels are 1s
#testname: test_mask_pbmask0
#
#39. Threshold PB Mask: create .mask from .pb with a pb threshold for masking
#rm 'tst.mask', usemask='pb', pbmask=0.995
#verify the mask file exists and only certain pixels are 1s
#testname: test_mask_pbmask9950
#
#40. Auto Mask: Should produce the same results as the test by the same name for tclean.
#imsize:100,cell:'8.0arcsec',deconvolver:'hogbom',interactive:0,usemask:'auto-multithresh',niter=10
#testname: test_mask_autobox_multithresh
#
#41. Auto Mask with New Noise Calc: Should produce the same results as the test by the same name for tclean.
#imsize:100,cell:'8.0arcsec',deconvolver:'hogbom',interactive:0,usemask:'auto-multithresh',niter=10,fastnoise=False
#testname: test_mask_autobox_multithresh_newnoise
#
#42. Auto Mask with Nsigma 3.0: Should produce the same results as the test by the same name for tclean.
#imsize:100,cell:'8.0arcsec',deconvolver:'hogbom',interactive:0,usemask:'auto-multithresh',niter=10,nsigma=3.0
#testname: test_mask_autobox_multithresh_with_nsigma
#
#43. Auto Mask with New Noise Calc and Nsigma 3.0: Should produce the same results as the test by the same name for tclean.
#imsize:100,cell:'8.0arcsec',deconvolver:'hogbom',interactive:0,usemask:'auto-multithresh',niter=10,fastnoise=False,nsigma=3.0
#testname: test_mask_autobox_multithresh_with_nsigma_newnoise
#
#44. Auto Mask with Pruning: Should produce the same results as the test by the same name for tclean.
#imsize:100,cell:'8.0arcsec',deconvolver:'hogbom',interactive:0,usemask:'auto-multithresh',niter=10,minbeamfrac=0.3
#testname: test_mask_autobox_multithresh_with_prune
#
#
#
#Show multiple executions of deconvolve get a probably-correct answer
#45. Deconvolve Hogbom + Hogbom: execute deconvolve once and compare the results to those of running deconvolve twice in a row (almost the same, as it is with tclean).
#imsize:100,cell:'8.0arcsec',deconvolver:'hogbom',threshold:'1mJy'
#first run to get expected value: niter=20
#second/third runs to get actual value: niter=10, niter=10
#testname: test_multirun_hogbomhogbom
#
#46. Deconvolve Clark + Clark: execute deconvolve once and compare the results to those of running deconvolve twice in a row.
#imsize:100,cell:'8.0arcsec',deconvolver:'clark',threshold:'1mJy'
#first run to get expected value: niter=20
#second/third runs to get actual value: niter=10, niter=10
#testname: test_multirun_clarkclark
#
#47. Deconvolve Clarkstokes + Clarkstokes: execute deconvolve once and compare the results to those of running deconvolve twice in a row.
#imsize:100,cell:'8.0arcsec',deconvolver:'clarkstokes','cyclefactor':0.1,'stokes':'I'
#first run to get expected value: niter=20
#second/third runs to get actual value: niter=10, niter=10
#testname: test_multirun_clarkstokesclarkstokes
#
#48. Deconvolve Multiscale + Multiscale: execute deconvolve once and compare the results to those of running deconvolve twice in a row.
#imsize:100,cell:'8.0arcsec',deconvolver:'multiscale',threshold:'1mJy',scales:[0,20,40,100]
#first run to get expected value: niter=20
#second/third runs to get actual value: niter=10, niter=10
#testname: test_multirun_multiscalemultiscale
#
#49. Deconvolve MTMFS + MTMFS: execute deconvolve once and compare the results to those of running deconvolve twice in a row.
#imsize:100,cell:'8.0arcsec',deconvolver:'mtmfs',threshold:'1mJy',scales:[0,20,40,100]
#first run to get expected value: niter=20
#second/third runs to get actual value: niter=10, niter=10
#testname: test_multirun_mtmfsmtmfs
#
#50. Deconvolve Multiscale + Hogbom: Tests the example use case of using hogbom to speed up cleaning after the benefits of multiscale have worn off.
#Note: we only test for task completion, don't know what the value should be at the end. (TODO needs validation)
#imsize:200, cell:'8.0arcsec', deconvolver:'multiscale', scales:[0,20,40,100],niter=10
#imsize:200, cell:'8.0arcsec', deconvolver:'hogbom', niter=10
#testname: test_multirun_multiscalehog
#
#51. Run then Restore: Tests the example use case of using task_deconvolve for just the restoration step.
#imsize:100, cell:['10.0arcsec','30.0arcsec'], restoration=False, niter=10
#imsize:100, cell:['10.0arcsec','30.0arcsec'], restoration=True, niter=0
#testname: test_multirun_norestore_restore
#
#
#
#'img' value input checking using the hogbom deconvolver
#52. Missing Table: tst.residual is always required
#rm 'tst.residual', niter=10
#testname: test_imgval_missingimgs_residual
#
#53. Missing Table: tst.psf is always required
#rm 'tst.psf', niter=10
#testname: test_imgval_missingimgs_psf
#
#54. Missing Table: tst.model is used to continue deconvolution, but is not required.
#rm 'tst.model', niter=10
#testname: test_imgval_missingimgs_model
#
#55. Missing Table: tst.sumwt is never required
#rm 'tst.sumwt', niter=10
#testname: test_imgval_missingimgs_sumwt
#
#56. Reorder Image Axes: tst.residual must have the axes as is given in tclean
#imtrans("tst.residual", order="3012"), niter=10
#testname: test_imgval_axesmismatch_residual
#
#57. Reorder Image Axes: tst.psf must have the axes as is given in tclean
#imtrans("tst.psf", order="3012"), niter=10
#testname: test_imgval_axesmismatch_psf
#
#58. Reorder Image Axes: tst.model must have the axes as is given in tclean
#imtrans("tst.model", order="3012"), niter=10
#testname: test_imgval_axesmismatch_model
#
#59. Reorder Image Axes: tst.pb must have the axes as is given in tclean
#imtrans("tst.pb", order="3012"), niter=10, usemask='pb', pbmask=0.2
#testname: test_imgval_axesmismatch_pb
#
#60. Rebin To Smaller Image: everything else must have the same shape as tst.residual
#imrebin("tst.residual", factor=[50,50]), niter=10
#testname: test_imgval_shapemismatch_residual
#
#61. Rebin To Smaller Image: tst.psf must have the same shape as tst.residual
#imrebin("tst.psf", factor=[50,50]), niter=10
#testname: test_imgval_shapemismatch_psf
#
#62. Rebin To Smaller Image: tst.model must have the same shape as tst.residual
#imrebin("tst.model", factor=[50,50]), niter=10
#testname: test_imgval_shapemismatch_model
#
#63. Rebin To Smaller Image: tst.pb must have the same shape as tst.residual
#imrebin("tst.pb", factor=[50,50]), niter=10, usemask='pb', pbmask=0.2
#testname: test_imgval_shapemismatch_pb
#
#64. Empty 'startmodel' Parameter String: Deconvolve should ignore all empty strings entered for the startmodel
#startmodel='', startmodel=['', '', ''], startmodel=['', '', self.mname2, '', '']
#testname: test_imgval_startmodel_empty
#
#65. Parameter 'startmodel' Does Not Exist: Throws an error if startmodel is set but does not exist
#startmodel='doesnotexists.model'
#testname: test_imgval_startmodel_dne
#
#66. Both 'startmodel' And tst.model Exist: Throws an error if startmodel is set and tst.model exists (must be one or the other, not both)
#startmodel='tst_2.model'
#testname: test_imgval_startmodel_model_exists
#
#67. Parameter 'startmodel' Set: Tests ability of deconvolve to copy startmodel to tst.model before starting deconvolution
#startmodel='tst_2.model'
#testname: test_imgval_startmodel_basic_copy
#
#68. Reordered Startmodel Axes: Try to deconvolve with mismatched axes between startmodel and psf (should fail).
#startmodel='tst_2.model', imtrans(order="3012")
#testname:test_imgval_startmodel_axesmismatch
#
#69. Coordinate System Change: Task deconvolve should regrid the csys of the startmodel to that of tst.residual
#set_crval0(51), niter=10
#testname:test_imgval_startmodel_csysmismatch
#
#70. Image Shape Change: Task deconvolve should regrid the shape of the startmodel to that of tst.residual
#imrebin(factor=[2,2]), niter=10
#testname:test_imgval_startmodel_shapemismatch
#
#
#
#'img' value input checking using the mtmfs deconvolver
#71. Missing Table: tst.residual.tt1 is always required
#rm 'tst.residual.tt1', niter=10
#testname: test_mtmfsimgval_missingimgs_residual
#
#72. Missing Table: tst.psf.tt1 is always required
#rm 'tst.psf.tt1', niter=10
#testname: test_mtmfsimgval_missingimgs_psf
#
#73. Missing Table: tst.model.tt1 is used to continue deconvolution, but is not required.
#rm 'tst.model.tt1', niter=10
#testname: test_mtmfsimgval_missingimgs_model
#
#74. Missing Table: tst.sumwt.tt1 is never required
#rm 'tst.sumwt.tt1', niter=10
#testname: test_mtmfsimgval_missingimgs_sumwt
#
#75. Reorder Image Axes: tst.residual.tt1 must have the axes as is given in tclean
#imtrans("tst.residual.tt1", order="3012"), niter=10
#testname: test_mtmfsimgval_axesmismatch_residual
#
#76. Reorder Image Axes: tst.psf.tt1 must have the axes as is given in tclean
#imtrans("tst.psf.tt1", order="3012"), niter=10
#testname: test_mtmfsimgval_axesmismatch_psf
#
#77. Reorder Image Axes: tst.model.tt1 must have the axes as is given in tclean
#imtrans("tst.model.tt1", order="3012"), niter=10
#testname: test_mtmfsimgval_axesmismatch_model
#
#78. Reorder Image Axes: tst.pb.tt1 must have the axes as is given in tclean
#imtrans("tst.pb.tt1", order="3012"), niter=10, usemask='pb', pbmask=0.2
#testname: test_mtmfsimgval_axesmismatch_pb
#
#79. Rebin To Smaller Image: everything else must have the same shape as tst.residual.tt1
#imrebin("tst.residual.tt1", factor=[2,2]), niter=10
#testname: test_mtmfsimgval_shapemismatch_residual
#
#80. Rebin To Smaller Image: tst.psf.tt1 must have the same shape as tst.residual.tt1
#imrebin("tst.psf.tt1", factor=[2,2]), niter=10
#testname: test_mtmfsimgval_shapemismatch_psf
#
#81. Rebin To Smaller Image: tst.model.tt1 must have the same shape as tst.residual.tt1
#imrebin("tst.model.tt1", factor=[2,2]), niter=10
#testname: test_mtmfsimgval_shapemismatch_model
#
#82. Rebin To Smaller Image: tst.pb.tt1 must have the same shape as tst.residual.tt1
#imrebin("tst.pb.tt1", factor=[2,2]), niter=10, usemask='pb', pbmask=0.2
#testname: test_mtmfsimgval_shapemismatch_pb
#
#83. Empty 'startmodel' Parameter String: Deconvolve should ignore all empty strings entered for the startmodel
#startmodel='', startmodel=['', '', ''], startmodel=['', '', self.mname2, '', '']
#testname: test_mtmfsimgval_startmodel_empty
#
#84. Parameter 'startmodel' Does Not Exist: Throws an error if startmodel is set but does not exist
#startmodel='doesnotexists.model'
#testname: test_mtmfsimgval_startmodel_dne
#
#85. Both 'startmodel' And tst.model.tt1 Exist: Throws an error if startmodel is set and tst.model.tt1 exists (must be one or the other, not both)
#startmodel='tst_2.model.tt1'
#testname: test_mtmfsimgval_startmodel_model_exists
#
#86. Parameter 'startmodel' Set: Tests ability of deconvolve to copy startmodel to tst.model.tt1 before starting deconvolution
#startmodel='tst_2.model.tt1'
#testname: test_mtmfsimgval_startmodel_basic_copy
#
#87. Reordered Startmodel Axes: Try to deconvolve with mismatched axes between startmodel and psf (should fail).
#startmodel='tst_2.model.tt1', imtrans(order="3012")
#testname:test_mtmfsimgval_startmodel_axesmismatch
#
#88. Coordinate System Change: Task deconvolve should regrid the csys of the startmodel to that of tst.residual.tt1
#set_crval0(51), niter=10
#testname:test_mtmfsimgval_startmodel_csysmismatch
#
#89. Image Shape Change: Task deconvolve should regrid the shape of the startmodel to that of tst.residual.tt1
#imrebin(factor=[2,2]), niter=10
#testname:test_mtmfsimgval_startmodel_shapemismatch
#
#
#
#Multiple deconvolves update the .residual
#90. Hogbom Updates Residual: Task deconvolve should update the .residual with every execution.
#This behavior is left up to each deconvolver. Test that hogbom does this correctly.
#'vis':'refim_eptwochan.ms', 'imsize':10, 'cell':'8.0arcsec', 'deconvolver':hogbom, 'niter':10
#testname: test_residual_update_hogbom
#
#91. Clark Updates Residual: Task deconvolve should update the .residual with every execution.
#This behavior is left up to each deconvolver. Test that clark does this correctly.
#'vis':'refim_eptwochan.ms', 'imsize':10, 'cell':'8.0arcsec', 'deconvolver':clark, 'niter':10
#testname: test_residual_update_clark
#
#92. Clarkstokes Updates Residual: Task deconvolve should update the .residual with every execution.
#This behavior is left up to each deconvolver. Test that clarkstokes does this correctly.
#'vis':'refim_point_linRL.ms', 'imsize':10, 'cell':'8.0arcsec', 'deconvolver':clark, 'niter':10, 'stokes':'I'
#testname: test_residual_update_clarkstokes
#
#93. Multiscale Updates Residual: Task deconvolve should update the .residual with every execution.
#This behavior is left up to each deconvolver. Test that hogbom does this correctly.
#'vis':'refim_eptwochan.ms', 'imsize':10, 'cell':'8.0arcsec', 'deconvolver':multiscale, 'niter':10
#testname: test_residual_update_multiscale
#
#94. MTMFS Updates Residual: Task deconvolve should update the .residual with every execution.
#This behavior is left up to each deconvolver. Test that hogbom does this correctly.
#'vis':'refim_eptwochan.ms', 'imsize':10, 'cell':'8.0arcsec', 'deconvolver':mtmfs, 'niter':10
#testname: test_residual_update_mtmfs
#
#
#
#Control .image image restoration
#Most deconvolvers do restoration the same with (mtmfs being the exception). These tests show that hogbom deconvolve restores correctly.
#95. Deconvolve but don't restore: should not create a .image image
#'vis':'refim_eptwochan.ms', 'imsize':10, 'cell':'8.0arcsec', 'niter':10, restoration=False
#testname: test_restoration_none
#
#96. Deconvolve and restore, and compare results with those from a tclean run
#'vis':'refim_eptwochan.ms', 'imsize':10, 'cell':'8.0arcsec', 'niter':10, restoration=True
#testname: test_restoration_basic(self):
#
#97. Deconvolve and don't restore, then restore and compare results with those from a tclean run
#run 1: 'vis':'refim_eptwochan.ms', 'imsize':10, 'cell':'8.0arcsec', 'niter':10, restoration=False
#run 2: 'niter':0, restoration=True
#testname: test_restoration_onlyrestore(self):
#
#98. Deconvolve and restore with a gigantic beam, then restore and compare results with those from a tclean run with a gigantic beam
#This test is here mainly to verify that "restoringbeam" is picked up by Deconvolve.
#'vis':'refim_eptwochan.ms', 'imsize':10, 'cell':'8.0arcsec', 'niter':10, restoringbeam='100.0arcsec'
#testname: test_restoration_bigbeam(self):
#
#
#
#N Iter Params tests: verify that we perform the same number of iterations as tclean for the same iteration parameters
#99. Deconvolve should execute 10 iterations for gain=0.2, just like the first major-minor cycle of tclean.
#'gain':0.2
#testname: test_niterparms_gain_1
#
#100. Deconvolve should execute 7 iterations for gain=0.3, just like the first major-minor cycle of tclean.
#'gain':0.3
#testname: test_niterparms_gain_2
#
#101. Deconvolve should execute 16 iterations for threshold=0.22, just like the first major-minor cycle of tclean.
#'threshold':0.22
#testname: test_niterparms_threshold_1
#
#102. Deconvolve should execute 19 iterations for threshold=0.18, just like the first major-minor cycle of tclean.
#'threshold':0.18
#testname: test_niterparms_threshold_2
#
#103. Deconvolve should execute 72 iterations for nsigma=0.9, just like the first major-minor cycle of tclean.
#'nsigma':0.9, 'maxpsffraction':0
#testname: test_niterparms_nsigma_1
#
#104. Deconvolve should execute 60 iterations for nsigma=1.5, just like the first major-minor cycle of tclean.
#'nsigma':1.5, 'maxpsffraction':0
#testname: test_niterparms_nsigma_2
#
#105. Deconvolve should execute 40 iterations for cyclefactor=0.1, just like the first major-minor cycle of tclean.
#'cyclefactor':0.1
#testname: test_niterparms_cyclefactor_1
#
#106. Deconvolve should execute 13 iterations for cyclefactor=2.0, just like the first major-minor cycle of tclean.
#'cyclefactor':2.0
#testname: test_niterparms_cyclefactor_2
#
#107. Deconvolve should execute 7 iterations for minpsffraction=0.5, just like the first major-minor cycle of tclean.
#'minpsffraction':0.5
#testname: test_niterparms_minpsffraction_1
#
#108. Deconvolve should execute 16 iterations for minpsffraction=0.2, just like the first major-minor cycle of tclean.
#'minpsffraction':0.2
#testname: test_niterparms_minpsffraction_2
#
#109. Deconvolve should execute 78 iterations for maxpsffraction=0.01, just like the first major-minor cycle of tclean.
#'maxpsffraction':0.01
#testname: test_niterparms_maxpsffraction_1
#
#110. Deconvolve should execute 40 iterations for maxpsffraction=0.05, just like the first major-minor cycle of tclean.
#'maxpsffraction':0.05
#testname: test_niterparms_maxpsffraction_2
#
##########################################################################

from __future__ import absolute_import
import os
import sys
import shutil
import unittest
import inspect
import numpy as np
import re
import glob
from casatestutils.imagerhelpers import TestHelpers

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import ctsys, quanta, measures, image, vpmanager, calibrater
    from casatasks import casalog, deconvolve, tclean, imtrans, imrebin, imregrid, imval
    from casatasks.private.parallel.parallel_task_helper import ParallelTaskHelper
    from casatasks.private.imagerhelpers.parallel_imager_helper import PyParallelImagerHelper

    _ia = image( )
    _vp = vpmanager( )
    _cb = calibrater( )
    _qa = quanta( )
    _me = measures( )

    refdatapath = ctsys.resolve('unittest/deconvolve/')
else:
    from __main__ import default
    from tasks import *
    from taskinit import *
    from parallel.parallel_task_helper import ParallelTaskHelper
    from imagerhelpers.parallel_imager_helper import PyParallelImagerHelper

    _ia = iatool( )
    _vp = vptool( )
    _cb = cbtool( )
    # not local tools
    _qa = qa
    _me = me

    refdatapath = os.environ.get('CASAPATH').split()[0] + '/casatestdata/unittest/deconvolve/'

th = TestHelpers()

## List to be run
def suite():
    return [test_onefield, test_iterbot, test_multifield, test_stokes, test_cube, test_mask, test_multirun, test_imgval, test_mtmfsimgval, test_residual_update, test_restoration, test_niterparms]

## Base Test class with Utility functions
class testref_base(unittest.TestCase):
    inptbls = [".residual", ".psf", ".pb", ".model"]
    img = "tst"

    @classmethod
    def setUpClass(cls):
        casalog.setlogfile('testlog.log')
        cls.msfile = ""           # the filename of the referenced ".ms" measurement set, with the ".ms" extension (note this variable is shadowed by and instance variable by the same name)
        cls.cachedir = ""         # name of the local directory used in staticCopyToCache() and staticCopyFromCache()

    def setUp(self):
        self.msfile = ""          # the filename of the referenced ".ms" measurement set, with the ".ms" extension (note this variable shadows the class variable by the same name)
        self.imgsrc = ""          # the filename of the referenced ".ms" measurement set, without the ".ms" extension
        self.img = type(self).img # prefix of the image names to generate
        self.parallel = False     # whether to run tclean in parallel
        if ParallelTaskHelper.isMPIEnabled():
            self.parallel = True

    def tearDown(self):
        # Default: delete all (input and output data)
        self.delData()
        # self.delData(delinput=False, deloutput=False) # Debugging: leave input and output data

    @classmethod
    def tearDownClass(cls):
        cls.staticDelData(cls.msfile, cls.img)
        # cls.staticDelData(cls.msfile, cls.img, delinput=False, deloutput=False) # Debugging: leave input and output data
        if cls.cachedir != "":
            os.system('rm -rf '+cls.cachedir)

    def _setImgsrcField(self,msfile):
        """ Set the self.msfile and self.imgsrc values based off of the given msfile """
        imgsrc=msfile
        if msfile.endswith(".ms"):
            self.msfile = msfile
            imgsrc = imgsrc[:-3]
        if imgsrc != "":
            self.imgsrc = imgsrc

    def copyMS(self,msfile):
        """ Copies the given measurement set file to the working directory """
        self._setImgsrcField(msfile)
        type(self).staticCopyMS(self.msfile)

    @classmethod
    def staticCopyMS(cls,msfile):
        """ Copies the given measurement set file to the working directory """
        if (cls.msfile == ""):
            cls.msfile = msfile
        shutil.copytree(os.path.join(refdatapath,msfile), msfile)

    @classmethod
    def staticCopyToCache(cls, msfile, imagename, cachedir):
        """ Copy all "self.img".* generated image files (eg tst.model/tst.model.tt0) from the working directory to a local directory "cachedir".
            Also creates the cachedir directory as necessary. """
        cls.msfile = msfile
        cls.img = imagename
        cls.cachedir = cachedir
        for tbl in cls.inptbls:
            for ttn in ["", ".tt0", ".tt1", ".tt2", ".tt3", ".tt4", ".tt5"]:
                ext = tbl+ttn
                srctbl = imagename + ext
                dsttbl = os.path.join(cachedir, msfile + ext)
                if os.path.exists(srctbl):
                    shutil.copytree(srctbl, dsttbl)

    @classmethod
    def staticCopyFromCache(cls, msfile="", imagename="", cachedir=""):
        """ Copy all "self.img".* generated image files (eg tst.model/tst.model.tt0) from the local directory "cachedir" to the working directory.
            If staticCopyToCache was previously used by this class (does not need to be the same instance of the class), then the msfile, imagename, and cachedir
            from that staticCopyToCache call will be used. """
        msfile    = msfile    if msfile!=""    else cls.msfile
        imagename = imagename if imagename!="" else cls.img
        cachedir  = cachedir  if cachedir!=""  else cls.cachedir
        for tbl in cls.inptbls:
            for ttn in ["", ".tt0", ".tt1", ".tt2", ".tt3", ".tt4", ".tt5"]:
                ext = tbl+ttn
                if os.path.exists(os.path.join(cachedir, msfile+ext)):
                    shutil.copytree(os.path.join(cachedir, msfile+ext), imagename+ext)

    @classmethod
    def staticPrepData(cls,msfile="",tclean_args={},parallel=False):
        """
        Copies the msfile ".ms" file to the local directory and
        prepares the tclean(niter=0) tst.* tables for use with deconvolve.

        1 supplement the given tclean_args with some default values, if those values aren't present
        2 copy the msfile ms file to the working directory
        3 execute tclean with the given tclean_args

        @param msfile A measurement set ".ms" file (include ".ms" extension, don't include directory)
        @param tclean_args Arguments to use when running tclean, for generating the input images to deconvolve.
        @return the summation of the given tclean_args and the default tclean_args.
        """
        # (1) set some default values for tclean
        defs = { 'vis'          : msfile,
                 'imagename'    : cls.img,
                 'niter'        : 0,
                 'restoration'  : False,
                 'calcres'      : True,
                 'pbcor'        : True,
                 'parallel'     : parallel }
        for k in defs.keys():
            tclean_args[k] = defs[k] if k not in tclean_args else tclean_args[k]
        
        # (2) copy msfile
        cls.staticCopyMS(msfile)

        # (3) create the residuals and other feeder tables
        tclean(**tclean_args)

        return tclean_args

    # Separate functions here, for special-case tests that need their own MS.
    def prepData(self,msfile="",tclean_args={},delold=True):
        """
        Deletes old data and uses staticPrepData to tclean new input images.

        1 update self.imgsrc and self.msfile with the value from msfile
        2 delete old data
        3 copy the ms and evaluate tclean using staticPrepData

        @param msfile A measurement set ".ms" file (include ".ms" extension, don't include directory)
        @param tclean_args Arguments to use when running tclean, for generating the input images to deconvolve.
        @param delold True to delete all data with self.delData()
                      False is mostly useful when creating mask files, so that the newly created mask isn't removed
        @return the summation of the given tclean_args and the default tclean_args.
        """
        # (1) update self.imgsrc and self.msfile
        self._setImgsrcField(msfile)

        # (2) delete old data
        if delold:
            self.delData()

        # print some debugging information
        print(refdatapath)
        print(self.imgsrc)

        # (3) evaluted with class method
        tclean_args = type(self).staticPrepData(self.imgsrc+'.ms', tclean_args, self.parallel)

        return tclean_args

    def copy_products(self, old_pname, new_pname, ignore=None):
        """ function to copy iter0 images to iter1 images
            (taken from test_stk_alma_pipeline_imaging.py, which in turn was taken from the pipeline)
        """
        imlist = glob.glob('%s.*' % old_pname)
        imlist = [xx for xx in imlist if ignore is None or ignore not in xx]
        for image_name in imlist:
            newname = image_name.replace(old_pname, new_pname)
            if image_name == old_pname + '.workdirectory':
                mkcmd = 'mkdir '+ newname
                os.system(mkcmd)
                self.copy_products(os.path.join(image_name, old_pname), \
                    os.path.join(newname, new_pname))
            else:
                shutil.copytree(image_name, newname, symlinks=True)

    @classmethod
    def staticDelData(cls,msfile="",imagename="",delinput=True,deloutput=True):
        """ Delete the .ms, imagename.*, and usermask.mask.txt files/directories from the working directory.
            If previously set from another static method, msfile and imagename will be provided default values. """
        if (os.path.exists(msfile) and delinput):
            os.system('rm -rf ' + msfile)
        if (os.path.exists('usermask.mask.txt') and delinput):
            os.system('rm -rf usermask.mask.txt')
        if deloutput and imagename != "":
            os.system('rm -rf ' + imagename+'*')

    def delData(self,msfile="",delinput=True,deloutput=True):
        """ Calls staticDelData with the instance values for self.img, and a default self.msfile is msfile isn't provided. """
        if msfile != "":
            self.msfile = msfile
        elif msfile == "":
            msfile = self.msfile
        type(self).staticDelData(msfile, self.img, delinput, deloutput)

    def checkfinal(self,pstr=""):
        """ Fails the test via self.fail() if either of the strings '(Fail' or '( Fail' are present in the given pstr. """
        # pstr += "["+inspect.stack()[1][3]+"] : To re-run this test :  runUnitTest.main(['test_req_task_deconvolve["+ inspect.stack()[1][3] +"]'])"
        casalog.post(pstr,'INFO')
        if ( re.search('\( ?Fail', pstr) != None ):
            self.fail("\n"+pstr)

## Python27 backports
if not hasattr(testref_base, 'assertRaisesRegex'):
    testref_base.assertRaisesRegex = testref_base.assertRaisesRegexp

##############################################
##############################################

##Task level tests : one field, 2chan.
class test_onefield(testref_base):

    # Test 1
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_onefield_hogbom(self):
        """ [onefield] test_onefield_hogbom """
        ######################################################################################
        # Test mfs with hogbom minor cycle. Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_twochan.ms', tclean_args={'imsize':100, 'cell':'8.0arcsec', 'deconvolver':'hogbom'})
        results = deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0)
        report=th.checkall(ret=results['retrec'], peakres=0.353, modflux=0.772, iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image',self.img+'.model'],
                           imgval=[(self.img+'.model',0.772,[50,50,0,0])])
        self.checkfinal(pstr=report)

    # Test 2
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_onefield_clark(self):
        """ [onefield] test_onefield_clark """
        ######################################################################################
        # Test mfs with clark minor cycle. Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_twochan.ms', tclean_args={'imsize':100, 'cell':'8.0arcsec', 'deconvolver':'clark'})
        results = deconvolve(imagename=self.img, niter=10, deconvolver='clark', interactive=0)
        report=th.checkall(ret=results['retrec'], peakres=0.392, modflux=0.733, iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image',self.img+'.model'],
                           imgval=[(self.img+'.model',0.733,[50,50,0,0])])
        self.checkfinal(pstr=report)

    # Test 3
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_onefield_multiscale(self):
        """ [onefield] test_onefield_multiscale """
        ######################################################################################
        # Test mfs with multiscale minor cycle. Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_eptwochan.ms', tclean_args={'imsize':200, 'cell':'8.0arcsec', 'deconvolver':'multiscale', 'scales':[0,20,40,100]})
        results = deconvolve(imagename=self.img, niter=10, deconvolver='multiscale', scales=[0,20,40,100], interactive=0)
        report=th.checkall(ret=results['retrec'], peakres=0.823, modflux=3.816, iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image',self.img+'.model'],
                           imgval=[(self.img+'.model',0.234,[94,107,0,0])])
        self.checkfinal(pstr=report)

    # Test 4
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_onefield_mtmfs(self):
        """ [onefield] test_onefield_mtmfs """
        ######################################################################################
        # Test mt-mfs with minor cycle iterations . Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_twochan.ms', tclean_args={'imsize':100, 'cell':'8.0arcsec', 'deconvolver':'mtmfs'})
        results = deconvolve(imagename=self.img, niter=10, deconvolver='mtmfs', interactive=0)
        report=th.checkall(ret=results['retrec'], peakres=0.392, modflux=0.732, iterdone=10, imgexist=[self.img+'.psf.tt0', self.img+'.residual.tt0', self.img+'.image.tt0', self.img+'.model.tt0',self.img+'.model.tt1',self.img+'.alpha'],
                           imgval=[(self.img+'.model.tt0',0.733,[50,50,0,0]),(self.img+'.image.tt1',0.019,[2,94,0,0])])
        ## iterdone=11 only because of the return (iterdone_p+1) in MultiTermMatrixCleaner::mtclean() !
        self.checkfinal(pstr=report)

    # Test 5
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_onefield_rectangular_pixels(self):
        """ [onefield] test_onefield_rectangular_pixels """
        ######################################################################################
        # Test restoration with rectangular pixels (cas-7171). Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_point.ms', tclean_args={'imsize':100, 'cell':['10.0arcsec','30.0arcsec']})
        deconvolve(imagename=self.img, niter=10)
        report=th.checkall(imgexist=[self.img+'.image'], imgval=[(self.img+'.image',0.482,[50,49,0,0])] )
        self.checkfinal(report)

##############################################
##############################################

##Task level tests : iteration controls
class test_iterbot(testref_base):

    # Test 6
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_iterbot_mfs_4(self):
        """ [iterbot] test_iterbot_mfs_4 """
        ######################################################################################
        # Test Iterations with high gain . Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_twochan.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'clark','gain':0.15})
        results = deconvolve(imagename=self.img, deconvolver='clark', niter=20, gain=0.15,interactive=0)
        report=th.checkall(ret=results['retrec'], iterdone=14, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'],
                           imgval=[(self.img+'.model',0.937,[50,50,0,0])])

        self.checkfinal(report)

    # Test 7
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_iterbot_mfs_5(self):
        """ [iterbot] test_iterbot_mfs_5 """
        ######################################################################################
        # Threshold test . Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_twochan.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'clark','threshold':'0.5Jy','gain':0.15})
        results = deconvolve(imagename=self.img, deconvolver='clark', niter=10, threshold='0.5Jy', gain=0.15, interactive=0)
        report=th.checkall(ret=results['retrec'], peakres=0.499, modflux=0.626, iterdone=5, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'],
                           imgval=[(self.img+'.model',0.626,[50,50,0,0])])

        self.checkfinal(report)

##############################################
##############################################

##Task level tests : multi-field, 2chan.
class test_multifield(testref_base):
    # n/a, outlier tests are for running a combined major cycle on many images
    # just run deconvolve on each of the individual images
    pass

##############################################
##############################################

##Task level tests : Stokes imaging options
class test_stokes(testref_base):

    # Test 8
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_stokes_mfs_I(self):
        """ [stokes] test_stokes_mfs_I """
        ######################################################################################
        # Test_Stokes_I_mfs mfs with stokes I. Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_point_linRL.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','stokes':'I'})
        deconvolve(imagename=self.img, niter=10)
        report=th.checkall(imgexist=[self.img+'.image'],
                           imgval=[(self.img+'.image',1.0,[50,50,0,0])])
        self.checkfinal(report)

    # Test 9
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_stokes_mtmfs_IQUV(self):
        """ [stokes] test_stokes_mtmfs_IQUV """
        ######################################################################################
        # Test mtmfs with stokes IQUV. Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_point_linRL.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','stokes':'IQUV','deconvolver':'mtmfs','nterms':2})
        deconvolve(imagename=self.img, niter=10, deconvolver='mtmfs', nterms=2)
        report=th.checkall(imgexist=[self.img+'.image.tt0'],imgexistnot=[self.img+'.image.alpha'],
                           imgval=[(self.img+'.image.tt0',1.0,[50,50,0,0]),(self.img+'.image.tt0',2.0,[50,50,1,0]), (self.img+'.image.tt0',3.0,[50,50,2,0]),(self.img+'.image.tt0',4.0,[50,50,3,0]) ])
        self.checkfinal(report)

##############################################
##############################################

##Task level tests : cube.
class test_cube(testref_base):

    def setUp(self):
        super(test_cube, self).setUp()

        ## Setup some variables to use in all the tests

        ## chan 4 (TOPO)
        qfstart=_qa.quantity("1.2GHz")
        #qvstart=_qa.quantity("-59958.5km/s")
        # for restf=1.25GHz
        qvstart=_qa.quantity("11991.7km/s")
        # ch10
        #qvstart=_qa.quantity("16788.4km/s")

        #mfstart=_me.frequency('LSRK',_qa.quantity("1.09999GHz"))
        # ch4 (for rest 1.25GHz)
        mfstart=_me.frequency('LSRK',_qa.quantity("1.199989GHz"))
        mvstart=_me.radialvelocity('BARY',_qa.quantity("11977.6km/s"))
        #dop = _me.todoppler('radio',mfstart,_qa.quantity('1.0GHz'))
        mfstart10=_me.frequency('LSRK',_qa.quantity(" 1.17999GHz"))                                                        
        # doppler with ch4 freq
        dop = _me.todoppler('radio',mfstart,_qa.quantity('1.25GHz'))                                              

        #1chan width 
        #qvwidth = _qa.quantity("11991.700km/s")
        #qvwidth = _qa.quantity("4796.7km/s")
        qvwidth = _qa.quantity("11991.7km/s")
        mvwidth = _me.radialvelocity('TOPO',qvwidth)

        # restf = 1.25GHz
        # vel range: 59961.1 -  -31174.7 km/s (lsrk/radio)
        #            74952.3 -  -28238.3 km/s (lsrk/optical)  

        self.testList = {
                    0:{'imagename':'Cubetest_chandefstdefwidth','spw':'0','start':0,'width':1,'outframe':'LSRK','veltype':'radio',
                      'desc':'channel, default start and width, LSRK'},
                    1:{'imagename':'Cubetest_chandefstdefwidthtopo','spw':'0','start':0,'width':1, 'outframe':'TOPO','veltype':'radio',
                      'desc':'channel, default start and width, TOPO'},
                    2:{'imagename':'Cubetest_chandefstwidth2','spw':'0','start':0,'width':2, 'outframe':'LSRK','veltype':'radio',
                      'desc':'channel, default start, width=2, LSRK'},
                    3:{'imagename':'Cubetest_chanst5wd1','spw':'0','start':5,'width':1, 'outframe':'LSRK','veltype':'radio',
                      'desc':'channel, start=5, default width, LSRK'},
                    # this will result in blank channnel images (calcChanFreqs requires start and width in channel       
                    # mode to be given in chan index                                                                 
                    4:{'imagename':'Cubetest_chandefstwd1spwsel','spw':'0:5~19','start':0,'width':1, 'outframe':'LSRK','veltype':'radio',
                      'desc':'channel, spw=0:5~19, LSRK'},
                    #5:{'imagename':'Cubetest_freqdefstwd2','spw':'0','start':'','width':'40MHz','outframe':'TOPO',
                    #  'desc':'frequency, default start, width=\'40MHz\', TOPO'},
                    # data set changed!
                    5:{'imagename':'Cubetest_freqdefstwd2','spw':'0','start':'','width':'100MHz','outframe':'TOPO','veltype':'radio',
                      'desc':'frequency, default start, width=\'100MHz\'(2 x chanwidth), TOPO'},
                    6:{'imagename':'Cubetest_freqst2defwd','spw':'0','start':'1.1GHz','width':'','outframe':'TOPO','veltype':'radio',
                      'desc':'frequency, start=\'1.1GHz\', default width, TOPO'},
                    7:{'imagename':'Cubetest_freqst2defwdspwsel','spw':'0:4~19','start':'1.1GHz','width':'','outframe':'TOPO','veltype':'radio',
                      'desc':'frequency, start=\'1.1GHz\', default width, spw=0:4~19, TOPO'},
                    8:{'imagename':'Cubetest_freqst10wdm','spw':'0','start':'1.5GHz','width':'-50MHz','outframe':'TOPO','veltype':'radio',
                      'desc':'frequency, start=\'1.5GHz\', width=\'-50MHz\', TOPO'},
                    9:{'imagename':'Cubetest_veldefstwd2','spw':'0','start':'','width':'23983.4km/s','outframe':'TOPO','veltype':'radio',
                      'desc':'frequency, default start, width=\'23983.4km/s\', TOPO'},
                   10:{'imagename':'Cubetest_veldefstwd2m','spw':'0','start':'','width':'-23983.4km/s','outframe':'TOPO','veltype':'radio',
                      'desc':'velocity, default start, width=\'-23983.4m/s\', TOPO'},
                   11:{'imagename':'Cubetest_velst4defwd','spw':'0','start':'11991.7km/s','width':'','outframe':'TOPO','veltype':'radio',
                      'desc':'velocity, start=\'11991.7km/s\', default width, TOPO'},
                   12:{'imagename':'Cubetest_velst4defwdbary','spw':'0','start':'11977.6km/s','width':'','outframe':'BARY','veltype':'radio',
                      'desc':'velocity, start=\'11977.6km/s\', default width, BARY'},
                   # currently 13 is not quite properly working, investigating - 2014.08.27 TT 
                   # for refim_point.ms ch9=-41347.8km/s (opt)
                   #13:{'imagename':'Cubetest_optvelst10wdeflsrk','spw':'0','start':'-49962.6km/s','width':'',
                   13:{'imagename':'Cubetest_optvelst19wdlsrk','spw':'0','start':'-41347.8km/s','width':'20000km/s',
                      'veltype':'optical','outframe':'LSRK',
                   ##   'desc':'velocity, start=\'74952.3km/s\', default width, veltype=optical LSRK'},
                   #   'desc':'velocity, start=\'-49962.6km/s\', default width, veltype=optical LSRK'},
                      'desc':'velocity, start=\'-41347.5km/s\', default width , veltype=optical LSRK'},
                   14:{'imagename':'Cubetest_stqfreqdefwd','spw':'0','start':qfstart,'width':'', 'veltype':'radio','outframe':'',
                      'desc':'frequency, start(quanity)=%s, default width, veltype=radio LSRK' % qfstart},
                   15:{'imagename':'Cubetest_stmfreqdefwd','spw':'0','start':mfstart,'width':'', 'veltype':'radio','outframe':'',
                      'desc':'frequency, start=%s, default width, veltype=radio LSRK' % mfstart},
                   16:{'imagename':'Cubetest_stqveldefwd','spw':'0','start':qvstart,'width':'','outframe':'TOPO','veltype':'radio',
                      'desc':'velocity(quantity), start=%s, default width, TOPO ' % qvstart},
                   17:{'imagename':'Cubetest_stmveldefwd','spw':'0','start':mvstart,'width':'','outframe':'TOPO','veltype':'radio',
                      'desc':'velocity(measure), start=%s, default width(outframe=TOPO will be overridden)' % mvstart},
                   18:{'imagename':'Cubetest_veldefstqvwidth','spw':'0','start':'','width':qvwidth,'outframe':'TOPO','veltype':'radio',
                      'desc':'velocity, default start, width(quantity)=%s' % qvwidth},
                   19:{'imagename':'Cubetest_veldefstmvwidth','spw':'0','start':'','width':mvwidth,'outframe':'TOPO','veltype':'radio',
                      'desc':'velocity, default start, width(measure)=%s, TOPO' % mvwidth},
                   20:{'imagename':'Cubetest_stdopdefwd','spw':'0','start':dop,'width':'','outframe':'LSRK','veltype':'radio',
                      'desc':'doppler, start=%s, default width, LSRK' % dop},
                   # with a gap in spw channel sel
                   21:{'imagename':'Cubetest_st4gap','spw':'0:4~9;12~14','start':4,'width':'','outframe':'LSRK','veltype':'radio',
                      'desc':'channel, start=%s, default width, channel gap (10-11) LSRK' % 4},
                   # stride > 1
                   22:{'imagename':'Cubetest_st4stride2','spw':'0:0~10^2','start':0,'width':'','outframe':'LSRK','veltype':'radio', 'interpolation':'nearest',
                      'desc':'channel, start=%s, default width, step=2 LSRK nearest' % 0},
                   23:{'imagename':'Cubetest_defstspwchansel4','spw':'0:4~13','start':'','width':'','outframe':'TOPO','veltype':'radio',
                      'desc':'spw with channel selection( 0:4~13 ), default start, LSRK nearest'}
                  }
        
        # self.test_cube_0.__func__.__doc__ %=self.testList[0]['desc']
    
    def get_cubetclean_args(self, testid):
        """ core function to execute a cube tclean """
        if 'interpolation' in self.testList[testid]:
            interpolation = self.testList[testid]['interpolation']
        else:
            interpolation = 'linear'

        tclean_args = { 'field':'0', 'imsize':100, 'cell':'8.0arcsec',
                        'specmode':'cube', 'nchan':10, 'restfreq':['1.25GHz'],
                        'phasecenter':"J2000 19:59:28.500 +40.44.01.50", 'deconvolver':'hogbom',
                        'spw':self.testList[testid]['spw'],
                        #'imagename':self.img+self.testList[testid]['imagename'],
                        'start':self.testList[testid]['start'],
                        'width':self.testList[testid]['width'],
                        'veltype':self.testList[testid]['veltype'],
                        'outframe':self.testList[testid]['outframe'],
                        'interpolation':interpolation }
        return tclean_args

    # Test 10
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_0(self):
        """ [cube] test_cube_0 """
        ######################################################################################
        # Test_Cube_0 new . Should produce the same results as tclean.
        ######################################################################################
        testid=0
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.50002,[50,50,0,0])])
        report2 = th.check_spec_frame(self.img+'.image','LSRK',999988750)
        self.checkfinal(report+report2)

    # Test 11
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_1(self):
        """ [cube] test_cube_1 """
        ######################################################################################
        # test_cube_1. Should produce the same results as tclean.
        ######################################################################################
        testid=1
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.50002,[50,50,0,0])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO', 9.9999999e8)

        self.checkfinal(report+report2)

    # Test 12
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_2(self):
        """ [cube] test_cube_2 """
        ######################################################################################
        # test_cube_2. Should produce the same results as tclean.
        ######################################################################################
        testid=2
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.4643,[50,50,0,0])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.02498846e9)
        self.checkfinal(report+report2)

    # Test 13
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_3(self):
        """ [cube] test_cube_3 """
        ######################################################################################
        # test_cube_3. Should produce the same results as tclean.
        ######################################################################################
        # start = 5 (1.25GHZ IN TOPO)
        testid=3
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.2000,[50,50,0,0])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.249985937e9)
        self.checkfinal(report+report2)

    # Test 14
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_4(self):
        """ [cube] test_cube_4 """
        ######################################################################################
        # test_cube_4. Should produce the same results as tclean.
        ######################################################################################
        testid=4
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.5000,[50,50,0,0])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.23998593e9)

    # Test 15
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_5(self):
        """ [cube] test_cube_5 """
        ######################################################################################
        # test_cube_5. Should produce the same results as tclean.
        ######################################################################################
        # width by freq (2x chanw) result should be the same as #2
        testid=5
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.4643,[50,50,0,0])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.025e9)
        self.checkfinal(report+report2)

    # Test 16
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_6(self):
        """ [cube] test_cube_6 """
        ######################################################################################
        # test_cube_6. Should produce the same results as tclean.
        ######################################################################################
        # start in freq=1.1GHz (=chan5)
        testid=6
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.36365354,[50,50,0,0])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.1e9)
        self.checkfinal(report+report2)

    # Test 17
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_7(self):
        """ [cube] test_cube_7 """
        ######################################################################################
        # test_cube_7. Should produce the same results as tclean.
        ######################################################################################
        # start 1.1GHz(TOPO)=chan2 spw=4~19
        # Currently different behaviors between serial and parallel on non-overlapping data and image
        # parameter selections.
        # serial: result in chan 0&1 psf blanked  
        # parallel: spw channel selection will be ignored and tuneselectdata will 
        # select overlapping data and image selections (this seems to me more correct? behavior)
        # as of 2019.01.08, this is no longer true, psf blanked for chan 0 and 1 for parallel
        testid=7
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',0.0,[50,50,0,0]),
                                     (self.img+'.image',1.2000,[50,50,0,3])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.1e9)
        self.checkfinal(report+report2)

    # Test 18
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_8(self):
        """ [cube] test_cube_8 """
        ######################################################################################
        # test_cube_8. Should produce the same results as tclean.
        ######################################################################################
        # start =1.5GHz(chan10)  width=-50MHz TOPO (descending freq)
        testid=8
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.42858946,[50,50,0,9])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.5e9)
        self.checkfinal(report+report2)

    # Test 19
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_9(self):
        """ [cube] test_cube_9 """
        ######################################################################################
        # test_cube_9. Should produce the same results as tclean.
        ######################################################################################
        # width in vel (=23983.4km/s=2xChanW) def start (=cube will be ascending order in vel)
        testid=9
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.46184647,[50,50,0,9])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.925e9)
        self.checkfinal(report+report2)

    # Test 20
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_10(self):
        """ [cube] test_cube_10 """
        ######################################################################################
        # test_cube_10. Should produce the same results as tclean.
        ######################################################################################
        # width in vel = -23983.4m/s def start (cube will be in descending order in vel)
        testid=10
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.46184647,[50,50,0,0])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.025e9)
        self.checkfinal(report+report2)

    # Test 21
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_11(self):
        """ [cube] test_cube_11 """
        ######################################################################################
        # test_cube_11. Should produce the same results as tclean.
        ######################################################################################
        # start 11991.7km/s (chan4)
        testid=11
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.50001776,[50,50,0,4])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.2e9)
        self.checkfinal(report+report2)

    # Test 22
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_12(self):
        """ [cube] test_cube_12 """
        ######################################################################################
        # test_cube_12. Should produce the same results as tclean.
        ######################################################################################
        # start 11977.6km/s (BARY) = chan4
        testid=12
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.50001931,[50,50,0,4])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','BARY',1.200058783e9)
        self.checkfinal(report+report2)

    # Test 23
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_13(self):
        """ [cube] test_cube_13 """
        ######################################################################################
        # test_cube_13. Should produce the same results as tclean.
        ######################################################################################
        # 
        testid=13
        print(" : " , self.testList[testid]['desc'])
        # use own tclean command as nchan need to modify
        t = {"vis":"refim_point.ms", "field":'0', "imsize":100, "cell":'8.0arcsec', "niter":0, "specmode":'cube', "nchan":8, "restfreq":['1.25GHz'], 
             "phasecenter":"J2000 19:59:28.500 +40.44.01.50", "deconvolver":'hogbom', "spw":self.testList[testid]['spw'], 
             "imagename":self.img+self.testList[testid]['imagename'], "start":self.testList[testid]['start'], 
             "width":self.testList[testid]['width'], "veltype":self.testList[testid]['veltype'], 
             "outframe":self.testList[testid]['outframe'], "parallel":self.parallel}
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img+self.testList[testid]['imagename'], niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+self.testList[testid]['imagename']+'.psf') and os.path.exists(self.img+self.testList[testid]['imagename']+'.residual') )
        #report=self.th.checkall(imgexist=[self.img+self.testList[testid]['imagename']+'.image'],
        #imgval=[(self.img+self.testList[testid]['imagename']+'.image',1.50001931,
        #[50,50,0,4])])
        # report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.2000e9)

    # Test 24
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_14(self):
        """ [cube] test_cube_14 """
        ######################################################################################
        # test_cube_14. Should produce the same results as tclean.
        ######################################################################################
        # start = quantity ('1.2GHz') frame default(LSRK)
        testid=14
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.25000215,[50,50,0,0])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.2e9)
        self.checkfinal(report+report2)

    # Test 25
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_15(self):
        """ [cube] test_cube_15 """
        ######################################################################################
        # test_cube_15. Should produce the same results as tclean.
        ######################################################################################
        # measure freq in LSRK ch4
        testid=15
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.25001216,[50,50,0,0])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.199989e9)
        self.checkfinal(report+report2)

    # Test 26
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_16(self):
        """ [cube] test_cube_16 """
        ######################################################################################
        # test_cube_16. Should produce the same results as tclean.
        ######################################################################################
        # start quantity vel=11991.7km/s outframe=topo (ascending vel order)
        testid=16
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.50001776,[50,50,0,4])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.2000e9)
        self.checkfinal(report+report2)

    # Test 27
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_17(self):
        """ [cube] test_cube_17 """
        ######################################################################################
        # test_cube_17. Should produce the same results as tclean.
        ######################################################################################
        # start measure vel=11977.6km/s BARY, outframe=TOPO will be overridedden (ascending vel order)
        testid=17
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.50001931,[50,50,0,4])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','BARY',1.200058783e9)
        self.checkfinal(report+report2)

    # Test 28
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_18(self):
        """ [cube] test_cube_18 """
        ######################################################################################
        # test_cube_18. Should produce the same results as tclean.
        ######################################################################################
        # defaut start, width in vel (quantity) +11991.7km/s (TOPO, radio)=datachan width, will be
        # ascending order in vel so highet DATA channel will be chan 0 in the image (image chan0=1.45GHz)
        testid=18
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.50001764,[50,50,0,9])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.45e9)
        self.checkfinal(report+report2)

    # Test 29
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_19(self):
        """ [cube] test_cube_19 """
        ######################################################################################
        # test_cube_19. Should produce the same results as tclean.
        ######################################################################################
        # default start, width in vel (measure) +11991.7km/s (TOPO, radio)
        testid=19
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.50001764,[50,50,0,9])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.45e9)
        self.checkfinal(report+report2)

    # Test 30
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_20(self):
        """ [cube] test_cube_20 """
        ######################################################################################
        # test_cube_20. Should produce the same results as tclean.
        ######################################################################################
        # doppler (with ch4 LSRK freq, rest freq=1.25GHz)
        testid=20
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.5000546,[50,50,0,4])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.199989152e9)
        self.checkfinal(report+report2)

    # Test 31
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_21(self):
        """ [cube] test_cube_21 """
        ######################################################################################
        # test_cube_21. Should produce the same results as tclean.
        ######################################################################################
        # data sel with channel gap (10,11 excluded) 4~9, 12~14
        testid=21
        self.testList[testid]['interpolation']='nearest'
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.250001562,[50,50,0,0]),
                                     (self.img+'.image',0.0,[50,50,0,6]),
                                     (self.img+'.image',0.0,[50,50,0,7])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',1.199986500e9)
        self.checkfinal(report+report2)

    # Test 32
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_22(self):
        """ [cube] test_cube_22 """
        ######################################################################################
        # test_cube_22. Should produce the same results as tclean.
        ######################################################################################
        # stride (step=2) use nearest interpolation (other interpotion methods
        # may not work well...)
        testid=22
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.5000546,[50,50,0,0])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','LSRK',0.999988750387e9)
        self.checkfinal(report+report2)

    # Test 33
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    # @unittest.skipIf(not hasattr(th, 'checkspecframe'), "Skip this test if checkspecframe hasn't been carried over to the new testing scripts yet")
    def test_cube_23(self):
        """ [cube] test_cube_23 """
        ######################################################################################
        # test_cube_23. Should produce the same results as tclean.
        ######################################################################################
        testid=23
        print(" : " , self.testList[testid]['desc'])
        t = self.get_cubetclean_args(testid)
        self.prepData('refim_point.ms', tclean_args=t)
        deconvolve(imagename=self.img, niter=10, deconvolver=t['deconvolver'])

        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.residual') )
        report = th.checkall(imgexist=[self.img+'.image'],
                             imgval=[(self.img+'.image',1.2500156,[50,50,0,0])])
        report2 = th.check_spec_frame(self.img+self.testList[testid]['imagename']+'.image','TOPO',1.20e9)
        self.checkfinal(report+report2)

    # Test 34
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_cube_chanchunks_auto(self):
        """ [cube] test_cube_chanchunks_auto """
        ######################################################################################
        # Test channel chunking for large cubes : automatic calc of nchanchunks . Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_point.ms', tclean_args={'specmode':'cube','imsize':100,'cell':'10.0arcsec','deconvolver':'hogbom'})
        deconvolve(imagename=self.img, niter=10, deconvolver='hogbom')
        self.assertTrue(os.path.exists(self.img+'.psf') and os.path.exists(self.img+'.image') )
        report=th.checkall(imgexist=[self.img+'.image'],imgval=[(self.img+'.image',1.5002,[50,50,0,0]) , (self.img+'.image',0.769,[50,50,0,19]) ])
        self.checkfinal(report)

##############################################
##############################################

##Task level tests : masks and clean boxes.
class test_mask(testref_base):

    # Test 35
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mask_1(self):
        """ [mask] test_mask_1 """
        ######################################################################################
        # Test Input mask as file and string : mfs . Should produce the same results as tclean.
        ######################################################################################
        mstr = 'circle[[50pix,80pix],10pix]'
        self.delData('refim_twochan.ms') # delete data here, since we're not doing that in prepData
        
        th.write_file(self.img+'.mask.txt', '#CRTFv0 CASA Region Text Format version 0\n'+mstr+'\n')
        # delold=False -> don't delete the mask file
        self.prepData('refim_twochan.ms', delold=False, tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'hogbom','usemask':'user','mask':self.img+'.mask.txt'})
        deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0, usemask='user', mask=self.img+'.mask.txt')
        report1=th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',0.0,[50,50,0,0]),(self.img+'.mask',1.0,[50,80,0,0])])

        self.prepData('refim_twochan.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'hogbom','usemask':'user','mask':mstr})
        deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0, usemask='user', mask=mstr)
        report2=th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',0.0,[50,50,0,0]),(self.img+'.mask',1.0,[50,80,0,0])])

        self.checkfinal(report1+report2)

    # Test 36
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mask_2(self):
        """ [mask] test_mask_2 """
        ######################################################################################
        # Test  Input mask as file and string : cube (few channels) . Should produce the same results as tclean.
        ######################################################################################
        mstr =  'circle[[50pix,50pix],10pix],range=[1.1GHz,1.5GHz]'
        self.delData('refim_point.ms') # delete data here, since we're not doing that in prepData

        th.write_file(self.img+'.mask.txt', '#CRTFv0 CASA Region Text Format version 0\n'+mstr+'\n')
        # delold=False -> don't delete the mask file
        self.prepData('refim_point.ms', delold=False, tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'hogbom','specmode':'cube','interactive':0,'usemask':'user','mask':self.img+'.mask.txt'})
        deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0, usemask='user', mask=self.img+'.mask.txt')
        report1=th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',0.0,[50,50,0,1]),(self.img+'.mask',1.0,[50,50,0,2]),(self.img+'.mask',1.0,[50,50,0,10]),(self.img+'.mask',0.0,[50,50,0,11])])

        self.prepData('refim_point.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'hogbom','specmode':'cube','interactive':0,'usemask':'user','mask':mstr})
        deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0, usemask='user', mask=mstr)
        report2=th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',0.0,[50,50,0,1]),(self.img+'.mask',1.0,[50,50,0,2]),(self.img+'.mask',1.0,[50,50,0,10]),(self.img+'.mask',0.0,[50,50,0,11])])

        self.checkfinal(report1+report2)

    # Test 37
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mask_missingfile(self):
        """ [mask] test_mask_missingfile """
        ######################################################################################
        # tst.mask is sometimes required
        ######################################################################################
        mstr = 'circle[[50pix,50pix],10pix],range=[1.1GHz,1.5GHz]'
        mstr = '#CRTFv0 CASA Region Text Format version 0\n'+mstr+'\n'
        mname = self.img+'.mask.txt'
        th.write_file(mname, mstr)
        self.prepData('refim_twochan.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'hogbom','usemask':'user','mask':mname})
        deconvolve_args={'imagename':self.img, 'niter':10, 'deconvolver':'hogbom', 'usemask':'user', 'mask':mname}

        # because tclean apparently deletes the mask file... ugh
        self.assertFalse(os.path.exists(mname))
        th.write_file(mname, mstr)

        # this should be fine ('tst.mask.txt' is present and 'tst.mask' is not)
        os.system("rm -rf "+self.img+".mask")
        deconvolve(**deconvolve_args)

        # this should error out ('tst.mask.txt' has been deleted)
        os.system("rm -rf "+self.img+".mask")
        os.system("rm -rf "+mname)
        strcheck = "'mask' parameter specified as a filename '"+mname+"', but no such file exists"
        with self.assertRaisesRegex(RuntimeError, strcheck):
            deconvolve(**deconvolve_args)

    # Test 38
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mask_pbmask0(self):
        """ [mask] test_mask_pbmask0 """
        ######################################################################################
        # Have deconvolve create the .mask image using the full pb
        ######################################################################################
        mname=self.img+".mask"
        self.prepData('refim_point.ms', tclean_args={'imsize':100, 'cell':'8.0arcsec', 'deconvolver':'hogbom', 'usemask':'pb'})

        # delete the tclean-created .mask image
        os.system("rm -rf "+mname)
        self.assertFalse(os.path.exists(mname))

        # run deconvolve to have it create the .mask image
        ret = deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0, usemask='pb')
        self.assertTrue(os.path.exists(mname), "Mask image was not created by deconvolve task!")

        # verify the mask has the right pixels highlighted
        _ia.open(mname)
        try:
            stats = _ia.statistics()
        finally:
            _ia.close()
        self.assertEqual(10000, stats['sum'][0], "Mask image does not contain the right number of masked pixels (should be 10000 but is {})!".format(stats['sum'][0]))

    # Test 39
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mask_pbmask995(self):
        """ [mask] test_mask_pbmask995 """
        ######################################################################################
        # Have deconvolve create the .mask image using the pb > 0.995
        ######################################################################################
        mname=self.img+".mask"
        self.prepData('refim_point.ms', tclean_args={'imsize':100, 'cell':'8.0arcsec', 'deconvolver':'hogbom', 'usemask':'pb'})

        # delete the tclean-created .mask image
        os.system("rm -rf "+mname)
        self.assertFalse(os.path.exists(mname))

        # run deconvolve to have it create the .mask image
        ret = deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0, usemask='pb', pbmask=0.995)
        self.assertTrue(os.path.exists(mname), "Mask image was not created by deconvolve task!")

        # verify the mask has the right pixels highlighted
        _ia.open(mname)
        try:
            stats = _ia.statistics()
        finally:
            _ia.close()
        self.assertEqual(989, stats['sum'][0], "Mask image does not contain the right number of masked pixels (should be 989 but is {})!".format(stats['sum'][0]))

    # AUTOMASK TESTS
    # Test 40
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mask_autobox_multithresh(self):
        """ [mask] test_mask_autobox_multithresh """
        ######################################################################################
        # Test multi-threshold Autobox (default). Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_twochan.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'hogbom','interactive':0,'usemask':'auto-multithresh'})
        deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0, usemask='auto-multithresh')
        report=th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
        self.checkfinal(report)

    # Test 41
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mask_autobox_multithresh_newnoise(self):
        """ [mask] test_mask_autobox_multithresh_newnoise """
        ######################################################################################
        # Test multi-threshold Autobox (new noise calculation). Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_twochan.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'hogbom','interactive':0,'usemask':'auto-multithresh','fastnoise':False})
        deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0, usemask='auto-multithresh', fastnoise=False)
        report=th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
        self.checkfinal(report)

    # Test 42
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mask_autobox_multithresh_with_nsigma(self):
        """ [mask] test_mask_autobox_multithresh_with_nsigma """
        ######################################################################################
        # Test multi-threshold Autobox (non-default nsigma). Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_twochan.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'hogbom','interactive':0,'usemask':'auto-multithresh'})
        deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0, usemask='auto-multithresh', nsigma=3.0)
        report=th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
        self.checkfinal(report)

    # Test 43
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mask_autobox_multithresh_with_nsigma_newnoise(self):
        """ [mask] test_mask_autobox_multithresh_with_nsigma_newnoise """
        ######################################################################################
        # Test multi-threshold Autobox (new noise calculation & non-default nsigma). Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_twochan.ms', tclean_args={'imsize':100,'cell':'8.0arcsec','deconvolver':'hogbom','interactive':0,'usemask':'auto-multithresh'})
        deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0, usemask='auto-multithresh', nsigma=3.0, fastnoise=False)
        report=th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[50,50,0,0]),(self.img+'.mask',0.0,[50,85,0,0])])
        self.checkfinal(report)

    # Test 44
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mask_autobox_multithresh_with_prune(self):
        """ [mask] test_mask_autobox_multithresh_with_prune """
        ######################################################################################
        # Test multi-threshold Autobox (with pruning). Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_twochan.ms', tclean_args={'imsize':1000,'cell':'8.0arcsec','deconvolver':'hogbom','interactive':0,'usemask':'auto-multithresh'})
        deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0, usemask='auto-multithresh', minbeamfrac=0.3)
        report=th.checkall(imgexist=[self.img+'.mask'], imgval=[(self.img+'.mask',1.0,[500,500,0,0]),(self.img+'.mask',0.0,[500,510,0,0])])
        self.checkfinal(report)

##############################################
##############################################

##Task level tests : run deconvolve multiple times in a row
class test_multirun(testref_base):

    # Test 45
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_multirun_hogbomhogbom(self):
        """ [multirun] test_multirun_hogbomhogbom """
        ######################################################################################
        # Test running hogbom twice in a row and show that it gets the same value as one run with twice the iterations. Should produce the same results as tclean.
        # Note: tclean produced different results when running once with 20 iterations vs twice with 10 iterations, so that is what we are comparing deconvolve against.
        #       We think this is because of a known off-by-one error for hogbom, where iterations aren't counted correctly (niter=10 really means niter=11). See CAS-13200.
        ######################################################################################
        tca = {'imsize':100,'cell':'8.0arcsec','deconvolver':'hogbom','threshold':'1mJy'}
        self.prepData('refim_twochan.ms', tclean_args=tca)
        results1 = deconvolve(imagename=self.img, deconvolver='hogbom', niter=20, threshold='1mJy', interactive=0)
        report1  = th.checkall(ret=results1['retrec'], peakres=0.148, modflux=1.008, iterdone=20, epsilon=0.005, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'],
                               imgval=[(self.img+'.model',0.917,[50,50,0,0])])

        self.prepData('refim_twochan.ms', tclean_args=tca)
        results2 = deconvolve(imagename=self.img, deconvolver='hogbom', niter=10, threshold='1mJy', interactive=0, restoration=False)
        report2  = th.checkall(ret=results2['retrec'], peakres=0.353, modflux=0.772, iterdone=10, epsilon=0.005, imgexist=[self.img+'.psf', self.img+'.residual'], imgexistnot=[self.img+'.image'],
                               imgval=[(self.img+'.model',0.772,[50,50,0,0])])
        results3 = deconvolve(imagename=self.img, deconvolver='hogbom', niter=10, threshold='1mJy', interactive=0)
        report3  = th.checkall(ret=results3['retrec'], peakres=0.142, modflux=1.023, iterdone=10, epsilon=0.005, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'],
                               imgval=[(self.img+'.model',0.917,[50,50,0,0])])

        self.checkfinal(report1 + report2 + report3)

    # Test 46
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_multirun_clarkclark(self):
        """ [multirun] test_multirun_clarkclark """
        ######################################################################################
        # Test running clark twice in a row and show that it gets the same value as one run with twice the iterations. Should produce the same results as tclean.
        ######################################################################################
        tca = {'imsize':100,'cell':'8.0arcsec','deconvolver':'clark','threshold':'1mJy'}
        self.prepData('refim_twochan.ms', tclean_args=tca)
        results1 = deconvolve(imagename=self.img, deconvolver='clark', niter=20, threshold='1mJy', interactive=0)
        report1  = th.checkall(ret=results1['retrec'], peakres=0.161, modflux=0.992, iterdone=20, epsilon=0.005, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'],
                               imgval=[(self.img+'.model',0.917,[50,50,0,0])])

        self.prepData('refim_twochan.ms', tclean_args=tca)
        results2 = deconvolve(imagename=self.img, deconvolver='clark', niter=10, threshold='1mJy', interactive=0, restoration=False)
        report2  = th.checkall(ret=results2['retrec'], peakres=0.392, modflux=0.733, iterdone=10, epsilon=0.005, imgexist=[self.img+'.psf', self.img+'.residual'], imgexistnot=[self.img+'.image'],
                               imgval=[(self.img+'.model',0.733,[50,50,0,0])])
        results3 = deconvolve(imagename=self.img, deconvolver='clark', niter=10, threshold='1mJy', interactive=0)
        report3  = th.checkall(ret=results3['retrec'], peakres=0.161, modflux=0.992, iterdone=10, epsilon=0.005, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'],
                               imgval=[(self.img+'.model',0.917,[50,50,0,0])])

        self.checkfinal(report1 + report2 + report3)

    # Test 47
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_multirun_clarkstokesclarkstokes(self):
        """ [multirun] test_multirun_clarkstokesclarkstokes """
        ######################################################################################
        # Test running clarkstokes twice in a row and show that it gets the same value as one run with twice the iterations. Should produce the same results as tclean.
        ######################################################################################
        tca = {'imsize':100,'cell':'8.0arcsec','deconvolver':'clarkstokes','cyclefactor':0.1,'stokes':'I'}
        tca=self.prepData('refim_point_linRL.ms', tclean_args=tca)
        results1 = deconvolve(imagename=self.img, deconvolver='clarkstokes', niter=20, cyclefactor=0.1, interactive=0)
        report1  = th.checkall(ret=results1['retrec'], peakres=0.122, modflux=0.878, iterdone=20, epsilon=0.005, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'],
                               imgval=[(self.img+'.model',0.878,[50,50,0,0])])

        self.prepData('refim_point_linRL.ms', tclean_args=tca)
        results2 = deconvolve(imagename=self.img, deconvolver='clarkstokes', niter=10, cyclefactor=0.1, interactive=0, restoration=False)
        report2  = th.checkall(ret=results2['retrec'], peakres=0.349, modflux=0.651, iterdone=10, epsilon=0.005, imgexist=[self.img+'.psf', self.img+'.residual'], imgexistnot=[self.img+'.image'],
                               imgval=[(self.img+'.model',0.651,[50,50,0,0])])
        results3 = deconvolve(imagename=self.img, deconvolver='clarkstokes', niter=10, cyclefactor=0.1, interactive=0)
        report3  = th.checkall(ret=results3['retrec'], peakres=0.122, modflux=0.878, iterdone=10, epsilon=0.005, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'],
                               imgval=[(self.img+'.model',0.878,[50,50,0,0])])

        self.checkfinal(report1 + report2 + report3)

    # Test 48
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_multirun_multiscalemultiscale(self):
        """" [multirun] test_multirun_multiscalemultiscale """
        ######################################################################################
        # Test running multiscale twice in a row and show that it gets the same value as one run with twice the iterations. Should produce the same results as tclean.
        ######################################################################################
        tca={'imsize':100,'cell':'8.0arcsec','deconvolver':'multiscale','scales':[0,20,40,100]}
        self.prepData('refim_twochan.ms', tclean_args=tca)
        results1 = deconvolve(imagename=self.img, deconvolver='multiscale', scales=[0,20,40,100], niter=20, threshold='1mJy', interactive=0)
        report1  = th.checkall(ret=results1['retrec'], peakres=0.158, modflux=1.024, iterdone=20, epsilon=0.005, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'],
                               imgval=[(self.img+'.model',0.917,[50,50,0,0])])

        self.prepData('refim_twochan.ms', tclean_args=tca)
        results2 = deconvolve(imagename=self.img, deconvolver='multiscale', scales=[0,20,40,100], niter=10, threshold='1mJy', interactive=0, restoration=False)
        report2  = th.checkall(ret=results2['retrec'], peakres=0.392, modflux=0.733, iterdone=10, epsilon=0.005, imgexist=[self.img+'.psf', self.img+'.residual'], imgexistnot=[self.img+'.image'],
                               imgval=[(self.img+'.model',0.733,[50,50,0,0])])
        results3 = deconvolve(imagename=self.img, deconvolver='multiscale', scales=[0,20,40,100], niter=10, threshold='1mJy', interactive=0)
        report3  = th.checkall(ret=results3['retrec'], peakres=0.158, modflux=1.024, iterdone=10, epsilon=0.005, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image'],
                               imgval=[(self.img+'.model',0.917,[50,50,0,0])])

        self.checkfinal(report1 + report2 + report3)

    # Test 49
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_multirun_mtmfsmtmfs(self):
        """" [multirun] test_multirun_mtmfsmtmfs """
        ######################################################################################
        # Test running mtmfs twice in a row and show that it gets the same value as one run with twice the iterations. Should produce the same results as tclean.
        ######################################################################################
        tca={'imsize':100,'cell':'8.0arcsec','deconvolver':'mtmfs','scales':[0,20,40,100]}
        self.prepData('refim_eptwochan.ms', tclean_args=tca)
        results1 = deconvolve(imagename=self.img, deconvolver='mtmfs', scales=[0,20,40,100], niter=20, threshold='1mJy', interactive=0)
        report11 = th.checkall( ret=results1['retrec'], peakres=0.444, modflux=3.450, iterdone=20, epsilon=0.005,
                                imgexist=[self.img+'.psf.tt0', self.img+'.psf.tt1', self.img+'.psf.tt2', self.img+'.residual.tt0', self.img+'.residual.tt1', self.img+'.image.tt0', self.img+'.image.tt1'] )
        report12 = th.checkall( imgval=[(self.img+'.model.tt0',0.008,[50,50,0,0]), (self.img+'.model.tt1',-0.007,[50,50,0,0])], epsilon=0.05 )

        self.prepData('refim_eptwochan.ms', tclean_args=tca)
        results2 = deconvolve(imagename=self.img, deconvolver='mtmfs', scales=[0,20,40,100], niter=10, threshold='1mJy', interactive=0, restoration=False)
        report21 = th.checkall( ret=results2['retrec'], peakres=0.782, modflux=2.442, iterdone=10, epsilon=0.005,
                                imgexist=[self.img+'.psf.tt0', self.img+'.psf.tt1', self.img+'.psf.tt2', self.img+'.residual.tt0', self.img+'.residual.tt1'], imgexistnot=[self.img+'.image.tt0', self.img+'.image.tt1'] )
        report22 = th.checkall( imgval=[(self.img+'.model.tt0',0.006,[50,50,0,0]), (self.img+'.model.tt1',-0.005,[50,50,0,0])], epsilon=0.05 )
        results3 = deconvolve(imagename=self.img, deconvolver='mtmfs', scales=[0,20,40,100], niter=10, threshold='1mJy', interactive=0)
        report31 = th.checkall( ret=results3['retrec'], peakres=0.444, modflux=3.450, iterdone=10, epsilon=0.005,
                                imgexist=[self.img+'.psf.tt0', self.img+'.psf.tt1', self.img+'.psf.tt2', self.img+'.residual.tt0', self.img+'.residual.tt1', self.img+'.image.tt0', self.img+'.image.tt1'] )
        report32 = th.checkall( imgval=[(self.img+'.model.tt0',0.008,[50,50,0,0]), (self.img+'.model.tt1',-0.007,[50,50,0,0])], epsilon=0.05 )

        self.checkfinal(report11 + report12 + report21 + report22 + report31 + report32)

    # Test 50
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_multirun_multiscalehog(self):
        """ [multirun] test_multirun_multiscalehog """
        ######################################################################################
        # Test running multiscale clean followed by hogbom. Should produce the same results as tclean.
        # Note: aren't completely sure of what the value should be at the end. (TODO needs validation)
        ######################################################################################
        self.prepData('refim_eptwochan.ms', tclean_args={'imsize':200, 'cell':'8.0arcsec', 'deconvolver':'multiscale', 'scales':[0,20,40,100], 'threshold':'1mJy'})
        results1 = deconvolve(imagename=self.img, niter=10, deconvolver='multiscale', scales=[0,20,40,100], interactive=0, restoration=False, threshold='1mJy')
        report1 = th.checkall(ret=results1['retrec'], peakres=0.822, modflux=3.816, iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.model'], imgexistnot=[self.img+'.image'])
        results2 = deconvolve(imagename=self.img, niter=10, deconvolver='hogbom', interactive=0)
        report2 = th.checkall(ret=results2['retrec'], peakres=0.283, modflux=4.395, iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.image', self.img+'.model'],
                              imgval=[(self.img+'.model',0.453,[94,107,0,0])])

        self.checkfinal(report1 + report2)

    # Test 51
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_multirun_norestore_restore(self):
        """ [multirun] test_multirun_norestore_restore """
        ######################################################################################
        # Test to test the retore-only feature . Should produce the same results as tclean.
        ######################################################################################
        self.prepData('refim_point.ms', tclean_args={'imsize':100, 'cell':['10.0arcsec','30.0arcsec']})
        results1 = deconvolve(imagename=self.img, niter=10, restoration=False)
        report1=th.checkall(ret=results1['retrec'], iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.model'], imgexistnot=[self.img+'.image'])
        results2 = deconvolve(imagename=self.img, niter=0, restoration=True)
        report2=th.checkall(ret=results2['retrec'], iterdone=10, imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.model', self.img+'.image'],
                            imgval=[(self.img+'.image',0.482,[50,49,0,0])] )
        self.checkfinal(report1 + report2)

##############################################
##############################################

##Task level tests : test validation of input images and copying/converting startmodel images for the hogbom deconvolver
# Note that all of these tests utilize using the same ms and same tclean args in order to be able to cache the results
# from tclean in order to execute faster.
class test_imgval(testref_base):
    @classmethod
    def setUpClass(cls):
        super(test_imgval, cls).setUpClass()
        msfile='refim_point.ms'
        cls.staticDelData(msfile)
        cls.staticPrepData(msfile, tclean_args={'imsize':100, 'cell':['10.0arcsec','30.0arcsec']})
        cls.staticCopyToCache(msfile, imagename=cls.img, cachedir='imgval_cache')

    def ivsetup(self):
        # we can use a cache here because tclean was run only run once, during setUpClass
        self.delData()
        type(self).staticCopyFromCache()
        self.mname = self.img + ".model"
        self.mname2 = self.img + "_2.model"

    def tearDown(self):
        super(test_imgval, self).tearDown()

        # remove all tst_2.* and tst_bak.* images
        for uname in ["_bak", "_2"]:
            fn = self.img + uname + ".*"
            if os.path.exists(fn):
                os.system('rm -rf '+fn)

    def get_csys_crval0(self, img):
        _ia.open(img)
        csys = _ia.coordsys()
        # pnt = csys.torecord()['direction0']['crval'][0]
        pnt = csys.torecord()['direction0']['crpix'][0]
        _ia.close()
        return csys, pnt

    def set_crval0(self, img, crval0):
        _ia.open(img)
        csys = _ia.coordsys()
        rec = csys.torecord()
        # rec['direction0']['crval'][0] = crval0
        rec['direction0']['crpix'][0] = crval0
        csys.fromrecord(rec)
        _ia.setcoordsys(csys.torecord())
        _ia.close()

    def get_shape(self, img):
        _ia.open(img)
        ret = _ia.shape()
        _ia.close()
        return ret

    def helper_imgval_missingimgs(self, ext):
        # helper method for test_imgval_missingimgs_*
        self.ivsetup()
        os.system("mv {0}{1} {0}_bak{1}".format(self.img, ext))

        strcheck = "missing one or more of the required images"
        with self.assertRaisesRegex(RuntimeError, strcheck):
            deconvolve(imagename=self.img, niter=10)

    def helper_imgval_axesmismatch(self, ext, deconvolve_args={}):
        # helper method for test_imgval_axesmismatch_*
        self.ivsetup()

        # cause the tst.model image to exist prior to looking for it
        deconvolve(imagename=self.img, niter=1, restoration=False, **deconvolve_args)

        # verify the images have the wrong format and an exception is thrown
        # Note: cpp code does not mind that the .pb image has a weird axes order
        fn1 = self.img + ext
        fn2 = self.img + "_bak" + ext
        os.system("mv {0} {1}".format(fn1, fn2))

        imtrans(imagename=self.img+"_bak"+ext, outfile=self.img+ext, order="3012")
        strcheck = "(There is a shape mismatch between existing images|There is a coordinate system mismatch between existing images on disk and current parameters)"
        with self.assertRaisesRegex(RuntimeError, strcheck):
            deconvolve(imagename=self.img, niter=10, **deconvolve_args)

    def helper_imgval_shapemismatch(self, ext, deconvolve_args={}):
        """All input images should have the same shape as the tst.residual image"""
        self.ivsetup()

        # cause the tst.model image to exist prior to looking for it
        deconvolve(imagename=self.img, niter=1, restoration=False, **deconvolve_args)

        # verify the images have the wrong format and an exception is thrown
        # Note: cpp code does not mind that the .pb image has a weird shape
        os.system("mv {0}{1} {0}_bak{1}".format(self.img, ext))

        imrebin(imagename=self.img+"_bak"+ext, outfile=self.img+ext, factor=[50,50])
        strcheck = "(There is a shape mismatch between existing images|There is a coordinate system mismatch between existing images on disk and current parameters)"
        with self.assertRaisesRegex(RuntimeError, strcheck):
            deconvolve(imagename=self.img, niter=10, **deconvolve_args)

    # Test 52
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_imgval_missingimgs_residual(self):
        """ [imgval] test_imgval_missingimgs_residual """
        ######################################################################################
        # tst.residual and tst.psf are always required
        ######################################################################################
        self.helper_imgval_missingimgs(".residual")

    # Test 53
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_imgval_missingimgs_psf(self):
        """ [imgval] test_imgval_missingimgs_psf """
        ######################################################################################
        # tst.residual and tst.psf are always required
        ######################################################################################
        # Note: cpp code doesn't throw an exception when psf is missing, just prints a warning
        self.helper_imgval_missingimgs(".psf")

    # Test 54
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_imgval_missingimgs_model(self):
        """ [imgval] test_imgval_missingimgs_model """
        ######################################################################################
        # tst.model is used to continue deconvolution, but is not required.
        ######################################################################################
        self.ivsetup()
        if (os.path.exists(self.img+'.model')):
            os.system("rm -rf "+self.img+".model")
        deconvolve(imagename=self.img, niter=10)

    # Test 55
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_imgval_missingimgs_sumwt(self):
        """ [imgval] test_imgval_missingimgs_sumwt """
        ######################################################################################
        # tst.sumwt is never required
        ######################################################################################
        self.ivsetup()
        if (os.path.exists(self.img+'.sumwt')):
            os.system("rm -rf "+self.img+".sumwt")

        # Should be fine. Sumwt should not be required for task deconvolve.
        deconvolve(imagename=self.img, niter=10)

    # Test 56
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_imgval_axesmismatch_residual(self):
        """ [imgval] test_imgval_axesmismatch_residual """
        ######################################################################################
        # all input images must have the axes as they are given in tclean
        ######################################################################################
        self.helper_imgval_axesmismatch(".residual")

    # Test 57
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_imgval_axesmismatch_psf(self):
        """ [imgval] test_imgval_axesmismatch_psf """
        ######################################################################################
        # all input images must have the axes as they are given in tclean
        ######################################################################################
        self.helper_imgval_axesmismatch(".psf")

    # Test 58
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_imgval_axesmismatch_model(self):
        """ [imgval] test_imgval_axesmismatch_model """
        ######################################################################################
        # all input images must have the axes as they are given in tclean
        ######################################################################################
        self.helper_imgval_axesmismatch(".model")

    # Test 59
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_imgval_axesmismatch_pb(self):
        """ [imgval] test_imgval_axesmismatch_pb """
        ######################################################################################
        # all input images must have the axes as they are given in tclean
        ######################################################################################
        self.helper_imgval_axesmismatch(".pb", deconvolve_args={'usemask':'pb', 'pbmask':0.2})

    # Test 60
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_imgval_shapemismatch_residual(self):
        """ [imgval] test_imgval_shapemismatch_residual """
        ######################################################################################
        # All input images should have the same shape as the tst.residual image
        ######################################################################################
        self.helper_imgval_shapemismatch(".residual")

    # Test 61
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_imgval_shapemismatch_psf(self):
        """ [imgval] test_imgval_shapemismatch_psf """
        ######################################################################################
        # All input images should have the same shape as the tst.residual image
        ######################################################################################
        self.helper_imgval_shapemismatch(".psf")

    # Test 62
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_imgval_shapemismatch_model(self):
        """ [imgval] test_imgval_shapemismatch_model """
        ######################################################################################
        # All input images should have the same shape as the tst.residual image
        ######################################################################################
        self.helper_imgval_shapemismatch(".model")

    # Test 63
    # @unittest.skip("The cpp code does not mind that the .pb image has a weird shape; no exception is thrown.")
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_imgval_shapemismatch_pb(self):
        """ [imgval] test_imgval_shapemismatch_pb """
        ######################################################################################
        # All input images should have the same shape as the tst.residual image
        ######################################################################################
        self.helper_imgval_shapemismatch(".pb", deconvolve_args={'usemask':'pb', 'pbmask':0.2})
    
    # TODO figure out why running the startmodel_axesmismatch test immediately before this test causes an exception to be thrown
    # Test 64
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    @unittest.skip("if test_imgval_startmodel_axesmismatch executes immediately before this test then this test fails")
    def test_imgval_startmodel_empty(self):
        """ [imgval] test_imgval_startmodel_empty """
        ######################################################################################
        # Deconvolve should ignore all empty strings entered for the startmodel
        ######################################################################################
        # self.test_imgval_startmodel_axesmismatch()
        self.ivsetup()

        # basic test with empty string as startmodel
        deconvolve(imagename=self.img, niter=10, startmodel='')

        # basic test with list of empty strings as startmodel
        deconvolve(imagename=self.img, niter=10, startmodel=['', '', ''])

        # basic copy test where empty string is discarded from list
        self.mname = self.img + ".model"
        self.mname2 = self.img + "_2.model"
        os.system("mv {0} {1}".format(self.mname, self.mname2))
        deconvolve(imagename=self.img, niter=10, startmodel=['', '', self.mname2, '', ''])

    # Test 65
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_imgval_startmodel_dne(self):
        """ [imgval] test_imgval_startmodel_dne """
        ######################################################################################
        # Throws an error if startmodel is set but does not exist
        ######################################################################################
        self.ivsetup()
        strcheck = "does not exist"
        with self.assertRaisesRegex(RuntimeError, strcheck):
            deconvolve(imagename=self.img, niter=10, startmodel='doesnotexists.model')

    # Test 66
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_imgval_startmodel_model_exists(self):
        """ [imgval] test_imgval_startmodel_model_exists """
        ######################################################################################
        # Throws an error if startmodel is set and tst.model exists (must be one or the other, not both)
        ######################################################################################
        self.ivsetup()
        deconvolve(imagename=self.img, niter=10)
        shutil.copytree(self.mname, self.mname2)
        
        strcheck = "exists"
        with self.assertRaisesRegex(RuntimeError, strcheck):
            deconvolve(imagename=self.img, niter=10, startmodel=self.mname2)

    # TODO figure out why running the startmodel_axesmismatch test immediately before this test causes an exception to be thrown
    # Test 67
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    @unittest.skip("if test_imgval_startmodel_axesmismatch executes immediately before this test then this test fails")
    def test_imgval_startmodel_basic_copy(self):
        """ [imgval] test_imgval_startmodel_basic_copy """
        ######################################################################################
        # Tests ability of deconvolve to copy startmodel to tst.model before starting deconvolution
        ######################################################################################
        # self.test_imgval_startmodel_axesmismatch()
        self.ivsetup()
        deconvolve(imagename=self.img, niter=10, restoration=False) # generate the first model to work off of
        os.system("mv {0} {1}".format(self.mname, self.mname2))
            
        self.assertTrue(os.path.exists(self.mname2) and not os.path.exists(self.mname), "File {0} did not get moved to {1}!".format(self.mname, self.mname2))
        deconvolve(imagename=self.img, niter=10, startmodel=self.mname2)
        self.assertTrue(os.path.exists(self.mname), "File {0} did not get copied!".format(self.mname))

    # Test 68
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_imgval_startmodel_axesmismatch(self):
        """ [imgval] test_imgval_startmodel_axesmismatch """
        ######################################################################################
        # Tests the existing functionality. If in the future the logic is added to auto-translate images, this test can be removed.
        ######################################################################################
        self.ivsetup()
        deconvolve(imagename=self.img, niter=10, restoration=False) # generate the first model to work off of
        imtrans(imagename=self.mname, outfile=self.mname2, order="3012")
            
        os.system("rm -rf "+self.mname)
        self.assertTrue(os.path.exists(self.mname2) and not os.path.exists(self.mname), "File {0} did not get translated to {1}!".format(self.mname, self.mname2))
        strcheck = "Error in setting"
        with self.assertRaisesRegex(RuntimeError, strcheck):
            deconvolve(imagename=self.img, niter=10, startmodel=self.mname2)

    # Test 69
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_imgval_startmodel_csysmismatch(self):
        """ [imgval] test_imgval_startmodel_csysmismatch """
        ######################################################################################
        # Tests the ability of the deconvolve regrid the csys of the startmodel to that of tst.residual
        ######################################################################################
        self.ivsetup()
        deconvolve(imagename=self.img, niter=10, restoration=False) # generate the first model to work off of
        os.system("mv {0} {1}".format(self.mname, self.mname2))

        # change csys
        csys, oldpnt = self.get_csys_crval0(self.mname2)
        newpnt = 51#oldpnt * 0.9999
        self.assertNotEqual(oldpnt, newpnt, "Change amount not big enough")
        self.set_crval0(self.mname2, newpnt)
        csys2, newpnt2 = self.get_csys_crval0(self.mname2)
        self.assertEqual(newpnt2, newpnt, "Image {0} did not get its csys.direction0.crval[0] value updated properly from {1} to the expected {2}! (actual value is {3})".format(self.mname2, oldpnt, newpnt, newpnt2))

        # test that deconvolve regrids the image
        os.system("rm -rf "+self.mname)
        self.assertTrue(os.path.exists(self.mname2) and not os.path.exists(self.mname), "File {0} did not get rebinned to {1}!".format(self.mname, self.mname2))
        deconvolve(imagename=self.img, niter=10, startmodel=self.mname2)
        self.assertTrue(os.path.exists(self.mname), "File {0} did not get regridded back!".format(self.mname))
        csys3, regridpnt = self.get_csys_crval0(self.mname)
        self.assertAlmostEqual(regridpnt, oldpnt, "Image {0} did not get its csys.direction0.crval[0] value regridded properly from {1} to {2}! (actual value is {3})".format(self.mname2, newpnt, oldpnt, regridpnt))

    # Test 70
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_imgval_startmodel_shapemismatch(self):
        """ [imgval] test_imgval_startmodel_shapemismatch """
        ######################################################################################
        # Tests the ability of the deconvolve regrid the shape of the startmodel to that of tst.residual
        ######################################################################################
        self.ivsetup()
        deconvolve(imagename=self.img, niter=10, restoration=False) # generate the first model to work off of

        # change shape
        s = self.get_shape(self.mname)
        self.assertEqual([s[0], s[1]], [100, 100], "Image shape ({0}) didn't start as expected ({1})".format(s, [100, 100]))
        imrebin(imagename=self.mname, outfile=self.mname2, factor=[2,2])

        # sanity: make sure the shape changed
        s = self.get_shape(self.mname2)
        self.assertEqual([s[0], s[1]], [50, 50], "Image shape ({0}) didn't get rebinned as expected ({1})".format(s, [50, 50]))

        # run deconvolve and make sure the shape gets regridded back in
        os.system("rm -rf "+self.mname)
        self.assertTrue(os.path.exists(self.mname2) and not os.path.exists(self.mname), "File {0} did not get rebinned to {1}!".format(self.mname, self.mname2))
        deconvolve(imagename=self.img, niter=10, startmodel=self.mname2)
        self.assertTrue(os.path.exists(self.mname), "File {0} did not get regridded back!".format(self.mname))
        s = self.get_shape(self.mname)
        self.assertEqual([s[0], s[1]], [100, 100], "Image shape ({0}) didn't get regridded as expected ({1})".format(s, [100, 100]))

##############################################
##############################################

##Task level tests : test validation of input images and copying/converting startmodel images for the mtmfs deconvolver
# Note that all of these tests utilize using the same ms and same tclean args in order to be able to cache the results
# from tclean in order to execute faster.
class test_mtmfsimgval(testref_base):
    @classmethod
    def setUpClass(cls):
        super(test_mtmfsimgval, cls).setUpClass()
        msfile='refim_eptwochan.ms'
        cls.staticDelData(msfile)
        cls.staticPrepData(msfile, tclean_args={'imsize':10, 'cell':'8.0arcsec', 'deconvolver':'mtmfs'})
        cls.staticCopyToCache(msfile, imagename=cls.img, cachedir='mtmfsimgval_cache')

    def ivsetup(self):
        # we can use a cache here because tclean was run only run once, during setUpClass
        self.delData()
        type(self).staticCopyFromCache()
        self.mname = self.img + ".model.tt1"
        self.mname2 = self.img + "_2.model.tt1"

    def tearDown(self):
        super(test_mtmfsimgval, self).tearDown()

        # remove all tst_2.* and tst_bak.* images
        for uname in ["_bak", "_2"]:
            fn = self.img + uname + ".*"
            if os.path.exists(fn):
                os.system('rm -rf '+fn)

    def get_csys_crval0(self, img):
        _ia.open(img)
        csys = _ia.coordsys()
        pnt = csys.torecord()['direction0']['crpix'][0]
        _ia.close()
        return csys, pnt

    def set_crval0(self, img, crval0):
        _ia.open(img)
        csys = _ia.coordsys()
        rec = csys.torecord()
        # rec['direction0']['crval'][0] = crval0
        rec['direction0']['crpix'][0] = crval0
        csys.fromrecord(rec)
        _ia.setcoordsys(csys.torecord())
        _ia.close()

    def get_shape(self, img):
        _ia.open(img)
        ret = _ia.shape()
        _ia.close()
        return ret

    def helper_mtmfsimgval_missingimgs(self, ext):
        # helper method for test_mtmfsimgval_missingimgs_*
        self.ivsetup()
        os.system("mv {0}{1}.tt1 {0}_bak{1}.tt1".format(self.img, ext))

        strcheck = "missing one or more of the required images"
        with self.assertRaisesRegex(RuntimeError, strcheck):
            deconvolve(imagename=self.img, niter=10, deconvolver='mtmfs')

    def helper_mtmfsimgval_axesmismatch(self, ext, ttn=".tt1", deconvolve_args={}):
        # helper method for test_mtmfsimgval_axesmismatch_*
        self.ivsetup()

        # cause the tst.model image to exist prior to looking for it
        deconvolve(imagename=self.img, niter=1, restoration=False, deconvolver='mtmfs', **deconvolve_args)

        # verify the images have the wrong format and an exception is thrown
        # Note: cpp code does not mind that the .pb image has a weird axes order
        fn1 = self.img + ext + ttn
        fn2 = self.img + "_bak" + ext + ttn
        os.system("mv {0} {1}".format(fn1, fn2))

        imtrans(imagename=self.img+"_bak"+ext+ttn, outfile=self.img+ext+ttn, order="3012")
        strcheck = "(There is a shape mismatch between existing images|There is a coordinate system mismatch between existing images on disk and current parameters)"
        with self.assertRaisesRegex(RuntimeError, strcheck):
            deconvolve(imagename=self.img, niter=10, deconvolver='mtmfs', **deconvolve_args)

    def helper_mtmfsimgval_shapemismatch(self, ext, ttn=".tt1", deconvolve_args={}):
        """All input images should have the same shape as the tst.residual image"""
        self.ivsetup()

        # cause the tst.model image to exist prior to looking for it
        deconvolve(imagename=self.img, niter=1, restoration=False, deconvolver='mtmfs', **deconvolve_args)

        # verify the images have the wrong format and an exception is thrown
        # Note: cpp code does not mind that the .pb image has a weird shape
        os.system("mv {0}{1}{2} {0}_bak{1}{2}".format(self.img, ext, ttn))

        imrebin(imagename=self.img+"_bak"+ext+ttn, outfile=self.img+ext+ttn, factor=[2,2])
        strcheck = "(There is a shape mismatch between existing images|There is a coordinate system mismatch between existing images on disk and current parameters)"
        with self.assertRaisesRegex(RuntimeError, strcheck):
            deconvolve(imagename=self.img, niter=10, deconvolver='mtmfs', **deconvolve_args)

    # Test 71
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mtmfsimgval_missingimgs_residual(self):
        """ [mtmfsimgval] test_mtmfsimgval_missingimgs_residual """
        ######################################################################################
        # tst.residual and tst.psf are always required
        ######################################################################################
        self.helper_mtmfsimgval_missingimgs(".residual")

    # Test 72
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mtmfsimgval_missingimgs_psf(self):
        """ [mtmfsimgval] test_mtmfsimgval_missingimgs_psf """
        ######################################################################################
        # tst.residual and tst.psf are always required
        ######################################################################################
        # Note: cpp code doesn't throw an exception when psf is missing, just prints a warning
        self.helper_mtmfsimgval_missingimgs(".psf")

    # Test 73
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mtmfsimgval_missingimgs_model(self):
        """ [mtmfsimgval] test_mtmfsimgval_missingimgs_model """
        ######################################################################################
        # tst.model is used to continue deconvolution, but is not required.
        ######################################################################################
        self.ivsetup()
        if (os.path.exists(self.img+'.model')):
            os.system("rm -rf "+self.img+".model")

        # Should be fine. Model should not be required for the first run of task deconvolve.
        deconvolve(imagename=self.img, niter=10, deconvolver='mtmfs')

    # Test 74
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mtmfsimgval_missingimgs_sumwt(self):
        """ [mtmfsimgval] test_mtmfsimgval_missingimgs_sumwt """
        ######################################################################################
        # tst.sumwt is never required
        ######################################################################################
        self.ivsetup()
        if (os.path.exists(self.img+'.sumwt')):
            os.system("rm -rf "+self.img+".sumwt")

        # Should be fine. Sumwt should not be required for task deconvolve.
        deconvolve(imagename=self.img, niter=10, deconvolver='mtmfs')

    # Test 75
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mtmfsimgval_axesmismatch_residual(self):
        """ [mtmfsimgval] test_mtmfsimgval_axesmismatch_residual """
        ######################################################################################
        # all input images must have the axes as they are given in tclean
        ######################################################################################
        self.helper_mtmfsimgval_axesmismatch(".residual")

    # Test 76
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mtmfsimgval_axesmismatch_psf(self):
        """ [mtmfsimgval] test_mtmfsimgval_axesmismatch_psf """
        ######################################################################################
        # all input images must have the axes as they are given in tclean
        ######################################################################################
        self.helper_mtmfsimgval_axesmismatch(".psf")

    # Test 77
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mtmfsimgval_axesmismatch_model(self):
        """ [mtmfsimgval] test_mtmfsimgval_axesmismatch_model """
        ######################################################################################
        # all input images must have the axes as they are given in tclean
        ######################################################################################
        self.helper_mtmfsimgval_axesmismatch(".model")

    # Test 78
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mtmfsimgval_axesmismatch_pb(self):
        """ [mtmfsimgval] test_mtmfsimgval_axesmismatch_pb """
        ######################################################################################
        # all input images must have the axes as they are given in tclean
        ######################################################################################
        self.helper_mtmfsimgval_axesmismatch(".pb", ttn=".tt0", deconvolve_args={'usemask':'pb', 'pbmask':0.2})

    # Test 79
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mtmfsimgval_shapemismatch_residual(self):
        """ [mtmfsimgval] test_mtmfsimgval_shapemismatch_residual """
        ######################################################################################
        # All input images should have the same shape as the tst.residual image
        ######################################################################################
        self.helper_mtmfsimgval_shapemismatch(".residual")

    # Test 80
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mtmfsimgval_shapemismatch_psf(self):
        """ [mtmfsimgval] test_mtmfsimgval_shapemismatch_psf """
        ######################################################################################
        # All input images should have the same shape as the tst.residual image
        ######################################################################################
        self.helper_mtmfsimgval_shapemismatch(".psf")

    # Test 81
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mtmfsimgval_shapemismatch_model(self):
        """ [mtmfsimgval] test_mtmfsimgval_shapemismatch_model """
        ######################################################################################
        # All input images should have the same shape as the tst.residual image
        ######################################################################################
        self.helper_mtmfsimgval_shapemismatch(".model")

    # Test 82
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mtmfsimgval_shapemismatch_pb(self):
        """ [mtmfsimgval] test_mtmfsimgval_shapemismatch_pb """
        ######################################################################################
        # All input images should have the same shape as the tst.residual image
        ######################################################################################
        self.helper_mtmfsimgval_shapemismatch(".pb", ttn=".tt0", deconvolve_args={'usemask':'pb', 'pbmask':0.2})
    
    # TODO figure out why running the startmodel_axesmismatch test immediately before this test causes an exception to be thrown
    # Test 83
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    @unittest.skip("if test_mtmfsimgval_startmodel_axesmismatch executes immediately before this test then this test fails")
    def test_mtmfsimgval_startmodel_empty(self):
        """ [mtmfsimgval] test_mtmfsimgval_startmodel_empty """
        ######################################################################################
        # Deconvolve should ignore all empty strings entered for the startmodel
        ######################################################################################
        # self.test_mtmfsimgval_startmodel_axesmismatch()
        self.ivsetup()

        # basic test with empty string as startmodel
        deconvolve(imagename=self.img, niter=10, startmodel='', deconvolver='mtmfs')

        # basic test with list of empty strings as startmodel
        deconvolve(imagename=self.img, niter=10, startmodel=['', '', ''], deconvolver='mtmfs')

        # basic copy test where empty string is discarded from list
        self.mname = self.img + ".model.tt1"
        self.mname2 = self.img + "_2.model.tt1"
        os.system("mv {0} {1}".format(self.mname, self.mname2))
        deconvolve(imagename=self.img, niter=10, startmodel=['', '', self.mname2, '', ''], deconvolver='mtmfs')

    # Test 84
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mtmfsimgval_startmodel_dne(self):
        """ [mtmfsimgval] test_mtmfsimgval_startmodel_dne """
        ######################################################################################
        # Throws an error if startmodel is set but does not exist
        ######################################################################################
        self.ivsetup()
        strcheck = "does not exist"
        with self.assertRaisesRegex(RuntimeError, strcheck):
            deconvolve(imagename=self.img, niter=10, startmodel='doesnotexists.model', deconvolver='mtmfs')

    # Test 85
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mtmfsimgval_startmodel_model_exists(self):
        """ [mtmfsimgval] test_mtmfsimgval_startmodel_model_exists """
        ######################################################################################
        # Throws an error if startmodel is set and tst.model exists (must be one or the other, not both)
        ######################################################################################
        self.ivsetup()
        deconvolve(imagename=self.img, niter=10, deconvolver='mtmfs')
        shutil.copytree(self.mname, self.mname2)
        
        strcheck = "exists"
        with self.assertRaisesRegex(RuntimeError, strcheck):
            deconvolve(imagename=self.img, niter=10, startmodel=self.mname2, deconvolver='mtmfs')

    # TODO figure out why running the startmodel_axesmismatch test immediately before this test causes an exception to be thrown
    # Test 86
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    @unittest.skip("if test_mtmfsimgval_startmodel_axesmismatch executes immediately before this test then this test fails")
    def test_mtmfsimgval_startmodel_basic_copy(self):
        """ [mtmfsimgval] test_mtmfsimgval_startmodel_basic_copy """
        ######################################################################################
        # Tests ability of deconvolve to copy startmodel to tst.model before starting deconvolution
        ######################################################################################
        # self.test_mtmfsimgval_startmodel_axesmismatch()
        self.ivsetup()
        deconvolve(imagename=self.img, niter=10, restoration=False, deconvolver='mtmfs') # generate the first model to work off of
        os.system("mv {0} {1}".format(self.mname, self.mname2))
            
        self.assertTrue(os.path.exists(self.mname2) and not os.path.exists(self.mname), "File {0} did not get moved to {1}!".format(self.mname, self.mname2))
        deconvolve(imagename=self.img, niter=10, startmodel=self.mname2, deconvolver='mtmfs')
        self.assertTrue(os.path.exists(self.mname), "File {0} did not get copied!".format(self.mname))

    # Test 87
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mtmfsimgval_startmodel_axesmismatch(self):
        """ [mtmfsimgval] test_mtmfsimgval_startmodel_axesmismatch """
        ######################################################################################
        # Tests the existing functionality. If in the future the logic is added to auto-translate images, this test can be removed.
        ######################################################################################
        self.ivsetup()
        deconvolve(imagename=self.img, niter=10, restoration=False, deconvolver='mtmfs') # generate the first model to work off of
        imtrans(imagename=self.mname, outfile=self.mname2, order="3012")
            
        os.system("rm -rf "+self.mname)
        self.assertTrue(os.path.exists(self.mname2) and not os.path.exists(self.mname), "File {0} did not get translated to {1}!".format(self.mname, self.mname2))
        strcheck = "Error in setting"
        with self.assertRaisesRegex(RuntimeError, strcheck):
            deconvolve(imagename=self.img, niter=10, startmodel=self.mname2, deconvolver='mtmfs')

    # Test 88
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mtmfsimgval_startmodel_csysmismatch(self):
        """ [mtmfsimgval] test_mtmfsimgval_startmodel_csysmismatch """
        ######################################################################################
        # Tests the ability of the deconvolve regrid the csys of the startmodel to that of tst.residual
        ######################################################################################
        self.ivsetup()
        deconvolve(imagename=self.img, niter=10, restoration=False, deconvolver='mtmfs') # generate the first model to work off of
        os.system("mv {0} {1}".format(self.mname, self.mname2))

        # change csys
        csys, oldpnt = self.get_csys_crval0(self.mname2)
        newpnt = 51#oldpnt * 0.9999
        self.assertNotEqual(oldpnt, newpnt, "Change amount not big enough")
        self.set_crval0(self.mname2, newpnt)
        csys2, newpnt2 = self.get_csys_crval0(self.mname2)
        self.assertEqual(newpnt2, newpnt, "Image {0} did not get its csys.direction0.crval[0] value updated properly from {1} to the expected {2}! (actual value is {3})".format(self.mname2, oldpnt, newpnt, newpnt2))

        # test that deconvolve regrids the image
        os.system("rm -rf "+self.mname)
        self.assertTrue(os.path.exists(self.mname2) and not os.path.exists(self.mname), "File {0} did not get rebinned to {1}!".format(self.mname, self.mname2))
        deconvolve(imagename=self.img, niter=10, startmodel=self.mname2, deconvolver='mtmfs')
        self.assertTrue(os.path.exists(self.mname), "File {0} did not get regridded back!".format(self.mname))
        csys3, regridpnt = self.get_csys_crval0(self.mname)
        self.assertAlmostEqual(regridpnt, oldpnt, "Image {0} did not get its csys.direction0.crval[0] value regridded properly from {1} to {2}! (actual value is {3})".format(self.mname2, newpnt, oldpnt, regridpnt))

    # Test 89
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_mtmfsimgval_startmodel_shapemismatch(self):
        """ [mtmfsimgval] test_mtmfsimgval_startmodel_shapemismatch """
        ######################################################################################
        # Tests the ability of the deconvolve regrid the shape of the startmodel to that of tst.residual
        ######################################################################################
        self.ivsetup()
        deconvolve(imagename=self.img, niter=10, restoration=False, deconvolver='mtmfs') # generate the first model to work off of

        # change shape
        s = self.get_shape(self.mname)
        self.assertEqual([s[0], s[1]], [10, 10], "Image shape ({0}) didn't start as expected ({1})".format(s, [10, 10]))
        imrebin(imagename=self.mname, outfile=self.mname2, factor=[2,2])

        # sanity: make sure the shape changed
        s = self.get_shape(self.mname2)
        self.assertEqual([s[0], s[1]], [5, 5], "Image shape ({0}) didn't get rebinned as expected ({1})".format(s, [5, 5]))

        # run deconvolve and make sure the shape gets regridded back in
        os.system("rm -rf "+self.mname)
        self.assertTrue(os.path.exists(self.mname2) and not os.path.exists(self.mname), "File {0} did not get rebinned to {1}!".format(self.mname, self.mname2))
        deconvolve(imagename=self.img, niter=10, startmodel=self.mname2, deconvolver='mtmfs')
        self.assertTrue(os.path.exists(self.mname), "File {0} did not get regridded back!".format(self.mname))
        s = self.get_shape(self.mname)
        self.assertEqual([s[0], s[1]], [10, 10], "Image shape ({0}) didn't get regridded as expected ({1})".format(s, [10, 10]))

##############################################
##############################################

##Task level tests : verify that the .residual image gets updated between consecutive runs of deconvolve
class test_residual_update(testref_base):
    def get_stats(self, statsmap, ismtmfs, ttnrange, ttnext):
        for ext in statsmap.keys():
            for ttn in range(ttnrange):
                imgname = self.img+"."+ext+( ttnext.format(ttn) )
                _ia.open(imgname)
                stats = _ia.statistics()
                statsmap[ext].append({
                    'max':  stats['max'][0],
                    'mean': stats['mean'][0],
                    'min':  stats['min'][0]
                })
                _ia.close()
        return statsmap

    def helper_residual_update(self, deconvolver, tclean_args={}):
        """Helper method to execute the residual update tests for non-multiterm images
           deconvolver (string):  name of the deconvolver to execute
        """
        # these all have to be here for the mtmfs deconvolver
        ismtmfs = False
        ttnrange = 1
        ttnext = ""
        ttnstr = ""
        if deconvolver is 'mtmfs':
            ismtmfs = True
            ttnrange = 2
            ttnext = ".tt{0}"
            ttnstr = ", tt{0}"

        # copy the measurement set and execute the first part of tclean (major cycle only)
        tca={'vis':'refim_eptwochan.ms', 'imagename':self.img, 'imsize':10, 'cell':'8.0arcsec', 'niter':0, 'deconvolver':deconvolver}
        for k in tclean_args.keys():
            tca[k] = tclean_args[k]
        tca=self.prepData(tca['vis'], tclean_args=tca)
        tclean_stats = self.get_stats({'residual':[]}, ismtmfs, ttnrange, ttnext, )
        delta=0.00001

        # execute a single deconvolve and compare to tclean
        deconvolve(imagename=self.img, deconvolver=deconvolver, niter=10)
        decon_stats1 = self.get_stats({'residual':[], 'model':[]}, ismtmfs, ttnrange, ttnext)
        for ttn in range(ttnrange):
            ts  = tclean_stats['residual'][ttn]
            ds1 = decon_stats1['residual'][ttn]
            msg = "No difference exists between initial tclean residual (dirty image) statistics, and residual for deconvolve:\n" +\
                  "Tclean residual stats:     {0}\nDeconvolve residual stats: {1}".format(ts, ds1)
            self.assertNotAlmostEqual(ds1['max'],  ts['max'],  delta=delta, msg=msg)
            self.assertNotAlmostEqual(ds1['min'],  ts['min'],  delta=delta, msg=msg)
            self.assertNotAlmostEqual(ds1['mean'], ts['mean'], delta=delta, msg=msg)

        # execute deconvolve again and compare to first results
        deconvolve(imagename=self.img, deconvolver=deconvolver, niter=10)
        decon_stats2 = self.get_stats({'residual':[], 'model':[]}, ismtmfs, ttnrange, ttnext)
        for ttn in range(ttnrange):
            # residual update check
            ds1 = decon_stats1['residual'][ttn]
            ds2 = decon_stats2['residual'][ttn]
            msg = "No difference exists between first deconvolve residual statistics, and residual for the second deconvolve:\n" +\
                  "First deconvolve stats:    {0}\nSecond deconvolve stats:   {1}".format(ds1, ds2)
            self.assertNotAlmostEqual(ds2['max'],  ds1['max'],  delta=delta, msg=msg)
            self.assertNotAlmostEqual(ds2['min'],  ds1['min'],  delta=delta, msg=msg)
            self.assertNotAlmostEqual(ds2['mean'], ds1['mean'], delta=delta, msg=msg)
            # mdoel update check
            ds1 = decon_stats1['model'][ttn]
            ds2 = decon_stats2['model'][ttn]
            msg = "No difference exists between first deconvolve model statistics, and model for the second deconvolve:\n" +\
                  "First deconvolve stats:    {0}\nSecond deconvolve stats:   {1}".format(ds1, ds2)
            if (ds1['max'] != 0 or ds2['max'] != 0):
                self.assertNotAlmostEqual(ds2['max'],  ds1['max'],  delta=delta, msg=msg)
            if (ds1['min'] != 0 or ds2['min'] != 0):
                self.assertNotAlmostEqual(ds2['min'],  ds1['min'],  delta=delta, msg=msg)
            self.assertNotAlmostEqual(ds2['mean'], ds1['mean'], delta=delta, msg=msg)

    # Test 90
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_residual_update_hogbom(self):
        """ [residual_update] test_residual_update_hogbom """
        ######################################################################################
        # Task deconvolve should update the .residual with every execution. This behavior is
        # left up to each deconvolver. Test that hogbom does this correctly.
        ######################################################################################
        self.helper_residual_update('hogbom')

    # Test 91
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_residual_update_clark(self):
        """ [residual_update] test_residual_update_clark """
        ######################################################################################
        # Task deconvolve should update the .residual with every execution. This behavior is
        # left up to each deconvolver. Test that clark does this correctly.
        ######################################################################################
        self.helper_residual_update('clark')

    # Test 92
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_residual_update_clarkstokes(self):
        """ [residual_update] test_residual_update_clarkstokes """
        ######################################################################################
        # Task deconvolve should update the .residual with every execution. This behavior is
        # left up to each deconvolver. Test that clarkstokes does this correctly.
        ######################################################################################
        self.helper_residual_update('clarkstokes', tclean_args={'vis':'refim_point_linRL.ms','stokes':'I'})

    # Test 93
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_residual_update_multiscale(self):
        """ [residual_update] test_residual_update_multiscale """
        ######################################################################################
        # Task deconvolve should update the .residual with every execution. This behavior is
        # left up to each deconvolver. Test that multiscale does this correctly.
        ######################################################################################
        self.helper_residual_update('multiscale')

    # Test 94
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_residual_update_mtmfs(self):
        """ [residual_update] test_residual_update_mtmfs """
        ######################################################################################
        # Task deconvolve should update the .residual with every execution. This behavior is
        # left up to each deconvolver. Test that mtmfs does this correctly.
        ######################################################################################
        self.helper_residual_update('mtmfs')

##############################################
##############################################

##Task level tests : verify that image restoration can be controlled with the restoration parameters
class test_restoration(testref_base):
    @classmethod
    def setUpClass(cls):
        super(test_restoration, cls).setUpClass()
        msfile='refim_point.ms'
        cls.staticDelData(msfile)
        cls.staticPrepData(msfile, tclean_args={'imsize':10, 'cell':'8.0arcsec', 'deconvolver':'hogbom'})
        cls.staticCopyToCache(msfile, imagename=cls.img, cachedir='restoration_cache')

    def rsetup(self):
        # we can use a cache here because tclean was run only run once, during setUpClass
        self.delData()
        type(self).staticCopyFromCache()
        os.system("rm -rf {}.image".format(self.img))

    # Test 95
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_restoration_none(self):
        """ [restoration] test_restoration_none """
        ######################################################################################
        # Deconvolve but don't restore: should not create a .image image
        ######################################################################################
        self.rsetup()
        results = deconvolve(imagename=self.img, niter=10, restoration=False)
        report = th.checkall(ret=results['retrec'],
                             imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.mask',self.img+'.model'],
                             imgexistnot=[self.img+'.image'])
        self.checkfinal(pstr=report)

    # Test 96
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_restoration_basic(self):
        """ [restoration] test_restoration_basic """
        ######################################################################################
        # Deconvolve and restore, and compare results with those from a tclean run
        ######################################################################################
        self.rsetup()
        results = deconvolve(imagename=self.img, niter=10, restoration=True)
        report = th.checkall(ret=results['retrec'], epsilon=0.005,
                             imgval=[(self.img+'.image',1.060,[5,5,0,0])],
                             imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.mask',self.img+'.model', self.img+'.image'])
        self.checkfinal(pstr=report)

    # Test 97
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_restoration_onlyrestore(self):
        """ [restoration] test_restoration_onlyrestore """
        ######################################################################################
        # Deconvolve and don't restore, then restore and compare results with those from a tclean run
        ######################################################################################
        self.rsetup()
        results = deconvolve(imagename=self.img, niter=10, restoration=False)
        report = th.checkall(ret=results['retrec'], peakres=0.333, imgval=[(self.img+'.model',0.727,[5,5,0,0])], epsilon=0.005,
                             imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.mask',self.img+'.model'],
                             imgexistnot=[self.img+'.image'])
        results = deconvolve(imagename=self.img, niter=0, restoration=True)
        report = th.checkall(ret=results['retrec'], peakres=0.333, imgval=[(self.img+'.model',0.727,[5,5,0,0])], epsilon=0.005,
                             imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.mask',self.img+'.model', self.img+'.image'])
        self.checkfinal(pstr=report)

    # Test 98
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Cube Parallel Output Can't be used. Revisit after CAS-9386")
    def test_restoration_bigbeam(self):
        """ [restoration] test_restoration_bigbeam """
        ######################################################################################
        # Deconvolve and restore with a gigantic beam, then restore and compare results with those from a tclean run with a gigantic beam
        ######################################################################################
        self.rsetup()
        results = deconvolve(imagename=self.img, niter=10, restoration=True, restoringbeam='100.0arcsec')
        report = th.checkall(ret=results['retrec'], epsilon=0.005,
                             imgval=[(self.img+'.image',1.436,[5,5,0,0])],
                             imgexist=[self.img+'.psf', self.img+'.residual', self.img+'.mask',self.img+'.model', self.img+'.image'])
        self.checkfinal(pstr=report)

##############################################
##############################################

##Task level tests : verify that we perform the same number of iterations as tclean for the same iteration parameters
class test_niterparms(testref_base):

    def helper_deconvolve_check_iterdone(self, param_name, param_val, expected_iter, extra_params=None):
        # casalog.post("Executing deconvolve with {}={}, expected iterations: {}".format(param_name, param_val, expected_iter), "WARN")

        # prepare arguments dicts
        ta = {'imsize':100, 'cell':'8.0arcsec', 'deconvolver':'clark', 'threshold':'1mJy', 'gain':0.1}
        da = {'threshold':'1mJy', 'gain':0.1}
        ta[param_name] = param_val
        da[param_name] = param_val
        if extra_params != None:
            for ep_name in extra_params:
                ta[ep_name] = extra_params[ep_name]
                da[ep_name] = extra_params[ep_name]

        # run tclean(niter=0) to get the image ready to work with
        self.prepData('refim_twochan.ms', tclean_args=ta)

        # run deconvolve
        results2 = deconvolve(imagename=self.img, deconvolver=ta['deconvolver'], niter=100, interactive=0, restoration=False, **da)

        # verify results
        report = th.checkall(ret=results2['retrec'], iterdone=expected_iter)
        return report

    # Test 99
    def test_niterparms_gain_1(self):
        """ [niterparms] test_niterparms_gain_1 """
        ######################################################################################
        # Deconvolve should execute 10 iterations for gain=0.2, just like the first major-minor cycle of tclean.
        ######################################################################################
        report = self.helper_deconvolve_check_iterdone('gain', param_val=0.2, expected_iter=10)
        self.checkfinal(report)

    # Test 100
    def test_niterparms_gain_2(self):
        """ [niterparms] test_niterparms_gain_2 """
        ######################################################################################
        # Deconvolve should execute 7 iterations for gain=0.3, just like the first major-minor cycle of tclean.
        ######################################################################################
        report = self.helper_deconvolve_check_iterdone('gain', param_val=0.3, expected_iter=7)
        self.checkfinal(report)

    # Test 101
    def test_niterparms_threshold_1(self):
        """ [niterparms] test_niterparms_threshold_1 """
        ######################################################################################
        # Deconvolve should execute 16 iterations for threshold=0.22, just like the first major-minor cycle of tclean.
        ######################################################################################
        report = self.helper_deconvolve_check_iterdone('threshold', param_val=0.22, expected_iter=16)
        self.checkfinal(report)

    # Test 102
    def test_niterparms_threshold_2(self):
        """ [niterparms] test_niterparms_threshold_2 """
        ######################################################################################
        # Deconvolve should execute 19 iterations for threshold=0.18, just like the first major-minor cycle of tclean.
        ######################################################################################
        report = self.helper_deconvolve_check_iterdone('threshold', param_val=0.18, expected_iter=19)
        self.checkfinal(report)

    # Test 103
    def test_niterparms_nsigma_1(self):
        """ [niterparms] test_niterparms_nsigma_1 """
        ######################################################################################
        # Deconvolve should execute 72 iterations for nsigma=0.9, just like the first major-minor cycle of tclean.
        ######################################################################################
        report = self.helper_deconvolve_check_iterdone('nsigma', param_val=0.9, expected_iter=72, extra_params={'maxpsffraction':0})
        self.checkfinal(report)

    # Test 104
    def test_niterparms_nsigma_2(self):
        """ [niterparms] test_niterparms_nsigma_2 """
        ######################################################################################
        # Deconvolve should execute 60 iterations for nsigma=1.5, just like the first major-minor cycle of tclean.
        ######################################################################################
        report = self.helper_deconvolve_check_iterdone('nsigma', param_val=1.5, expected_iter=60, extra_params={'maxpsffraction':0})
        self.checkfinal(report)

    # Test 105
    def test_niterparms_cyclefactor_1(self):
        """ [niterparms] test_niterparms_cyclefactor_1 """
        ######################################################################################
        # Deconvolve should execute 40 iterations for cyclefactor=0.1, just like the first major-minor cycle of tclean.
        ######################################################################################
        report = self.helper_deconvolve_check_iterdone('cyclefactor', param_val=0.1, expected_iter=40)
        self.checkfinal(report)

    # Test 106
    def test_niterparms_cyclefactor_2(self):
        """ [niterparms] test_niterparms_cyclefactor_2 """
        ######################################################################################
        # Deconvolve should execute 13 iterations for cyclefactor=2.0, just like the first major-minor cycle of tclean.
        ######################################################################################
        report = self.helper_deconvolve_check_iterdone('cyclefactor', param_val=2.0, expected_iter=12)
        self.checkfinal(report)

    # Test 107
    def test_niterparms_minpsffraction_1(self):
        """ [niterparms] test_niterparms_minpsffraction_1 """
        ######################################################################################
        # Deconvolve should execute 7 iterations for minpsffraction=0.5, just like the first major-minor cycle of tclean.
        ######################################################################################
        report = self.helper_deconvolve_check_iterdone('minpsffraction', param_val=0.5, expected_iter=7)
        self.checkfinal(report)

    # Test 108
    def test_niterparms_minpsffraction_2(self):
        """ [niterparms] test_niterparms_minpsffraction_2 """
        ######################################################################################
        # Deconvolve should execute 16 iterations for minpsffraction=0.2, just like the first major-minor cycle of tclean.
        ######################################################################################
        report = self.helper_deconvolve_check_iterdone('minpsffraction', param_val=0.2, expected_iter=16)
        self.checkfinal(report)

    # Test 109
    def test_niterparms_maxpsffraction_1(self):
        """ [niterparms] test_niterparms_maxpsffraction_1 """
        ######################################################################################
        # Deconvolve should execute 78 iterations for maxpsffraction=0.01, just like the first major-minor cycle of tclean.
        ######################################################################################
        report = self.helper_deconvolve_check_iterdone('maxpsffraction', param_val=0.01, expected_iter=78)
        self.checkfinal(report)

    # Test 110
    def test_niterparms_maxpsffraction_2(self):
        """ [niterparms] test_niterparms_maxpsffraction_2 """
        ######################################################################################
        # Deconvolve should execute 40 iterations for maxpsffraction=0.05, just like the first major-minor cycle of tclean.
        ######################################################################################
        report = self.helper_deconvolve_check_iterdone('maxpsffraction', param_val=0.05, expected_iter=40)
        self.checkfinal(report)

if __name__ == '__main__':
    unittest.main()