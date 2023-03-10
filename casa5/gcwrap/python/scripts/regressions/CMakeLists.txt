#
# CASA - Common Astronomy Software Applications
# Copyright (C) 2010 by ESO (in the framework of the ALMA collaboration)
#
# This file is part of CASA.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

casa_add_python( tools tools_regressions_python lib/python${PYTHONV}/regressions
  3C129_regression.py
  3C219D_test_regression.py
  asdm-imexport+analysis_regression.py
  equinox_vis_regression.py
  ic2233_regression.py
  ngc2403_regression.py
  ngc5921redux_regression.py
  pylabmem_regression.py
  wproject_regression.py
  calstat_test_regression.py
  3c129_tutorial_regression.py
  3C219D_regression.py
  3C391polim_regression.py
#  accum_regression.py
  alma-m100-analysis-regression.py
  alma-m100-analysis-hpc-regression.py
  alma-m100-analysis-hpc-regression-tclean.py
  alma-m100-sd-regression.py
  alma-ngc3256-analysis-regression.py
  alma-titan-ephemeris-regression.py
  alma_asdm+reimported-asdm_sf.py
  async-segfault_regression.py
  b0319_regression.py
  cleanchan_regression.py
  cleanhelper_regression.py
  coordsystest_regression.py
  cvel+clean_test_regression.py
  cvel_regression2.py
  fits-import-export_regression.py
  g192_regression.py
  h121_regression.py
  imagepoltest_regression.py
  imagetest_regression.py
  l02d_regression.py
  leo2pt_regression.py
  listcal_regression.py
  listvis_regression.py
  maclog_regression.py
  nf_regression.py
  ngc1333_regression.py
  ngc2403_tutorial_regression.py
  ngc4826c_regression.py
  ngc4826_regression.py
  ngc4826_tutorial_regression.py
  ngc4826-mms_tutorial_regression.py
  ngc5921_regression.py
  ngc7538_regression.py
  orion_regression.py
#  plotxy_regression.py
  pointing_regression.py
  polcal_20080224_cband_regression.py
  stokesimagingtest_regression.py
  unit_tests.py
  vla_uvfits_line_sf.py
  wf_regression.py
  wideband_regression.py
  widebandmosaic_regression.py
  wproj3ddat_regression.py
  run_orionmos4sim.py
  m51_3sim_regression.py
  m51_tpsim_regression.py
  ppdisk2_regression.py
  testcube2_regression.py
  testcompsim.py
  pipelineTest_almaif.py
  pipelineTest_regression_almaif.py
  pipelineTest_almasd.py
  pipelineTest_regression_almasd.py
  pipelineTest_vla.py
  pipelineTest_regression_vla.py
  time_then_chan_avg.py
  calstat_test_regression.py
  wideband.py
  wproject.py
  polcal_20080224_cband.py
  alma-m100-analysis.py
  alma-m100-sd.py
  alma-m100-analysis-hpc-tclean.py
  alma-m100-analysis-hpc.py
  alma-ngc3256-analysis.py
  alma-titan-ephemeris.py

)
