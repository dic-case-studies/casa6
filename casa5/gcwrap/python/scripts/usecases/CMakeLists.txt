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

casa_add_python( tools tools_usecases_python lib/python${PYTHONV}/regressions
  #jupiter6cm_usecase.py
  ngc5921_usecase.py
  run_polcal_20080224_cband_usecase.py
  run_usecase_orions.py
  sdnewfeatures_b3_tutorial.py
  test_task_cvel.py
  test_task_exportasdm.py
)
