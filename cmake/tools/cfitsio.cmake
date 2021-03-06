#
#  $Id: cfitsio.cmake 170121 2015-09-10 21:35:14Z nega $
#  
#  Copyright (C) 2007   Troy D. Straszheim  <troy@icecube.umd.edu>
#  and the IceCube Collaboration <http://www.icecube.wisc.edu>
#  
#  This file is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>
#  
TOOLDEF (cfitsio
  include
  fitsio.h
  lib
  NONE
  cfitsio
  )

## for Fedora...
if(NOT CFITSIO_FOUND)
TOOLDEF (cfitsio
  include
  cfitsio/fitsio.h
  lib
  NONE
  cfitsio
  )
if(CFITSIO_FOUND)
  set(CFITSIO_INCLUDE_DIR ${CFITSIO_INCLUDE_DIR}/cfitsio CACHE FILEPATH "cfitsio header path" FORCE)
endif(CFITSIO_FOUND)
endif(NOT CFITSIO_FOUND)
