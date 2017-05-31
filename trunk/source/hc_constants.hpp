//! \file hc_constants.hpp
//
// Copyright (c) 2017 by Guillaume Oliviéro <goliviero@lpccaen.in2p3.fr>
//
// Half Commissioning constants
// (to be include in Falaise/src at term)
//

#ifndef HC_CONSTANTS_HPP
#define HC_CONSTANTS_HPP

// Standard library:
#include <cstdint>

// Third party:
#include <bayeux/datatools/units.h>
#include <bayeux/datatools/clhep_units.h>

struct hc_constants
{
  static const uint16_t NUMBERS_OF_CALO_PER_COLUMN = 13;
	static const uint16_t NUMBERS_OF_COLUMNS = 20;



};

#endif // HC_CONSTANTS_HPP

// Local Variables: --
// Mode: c++ --
// c-file-style: "gnu" --
// tab-width: 2 --
// End: --
