
/***  SMARTPOP, Simulating Mating Alliances as a Reproductive Tactic for Popuplations
      Copyright (C) 2015 E.G.Guillot, M.P.Cox

      This program is free software: you can redistribute it and/or modify
      it under the terms of the GNU General Public License as published by
      the Free Software Foundation, either version 3 of the License, or
      (at your option) any later version.

      This program is distributed in the hope that it will be useful,
      but WITHOUT ANY WARRANTY; without even the implied warranty of
      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
      GNU General Public License for more details.

      You should have received a copy of the GNU General Public License
      along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
      Contact: elza.guillot@gmail.com
***/
#ifndef PROGRAM_OPTIONS2_H
#define PROGRAM_OPTIONS2_H

#include <boost/program_options.hpp>
#include "boost/bind.hpp"

#include  "simulation.h"
#include "dna.h"

namespace po = boost::program_options;

DNA parse_mycommand_line(int argc,char* argv[],simuParam  & param);

#endif
