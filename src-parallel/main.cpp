/***  SMARTPOP, Simulating Mating Alliances as a Reproductive Tactic for POPuplations
      Copyright (C) 2015 E.G.Guillot

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

#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include "simulation.h"
#include "io_interface.h"// smartpop -h for help
#include "myprogram_options.h"
#include "boost/uuid/seed_rng.hpp"

using namespace std;

// NB : mu must be /generation


int main(int argc, char* argv[])
{
  unsigned seed = boost::uuids::detail::seed_rng()();// time(0)*(id+1); PB parallel seed cant be time(0)
  Simulation sim(seed);
  parse_mycommand_line(argc,argv,sim,1,0);
  srand(sim.get_seed());
  sim.check_parameters();
  sim.run();
  return 0;
}

