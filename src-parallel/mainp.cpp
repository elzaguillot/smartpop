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
    
    Contact: e.guillot@massey.ac.nz
***/

#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include "simulation.h"
#include "mpi.h" // must have mpi.h in the path
#include "myprogram_options.h"
#include "boost/uuid/seed_rng.hpp"

using namespace std;

// -h help -v verbose


int main(int argc, char* argv[])
{
  int id, p;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  unsigned seed = boost::uuids::detail::seed_rng()();// time(0)*(id+1); PB parallel seed cant be time(0)
  // see http://stackoverflow.com/questions/12779235/how-to-properly-choose-rng-seed-for-parallel-processes
  std::string expFile="SMARTPOP_parameters.txt";
  Simulation sim(seed);
  parse_mycommand_line(argc,argv,sim,p,id);
  int nbSimu;
  nbSimu = sim.get_nSimu();
  int nbsimuAdjustement=0;
  if((nbSimu%p)/(id+1)>=1)
    nbsimuAdjustement=1;
  nbSimu=nbSimu/p+nbsimuAdjustement;// second part bigger than 1 for as many process as difference
  sim.set_nSimu(nbSimu);
  if(sim.get_verbose() && (id==0))
    std::cout << " Proc " << id << " will compute " << nbSimu << " simulations" << std::endl;
  std::cout << " Proc " << id << " will compute " << nbSimu << " simulations" << std::endl;
  sim.run();
  MPI_Finalize();
  return 0;
}
  
