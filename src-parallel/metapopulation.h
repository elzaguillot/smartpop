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

#ifndef METAPOP_H 
#define METAPOP_H
#include "utils.h"
#include "population.h"
#include <string.h>
#include <stdlib.h>
#include <cstdlib>
#include <vector>
class Metapop
{
 public:
  // constructor & destructor
  Metapop();
  Metapop(int,int); // nb populations nb individuals
  Metapop(int npop,int nindiv,std::vector<int> niList); // nb populations nb individuals
  Metapop(Metapop const & );
  Metapop(std::string & filename );
  ~Metapop();
  Metapop& operator = (Metapop const & );
  friend  std::ostream & operator<<(std::ostream & os, Metapop const & );
  std::ostream& print(std::ostream &) const;
  void set_nbPop(int npop, std::vector<int> nindiv);
  void set_nbSources(int);
  void Init(int nstep);
  void Init_succ(int nstep);
  void Init_copy(const Metapop & a);
  void reset(int, std::vector<int>);
  // Evolution function
  void evolve(int nbGene,bool verbose);
  void evolve_one();
  void burnin(double diversityThreshold,int gen,double facMu,int dnaType); // output the number of gen run under the burnin
  void mutation();
  void migration();
  void reproduction();
  void reproduction_burnin();
  void update_netmating();
  void demographic_update();
  // Get / set
  void set_mu(double,double,double,double,double,double,double,double); /** Set mutation rate : mtDNA transition, mtDNA transversion, X transition, X transversion, Y transition, Y transversion, austomomal transition, autosomal transversion */
  void set_pmate(double a);
  void set_pmig(double a);

  // Input Output
  std::vector<Human> sample_local(int nindiv,int idpop) const;
  std::vector<Human> sample_local_random(int nindiv) const;
  std::vector<Human> sample_scattered(int nindiv) const;
  std::vector<Human> sample_pooled(int nindiv,int npop) const;
  std::vector<Human> sample(int nindiv,int npop) const;// nb of sample + nb of pop
  std::vector<Human> sample(int nindiv,int npop,int popID) const;// nb of sample + nb of pop
  std::vector<Human> sample(int nindiv, int npop, std::vector<double> prop) const; // nb of sample + nb of pop + proportion in each pop
  std::vector<Human> sample(int nindiv, std::vector<int> npop,std::vector<double> prop) const; // nb of sample + IDs of pop + proportion in each pop
  std::vector<Human> sample(int nindiv, std::vector<int> npop,std::vector<double> prop,std::vector<double> propf) const; // nb of sample + IDs of pop + proportion in each pop + prop of females
  void output_pop_size() const; /** Output the population size in the terminal */
  void write_smartFormat(std::string const& filename) const; /** Output the population in a file that can be realoaded by SMARTPOP. This procedure is used to model complex scenario */ 

  //  int popsize;
  //  int nbFemale;
  int migrationType; // 0 endogamy (no migration), 1 using a migration rate or 2 using wife givers/wife takers
  int nbPop;
  int popsize;
  int nbFemale;
  int generation;
  double pmate; // probability of marrying someone random and not cousin
  double pmig; // probability  of changing patrilineal destination
  int nbSources;
  std::vector<int> popID;
  std::vector<int> maxpopsizes;
  std::vector< std::vector<double> > migrationRateF;
  std::vector< std::vector<double> > migrationRateM;
  std::vector< std::vector<int> > nbMigrantsF;
  std::vector< std::vector<int> > nbMigrantsM;
  MSNetwork netmating;
  Population ** poptab;
  std::vector<std::vector<Human*> > female;
  std::vector<std::vector<Human*> > male;
  //  int* humanIDNet;
  //  std::vector<int> xcoor;
  //  std::vector<int> ycoor;
  void create_circleMigration();
  void create_APA();
  void create_SPA();
  void split_migrate(int popID1, int popID2);
  void move_indiv(int pop1, int id1, int pop2, bool sex);
};


#endif
