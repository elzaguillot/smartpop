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

#ifndef DEF_POPULATION_H
#define DEF_POPULATION_H
#include <iterator>
#include <algorithm>
#include "human.h"
#include "utils.h"
#include "demographic_function.h"
#include "boost/utility/enable_if.hpp"
#include "boost/type_traits/is_base_of.hpp"
//#include "matingSystem.h"
#include "msnetwork.h"


/* Population is a class that contains  2 array of human: female and male

// this population evolve through non overlapping generation with cycles of Mutation - Reproduction*/

class Population //
{
  // private:
 public:
  // *** Public attributes ***
  double popsizeDemo; /** Size of the population to reach at generation + 1. The demography function update this number at each generation. The number of children to make is this number which may be different from the population size a the case of a growth or decrease */
  Human tmpmale; // handy and fast
  Human tmpfemale; // handy and fast
  std::vector<int> matingSystem; /** Number representing the mating system. **/
  int polygamy; // nb of wife allowed per male
  int polygamyf; // nb of husband allowed per female
  int inbreeding;  /** If 1, alliances are forbidden between siblings */ // 1 : no inbreeding between brothers and sisters
  double varChild;
  double growthFactor[3]; // exponential growth factor
  /** Arrays of 3 doubles [a,b,c] that represent the demography function. The function is p(t+1) = a+  b*p(t) + c*p(t)*p(t)*/

 public:
  int generation; /** Counter of generation. It is increased by one after each evolution cycle */
  int popsize; /** Number of individuals in the population */
  int nbFemale; /** Number of females in the population */
  int maxsize;
  int offsetNetwork;
  int migrationType;
  std::vector< std::vector<  Human* > > female; /** Vector containing all the female of the current generation, and of the previous (in phase before reproduction) or the next generation (in phases after reproduction)*/
  std::vector< std::vector< Human* > > male; /** Vector containing all the female of the current generation, and of the previous (in phase before reproduction) or the next generation (in phases after reproduction)*/

 public:
  // *** Constructor / destructor ***  
  Population();
  explicit Population(int N);
  Population(Population const& a);
  //  explicit Population(std::string const& filename); /** Constructor of a population from a file in SMARTPOP format */
  ~Population(); /** Destructor */
  /*** Init - shuffle  ***/
  void random_shuffle_pop(); // handwritten shuffling of indiv of current generation: fast handy and efficient
  /** Shuffling of the arrays containing the individuals in the population following Knuth's algorithm */
  void initiate_pop(int N);
  /** Initialize all the values of the populations to defaults values for a given population size */
  void reset(int N);
  
  // *** std get/set ***

  Population& operator = (Population const&);
  friend std::ostream& operator << (std::ostream& os, Population const&a)
    {
      a.print(os);
      return os;
    }
  void set_matingSystem(std::vector<int>); /** Change mating system */
  std::vector<int> get_matingSystem() const; /** Return mating system */
  void set_inbreeding(int); /** Change inbreeding */
  int get_inbreeding() const; /** Return inbreeding */
  void set_varChild(double); /** Change inbreeding */
  double get_varChild() const; /** Return inbreeding */
  void set_polygamy(int); /** Change polygamy */
  int get_polygamy() const; /** Return polygamy */
  void set_polygamyf(int); /** Change polygamy */
  int get_polygamyf() const; /** Return polygamy */
  void set_growthFactor(double*); /** Change the values a,b,c of the demogrphy function. p(t+1) = a + b*p(t) + c*p(t)*p(t) */
  void set_growthFactor(std::vector<double>); /** Change the values a,b,c of the demogrphy function. p(t+1) = a + b*p(t) + c*p(t)*p(t) */
  void resize(int); /** Resize the arrays containing the humans. WARNING to do once only if posisble as time consuming */
  int get_generation() const; /** Return the number of generation since the begining */
  int get_maxpop(int nstep) const;
  int get_popsize() const; /** Return the current population size */
  int get_nbFemale() const; /** Return the current female population size */
  int get_nbMale() const; /** Return the current male population size */
  double get_popsizeDemo() const; /** Return the population size for the children */
  std::vector<double> get_growthFactor() const; /** Return demographic parameters */
  std::vector<int> get_nb_mates() const; /** Return a vector containing the number of mates per indiv */
  std::vector<int> get_nb_children() const; /** Return a vector containing the number of children per indiv */
  // *** Sub sample ***
  Population subsample(int sampleSize, int nbWomen); /** Return a sample of the population (of type Population aswell) for a given sample size, with a given number of female in this sample*/
  Population subsample_female(); /** Return a sample of the population (of type Population aswell) with the female only */ // fast and easy, if want to subsample a # of women use subsample(#,#)
  Population subsample_male(); /** Return a sample of the population (of type Population aswell) with the male only */ // idem subsample(#,0)
  void analyse(); /** Output different estimators of diversity in a file + fasta file to load in arlequin */ // TODO

  // *** Outputs ***
  std::ostream& print(std::ostream&)const;
  void output_pop_size() const; /** Output the population size in the terminal */
  // private:
  //  void write_smartFormat(std::string const& filename) const; /** Output the population in a file that can be realoaded by SMARTPOP. This procedure is used to model complex scenario */ 
  void output_diversity();
  void output_nb_child(std::string filename,bool) const; /** Output the number of child for each human in a file */
  void output_nb_mates(std::string filename,bool) const; /** Output the number of mates for each human in the command line */
  void test(int); // dummy but <> can test on humanmt ror humanddna some stuff

  // *** Diversity estimators ***
  typeDNA polySites(int nuclear) const; /** Number of segregating (polymorphic) sites for a given chromosome */ 
  typeDNA polySites_mtdna() const; /** Number of segregating (polymorphic) sites on the mtDNA */ 
  double pairwise_distance_mtdna() const; /** Mean pairwise distance on the mtDNA in the population (used for burnin phase) */ 
  double pairwise_distance_Y() const; /** Mean pairwise distance on the mtDNA in the population (used for burnin phase) */ 
  double pairwise_distance_X() const; /** Mean pairwise distance on the mtDNA in the population (used for burnin phase) */ 
  double pairwise_distance_A() const; /** Mean pairwise distance on the mtDNA in the population (used for burnin phase) */
  double diversity_A() const; /** Proportion of polymorphic sites on the autosomes. */ 
  double diversity_Y() const; /** Proportion of polymorphic sites on the Y chromsome . */ 
  double diversity_X() const; /** Proportion of polymorphic sites on the X chromosome. */ 
  double diversity_mtdna() const; /** Proportion of polymorphic sites on the mtDNA. */ 

  // *** Evolution functions **
  void mutation(); /** Phase in which all the human mutates according to the mutation scheme */
  void reproduction(MSNetwork & netmating,int populationID,double pmate, double pmig); /** Phase in which inidividuals mate and produce the next generation */
  void reproduction_burnin(); /** Phase in which inidividuals mate and produce the next generation */
  void reproductionNew(MSNetwork & netmating,int populationID);
  void birth_boy(int indexfemale,int indexmale,int nbBabyMale,int indexchild,MSNetwork & netmating); /** Make a new baby male out of the DNA two  parents */  
  void birth_girl(int indexfemale,int indexmale,int nbBabyFemale,int indexchild,MSNetwork & netmating); /** Make a new baby female out of the DNA of the two parents */ 
  void birth_boy_burnin(int indexfemale,int indexmale,int nbBabyMale); /** Make a new baby male out of the DNA two  parents */  
  void birth_girl_burnin(int indexfemale,int indexmale,int nbBabyFemale);
  void demographic_update(); /** Update popsizedemo to follow the demographic scheme */
  int burnin(double diversityTreshold,int gen,double facMu,int dnaType, MSNetwork & netmating,int popID);  // output the number of gen run under the burnin
  void update_netmating(MSNetwork & netmating, int popID);
  void set_maleIDnetwork(MSNetwork & netmating,int popID);
  void set_femaleIDnetwork(MSNetwork &  netmating, int popID);
};



#endif
