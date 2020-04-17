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

#ifndef DEF_SIMULATION
#define DEF_SIMULATION
#include "dna.h"
#include "io_interface.h"

class Simulation
{
 public:
  Simulation(); /** Constructor by default. This should not be called in the normal set of simulations */
  explicit Simulation(unsigned randomSeed); /** Constructor from the random seed. This constructor is the one called in the normal set of simulations*/
  Simulation(Simulation const &); /** Constructor of copy */
  ~Simulation(); /** Destructor */
  Simulation& operator = (Simulation const&);
  // main operations
  int run(); /** Main function. It is called for every set of simulations after the parameters have been set.  */
  int print_help(); /** Output the help in command lines. It is called by default if the simulation set an not start */
  void load_file(std::string); /** Load a file of parameters to start a simulation. This method is an alternative to the inoput of parameters in command lines */
  void load_directory(std::string); /** Load a directory containing the data from a previous simulation set. WARNING : for parallel processes the new call must start the same number of process that the previous one! */
  void export_file(std::string filename); /** Output a file with all the parameters set for the simulations. By default it will be called smartpop_parameters.txt */
  //  get /set
  void set_seed(unsigned int); /** Set the random seed */
  void set_popsize(int psize); /*Set the population size */
  void set_popsize(int psize,int popID); /*Set the population size */
  void set_nbPop(int);
  void set_nbSources(int);
  void set_sample(int); /**  Set the sample size */
  void set_nbsamplepop(int); /**  Set the sample size */
  void set_sampleType(int); /**  Set the sample size */
  void set_sampleBetween_on();
  void set_step(int); /** Set the number of generation between two sampling */
  void set_nstep(int); /** Set the number of sampling */
  void set_matsys(int);
  void set_matsys(int a,int popID);
  void set_inbreeding(int); /** Set inbreeding scheme : \n 0 all accepted \n 1 alliances forbidden between siblings */
  void set_inbreeding(int a,int popID); /** same for a given population */
  void set_varChild(double); /** Set the variance of children per female*/
  void set_varChild(double a,int popID);
  void set_polygamy(int); /** Set inbreeding scheme : \n 0 all accepted \n 1 alliances forbidden between siblings */
  void set_polygamyall(int); /** Set inbreeding scheme : \n 0 all accepted \n 1 alliances forbidden between siblings */
  void set_polygamy(int a,int popID); /** same for a given population */
  void set_polygamyf(int); /** Set inbreeding scheme : \n 0 all accepted \n 1 alliances forbidden between siblings */
  void set_polygamyf(int a,int popID); /** same for a given population */
  void set_nSimu(int); /** Set the number of simulation to be run */
  void set_sizeMt(int);
  void set_sizeY(int);
  void set_sizePerX(int);
  void set_sizePerA(int);
  void set_nbLociA(int);
  void set_nbLociX(int);
  void set_mu(double a,double b, double c, double d);
  void set_mu(double,double,double,double,double,double,double,double); /** Set 8 mutation rates : mtDNA transition, mtDNA transversion, X transition, X transversion, Y transition, Y transversion, autosomal transition, autosomal transversion **/
  void set_growthFactor(double,double,double); /** Set the demography parameters. The demography is defined a as discrete process:\n p(t+1) = a + b*p(t) + c*p(t)*p(t)\n For a constant population: a=0, b=1 , c=0\n For aan exponentially growing (/decreasing) population : a=0, b>1 (b<1), c=0 \n For a constant growth (decrease): a>0, b=1, c=0*/
  void set_growthFactor(double,double,double,int popID); // set for one pop only
  void set_burninT(double); /** Set the diversity treshold for the burnin phae */
  void set_burninDna(int); /** Set the diversity type for the burnin phae */
  void set_diversity(int a); /** Set on which DNA type the diversity should be computed :\n 0 all DNA types (mtDNA, X, Y, autosomes) \n 1 mtDNA \n 2 X\n 3 Y \n 4 autosomes */
  void set_migrationType(int );
  void set_pmate(double);
  void set_pmig(double);
  void set_fileoutput(std::string); /** Set the name of the file for the outputs */
  void set_verbose_on(); /** Change the parameter verbose to on. The program will then output a lot of information in the terminal */
  void set_fasta_on(); /** Create a Fasta file for each population at the end of the simulation */
  void set_save_on(std::string dirName); /** Save the all the populations at the end of the simulatin in a SMARTPOP format */
  //  void set_save_xml();
  void set_arlequin_on(); /** Create a Arlequin file for each population at the end of the simulation */
  void set_haplotype_on();
  void set_sfs_on();
  void set_APA_mig();
  void set_header_on(); /** Set the parameter to write the header at the top of the diversity file */
  void set_migrationRate(int pop1,int pop2, bool sex,double rate);

  int get_nSimu(); /** Return the number of simulations */
  bool get_verbose(); /** Return the verbose value (true or false) */
  unsigned int get_seed();
  int get_popsize();
  int get_nbPop();
  int get_popsize(int popID);
  int get_sample();
  int get_sample(int popID);
  int get_inbreeding(int popID);
  double get_varChild(int popID);
  int get_polygamy(int popID);
  int get_polygamyf(int popID);
  int get_step();
  int get_nstep();
  std::vector<int> get_matsys(int popID);
  double* get_growthFactor(int popID);
  double get_burninT();
  int get_burninDna();
  std::vector<int> get_diversity();
  std::string get_fileoutput();
  bool get_fasta();
  bool get_save();
  bool get_arlequin();
  bool get_header();

  void header_write(); /** Write the header of the diversity file*/
  void check_parameters(); /** Check that parameters chosen are OK before starting the set of simulations */
 private: // ATTRIBUTES
  simuParam simuP;
  int counter; /** Increasing counter of the number of simulation that have run */  
  std::vector< struct popVal > popValues;
  struct metapopVal metapopValues;
  double maxpop; /** The maximum value that the population size will reach during the simulation */
  // hidden functions
  int classic(); /** Run for a set of simulations with the whole DNA set */
  int successive_classic(); /** Run for a set of simulations with the whole DNA set, starting from population files */
  int output_main(Metapop & epsilon); /** Output diversity in a file */
  void output_sfs(std::vector<Human> &,int,int);
  void output_afs(Metapop & epsilon,int);
  void output_haplotype(Metapop & epsilon, int );
  int verbose_start(); /** Output in the terminal lines all the parameters chosen in the case where verbose has been set on */
  int verbose_end(); /** Output in the terminal all the informations at the end of the set of simulations*/
  void save(Metapop & epsilon); /** Save a population in a file */
  void evolve(Metapop & epsilon); /** Make a population evolve and output its differents bits */
};



#endif
