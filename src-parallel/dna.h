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

#ifndef DEF_GENETICS
#define DEF_GENETICS

#include <cassert>
#include "metapopulation.h"

typedef std::vector< typeDNA* > array_type;

class DNA
{
 public:
  // *** constructor destructor *
  DNA(); /** constructor */
  /** create a DNA instance from a human instance. This constructor will produce an array containing the mtDNA of the whole population*/
  explicit DNA(Population &pop);
  DNA(Population & pop, int tdna); /** Create the DNA population for 1 chromosome \n WARNING we should always care about pairs of chromosome in the case nuclear chromosome */
  DNA(Population  & pop,int tdna,int loci); /* extract a loci from X chromosome or autosome */
  DNA(Metapop & pop, int tdna); /** Create the DNA population for 1 chromosome \n WARNING we should always care about pairs of chromosome in the case nuclear chromosome */
  DNA(Metapop  & pop,int tdna,int loci); /* extract a loci from X chromosome or autosome */
  DNA(std::vector<Human> & indivList, int nHuman, int nFemale,int tdna); // females are first in the list
  DNA(std::vector<Human> & indivList, int nHuman, int nFemale,int tdna,int loci); // females are first in the list
  //  DNA(typeDNA * input,int size);
  DNA(array_type const  & input,int sizeN, int sizeX);
  DNA(std::vector< std::vector< typeDNA > > input);
  DNA(DNA const & a); /** Constructor of copy */ 
  //  DNA(DNA const & a, DNA const & b); /** mMrge two DNA instances.\n For example in the case of autosome we need to consider both copy of the autosome for each human. SMARTPOP extracts the first copy, then the second and then merge both into a new instance. */
  //  DNA(typeDNA* input);
  ~DNA();
  friend std::ostream &operator<<(std::ostream&,DNA const&);
  //  friend DNA operator+(DNA const&,DNA const&); /** Another way to merge two instances */
  DNA& operator=(DNA const&);
  int get_N() const; /** Return the number of sequences */
  int get_nbXcell() const; /** Return the number of cell that represent one sequence (number of sites divided by 32)*/
  int get_nbS() const;   /** Return the number of sites per sequence */
  int get_nbSreal() const;
  void set_nbSreal(int a);
  long poly_sites() const; /** Number of polymorphic sites in the DNA array */
  long poly_sites2(std::vector<int> tsfs) const; /** Number of polymorphic sites in the DNA array */
  double prop_poly_sites() const; /** Proportion of polymorphic sites : number of polymorphic sites in the DNA array divived by the total number of polymorphic sites */
  double mean_pw_diff() const; /** Mean pairwise difference */
  double mean_pw_diff2(std::vector<int>) const; /** Mean pairwise difference */
  double mean_pw_diff_between(DNA & dna2) const; /** Mean pairwise difference */
  double Nei_D(DNA & dna2) const; /** Nei distance 1972*/
  double privateSNP(DNA & dna2) const;
  double Chord_D(DNA & dna2) const; /** Cavalli-Sforza distance 1967 */
  double RWC_D(DNA & dna2) const; /** Reynold's Weir Cochran */
  double theta_waterson() const; //** Return Watterson's Theta for the population */
  double theta_waterson(double polysites) const; //** Return Watterson's Theta for the population */
  double theta_hom() const; //** Return a theta computed from the homozygiosity */
  double theta_pi() const; /** Return theta Pi */
  int nb_of_haplotypes() const; /** Number of haplotypes. An haplotypes corresponds to a define sequences of sites, any polymorphism between two sequences will create a new haplotype. */
  DNA haplotypes() const; /** Array for the aplotypes found in the population */
  std::vector<double> hap_frequency() const; /** Array of the frquency of each haplotypes */
  std::vector<double> site_frequency(int site) const; /** Frequency of each nucleotide found at a given site */
  std::vector<int> site_frequency_int(int site) const; /** Frequency of each nucleotide found at a given site */
  double heterozigosity_hap() const; /** Return the allele diversity */
  double heterozigosity_sites() const; /** Return the sequence diversity */
  double heterozigosity_sites2(std::vector<int> tsfs) const; /** Return the sequence diversity */
  double shannon_index() const; /** Return the sequence diversity */
  double Neff(double theta, double mu) const; /** Return the effective population size, computed from a given theta and a given mutation rate */
  double tajima_D(long,long,long) const; /** Return Tajma's D*/
  std::vector<double> site_frequency_DNA(int site) const; /** Frequency of each nucleotide A C T G at a given site. A contrario de site_frequency it will return 0 for nucletide not found, therefore is always of size 4 */
  std::vector<double> site_frequency_DNA_crumb(int site) const; /** Frequency of each nucleotide A C T G at a given site. A contrario de site_frequency it will return 0 for nucletide not found, therefore is always of size 4 */
  std::vector<int> site_frequency_DNA_int(int site) const; /** Frequency of each nucleotide A C T G at a given site. A contrario de site_frequency it will return 0 for nucletide not found, therefore is always of size 4 */
  typeDNA dist_two(int a, int b) const; /** Distance between two given sequence. Number of sites that differ between both */
  typeDNA dist_two(int a, int b, DNA & dna2) const; /** Distance between two given sequence. Number of sites that differ between both */
  std::vector<double> output_diversity_with_time(std::string,int Ntheo, int timeG); /** Output diversity values in a file with time */
  void output_diversity_with_time(std::string,int Ntheo, int timeG,double pmate, double pmig) ; /** Output diversity values in a file with time */
  void output_diversity_with_time_nosimu(std::string,int Ntheo); /** Output diversity values in a file with time */
  //  void output_diversity_with_time_fast(std::string,int Ntheo, int timeG) const; /** Output diversity values in a file with time */
  //  void output_diversity_with_time(DNA & dna2,std::string,int Ntheo, int timeG) const; /** Output diversity values in a file with time */
  std::vector<double> output_diversity(DNA & dna2,std::string,int Ntheo, int timeG); /** Output diversity values in a file with time */
  void output_diversity(std::string  filename, int Ntheo) ; /** Output diversity in a file with the populaton size. In case of a computation of diversity from a sample, this will be the real population size */
  void output_arlequin(std::string filename) const; /** Output the sequences in a file corresponding to the format of arlequin version 3.5 */
  void output_haplotypes(std::string filename, int timeG);
  void output_fasta(std::string filename) const; /** Output the sequences in a file in FASTA format */
  void output_nexus(std::string filename) const; /** Output the sequences in a file in NEXUS format */
  void SFS(std::string filename) const;
  std::vector<int> SFS() const;
  void SFS_folded(std::string filename) const;
  std::vector<double> SFS_folded_dbl(std::vector<int> tsfs) const;
  std::vector<int> SFS_folded() const;
  std::vector<double> MFS(std::vector<int> msfs) const;
  //  std::vector<double> site_frequency_DNA(int site) const;
  void SFS_DNA(std::string filename, int generation) const;
  void AFS_DNA(std::string filename, int generation) const;
 private:
  std::ostream & print(std::ostream&) const;
  int N; /** Number of sequences */
  int nbXcell; /** Number of cells per sequence (number of sites divided by 32) */
  int nbS; /** Number of sites per sequence */
  int nbSreal; /** Number of sites per sequence */
  typeDNA ** tab; /** Array containing all the sequences of size nbXcell times N with each cell being a 64 bits number representing 32 sites */
  double hh; /* heterozygosity haplotype*/
  double pw; /* mean pairwise difference */
  double tw; /* mean pairwise difference */
  int ps; /* nb of polymorphic sites*/
  double hs; /* heterozygosity sites*/
  int nbh; /* nb of haplotypes */
};
std::vector< std::vector<typeDNA> >  input_fasta(std::string filename);
std::vector< std::vector<typeDNA> >  input_ped(std::string filename);
std::vector< std::vector<typeDNA> >  input_pedY(std::string filename);
int  input_ped_nbS(std::string filename);
int  input_ped_nbSY(std::string filename);
void output_diversity_headline(std::string); /** Output the headline in a file where the diversity will be written */
void output_diversity_headline_with_time(std::string); /** Output the headline in a file where the diversity will be written, in case we write the time with it */
void output_diversity_headline_with_time(std::string,std::string,std::string); /** Output the headline in a file where the diversity will be written, in case we write the time with it */
void output_diversity_headline_nosimu(std::string); /** Output the headline in a file where the diversity will be written, in case we write the time with it */
void output_diversity_headline_with_time2(std::string); /** Output the headline in a file where the diversity will be written, in case we write the time with it */
double ssq_sfs(  std::vector<double>, std::vector<double>);
double sabs_sfs(  std::vector<double>, std::vector<double>);
std::vector<double> resize_sfs(std::vector<double> tsfs) ;
std::vector<double> resize_simple_sfs(std::vector<double> tsfs);

#endif
