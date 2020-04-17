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

#include "population.h"

using namespace std;

Population::Population()
{
  initiate_pop(0);
}

Population::Population(int N) /** Constructor from the population size */
{
  initiate_pop(N);
} // assume mu=0

Population::Population(Population const& a) /** Constructor of copy */
{
  initiate_pop(a.popsize);
  generation = a.generation;
  nbFemale = a.nbFemale;
  popsizeDemo = a.popsizeDemo;
  matingSystem = a.matingSystem;
  inbreeding = a.inbreeding;
  varChild = a.varChild;
  polygamy = a.polygamy;
  polygamyf = a.polygamyf;
  female[0] = a.female[0];
  male[0] = a.male[0];  
  female[1] = a.female[1];
  male[1] = a.male[1]; 
  memcpy(growthFactor,a.growthFactor,3*sizeof(double));

  maxsize = a.maxsize;
  migrationType = a.migrationType;
  memcpy(growthFactor,a.growthFactor,3*sizeof(double));
}

Population::~Population()
{
  //  NOTHING
}


void Population::set_matingSystem(std::vector<int> a)
{
  matingSystem = a;
}

std::vector<int> Population::get_matingSystem() const
{
  return matingSystem;
}

void Population::set_inbreeding(int a)
{
  inbreeding = a;
}

int Population::get_inbreeding() const
{
  return inbreeding;
}

void Population::set_varChild(double a)
{
  varChild = a;
}

double Population::get_varChild() const
{
  return varChild;
}

void Population::set_polygamy(int a)
{
  polygamy = a;
}

void Population::set_polygamyf(int a)
{
  polygamyf = a;
}

int Population::get_polygamy() const
{
  return polygamy;
}

int Population::get_polygamyf() const
{
  return polygamyf;
}

void Population::set_growthFactor(double * a)
{
  memcpy(growthFactor,a,3*sizeof(double));
}

void Population::set_growthFactor(std::vector<double> a)
{
  if(a.size()==3)
    {
      growthFactor[0] = a[0];
      growthFactor[1] = a[1];
      growthFactor[2] = a[2];
    }
  else
    {
      cerr << "ERROR setting the demographic function" << endl;
      exit(2);
    }
}

int Population::get_generation() const
{
  return generation;
}

int Population::get_maxpop(int nstep) const
{
  if(!((growthFactor[0]==0)&&(growthFactor[1]==1)&&(growthFactor[2]==0)))
    return max_demo(popsize,growthFactor[0],growthFactor[1],growthFactor[2],nstep,0);
  else
    return popsize;
}

int Population::get_popsize() const
{
  return popsize;
}

int Population::get_nbFemale() const
{
  return nbFemale;
}

int Population::get_nbMale() const
{
  return popsize-nbFemale;
}

double Population::get_popsizeDemo() const
{
  return popsizeDemo;
}

std::vector<int> Population::get_nb_mates() const
{
  int parentG = !(generation&1);
  std::vector<int> mates(popsize);
  for(int i=0; i<nbFemale; ++i)
    mates[i] = female[parentG][i]->getNbMate(); // careful we have to look at the parent generation
  for(int i=0; i<(popsize-nbFemale); ++i)
    mates[i+nbFemale] = male[parentG][i]->getNbMate();
  return mates;
}


std::vector<int> Population::get_nb_children() const
{
  int parentG = (generation&1)^1;
  std::vector<int> children(popsize);
  for(int i=0; i<nbFemale; ++i)
    children[i] = female[parentG][i]->getNbChild(); // careful we have to look at the parent generation
  for(int i=0; i<(popsize-nbFemale); ++i)
    children[i+nbFemale] = male[parentG][i]->getNbChild();
  return children;
}

std::vector<double> Population::get_growthFactor() const
{
  std::vector<double> out;
  out.push_back(growthFactor[0]);
  out.push_back(growthFactor[1]);
  out.push_back(growthFactor[2]);
  return out;
}


Population Population::subsample(int sampleSize, int nbWomen)
{
  // Subsample the Population for analysis
  // Create a new object of he class Population
  int currentG = generation&1;
  if(sampleSize>popsize)
    {
      sampleSize = popsize;
    }
  if(nbWomen>nbFemale)
    {
      nbWomen = nbFemale;
    }
  Population sample(sampleSize);
  // NB : the Population sample will start at generation 0 -> must initiate at 0 and not currentG
  random_shuffle_pop();
  sample.nbFemale=nbWomen;
  sample.matingSystem=matingSystem;
  for(int i(0);i<nbWomen;++i)
    {
      sample.female[0][i]->reborn(female[currentG][i]); // pas referemce instead ? TODO
    }
  for(int i(0);i<sampleSize-nbWomen;++i)
    {
      sample.male[0][i]->reborn(male[currentG][i]);
    }
  return sample;
}

Population Population::subsample_female()
{
  // Subsample the Population for analysis
  // Create a new object of he class Population
  int currentG = generation&1;
  Population sample(2*nbFemale); // Default attribute for parameters - 2times female because we need female tab big enough
  sample.popsize = nbFemale; // then restore population size equals to female -> 'kill' the men
  for(int i(0);i<nbFemale;++i)
    {
      sample.female[0][i]->reborn(female[currentG][i]);
      // To not resample twice the same individual we swap the sampled ones and check somewhere else
    }
  return sample;
}

Population Population::subsample_male()
{
  int currentG = generation&1;
  int nbMale = popsize-nbFemale;
  Population sample(2*nbMale); // Default attribute for parameters -2times male because we need male tab big enough
  sample.popsize = nbMale; // then we 'kill' the female
  sample.nbFemale = 0;
  for(int i(0);i<(popsize-nbFemale);++i)
    {
      sample.male[0][i]->reborn(male[currentG][i]);
      // To not resample twice the same individual we swap the sampled ones and check somewhere else
    }
  return sample;
}



void Population::reproduction(MSNetwork & netmating,int populationID,double pmate,double pmig)
{
  if((nbFemale>0)&&(popsize>nbFemale)) // there are at least one male and one female
    {
      int currentG = (generation&1);
      int nbMating(0),nbBabyToMake(floor(popsizeDemo));
      int randomNbChildren,nbBabyFemale(0),nbBabyMale(0);
      int indexmale(-1),indexfemale(-1);// force the segfault if not initialized
      bool sexChild;
      //      std::vector<int> mating;
      // *********  Network setup *************
      GraphTraits::vertex_descriptor fatherNode = boost::vertex(netmating.nbNodes-1-2*populationID,netmating.myg);
      GraphTraits::vertex_descriptor motherNode = boost::vertex(netmating.nbNodes-2-2*populationID,netmating.myg);
      int gchild = (generation+1) % netmating.storedGeneration;
      //      int gparent = generation % netmating.storedGeneration;
      int indexchild = gchild * netmating.nnpg +offsetNetwork; // to sort out
      int chosen = -1;
      int counterPrescMating = 0;
      int apa=0;
      if(nbFemale>popsize/2)
	nbMating = popsize-nbFemale;
      else
	nbMating = nbFemale;
      GraphTraits::vertex_descriptor vdes;
      GraphTraits::vertex_descriptor vmom;
      GraphTraits::vertex_descriptor vchild;
      //	    GraphTraits::vertex_descriptor fatherNode = boost::vertex(netmating.nbNodes-1,netmating.myg);
      std::vector<GraphTraits::vertex_descriptor > possible;
      int npossible(0);
      //      cout << in_degree(motherNode,netmating.myg) << " " << nbFemale << " " << generation << endl;
      assert((int) in_degree(motherNode,netmating.myg)==nbFemale);
      assert((int) in_degree(fatherNode,netmating.myg)==(popsize-nbFemale));
      // std::vector<int>
      for(int i(0); i<nbMating; ++i)
	{
	  apa = false;
	  indexmale = -1;
	  indexfemale = i;
	  for(int pg(0);pg<polygamyf;pg++)
	    {
	      if(varChild == 0)
		randomNbChildren =2;
	      else
		{
		  randomNbChildren = rBinom(nbBabyToMake,1/(double)(nbMating-i));
		}
	      vmom = boost::vertex(female[currentG][indexfemale]->IDNetwork,netmating.myg);
	      double cheating = rBinom(1,pmate);
	      if(cheating)
		{
		  vdes = netmating.get_random_male(populationID);
		  indexmale = netmating.myg[vdes].ID;
		}
	      else
		{
		  possible = netmating.get_male(vmom,inbreeding,matingSystem,populationID);
		  possible = netmating.population_filter(possible,populationID);	      
		  npossible = possible.size();
		  if(npossible > 0)
		    {
		      chosen = myrand(npossible);
		      vdes = possible[chosen];
		      apa = true;
		    }
		  else // take a randon male via super node 
		    {
		      vdes = netmating.get_random_male(populationID);
		    }
		  indexmale = netmating.myg[vdes].ID;
		  while((male[currentG][indexmale]->getNbMate()>polygamy-1)&&npossible>-1) // monogamy
		    {
		      --npossible;
		      if(npossible>0)
			{
			  possible.erase(possible.begin()+chosen);
			  chosen = myrand(npossible);
			  vdes = possible[chosen];
			}
		      else
			{
			  vdes = netmating.get_random_male(populationID);
			  apa = false;
			}
		    }
		  indexmale = netmating.myg[vdes].ID;
		}
	      if(male[currentG][indexmale]->getNbMate()>polygamy-1)
		boost::remove_edge(vdes,fatherNode,netmating.myg);
	      if(TRACK_MATE_AND_CHILDREN)
		{
		  female[currentG][indexfemale]->hadChild(randomNbChildren);
		  male[currentG][indexmale]->hadChild(randomNbChildren);
		  female[currentG][indexfemale]->newMate();
		  male[currentG][indexmale]->newMate();
		}
	      if(apa)
		++counterPrescMating;
	      for(int k(0); k<randomNbChildren; ++k)
		{
		  if(varChild>0)
		    sexChild = rand2;
		  else
		    sexChild = k;
		  if(sexChild == 1)
		    {
		      birth_girl(indexfemale,indexmale,nbBabyFemale,indexchild,netmating);
		      ++nbBabyFemale;
		    }
		  else
		    {
		      birth_boy(indexfemale,indexmale,nbBabyMale,indexchild,netmating);
		      ++nbBabyMale;
		    }
		  ++indexchild;
		  --nbBabyToMake;
		}
	    }
	}
      nbFemale = nbBabyFemale;
      popsize = nbBabyMale+nbBabyFemale;      
      if(popsize > popsizeDemo)
	{
	  std::cerr << "ERROR : reproduction scheme is not working, too many baby were created " <<popsize << " / " << popsizeDemo << std::endl;
	}
    }
  else
    {
      cerr << "Population " << populationID <<  " is exctinct at generation " << generation << endl;
      popsize=0;
      nbFemale=0;
    }
}

void Population::reproduction_burnin()
{
  if((nbFemale>0)&&(popsize>nbFemale)) // there are at least one male and one female
    {
      int currentG = (generation&1);
      int nbMating(0),nbBabyToMake(floor(popsizeDemo));
      int randomNbChildren,nbBabyFemale(0),nbBabyMale(0);
      int indexmale(-1),indexfemale(-1);// force the segfault if not initialized
      bool sexChild;
      if(nbFemale>popsize/2)
	nbMating = popsize-nbFemale;
      else
	nbMating = nbFemale;
      for(int i(0); i<nbMating; ++i)
	{
	  indexmale = i;
	  indexfemale = i;
	  if(varChild == 0)
	    randomNbChildren =2;
	  else
	    {
	      randomNbChildren = rBinom(nbBabyToMake,1/(double)(nbMating-i));
	    }
	  for(int k(0); k<randomNbChildren; ++k)
	    {
	      if(varChild>0)
		{
		  sexChild = rand2;
		}
	      else
		sexChild = k;
	      if(sexChild == 1)
		{
		  birth_girl_burnin(indexfemale,indexmale,nbBabyFemale);
		  ++nbBabyFemale;
		}
	      else
		{
		  birth_boy_burnin(indexfemale,indexmale,nbBabyMale);
		  ++nbBabyMale;
		}
	      --nbBabyToMake;
	    }
	}
      nbFemale = nbBabyFemale;
      popsize = nbBabyMale+nbBabyFemale;
      //	}
      if(popsize > popsizeDemo)
	{
	  std::cerr << "ERROR : reproduction scheme is not working, too many baby were created " <<popsize << " / " << popsizeDemo << std::endl;
	}
    }
  else
    {
      popsize=0;
      nbFemale=0;
    }
}

void Population::demographic_update() // update function
{
  popsizeDemo = demographic_function(popsizeDemo,growthFactor[0],growthFactor[1],growthFactor[2]);
}

Population& Population::operator = (Population const& a)
{
  initiate_pop(a.popsize);
  female[0] = a.female[0];
  male[0] = a.male[0];
  female[1] = a.female[1];
  male[1] = a.male[1];
  popsizeDemo = a.popsizeDemo;
  popsize = a.popsize;
  nbFemale = a.nbFemale;
  generation = a.generation;
  memcpy(growthFactor,a.growthFactor,3*sizeof(double));
  matingSystem = a.matingSystem;
  inbreeding = a.inbreeding;
  polygamy = a.polygamy;
  polygamyf = a.polygamyf;
  varChild = a.varChild;
  maxsize = a.maxsize;
  migrationType = a.migrationType;
  return *this;
}

std::ostream& Population::print(std::ostream & os) const
{
  int currentG = (generation&1);
  os << "Pop with " << popsize << "individuals with a sex ratio of " << nbFemale/popsize << std::endl;
  for(int i(0);i<nbFemale;++i)
    {
      os << female[currentG][i];
    }
  for(int i(0);i<(popsize-nbFemale);++i)
    {
      os << male[currentG][i];
    }
  return os;
}


void Population::output_pop_size() const
{
  std::cout << "generation " << generation << " -- nb humans " << popsize << " -- Sex ratio " << (double)nbFemale/(double)popsize << std::endl;
}

void Population::reset(int N)
{
  popsize = N;
  generation = 0;
  nbFemale = N/2;
  popsizeDemo = N;
}

void Population::initiate_pop(int N)
{
  popsize = N;
  generation = 0;
  nbFemale = N/2;
  popsizeDemo = N;
  matingSystem.push_back(1);
  inbreeding = 0;
  polygamy = 1;
  polygamyf = 1;
  varChild = 1;
  growthFactor[0]=0;
  growthFactor[1]=1; // by default popsize = popsize
  growthFactor[2]=0;
  std::vector<double> prop;
  prop.push_back(0.5);
  prop.push_back(0.5);
  std::vector<  Human* > tmpSamp;
  female.push_back(tmpSamp);
  female.push_back(tmpSamp);
  male.push_back(tmpSamp);
  male.push_back(tmpSamp);
  female[0].resize(N); // alternatives resize or reserve
  female[1].resize(N);
  male[0].resize(N);
  male[1].resize(N);
  maxsize = N;
  offsetNetwork = 0;
  migrationType = 0;
}

void Population::resize(int nstep)
{
  maxsize = get_maxpop(nstep);
  female[0].resize(maxsize); // alternative resize
  female[1].resize(maxsize);
  male[0].resize(maxsize);
  male[1].resize(maxsize);
  maxsize = maxsize;
}

void Population::random_shuffle_pop() // knuth fisher yates algorithm taken from http://www.codinghorror.com/blog/2007/12/the-danger-of-naivete.html
{
  int currentG = (generation&1);
  int n;
  Human * tmp;
  for (int i = nbFemale - 1; i > 0; --i)
    {
      n = myrand(i-1);
      //swap(female[currentG][i],female[currentG][n]); // from algorithm
      tmp = female[currentG][i];
      female[currentG][i] = female[currentG][n];
      female[currentG][n] = tmp;
    }
  for (int i = popsize - nbFemale - 1; i > 0; --i)
    {
      n = myrand(i-1);
      tmp = male[currentG][i];
      male[currentG][i] = male[currentG][n];
      male[currentG][n] = tmp;
    }
}

void Population::output_diversity()
{
  std::cout <<  generation << "\t" << generation << "\t" << diversity_mtdna() <<  "\t"<< pairwise_distance_mtdna() << " \t " << diversity_A() << "\t" << diversity_X() << "\t" << diversity_Y() << std::endl;
}

void Population::output_nb_child(std::string filename,bool sex) const
{
  std::ofstream fs;
  std::string out = "";
  std::vector<int> nbChildren = get_nb_children();
  if(sex)
    for(int i(0);i<nbFemale;++i)
      out += to_string(nbChildren[i]) + "\n";
  else
    {
      for(int i(nbFemale);i<popsize;++i)
	out += to_string(nbChildren[i]) + "\n";
    }
  fs.open(filename.c_str(),std::ios::app);
  fs << out;
  fs.close();
}

void Population::output_nb_mates(std::string filename,bool sex) const
{
  std::ofstream fs;
  std::string out = "";
  std::vector<int> nbMates = get_nb_mates();
  if(sex)
    for(int i(0);i<nbFemale;++i)
      out += to_string(nbMates[i]) + "\n";
  else
    for(int i(nbFemale);i<popsize;++i)
      out += to_string(nbMates[i]) + "\n";
  fs.open(filename.c_str(),std::ios::app);
  fs << out;
  fs.close();
}

double Population::diversity_mtdna() const
{ // theta Pi
  typeDNA tmp(0),tmp2(0),x(0),c(0);
  int factoriel = 0;
  double tot = 0;
  int currentG = (generation&1);
  for(int k(0); k<DNAS.nbCrumbMt; ++k)
    {
      x = 0;
      factoriel = 0;
      for(int i(0); i<popsize; ++i)
	{
	  if(i<nbFemale)
	    tmp = female[currentG][i]->mtDna[k];
	  else
	    tmp = male[currentG][i]->mtDna[k];
	  for(int j(i+1); j<popsize; ++j)
	    {
	      if(j<nbFemale)
		tmp2 = female[currentG][j]->mtDna[k];
	      else
		tmp2 = male[currentG][j]->mtDna[k];
	      c = tmp ^ tmp2;
	      x += dna_distance(c);
	      ++factoriel;
	    }
	}
    }
  tot += (double) x /((double) factoriel);
  return tot;
}

void Population::mutation()
{
  // mtDNA
  int currentG = (generation&1);
  int nbBase = popsize*DNAS.sizeMt;
  int nbTransition = rBinom(nbBase,RATES.muMt1);
  int randomHuman, Xmutant, randomBase;
  typeDNA newdna, newBase;
  for(int i(0);i<nbTransition;++i)
    {
      randomHuman = myrand(popsize);
      Xmutant = myrand(DNAS.nbCrumbMt);
      randomBase = myrand(CRUMB);
      newBase = 1;
      newdna = newBase<<(randomBase<<1);
      if(randomHuman >= nbFemale){
	male[currentG][randomHuman-nbFemale]->mtDna[Xmutant] ^= newdna;
      }
      else
	female[currentG][randomHuman]->mtDna[Xmutant] ^= newdna;
    }
  int nbTransversion = rBinom(nbBase,RATES.muMt2);
  for(int i(0);i<nbTransversion;++i)
    {
      randomHuman = myrand(popsize);
      Xmutant = myrand(DNAS.nbCrumbMt);
      randomBase = myrand(CRUMB);
      newBase = rand2 + 2;
      newdna = newBase<<(randomBase<<1);
      if(randomHuman >= nbFemale)
	male[currentG][randomHuman-nbFemale]->mtDna[Xmutant] ^= newdna;
      else
	female[currentG][randomHuman]->mtDna[Xmutant] ^= newdna;
    }
  // Autosomal DNA
  if(DNAS.nbChromosomesA>0)
    { // if there is actually autosomal dna
      nbBase = popsize * DNAS.nbChromosomesA * DNAS.sizePerA*2;
      nbTransition = rBinom(nbBase, RATES.muA1);
      for(int i(0);i<nbTransition;++i)
	{
	  randomHuman =  myrand(popsize);
	  Xmutant =  myrand(DNAS.nbCrumbA);
	  randomBase = myrand(CRUMB);
	  newBase = 1; // always the opposite
	  typeDNA newdna = newBase<<(randomBase<<1);
	  if(randomHuman >= nbFemale)
	    male[currentG][randomHuman-nbFemale]->dna[Xmutant] ^= newdna;
	  else
	    female[currentG][randomHuman]->dna[Xmutant] ^= newdna;
	}
      nbTransversion = rBinom(nbBase, RATES.muA2);
      for(int i(0);i<nbTransversion;++i)
	{
	  randomHuman = myrand(popsize);
	  Xmutant = myrand(DNAS.nbCrumbA);
	  randomBase = myrand(CRUMB);
	  newBase = rand2 + 2;//  =  2 ou 3, with ^ obtain 0 1 for 2 or 3 and vice versa
	  newdna = newBase<<(randomBase<<1);
	  if(randomHuman >= nbFemale)
	    male[currentG][randomHuman-nbFemale]->dna[Xmutant] ^= newdna;
	  else
	    female[currentG][randomHuman]->dna[Xmutant] ^= newdna;
	}
    }
  // Y chromosome
  nbBase = (popsize -nbFemale) * DNAS.sizeY;
  nbTransition = rBinom(nbBase, RATES.muY1);
  for(int i(0);i<nbTransition;++i)
    {
      randomHuman =  myrand(popsize-nbFemale);
      Xmutant =  myrand(DNAS.nbCrumbY);
      randomBase = myrand(CRUMB);
      newBase = 1; // always the opposite
      newdna = newBase<<(randomBase<<1);
      male[currentG][randomHuman]->dna[Xmutant+DNAS.nbCrumbA+DNAS.nbCrumbX] ^= newdna;
    }
  nbTransversion = rBinom(nbBase, RATES.muY2);
  for(int i(0);i<nbTransversion;++i)
    {
      randomHuman =  myrand(nbFemale);
      Xmutant = myrand(DNAS.nbCrumbY);
      randomBase = myrand(CRUMB);
      newBase = rand2 + 2;//  =  2 ou 3, with ^ obtain 0 1 for 2 or 3 and vice versa
      newdna = newBase<<(randomBase<<1);
      male[currentG][randomHuman]->dna[Xmutant+DNAS.nbCrumbA+DNAS.nbCrumbX] ^= newdna;
    }
  // X chromosome
  nbBase = (popsize + nbFemale) * DNAS.sizePerX * DNAS.nbChromosomesX;
  nbTransition = rBinom(nbBase, RATES.muX1);
  for(int i(0);i<nbTransition;++i)
    {
      randomHuman =  myrand(popsize+nbFemale); // small trickif bigger than popsize is female -> twice as many X for female
      Xmutant =  myrand(DNAS.nbCrumbX);
      randomBase = myrand(CRUMB);
      newBase = 1; // always the opposite
      newdna = newBase<<(randomBase<<1);
      if(randomHuman >= popsize)
	female[currentG][randomHuman-popsize]->dna[Xmutant + DNAS.nbCrumbX + DNAS.nbCrumbA] ^= newdna;
      else if(randomHuman >= nbFemale)
	male[currentG][randomHuman-nbFemale]->dna[Xmutant + DNAS.nbCrumbA] ^= newdna;
      else
	female[currentG][randomHuman]->dna[Xmutant + DNAS.nbCrumbA] ^= newdna;
    }
  nbTransversion = rBinom(nbBase, RATES.muX2);
  for(int i(0);i<nbTransversion;++i)
    {
      randomHuman =  myrand(popsize+nbFemale); // small trickif bigger than popsize is female -> twice as many X for female
      Xmutant =  myrand(DNAS.nbCrumbX);
      randomBase = myrand(CRUMB);
      newBase = rand2 + 2;//  =  2 ou 3, with ^ obtain 0 1 for 2 or 3 and vice versa
      newdna = newBase<<(randomBase<<1);
      if(randomHuman >= popsize)
	female[currentG][randomHuman-popsize]->dna[Xmutant + DNAS.nbCrumbX + DNAS.nbCrumbA] ^= newdna;
      else if(randomHuman >= nbFemale)
	male[currentG][randomHuman-nbFemale]->dna[Xmutant + DNAS.nbCrumbA] ^= newdna;
      else
	female[currentG][randomHuman]->dna[Xmutant + DNAS.nbCrumbA] ^= newdna;
    }
}

void Population::birth_boy(int indexfemale,int indexmale,int nbBabyMale,int indexchild,MSNetwork & netmating)
{
  int index,randomX;
  typeDNA gamete[DNAS.nbCrumbNuc];
  int currentG  =  (generation&1);
  GraphTraits::vertex_descriptor mum,dad,kid;
  if(nbBabyMale<popsizeDemo){
    for(int k(0);k<DNAS.nbChromosomesA;++k) // all the chromosomes except the sexual ones , the last 2 ones
      {
	index = k * DNAS.nbCellPerA*2;
	randomX = rand2; // random on the kith mother chromosome
	memcpy(&(gamete[index]),&(female[currentG][indexfemale]->dna[index + randomX * DNAS.nbCellPerA]),DNAS.nbCellPerA * sizeof(typeDNA));
	randomX = rand2; // random on the father
	memcpy(&(gamete[index + DNAS.nbCellPerA]),&(male[currentG][indexmale]->dna[index + randomX * DNAS.nbCellPerA]),DNAS.nbCellPerA * sizeof(typeDNA));
      }
    // pick randomly an X chromosome from the mother, no recombination
    if(DNAS.nbChromosomesX>0)
      {
	index = DNAS.nbCrumbA;
	randomX = rand2; // random on the kith mother chromosome
	memcpy(&(gamete[index]),&(female[currentG][indexfemale]->dna[index + randomX * DNAS.nbCrumbX]),DNAS.nbCrumbX * sizeof(typeDNA));
      }
    // copy the Y chromosome from the father
    if(DNAS.sizeY>0)
      {
	index = DNAS.nbCrumbA + DNAS.nbCrumbX ;
	memcpy(&(gamete[index]),&(male[currentG][indexmale]->dna[index]),DNAS.nbCrumbY * sizeof(typeDNA));
      }
    male[!currentG][nbBabyMale]->reborn(gamete,female[currentG][indexfemale]->mtDna);
    if(((migrationType==3)||(migrationType==4))||(migrationType==6))
      male[!currentG][nbBabyMale]->destinationPop = male[currentG][indexmale]->destinationPop;
    else if(migrationType==5)
      male[!currentG][nbBabyMale]->destinationPop = female[currentG][indexfemale]->destinationPop;
    if((male[currentG][indexmale]->IDNetwork <0)||(female[currentG][indexfemale]->IDNetwork<0))
      {
	cerr << "ERROR IDNetwork <0 m index" << indexmale << " m idnetwork " << male[currentG][indexmale]->IDNetwork << " --  f index " <<  indexfemale << " f idnetwork " << female[currentG][indexfemale]->IDNetwork << " Popsize " << popsize << " nfemale " << nbFemale << endl;
	exit(2);
      }
    dad = boost::vertex(male[currentG][indexmale]->IDNetwork,netmating.myg);
    mum = boost::vertex(female[currentG][indexfemale]->IDNetwork,netmating.myg);
    kid = boost::vertex(indexchild,netmating.myg);
    assert( netmating.myg[dad].generation ==generation);
    assert( netmating.myg[mum].generation ==generation);
    boost::add_edge(dad,kid,netmating.myg);
    boost::add_edge(mum,kid,netmating.myg);
    male[!currentG][nbBabyMale]->IDNetwork = indexchild; ///toto
    if(netmating.myg[kid].population>-1)
      cout << " ERROR kid.pop already assigned " << netmating.myg[kid].ID << " " << netmating.myg[kid].sex << " " << netmating.myg[kid].generation << " " << netmating.myg[kid].population << " " << indexchild << " " << nbBabyMale << endl;
    netmating.myg[kid].ID = nbBabyMale;
    netmating.myg[kid].sex = false;
    netmating.myg[kid].generation = generation+1;
    netmating.myg[kid].population = netmating.myg[dad].population;
  }
  else
    {
      std::cerr<<"ERROR : attempts to create too many boys"<<std::endl;
      exit(2);
    }
}

void Population::birth_girl(int indexfemale,int indexmale,int nbBabyFemale,int indexchild,MSNetwork & netmating)
{
  int index,randomX,randomP;
  typeDNA gamete[DNAS.nbCrumbNuc];
  int currentG  =  (generation&1);
  GraphTraits::vertex_descriptor mum, dad, kid;
  if(nbBabyFemale<popsizeDemo){
    for(int k(0);k<DNAS.nbChromosomesA;++k) // all the chromosomes except the sexual ones , the last 2 ones
      {
	index = k * DNAS.nbCellPerA*2;
	randomX = rand2; // random on the kith mother chromosome
	memcpy(&(gamete[index]),&(female[currentG][indexfemale]->dna[index + randomX * DNAS.nbCellPerA]),DNAS.nbCellPerA * sizeof(typeDNA));
	randomX = rand2; // random on the father
	memcpy(&(gamete[index + DNAS.nbCellPerA]),&(male[currentG][indexmale]->dna[index + randomX * DNAS.nbCellPerA]),DNAS.nbCellPerA * sizeof(typeDNA));
      }
    // pick randomly an X chromosome from the mother and place as the nbChromosomes-2 chromosome
    index = DNAS.nbCrumbA;
    for(int k(0);k<DNAS.nbChromosomesX; ++k) // all the chromosomes except the sexual ones , the last 2 ones
      {
	index = k * DNAS.nbCellPerX + DNAS.nbCrumbA;
	randomX = rand2; // random on the kith mother chromosome
	randomP = rand2; // random on whether it recombines or not (on which allele it lands)
	memcpy(&(gamete[index+randomP*DNAS.nbCrumbX]),&(female[currentG][indexfemale]->dna[index + randomX * DNAS.nbCrumbX]),DNAS.nbCellPerX * sizeof(typeDNA));
	randomP ^= 1;
	memcpy(&(gamete[index+randomP*DNAS.nbCrumbX]),&(male[currentG][indexmale]->dna[index]),DNAS.nbCellPerX * sizeof(typeDNA));
      }
    // copy the X chromosome from the father, no recomb possible
    female[!currentG][nbBabyFemale]->reborn(gamete,female[currentG][indexfemale]->mtDna);
    if((male[currentG][indexmale]->IDNetwork <0)||(female[currentG][indexfemale]->IDNetwork<0))
      {
	cerr << "ERROR IDNetwork <0 m index" << indexmale << " m idnetwork " << male[currentG][indexmale]->IDNetwork << " --  f index " <<  indexfemale << " f idnetwork " << female[currentG][indexfemale]->IDNetwork << " Popsize " << popsize << " nfemale " << nbFemale << endl;
	exit(2);
      }
    if(((migrationType==3)||(migrationType==4))||(migrationType==6))
      female[!currentG][nbBabyFemale]->destinationPop=male[currentG][indexmale]->destinationPop;
    else if(migrationType==5)
      female[!currentG][nbBabyFemale]->destinationPop=female[currentG][indexfemale]->destinationPop;
    dad = boost::vertex(male[currentG][indexmale]->IDNetwork,netmating.myg);
    mum = boost::vertex(female[currentG][indexfemale]->IDNetwork,netmating.myg);
    kid = boost::vertex(indexchild,netmating.myg);
    boost::add_edge(dad,kid,netmating.myg);
    boost::add_edge(mum,kid,netmating.myg);
    if(netmating.myg[kid].population>-1)
      cout << " ERROR kid.pop already assigned " << netmating.myg[kid].ID << " " << netmating.myg[kid].sex << " " << netmating.myg[kid].generation << " " << netmating.myg[kid].population << " " << indexchild << " " << nbBabyFemale << endl;
    netmating.myg[kid].ID = nbBabyFemale;
    netmating.myg[kid].sex = true;
    netmating.myg[kid].generation = generation+1;
    netmating.myg[kid].population = netmating.myg[dad].population;
    female[!currentG][nbBabyFemale]->IDNetwork = indexchild;
  }
  else
    {
      std::cerr<<"ERROR : attemps to create too many girls"<<std::endl;
      exit(2);
    }
}

void Population::birth_boy_burnin(int indexfemale,int indexmale,int nbBabyMale)
{
  int index,randomX;
  typeDNA gamete[DNAS.nbCrumbNuc];
  int currentG  =  (generation&1);
  if(nbBabyMale<popsizeDemo){
    for(int k(0);k<DNAS.nbChromosomesA;++k) // all the chromosomes except the sexual ones , the last 2 ones
      {
	index = k * DNAS.nbCellPerA*2;
	randomX = rand2; // random on the kith mother chromosome
	memcpy(&(gamete[index]),&(female[currentG][indexfemale]->dna[index + randomX * DNAS.nbCellPerA]),DNAS.nbCellPerA * sizeof(typeDNA));
	randomX = rand2; // random on the father
	memcpy(&(gamete[index + DNAS.nbCellPerA]),&(male[currentG][indexmale]->dna[index + randomX * DNAS.nbCellPerA]),DNAS.nbCellPerA * sizeof(typeDNA));
      }
    // pick randomly an X chromosome from the mother, no recombination
    if(DNAS.nbChromosomesX>0)
      {
	index = DNAS.nbCrumbA;
	randomX = rand2; // random on the kith mother chromosome
	memcpy(&(gamete[index]),&(female[currentG][indexfemale]->dna[index + randomX * DNAS.nbCrumbX]),DNAS.nbCrumbX * sizeof(typeDNA));
      }
    // copy the Y chromosome from the father
    if(DNAS.sizeY>0)
      {
	index = DNAS.nbCrumbA + DNAS.nbCrumbX ;
	memcpy(&(gamete[index]),&(male[currentG][indexmale]->dna[index]),DNAS.nbCrumbY * sizeof(typeDNA));
      }
    male[!currentG][nbBabyMale]->reborn(gamete,female[currentG][indexfemale]->mtDna);
    male[!currentG][nbBabyMale]->destinationPop = male[currentG][indexmale]->destinationPop;
  }
  else
    {
      std::cerr<<"ERROR : attempts to create too many boys"<<std::endl;
      exit(2);
    }
}

void Population::birth_girl_burnin(int indexfemale,int indexmale,int nbBabyFemale)
{
  int index,randomX,randomP;
  typeDNA gamete[DNAS.nbCrumbNuc];
  int currentG  =  (generation&1);
  if(nbBabyFemale<popsizeDemo){
    for(int k(0);k<DNAS.nbChromosomesA;++k) // all the chromosomes except the sexual ones , the last 2 ones
      {
	index = k * DNAS.nbCellPerA*2;
	randomX = rand2; // random on the kith mother chromosome
	memcpy(&(gamete[index]),&(female[currentG][indexfemale]->dna[index + randomX * DNAS.nbCellPerA]),DNAS.nbCellPerA * sizeof(typeDNA));
	randomX = rand2; // random on the father
	memcpy(&(gamete[index + DNAS.nbCellPerA]),&(male[currentG][indexmale]->dna[index + randomX * DNAS.nbCellPerA]),DNAS.nbCellPerA * sizeof(typeDNA));
      }
    // pick randomly an X chromosome from the mother and place as the nbChromosomes-2 chromosome
    index = DNAS.nbCrumbA;
    for(int k(0);k<DNAS.nbChromosomesX; ++k) // all the chromosomes except the sexual ones , the last 2 ones
      {
	index = k * DNAS.nbCellPerX + DNAS.nbCrumbA;
	randomX = rand2; // random on the kith mother chromosome
	randomP = rand2; // random on whether it recombines or not (on which allele it lands)
	memcpy(&(gamete[index+randomP*DNAS.nbCrumbX]),&(female[currentG][indexfemale]->dna[index + randomX * DNAS.nbCrumbX]),DNAS.nbCellPerX * sizeof(typeDNA));
	randomP ^= 1;
	memcpy(&(gamete[index+randomP*DNAS.nbCrumbX]),&(male[currentG][indexmale]->dna[index]),DNAS.nbCellPerX * sizeof(typeDNA));
      }
    // copy the X chromosome from the father, no recomb possible
    female[!currentG][nbBabyFemale]->reborn(gamete,female[currentG][indexfemale]->mtDna);
    female[!currentG][nbBabyFemale]->destinationPop=male[currentG][indexmale]->destinationPop;
  }
  else
    {
      std::cerr<<"ERROR : attemps to create too many girls"<<std::endl;
      exit(2);
    }
}


typeDNA Population::polySites_mtdna() const
{
  typeDNA tot(0),c(0);
  int currentG  =  (generation&1);
  typeDNA tmp; // storing the number is faster
  for(int k(0);k<DNAS.nbCrumbMt;++k)
    {
      c = 0;
      tmp = female[currentG][0]->mtDna[k];
      for(int j(1);j<nbFemale;++j)
	{
	  c |= tmp^female[currentG][j]->mtDna[k];
	}
      for(int j(0);j<popsize-nbFemale;++j)
	{
	  c |= tmp^male[currentG][j]->mtDna[k];
	}
      tot+=dna_distance(c);
    }
  return tot;
}

double Population::pairwise_distance_mtdna() const
{
  typeDNA c(0),x(0);	//x(0),e(0),d(0),c(0)
  double tot(0);
  int factoriel(0);
  int currentG = (generation&1);
  typeDNA tmp; // storing the number is faster
  for(int k(0);k<DNAS.nbCrumbMt;++k)
    {
      x = 0;
      factoriel = 0;
      for(int i(0);i<nbFemale;++i)
	{
	  tmp = female[currentG][i]->mtDna[k];
	  for(int j(i+1);j<nbFemale;++j)
	    {
	      c = tmp ^ female[currentG][j]->mtDna[k];
	      x += dna_distance(c); // nbFemale!
	      ++factoriel;
	    }
	  for(int j(0);j<popsize-nbFemale;++j)
	    {
	      c = tmp ^ male[currentG][j]->mtDna[k];
	      x += dna_distance(c); // nbFemale x nbMale
	      ++factoriel;
	    }
	}
      for(int i(0);i<popsize-nbFemale;++i)
	{
	  tmp = male[currentG][i]->mtDna[k];
	  for(int j(i+1);j<popsize-nbFemale;++j)
	    {
	      c = tmp^male[currentG][j]->mtDna[k];
	      x += dna_distance(c); // nbMale!
	      ++factoriel;
	    }
	}
      tot += (double) x /((double) factoriel);
    }
  return tot;// / (double) SIZE_MTDNA;
}

double Population::pairwise_distance_Y() const
{
  typeDNA c(0),x(0);
  double tot(0);
  int factoriel(0);
  int currentG = (generation&1);
  typeDNA tmp; // storing the number is faster
  int offset = DNAS.nbCrumbX + DNAS.nbCrumbY;
  int nbMale = popsize - nbFemale;
  for(int k(0);k<DNAS.nbCrumbY;++k)
    {
      x = 0;
      factoriel = 0;
      for(int i(0);i<nbMale;++i)
	{
	  tmp = male[currentG][i]->dna[k+offset];
	  for(int j(i+1);j<nbMale;++j)
	    {
	      c = tmp ^ male[currentG][j]->dna[k+offset];
	      x += dna_distance(c); // nbFemale!
	      ++factoriel;
	    }
	}
      tot += (double) x /((double) factoriel);
    }
  return tot;// / (double) SIZE_MTDNA;
}

double Population::pairwise_distance_X() const
{
  typeDNA c(0),x(0);	//x(0),e(0),d(0),c(0)
  double tot(0);
  int factoriel(0);
  int currentG = (generation&1);
  typeDNA tmp; // storing the number is faster
  int offset = DNAS.nbCrumbA;
  for(int k(0);k<DNAS.nbCrumbX;++k)
    {
      x = 0;
      factoriel = 0;
      for(int i(0);i<nbFemale;++i)
	{
	  tmp = female[currentG][i]->dna[k+offset];
	  c = tmp ^ female[currentG][i]->dna[k+offset+DNAS.nbCrumbX];
	  x += dna_distance(c); // nbFemale!
	  ++factoriel;
	  for(int j(i+1);j<nbFemale;++j)
	    {
	      c = tmp ^ female[currentG][j]->dna[k+offset];
	      x += dna_distance(c); // nbFemale!
	      ++factoriel;
	      c = tmp ^ female[currentG][j]->dna[k+offset+DNAS.nbCrumbX];
	      x += dna_distance(c); // nbFemale!
	      ++factoriel;
	    }
	  for(int j(0);j<popsize-nbFemale;++j)
	    {
	      c = tmp ^ male[currentG][j]->dna[k+offset];
	      x += dna_distance(c);
	      ++factoriel;
	    }
	}
      for(int i(0);i<popsize-nbFemale;++i)
	{
	  tmp = male[currentG][i]->dna[k];
	  for(int j(i+1);j<popsize-nbFemale;++j)
	    {
	      c = tmp^male[currentG][j]->dna[k+offset];
	      x += dna_distance(c); // nbMale!
	      ++factoriel;
	    }
	}
      tot += (double) x /((double) factoriel);
    }
  return tot;// / (double) SIZE_MTDNA;
}

double Population::pairwise_distance_A() const
{
  typeDNA c(0),x(0);	//x(0),e(0),d(0),c(0)
  double tot(0);
  int factoriel(0);
  int currentG = (generation&1);
  typeDNA tmp; // storing the number is faster
  for(int k(0);k<DNAS.nbCrumbA; k += 2)
    {
      x = 0;
      factoriel = 0;
      for(int i(0);i<nbFemale;++i)
	{
	  tmp = female[currentG][i]->dna[k];
	  c = tmp ^ female[currentG][i]->dna[k+1];
	  x += dna_distance(c);
	  ++factoriel;
	  for(int j(i+1);j<nbFemale;++j)
	    {
	      c = tmp ^ female[currentG][j]->dna[k];
	      x += dna_distance(c); // nbFemale!
	      ++factoriel;
	      c = tmp ^ female[currentG][j]->dna[k+1];
	      x += dna_distance(c); // nbFemale!
	      ++factoriel;
	    }
	  for(int j(0);j<popsize-nbFemale;++j)
	    {
	      c = tmp ^ male[currentG][j]->dna[k];
	      x += dna_distance(c); // nbFemale x nbMale
	      ++factoriel;
	      c = tmp ^ male[currentG][j]->dna[k+1];
	      x += dna_distance(c); // nbFemale x nbMale
	      ++factoriel;
	    }
	}
      for(int i(0);i<popsize-nbFemale;++i)
	{
	  tmp = male[currentG][i]->dna[k];
	  for(int j(i+1);j<popsize-nbFemale;++j)
	    {
	      c = tmp^male[currentG][j]->dna[k];
	      x += dna_distance(c); // nbMale!
	      ++factoriel;
	      c = tmp^male[currentG][j]->dna[k+1];
	      x += dna_distance(c); // nbMale!
	      ++factoriel;
	    }
	}
      tot += (double) x /((double) factoriel);
    }
  return tot;
}

typeDNA Population::polySites(int nuclear) const
{ // WARNING : this work only if typeDNA is 64 bits integer
  typeDNA tot(0),c(0);
  int currentG = (generation&1);
  int allele1,allele2;
  typeDNA tmp,tmplong; // storing the number is faster
  switch(nuclear)
    {
    case 0 : // autosomes -> two similar chr per human
      for(int k(0);k<DNAS.nbChromosomesA;++k) // take not the sexual chromsome in this calculus
	for(int l(0);l<(DNAS.nbCellPerA);l++) // take not the sexual chromsome in this calculus
	  {
	    c = 0;
	    allele1=k*DNAS.nbCellPerA+l;
	    allele2=(k+1)*DNAS.nbCellPerA+l;
	    tmp = female[currentG][0]->dna[allele1];
	    tmplong = tmp^female[currentG][0]->dna[allele2];
	    c |= tmplong;
	    for(int j(1);j<nbFemale;++j)
	      {
		//	      c |= female[currentG][i]->dna[k]^female[currentG][j]->dna[k];
		tmplong = tmp^female[currentG][j]->dna[allele1];
		c |= tmplong;
		tmplong = tmp^female[currentG][j]->dna[allele2];
		c |=  tmplong;
	      }
	    for(int j(0);j<(popsize-nbFemale);++j)
	      {
		tmplong = tmp^male[currentG][j]->dna[allele1];
		c |= tmplong;
		tmplong = tmp^male[currentG][j]->dna[allele2];
		c |= tmplong;
	      }
	    tot+=dna_distance(c);
	  }
      break;
    case 1 : // X chromosome -> two per female, one per male
      for(int l(0);l<(DNAS.nbCrumbX);l++) // take not the sexual chromsome in this calculus
	{
	  c = 0;
	  allele1=DNAS.nbCrumbA+l;
	  allele2=DNAS.nbCrumbA+DNAS.nbCrumbX+l;
	  tmp = female[currentG][0]->dna[allele1];
	  tmplong = tmp^female[currentG][0]->dna[allele2];
	  c |= tmplong;
	  for(int j(1);j<nbFemale;++j)
	    {
	      tmplong = tmp^female[currentG][j]->dna[allele1];
	      c |= tmplong;
	      tmplong = tmp^female[currentG][j]->dna[allele2];
	      c |= tmplong;
	    }
	  for(int j(0);j<(popsize-nbFemale);++j)
	    {
	      tmplong = tmp^male[currentG][j]->dna[allele1];
	      c |= tmplong;
	    }
	  tot+=dna_distance(c);
	}
      break;
    case 2 : // Y chromosome -> one per male
      for(int l(0);l<(DNAS.nbCrumbY);l++) // take not the sexual chromsome in this calculus
	{
	  c = 0;
	  allele1=DNAS.nbCrumbA+DNAS.nbCrumbX+l;
	  tmp = male[currentG][0]->dna[allele1];
	  for(int j(1);j<(popsize-nbFemale);++j)
	    {
	      tmplong = tmp^male[currentG][j]->dna[allele1];
	      c |= tmplong;
	    }
	  tot+=dna_distance(c);
	}
      break;
    default:
      break;
    }
  return tot;
}

double Population::diversity_A() const
{
  typeDNA nbPolySites;
  nbPolySites = polySites(0);
  return (2 * (double) nbPolySites)/(double) (DNAS.sizePerA * (DNAS.nbChromosomesA));
}

double Population::diversity_Y() const
{
  typeDNA nbPolySites;
  nbPolySites = polySites(2);
  return ((double) nbPolySites)/(double) (DNAS.sizePerA);
}

double Population::diversity_X() const
{
  typeDNA nbPolySites;
  nbPolySites = polySites(1);
  return ((double) nbPolySites)/(double) (DNAS.sizePerA);
}

void Population::test(int i)
{
  int currentG = generation&1;
  for(int j(0);j<maxsize;++j)
    {
      if(male[currentG][j]==0)
	cout << "male " << j << " of population " << i << " is null " << endl;
      if(female[currentG][j]==0)
	cout << "female " << j << " of population " << i << " is null " << endl;
    }

}

int Population::burnin(double diversityTreshold, int genLimit, double facMu,int dnaType,MSNetwork & netmating,int popID)
{
  if(facMu>0)
    {
      std::vector<int> realMatingSystem;
      double diversity = 0;
      int returnGeneration;
      Population A;
      realMatingSystem = matingSystem;
      matingSystem = std::vector<int>(1,2);
      while ((diversityTreshold>diversity) && (generation<genLimit))
	{
	  int gparent = generation % netmating.storedGeneration;
	  random_shuffle_pop();
	  mutation();
	  demographic_update();
	  update_netmating(netmating,popID);
	  netmating.clear_premating(gparent);
	  reproduction(netmating,popID,1,1);	  
	  netmating.clear_postmating();
	  ++generation;
	  if((generation==floor(genLimit*3.0/4.0)))
	    { // accelerate the burnin with 10 times increase in the mutation rate
	      RATES.muMt1 *= 10;
	      RATES.muMt2 *= 10;
	      RATES.muX1 *= 10;
	      RATES.muX2 *= 10;
	      RATES.muY1 *= 10;
	      RATES.muY2 *= 10;
	      RATES.muA1 *= 10;
	      RATES.muA2 *= 10;
	    }
	  switch(dnaType)
	    {
	    case 1:
	      {
		diversity = A.pairwise_distance_mtdna();
		break;
	      }
	    case 2:
	      {
		diversity = A.pairwise_distance_Y();
		break;
	      }
	    case 3:
	      {
		diversity = A.pairwise_distance_X();
		break;
	      }
	    case 4:
	      {
		diversity = A.pairwise_distance_A();
		break;
	      }
	    default:
	      {
		diversity = A.pairwise_distance_mtdna();
		break;
	      }
	    }
	}
      matingSystem = realMatingSystem;
      returnGeneration = generation;
      generation = 0;
      return returnGeneration;
    }
  else
    {
      cerr << "ERROR in the burnin parameter, factor must be bigger than 0" << endl;
      exit(2);
    }
  return -1;
}

void Population::update_netmating(MSNetwork & netmating, int popID)
{
  int currentG(generation&1);
  for(int i(0); i<nbFemale; ++i)
    {
      if(female[currentG][i]->IDNetwork<0)
	cout << "ERROR pop burnin " << " f " << i<< " nbFemale " << nbFemale <<endl;
      else
	{
	  netmating.myg[boost::vertex(female[currentG][i]->IDNetwork,netmating.myg)].population = popID;
	  netmating.myg[boost::vertex(female[currentG][i]->IDNetwork,netmating.myg)].ID = i;
	}
    }
  for(int i(0); i<(popsize - nbFemale); ++i)
    {
      if(male[currentG][i]->IDNetwork<0)
	cout << "ERROR pop burnin  " << " m " << i<< " msize " << popsize-nbFemale << endl;
      else
	{
	  netmating.myg[boost::vertex(male[currentG][i]->IDNetwork,netmating.myg)].population = popID;
	  netmating.myg[boost::vertex(male[currentG][i]->IDNetwork,netmating.myg)].ID = i;
	}
    }
}

void Population::set_maleIDnetwork(MSNetwork & netmating, int popID)
{
  int currentG = generation & 1;
  GraphTraits::vertex_descriptor p;
  for(int i(0); i<popsize-nbFemale;++i)
    {
      p = boost::vertex(male[currentG][i]->IDNetwork,netmating.myg);
      netmating.myg[p].ID = i;
      netmating.myg[p].population = popID;
    }
}

void Population::set_femaleIDnetwork(MSNetwork & netmating, int popID)
{
  int currentG = generation & 1;
  GraphTraits::vertex_descriptor p;
  for(int i(0); i<nbFemale;++i)
    {
      if(female[currentG][i]->IDNetwork>-1)
	{
	  p = boost::vertex(female[currentG][i]->IDNetwork,netmating.myg);
	  netmating.myg[p].ID = i;
	  netmating.myg[p].population = popID;
	}
      else
	cerr << "ERROR female " << i << " of generation " << generation << "does not belong to the network"<<endl;
    }
}
