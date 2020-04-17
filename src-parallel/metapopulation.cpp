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

#include "metapopulation.h"

using namespace std;

Metapop::Metapop():netmating(0)
{
  migrationType = 0;
  nbPop = 0;
  popsize = 0;
  generation = 0;
  nbFemale = 0;
  pmate = 0;
  pmig = 0;
  nbSources = 0;
}

Metapop::Metapop(int a, int b):netmating(b)
{
  nbPop = a;
  migrationType = 0;
  poptab = new Population*[a];
  popsize = b;
  nbFemale = 0;
  for(int i=0;i<nbPop;++i)
    {
      poptab[i] = new Population(b/a);
      popID.push_back(i);
      nbFemale += (b/a)/2; // repeat the exact operations -> no bad surprises with rounding
      migrationRateF.push_back(std::vector<double>(nbPop,0));
      migrationRateM.push_back(std::vector<double>(nbPop,0));
      nbMigrantsF.push_back(std::vector<int>(nbPop,0));
      nbMigrantsM.push_back(std::vector<int>(nbPop,0));
    }
  generation = 0;
  pmate = 0;
  pmig = 0;
  nbSources = 0;
}

Metapop::Metapop(int npop, int nindiv, std::vector<int> niList):netmating(nindiv)
{
  nbPop = npop;
  migrationType = 0;
  poptab = new Population*[npop];
  popsize = nindiv;
  nbFemale = 0;
  for(int i=0;i<nbPop;++i)
    {
      poptab[i] = new Population(niList[i]);
      popID.push_back(i);
      nbFemale += (niList[i])/2; // repeat the exact operations -> no bad surprises with rounding
      migrationRateF.push_back(std::vector<double>(nbPop));
      migrationRateM.push_back(std::vector<double>(nbPop));
      nbMigrantsF.push_back(std::vector<int>(nbPop));
      nbMigrantsM.push_back(std::vector<int>(nbPop));
    }
  generation = 0;
  pmate = 0;
  pmig = 0;
  nbSources = 0;
}

Metapop::Metapop(Metapop const & a):netmating(a.netmating)
{
  for(unsigned i=0;i<female[0].size();++i)
    delete female[0][i];
  for(unsigned i=0;i<female[1].size();++i)
    delete female[1][i];
  for(unsigned i=0;i<male[0].size();++i)
    delete male[0][i];
  for(unsigned i=0;i<male[1].size();++i)
    delete male[1][i];
  male[0].clear();
  male[1].clear();
  female[0].clear();
  female[1].clear();
  delete [] poptab;
  popID = a.popID;
  nbPop = a.nbPop;
  popsize = a.popsize;
  nbFemale = a.nbFemale;
  migrationType = a.migrationType;
  nbMigrantsF = a.nbMigrantsF;
  nbMigrantsM = a.nbMigrantsM;
  generation = a.generation;
  poptab = new Population*[nbPop];
  for(int i=0;i<nbPop;++i)
    {
      poptab[i] = new Population(*(a.poptab[i]));
    }
  netmating = a.netmating;
  Init_copy(a);
  pmate = a.pmate;
  pmig = a.pmig;
  nbSources = a.nbSources;
}

Metapop::Metapop(std::string & filename)
{
}

Metapop::~Metapop()
{
  for(unsigned i=0;i<nbPop;++i)
    delete poptab[i];;
  if(nbPop >0)
    delete [] poptab;
  for(unsigned i=0;i<female[0].size();++i)
    delete female[0][i];
  for(unsigned i=0;i<female[1].size();++i)
    delete female[1][i];
  for(unsigned i=0;i<male[0].size();++i)
    delete male[0][i];
  for(unsigned i=0;i<male[1].size();++i)
    delete male[1][i];
}

Metapop& Metapop::operator = (Metapop const & a)
{
  if(this != &a)
    {
      for(unsigned i=0;i<female[0].size();++i)
	delete female[0][i];
      for(unsigned i=0;i<female[1].size();++i)
	delete female[1][i];
      for(unsigned i=0;i<male[0].size();++i)
	delete male[0][i];
      for(unsigned i=0;i<male[1].size();++i)
	delete male[1][i];
      male[0].clear();
      male[1].clear();
      female[0].clear();
      female[1].clear();
      delete [] poptab;
      popID = a.popID;
      nbPop = a.nbPop;
      popsize = a.popsize;
      nbFemale = a.nbFemale;
      migrationType = a.migrationType;
      nbMigrantsF = a.nbMigrantsF;
      nbMigrantsM = a.nbMigrantsM;
      generation = a.generation;
      poptab = new Population*[nbPop];
      for(int i=0;i<nbPop;++i)
	{
	  poptab[i] = new Population(*(a.poptab[i]));
	}
      netmating = a.netmating;
      Init_copy(a);
      pmate = a.pmate;
      pmig = a.pmig;
      nbSources = a.nbSources;
    }
  return *this;
}

void Metapop::set_nbPop(int npop, std::vector<int> nindiv)
{
  int nindivSum = 0;
  for(int i(0);i<npop;++i)
    nindivSum += nindiv[i];
  popsize = nindivSum;
  // TODO nbfemale
  if(nbPop != npop)
    {
      Population ** oldPop = poptab;
      poptab = new Population*[npop];
      if(nbPop<npop)
	{
	  for(int i(0);i<nbPop;++i)
	    poptab[i] = oldPop[i];
	  for(int i(nbPop);i<npop;++i)
	    {
	      Population A(nindiv[i]);
	      poptab[i] = &A;
	    }
	}
      else
	for(int i(0);i<npop;++i)
	  poptab[i] = oldPop[i];
      delete [] oldPop;
      nbPop = npop;
    }
}

void Metapop::set_nbSources(int a)
{
  nbSources = a;
}

void Metapop::Init(int nstep)
{
  int maxpop = 0;
  for(int i=0;i<nbPop;++i)
    {
      poptab[i]->resize(nstep);
      maxpop +=  poptab[i]->maxsize;
      maxpopsizes.push_back(poptab[i]->maxsize);
    }
  if(maxpop>1000000)
    {
      cerr<< " ERROR : growth rate and nb of step will lead you to a population size at the end of size " << maxpop << " This is too big ! " << endl;
      exit(2);
    }
  female.resize(2);
  male.resize(2);
  for(int i=0; i<maxpop; ++i)
    {
      female[0].push_back(new Human);
      female[1].push_back(new Human);
      male[0].push_back(new Human);
      male[1].push_back(new Human);
    }
  int offset = 0;
  int counterIDNetwork = 0;
  std::vector<int> nbFemales;
  std::vector<int> popsizes;
  for(int k=0;k<nbPop;++k)
    {
      maxpop = maxpopsizes[k];
      poptab[k]->offsetNetwork = offset;
      for(int j=0;j<maxpop;++j)
	{
	  poptab[k]->female[0][j]=female[0][j+offset];
	  poptab[k]->female[1][j]=female[1][j+offset];
	  poptab[k]->male[0][j]=male[0][j+offset];
	  poptab[k]->male[1][j]=male[1][j+offset];
	}
      nbFemales.push_back(poptab[k]->nbFemale);
      popsizes.push_back(poptab[k]->popsize);
      // init genealogy network
      for(int j=0;j<(poptab[k]->nbFemale);++j)
	{
	  female[0][j+offset]->IDNetwork = counterIDNetwork;
	  ++counterIDNetwork;
	}
      for(int j=0;j<(poptab[k]->popsize-poptab[k]->nbFemale);++j) // male
	{
	  male[0][j+offset]->IDNetwork = counterIDNetwork;
	  ++counterIDNetwork;
	}
      offset += maxpop;
    }
  netmating.set_nbPop(nbPop,maxpopsizes,popsizes,nbFemales);
  netmating.Init(nbFemales,popsizes,generation);
  if(migrationType>=3)
    if(migrationType==6)
      create_SPA();
    else      
      create_APA();
}

void Metapop::Init_succ(int nstep)
{
  int maxpop = 0;
  for(int i=0;i<nbPop;++i)
    {
      poptab[i]->resize(nstep);
    }
  if(maxpop>1000000)
    {
      cerr<< " ERROR : growth rate and nb of step will lead you to a population size at the end of size " << maxpop << " This is too big ! " << endl;
      exit(2);
    }
}

void Metapop::Init_copy(const Metapop & a)
{
  int maxpop = 0;
  maxpopsizes = a.maxpopsizes;
  for(int i(0);i<nbPop;++i)
    maxpop += a.maxpopsizes[i];
  if(maxpop>1000000)
    {
      cerr<< " ERROR : growth rate and nb of step will lead you to a population size at the end of size " << maxpop << " This is too big ! " << endl;
      exit(2);
    }
  female[0].clear();
  female[1].clear();
  male[0].clear();  
  male[1].clear();
  for(int i=0; i<maxpop; ++i)
    {
      female[0].push_back(new Human(*(a.female[0][i])));
      female[1].push_back(new Human(*(a.female[1][i])));
      male[0].push_back(new Human(*(a.male[0][i])));
      male[1].push_back(new Human(*(a.male[1][i])));
    }
  int offset = 0;
  int counterIDNetwork = 0;
  std::vector<int> nbFemales;
  std::vector<int> popsizes;
  for(int k=0;k<nbPop;++k)
    {
      //      maxpop = poptab[i]->maxsize;
      maxpop = a.maxpopsizes[k];
      poptab[k]->offsetNetwork = offset;
      for(int j=0;j<maxpop;++j)
	{
	  poptab[k]->female[0][j]=female[0][j+offset];
	  poptab[k]->female[1][j]=female[1][j+offset];
	  poptab[k]->male[0][j]=male[0][j+offset];
	  poptab[k]->male[1][j]=male[1][j+offset];
	}
      nbFemales.push_back(poptab[k]->nbFemale);
      popsizes.push_back(poptab[k]->popsize);
      // init genealogy network
      for(int j=0;j<(poptab[k]->nbFemale);++j)
	{
	  female[generation&1][j+offset]->IDNetwork = counterIDNetwork;
	  ++counterIDNetwork;
	}
      for(int j=0;j<(poptab[k]->popsize-poptab[k]->nbFemale);++j) // male
	{
	  male[generation&1][j+offset]->IDNetwork = counterIDNetwork;
	  ++counterIDNetwork;
	}
      offset += maxpop;
    }
  netmating.Init(nbFemales,popsizes,generation);
}


void Metapop::reset(int nindiv,std::vector<int> nindivL)
{
  generation = 0;
  std::vector<int> nbFemales;
  int maxpop = 0;
  for(int i(0);i<nbPop;++i)
    nbFemales.push_back((int) (nindivL[i] * 0.5));
  netmating.reset(nbFemales,nindivL);
  for(int i(0);i<nbPop;++i)
    {
      poptab[i]->reset(nindivL[i]);
    }
  for(unsigned i(0);i<female[0].size();++i)
    {
      female[0][i]->reset();
      female[1][i]->reset();
    }
  for(unsigned i(0);i<male[0].size();++i)
    {
      male[0][i]->reset();
      male[1][i]->reset();
    }
  int offset = 0;
  int counterIDNetwork = 0;
  std::vector<int> popsizes;
  for(int k=0;k<nbPop;++k)
    {
      maxpop = maxpopsizes[k];
      poptab[k]->offsetNetwork = offset;
      for(int j=0;j<maxpop;++j)
	{
	  poptab[k]->female[0][j]=female[0][j+offset];
	  poptab[k]->female[1][j]=female[1][j+offset];
	  poptab[k]->male[0][j]=male[0][j+offset];
	  poptab[k]->male[1][j]=male[1][j+offset];
	}
      popsizes.push_back(poptab[k]->popsize);
      // init genealogy network
      for(int j=0;j<(poptab[k]->nbFemale);++j)
	{
	  female[0][j+offset]->IDNetwork = counterIDNetwork;
	  ++counterIDNetwork;
	}
      for(int j=0;j<(poptab[k]->popsize-poptab[k]->nbFemale);++j) // male
	{
	  male[0][j+offset]->IDNetwork = counterIDNetwork;
	  ++counterIDNetwork;
	}
      offset += maxpop;
    }
  netmating.Init(nbFemales,popsizes,generation);
  popsize = nindiv;
  if(migrationType>=3)
    create_APA();
}

void Metapop::evolve(int nbGene,bool verbose)
{
    for(int i(0);i<nbGene;++i){
    if(popsize>0)
      {
	evolve_one();    
	if(verbose)
	  {
	    cout <<"start generation " << generation <<endl;
	  }
      }
    }
}

void Metapop::evolve_one()
{
  for(int i(0);i<nbPop;++i)
    {
      poptab[i]->random_shuffle_pop();
      poptab[i]->test(i);
    }
  mutation();
  demographic_update();  
  migration();
  for(int i(0);i<nbPop;++i)
    {
      if(nbPop>1)
	if((poptab[i]->popsize==poptab[i]->nbFemale)||(poptab[i]->nbFemale==0))
	  {
	    int r=myrand(nbPop);
	    if(r==i){
	      if(i<nbPop-1)
		++r;
	      else
		--r;}
	    split_migrate(i,r);		
	  }
    }
  for(int i(0);i<nbPop;++i)
    {
      poptab[i]->set_maleIDnetwork(netmating,i);
      poptab[i]->set_femaleIDnetwork(netmating,i);
    }
  update_netmating();
  reproduction();
  ++generation;
  for(int i(0);i<nbPop;++i)
    {
      ++ poptab[i]->generation;
    }
}

void Metapop::burnin(double diversityThreshold,int gen,double facMu,int dnaType)
{
  if(diversityThreshold!=0)
    {
      double  mt1,mt2,x1,x2,y1,y2,a1,a2;
      // save real paramaters
      mt1 = RATES.muMt1;
      mt2 = RATES.muMt2;
      x1 = RATES.muX1;
      x2 = RATES.muX2;
      y1 = RATES.muY1;
      y2 = RATES.muY2;
      a1 = RATES.muA1;
      a2 = RATES.muA2;
      // burnin parameters
      RATES.muMt1 *= facMu;
      RATES.muMt2 *= facMu;
      RATES.muX1 *= facMu;
      RATES.muX2 *= facMu;
      RATES.muY1 *= facMu;
      RATES.muY2 *= facMu;
      RATES.muA1 *= facMu;
      RATES.muA2 *= facMu;
      Population BurninPop(popsize);
      for(int i(0);i<female[0].size();++i)
	{
	  BurninPop.female[0][i]=female[0][i];
	  BurninPop.female[1][i]=female[1][i];
	}
      for(int i(0);i<male[0].size();++i)
	{
	  BurninPop.male[0][i]=male[0][i];
	  BurninPop.male[1][i]=male[1][i];
	}
      BurninPop.nbFemale = nbFemale;
      for(int g(0);g<gen;++g)
	{
	  BurninPop.random_shuffle_pop();
	  BurninPop.mutation();
	  BurninPop.reproduction_burnin();
	  for(int i(0);i<nbPop;++i)
	    {
	      ++(poptab[i]->generation);
	    }
	  ++(BurninPop.generation);
	  ++generation;
	  popsize = BurninPop.popsize;
	  nbFemale = BurninPop.nbFemale;
	}
      RATES.muMt1 = mt1;
      RATES.muMt2 = mt2;
      RATES.muX1 = x1;
      RATES.muX2 = x2;
      RATES.muY1 = y1;
      RATES.muY2 = y2;
      RATES.muA1 = a1;
      RATES.muA2 = a2;      
      for(int g(0);g<diversityThreshold;++g)
	{
	  BurninPop.random_shuffle_pop();
	  BurninPop.mutation();
	  BurninPop.reproduction_burnin();
	  for(int i(0);i<nbPop;++i)
	    {
	      ++(poptab[i]->generation);
	    }
	  ++(BurninPop.generation);
	  ++generation;
	  popsize = BurninPop.popsize;
	  nbFemale = BurninPop.nbFemale;
	}
      assert(BurninPop.generation == poptab[0]->generation);
      assert(BurninPop.generation == generation);
      int nbF=BurninPop.nbFemale;
      int nbFperP = floor(nbF/nbPop);
      int countercounter=BurninPop.maxsize - 1;
      for(int popid(0);popid<nbPop-1;++popid)
	{
	  for(int j(0);j<nbFperP;++j)
	    {
	      poptab[popid]->female[generation&1][j] = BurninPop.female[generation&1][j+popid*nbFperP];	    
	    }
	  for(int j(nbFperP);j<poptab[popid]->maxsize;++j)
	    {
	      poptab[popid]->female[generation&1][j] = BurninPop.female[generation&1][countercounter];
	      --countercounter;
	    }
	  poptab[popid]->nbFemale = nbFperP;
	  poptab[popid]->popsize = nbFperP;
	}
      int popid = nbPop-1;
      int nbFthisP = nbF-nbFperP*(nbPop-1);
      for(int j(0);j<nbFthisP;++j)
	{
	  poptab[popid]->female[generation&1][j] = BurninPop.female[generation&1][j+popid*nbFperP];
	}
      for(int j(nbFthisP);j<poptab[popid]->maxsize;++j)
	{
	  poptab[popid]->female[generation&1][j] = BurninPop.female[generation&1][countercounter];
	  --countercounter;
	}
      poptab[popid]->nbFemale = nbFthisP;
      poptab[popid]->popsize = nbFthisP;
      int nbM=BurninPop.popsize-BurninPop.nbFemale;
      int nbMperP = floor(nbM/nbPop);
      countercounter = BurninPop.maxsize - 1;
      for(int popid(0);popid<nbPop-1;++popid)
	{
	  for(int j(0);j<nbMperP;++j)
	    {
	      poptab[popid]->male[generation&1][j] = BurninPop.male[generation&1][j+popid*nbMperP];
	    }	 
	  for(int j(nbMperP);j<(poptab[popid]->maxsize);++j)
	    {
	      poptab[popid]->male[generation&1][j] = BurninPop.male[generation&1][countercounter];
	      --countercounter;
	    }
	  poptab[popid]->popsize += nbMperP;
	}
      popid = nbPop-1;
      int nbMthisP = nbM-nbMperP*(nbPop-1);
      for(int j(0);j<nbMthisP;++j)
	{
	  poptab[popid]->male[generation&1][j] = BurninPop.male[generation&1][j+popid*nbMperP];
	}
      for(int j(nbMthisP);j<(poptab[popid]->maxsize);++j)
	{
	  poptab[popid]->male[generation&1][j] = BurninPop.male[generation&1][countercounter];
	  --countercounter;
	}
      poptab[popid]->popsize += nbMthisP;      
      RATES.muMt1 = mt1;
      RATES.muMt2 = mt2;
      RATES.muX1 = x1;
      RATES.muX2 = x2;
      RATES.muY1 = y1;
      RATES.muY2 = y2;
      RATES.muA1 = a1;
      RATES.muA2 = a2;
      std::vector<int> nbFemales;
      std::vector<int> popsizes;
      for(int k=0;k<nbPop;++k)
	{
	  nbFemales.push_back(poptab[k]->nbFemale);
	  popsizes.push_back(poptab[k]->popsize);
	}
      netmating.Init2();
      int counterIDNetwork = (generation% netmating.storedGeneration)*netmating.nnpg;
      for(int k=0;k<nbPop;++k)
	{
	  for(int j=0;j<(poptab[k]->nbFemale);++j)
	    {
	      poptab[k]->female[generation&1][j]->IDNetwork = counterIDNetwork;
	      ++counterIDNetwork;
	      netmating.myg[poptab[k]->female[generation&1][j]->IDNetwork].generation=generation;
	      netmating.myg[poptab[k]->female[generation&1][j]->IDNetwork].ID=j;
	      netmating.myg[poptab[k]->female[generation&1][j]->IDNetwork].sex=1;
	      netmating.myg[poptab[k]->female[generation&1][j]->IDNetwork].population=k;
	    }
	  for(int j=0;j<(poptab[k]->popsize-poptab[k]->nbFemale);++j) // male
	    {
	      poptab[k]->male[generation&1][j]->IDNetwork = counterIDNetwork;
	      ++counterIDNetwork;
	      netmating.myg[poptab[k]->male[generation&1][j]->IDNetwork].generation=generation;
	      netmating.myg[poptab[k]->male[generation&1][j]->IDNetwork].ID=j;
	      netmating.myg[poptab[k]->male[generation&1][j]->IDNetwork].population=k;
	      netmating.myg[poptab[k]->male[generation&1][j]->IDNetwork].sex=0;
	    }
	  poptab[k]->set_maleIDnetwork(netmating,k);
	  poptab[k]->set_femaleIDnetwork(netmating,k);
	}
      if((migrationType==3)||(migrationType==4))
	create_APA();
    }
}

void Metapop::mutation()
{
  for(int i(0);i<nbPop;++i)
    {
      poptab[i]->mutation();
    }
}

void Metapop::reproduction()
{ // make the reproduction inside each population (i.e. must have migrated before or can't mate)
  // plus update the total population size
  int gparent = generation % netmating.storedGeneration;
  netmating.clear_premating(gparent);
  int popsizeNew = 0;
  int nbFemaleNew = 0;
  for(int i(0); i<nbPop; ++i)
    {
      poptab[i]->reproduction(netmating,i,pmate,pmig);
      popsizeNew += poptab[i]->popsize;
      nbFemaleNew += poptab[i]->nbFemale;
    }
  popsize = popsizeNew;
  nbFemale = nbFemaleNew;
  netmating.clear_postmating();
}

void Metapop::reproduction_burnin()
{ // make the reproduction inside each population (i.e. must have migrated before or can't mate)
  // plus update the total population size
  int popsizeNew = 0;
  int nbFemaleNew = 0;
  for(int i(0); i<nbPop; ++i)
    {
      poptab[i]->reproduction_burnin();
      popsizeNew += poptab[i]->popsize;
      nbFemaleNew += poptab[i]->nbFemale;
    }
  popsize = popsizeNew;
  nbFemale = nbFemaleNew;
}

void Metapop::update_netmating()
{
  int psize(0),fsize(0),currentG(generation&1);
  for(int n(0); n<nbPop; ++n)
    {
      psize = poptab[n]->popsize;
      fsize = poptab[n]->nbFemale;
      for(int i(0); i<fsize; ++i)
	{
	  if(poptab[n]->female[currentG][i]->IDNetwork<0)
	    cout << "ERROR ID < 0 pop " << n << " f " << i<< " fsize " << fsize <<endl;
	  else
	    {
	      netmating.myg[boost::vertex(poptab[n]->female[currentG][i]->IDNetwork,netmating.myg)].population = n;
	      netmating.myg[boost::vertex(poptab[n]->female[currentG][i]->IDNetwork,netmating.myg)].ID = i;}
	}
      for(int i(0); i<(psize - fsize); ++i)
	{
	  if(poptab[n]->male[currentG][i]->IDNetwork<0)
	    cout << "ERROR ID <0 pop " << n << " m " << i<< " msize " << psize-fsize << endl;
	  else
	    {
	      netmating.myg[boost::vertex(poptab[n]->male[currentG][i]->IDNetwork,netmating.myg)].population = n;
	      netmating.myg[boost::vertex(poptab[n]->male[currentG][i]->IDNetwork,netmating.myg)].ID = i;
	    }
	}
    }
}

void Metapop::demographic_update()
{
  for(int i(0); i<nbPop; ++i)
    {
      poptab[i]->demographic_update(); // update popsizedemo which will be used in the following generation as the new population size
    }
}

void Metapop::migration()
{
  //  std::vector< std::vector<int> > nbMigrantsM(nbPop); create as an attribute
  int r;
  int nbMF = 0;
  int nbMM = 0;
  int indexlastn = 0;
  int indexlastk = 0;
  int currentG = generation&1;
  Human * tmp;
  if(migrationType == 1)
    {
      int oldIndivF[nbPop]; 
      // keep track of who was already in the population before the migration started 
      // the old indiv are concentrated at the fron of the array
      int oldIndivM[nbPop];
      int indexoldlastn;
      for(int n(0);n<nbPop;++n)
	{
	  oldIndivF[n] = poptab[n]->nbFemale;
	  oldIndivM[n] = poptab[n]->popsize - poptab[n]->nbFemale;
	  for(int k(0);k<nbPop;++k)
	    {
	      nbMigrantsF[n][k] =  rBinom(poptab[n]->nbFemale,migrationRateF[n][k]);
	      nbMigrantsM[n][k] = rBinom(poptab[n]->popsize - poptab[n]->nbFemale,migrationRateM[n][k]);
	    }
	}
      for(int n(0);n<nbPop;++n)
	for(int k(0);k<nbPop;++k)
	  {
	    // FEMALE migration
	    nbMF = nbMigrantsF[n][k];
	    indexlastn = poptab[n]->nbFemale - 1;
	    indexoldlastn = oldIndivF[n]-1;
	    indexlastk = poptab[k]->nbFemale - 1;
	    for(int i(0);i<nbMF;++i)
	      {
		r = myrand(indexoldlastn);
		tmp = poptab[k]->female[currentG][indexlastk+1];
		poptab[k]->female[currentG][indexlastk+1] = poptab[n]->female[currentG][r];
		poptab[n]->female[currentG][r] = poptab[n]->female[currentG][indexoldlastn];
		poptab[n]->female[currentG][indexoldlastn] = poptab[n]->female[currentG][indexlastn];
		poptab[n]->female[currentG][indexlastn] = tmp;
		--indexlastn;
		--indexoldlastn;
		++indexlastk;
	      }
	    poptab[n]->nbFemale = poptab[n]->nbFemale - nbMF;
	    poptab[k]->nbFemale = poptab[k]->nbFemale + nbMF;
	    poptab[n]->popsize = poptab[n]->popsize - nbMF;
	    poptab[k]->popsize = poptab[k]->popsize + nbMF;
	    oldIndivF[n] -= nbMF;
	    // MALE migration
	    nbMM = nbMigrantsM[n][k];
	    indexlastn = poptab[n]->popsize - poptab[n]->nbFemale - 1;
	    indexoldlastn = oldIndivM[n] - 1;
	    indexlastk = poptab[k]->popsize - poptab[k]->nbFemale - 1;
	    for(int i(0);i<nbMM;++i)
	      {
		r = myrand(indexoldlastn);
		tmp = poptab[k]->male[currentG][indexlastk+1];
		poptab[k]->male[currentG][indexlastk+1] = poptab[n]->male[currentG][r];
		poptab[n]->male[currentG][r] = poptab[n]->male[currentG][indexoldlastn];
		poptab[n]->male[currentG][indexoldlastn] = poptab[n]->male[currentG][indexlastn];
		poptab[n]->male[currentG][indexlastn] = tmp;
		--indexlastn;
		--indexoldlastn;
		++indexlastk;
	      }
	    poptab[n]->popsize = poptab[n]->popsize - nbMM;
	    poptab[k]->popsize = poptab[k]->popsize + nbMM;
	    oldIndivM[n] -= nbMM;
	  }
      // do some kind of migration from a migration matrix
    }
  else
    if(migrationType == 2)
      { // must have a true circle
	std::vector< std::vector<Human* > > oldF(nbPop);
	std::vector< std::vector<Human* > > oldM(nbPop);
	int oldFn[nbPop];
	int oldMn[nbPop];
	for(int n(0);n<nbPop;++n)
	  {
	    oldF[n] = poptab[n]->female[currentG];
	    oldM[n] = poptab[n]->male[currentG];
	    oldFn[n] = poptab[n]->nbFemale;
	    oldMn[n] = poptab[n]->popsize - poptab[n]->nbFemale;
	  }	
	for(int n(0);n<nbPop;++n)
	  for(int k(0);k<nbPop;++k)
	    {
	      // FEMALE migration from k to n
	      if(migrationRateF[n][k] == 1)
		{
		  poptab[k]->female[currentG] = oldF[n];
		  int diff = oldFn[n] - oldFn[k];
		  poptab[k]->nbFemale = oldFn[n];
		  poptab[k]->popsize += diff;
		}
	      // MALE migration
	      if(migrationRateM[n][k] == 1)
		{
		  poptab[k]->male[currentG] = oldM[n];
		  int nbMaleN = oldMn[n]; // poptab[n]->popsize - poptab[n]->nbFemale ;
		  int nbMaleK = oldMn[k]; // poptab[k]->popsize - poptab[k]->nbFemale ;
		  int diff = nbMaleN - nbMaleK;
		  poptab[n]->popsize += diff;
		}
	    }
      }
  if(migrationType == 3)
    { // FEMALE MIGRATION FROM MALE FAMILY DEST
      std::vector< std::vector< Human* > > oldF(nbPop);
      int oldFn[nbPop];
      int oldMn[nbPop];
      int index1;
      int newpopID;
      vector<Human*> freeFemale;
      int totNBmale=0;
      for(int n(0);n<nbPop;++n)
	{
	  oldF[n] = poptab[n]->female[currentG];
	  oldFn[n] = poptab[n]->nbFemale;
	  oldMn[n] = poptab[n]->popsize - poptab[n]->nbFemale;
	  totNBmale += oldMn[n];
	  poptab[n]->nbFemale=0;
	  poptab[n]->popsize=oldMn[n];
	}
      int nbchanges = rBinom(totNBmale,pmig);
      for(int i(0);i<nbchanges;++i)
	{
	  int rpop = myrand(nbPop);
	  int rindiv = myrand(oldMn[rpop]);
	  int newD = myrand(nbPop);
	  poptab[rpop]->male[currentG][rindiv]->destinationPop=newD;
	}
      for(int n(0);n<nbPop;++n)
	{
	  // FEMALE migration from k to n
	  for(int p(0);p<oldFn[n];++p)
	    {
	      newpopID=oldF[n][p]->destinationPop;
	      index1=poptab[newpopID]->nbFemale;
	      if(index1<(poptab[newpopID]->maxsize))
		{
		  poptab[newpopID]->female[currentG][index1]=oldF[n][p];
		  poptab[newpopID]->nbFemale+=1;
		  poptab[newpopID]->popsize+=1;
		}
	      else
		{
		  for(int kn(0);kn<nbPop;++kn)
		    {
		      if(poptab[kn]->nbFemale <(poptab[kn]->maxsize))
			{
			  index1 = poptab[kn]->nbFemale;
			  poptab[kn]->female[currentG][index1]=oldF[n][p];
			  poptab[kn]->nbFemale+=1;
			  poptab[kn]->popsize+=1;
			  kn = nbPop; // breal loop
			}
		    }
		}
	    }
	  for(int p(oldFn[n]);p<oldF[n].size();++p)
	    {
	      freeFemale.push_back(oldF[n][p]);
	    }
	}
      for(int n(0);n<nbPop;++n)
	{	  
	  for(int i(poptab[n]->nbFemale);i<oldF[n].size();++i)
	    {
	      poptab[n]->female[currentG][i] = freeFemale.back();
	      freeFemale.pop_back();
	    }
	}
      assert(freeFemale.size()==0);
    }
  if((migrationType == 4) || (migrationType == 6))
    { // FEMALE MIGRATION FROM MALE FAMILY DEST
      std::vector< std::vector<Human* > > oldF(nbPop);
      std::vector<int> oldFn(nbPop);
      std::vector<int> oldMn(nbPop);
      int index1;
      int newpopID;
      vector<Human*> freeFemale;
      vector<int> destinationCount(nbPop,0);
      vector<int> source(nbPop,0);
      int biggestpop;
      int totNBmale=0;
      int totNBfemale=0;
      for(int n(0);n<nbPop;++n)
	{
	  oldF[n] = poptab[n]->female[currentG];
	  oldFn[n] = poptab[n]->nbFemale;
	  oldMn[n] = poptab[n]->popsize - poptab[n]->nbFemale;
	  totNBmale += oldMn[n];
	  totNBfemale += oldFn[n];
	  poptab[n]->nbFemale=0;
	  poptab[n]->popsize=oldMn[n];
	}
      int nbchanges = rBinom(totNBfemale,pmig);
      for(int i(0);i<nbchanges;++i)
	{
	  int rpop = myrand(nbPop);
	  int rindiv = myrand(oldFn[rpop]);
	  int newD = myrand(nbPop);
	  poptab[rpop]->female[currentG][rindiv]->destinationPop=newD;
	}
      for(int n(0);n<nbPop;++n)
	{
	  for(int i(0);i<oldMn[n];++i)
	    {
	      destinationCount[poptab[n]->male[currentG][i]->destinationPop]+=1;
	      source[poptab[n]->male[currentG][i]->destinationPop]=n;
	    }
	}
      for(int n(0);n<nbPop;++n)
	if(destinationCount[n]<20)
	  {
	    int rnb = 10 - destinationCount[n];
	    int rpop = source[n];
	    for(int k(0);k<rnb;++k)
	      {
		int rindiv = myrand(oldMn[rpop]);
		poptab[rpop]->male[currentG][rindiv]->destinationPop=n;
	      }
	  }
      for(int n(0);n<nbPop;++n)
	{
	  // FEMALE migration from k to n
	  for(int p(0);p<oldFn[n];++p)
	    {
	      newpopID=oldF[n][p]->destinationPop;
	      index1=poptab[newpopID]->nbFemale;
	      if(index1 < (poptab[newpopID]->maxsize))
		{
		  poptab[newpopID]->female[currentG][index1]=oldF[n][p];
		  poptab[newpopID]->nbFemale+=1;
		  poptab[newpopID]->popsize+=1;
		}
	      else
		for(int kn(0);kn<nbPop;++kn)
		  {
		    if(poptab[kn]->nbFemale <(poptab[kn]->maxsize))
		      {
			index1 = poptab[kn]->nbFemale;
			poptab[kn]->female[currentG][index1]=oldF[n][p];
			poptab[kn]->nbFemale+=1;
			poptab[kn]->popsize+=1;
			kn = nbPop; // breal loop
		      }
		  }
	    }
	  for(int p(oldFn[n]);p<oldF[n].size();++p)
	    {
	      freeFemale.push_back(oldF[n][p]);
	    }
	}
      for(int n(0);n<nbPop;++n)
	{	  
	  for(int i(poptab[n]->nbFemale);i<oldF[n].size();++i)
	    {
	      poptab[n]->female[currentG][i] = freeFemale.back();
	      freeFemale.pop_back();
	    }
	}
      assert(freeFemale.size()==0);
    }
  if(migrationType == 5)
    { // MALE MIGRATION FROM FEMALE FAMILY DEST
      std::vector< std::vector< Human* > > oldM(nbPop);
      std::vector<int > oldFn(nbPop);
      std::vector<int > oldMn(nbPop);
      int index1;
      int newpopID;
      vector<Human*> freeMale;
      vector<int> destinationCount(nbPop,0);
      vector<int> source(nbPop,0);
      int biggestpop;
      int totNBfemale=0;
      for(int n(0);n<nbPop;++n)
	{
	  oldM[n] = poptab[n]->male[currentG];
	  oldFn[n] = poptab[n]->nbFemale;
	  oldMn[n] = poptab[n]->popsize - poptab[n]->nbFemale;
	  totNBfemale += oldFn[n];
	  poptab[n]->popsize = oldFn[n];
	}
      int nbchanges = rBinom(totNBfemale,pmig);
      for(int i(0);i<nbchanges;++i)
	{
	  int rpop = myrand(nbPop);
	  int rindiv = myrand(oldFn[rpop]);
	  int newD = myrand(nbPop);
	  poptab[rpop]->female[currentG][rindiv]->destinationPop=newD;
	}
      for(int n(0);n<nbPop;++n)
	{
	  for(int i(0);i<oldFn[n];++i)
	    {
	      destinationCount[poptab[n]->female[currentG][i]->destinationPop]+=1;
	      source[poptab[n]->female[currentG][i]->destinationPop]=n;
	    }
	}
      // put a minimum for the number of male sending their women to a population otherwise this population will die, here minimum is 10
      for(int n(0);n<nbPop;++n)
	if(destinationCount[n]<20)
	  {
	    int rnb = 10 - destinationCount[n];
	    int rpop = source[n];
	    for(int k(0);k<rnb;++k)
	      {
		int rindiv = myrand(oldMn[rpop]);
		poptab[rpop]->female[currentG][rindiv]->destinationPop=n;
	      }
	  }
      for(int n(0);n<nbPop;++n)
	{
	  // MALE migration from k to n
	  for(int p(0);p<oldMn[n];++p)
	    {
	      newpopID=oldM[n][p]->destinationPop;
	      index1=poptab[newpopID]->popsize-poptab[newpopID]->nbFemale;
	      if(index1 < (poptab[newpopID]->maxsize))
		{
		  poptab[newpopID]->male[currentG][index1]=oldM[n][p];
		  poptab[newpopID]->popsize+=1;
		}
	      else
		for(int kn(0);kn<nbPop;++kn)
		  {
		    index1 = poptab[kn]->popsize - poptab[kn]->nbFemale;
		    if(index1 <(poptab[kn]->maxsize))
		      {			  
			poptab[kn]->male[currentG][index1]=oldM[n][p];
			poptab[kn]->popsize+=1;
			kn = nbPop; // breal loop
		      }
		    if(kn==(nbPop-1))
		      cerr<< "FAIL"<< endl;
		  }
	    }
	  for(int p(oldMn[n]);p<oldM[n].size();++p)
	    {
	      freeMale.push_back(oldM[n][p]);
	    }
	}
      for(int n(0);n<nbPop;++n)
	{
	  for(int i(poptab[n]->popsize-poptab[n]->nbFemale);i<oldM[n].size();++i)
	    {
	      poptab[n]->male[currentG][i] = freeMale.back();
	      freeMale.pop_back();
	    }
	}
      assert(freeMale.size()==0);
    }
  // total popsize and nb of females should not change
  // if migrtaionType == 0 : no migration
}

void Metapop::split_migrate(int popID1, int popID2)
{// if 1 population is empty then another population ( random) will split in half and recolonize it
  // here pop1 split and moves in pop2
  int currentG = generation&1;
  int nf1 = poptab[popID1]->nbFemale;
  int nm1 = poptab[popID1]->popsize-nf1;
  int nf2 = poptab[popID2]->nbFemale;
  int nm2 = poptab[popID2]->popsize-nf2;
  int migf = nf2/2;
  int migm = nm2/2;
  Human * tmp;
  for(int i(0);i<migf;++i)
    {
      tmp = poptab[popID1]->female[currentG][nf1];
      poptab[popID1]->female[currentG][nf1] = poptab[popID2]->female[currentG][nf2-1];
      poptab[popID2]->female[currentG][nf2-1] = tmp;
      ++nf1;
      --nf2;
    }
  poptab[popID1]->nbFemale += migf;
  poptab[popID2]->nbFemale -= migf;
  poptab[popID1]->popsize += migf;
  poptab[popID2]->popsize -= migf;
  for(int i(0);i<migm;++i)
    {
      tmp = poptab[popID1]->male[currentG][nm1];
      poptab[popID1]->male[currentG][nm1] = poptab[popID2]->male[currentG][nm2-1];
      poptab[popID2]->male[currentG][nm2-1] = tmp;
      ++nm1;
      --nm2;
    }
  poptab[popID1]->popsize += migm;
  poptab[popID2]->popsize -= migm;
}
      
void Metapop::create_circleMigration()
{
}

void Metapop::create_APA()
{
  if(nbSources>=nbPop)
    {
      cerr << " You can not have an APA or an SPA with only " << nbPop << " populations" << endl;
      exit(1);
    }
  assert(nbSources<nbPop);
  bool test = true;
  vector< vector<int> > sources;
  vector< vector<int> > sinks;
  // first circle = simple
  int nbEdges= 0;
  int iteration =0;
  migrationRateF.clear();
  std::vector<int> available_sources(nbPop);
  for(int n(0);n<nbPop;++n)
    {
      //      migrationRateF[n].clear();
      poptab[n]->migrationType = migrationType;
      vector< int > tmp2;
      vector< double > tmp3;
      sources.push_back(tmp2);
      sinks.push_back(tmp2);
      for(int k(0);k<nbPop;++k)
	tmp3.push_back(0);
      migrationRateF.push_back(tmp3);
      available_sources[n]=n;
    }
  // gives a npop * npop matrix
  while(test)
    {
      int avail_size = available_sources.size();
      int circlesize = myrand(avail_size-3)+3; // circle of minimum size 3
      std::vector<int> circle = sample_without_replacement(avail_size,circlesize);
      bool ok = true;
      ++iteration;
      for(int k(0);k<circlesize;++k)
	{
	  int source = available_sources[circle[k]];
	  int sink = available_sources[circle[(k+1)%circlesize]];
	  if((migrationRateF[source][sink]!=0)||(source==sink))
	    {
	      ok = false;
	    }
	  if((sources[sink].size()>=nbSources)||(sinks[source].size()>=nbSources))
	    {
	      ok = false;
	    }
	}
      if(ok)
	{
	  vector<int> toerase;
	  for(int k(0);k<circlesize;++k)
	    {
	      int source = available_sources[circle[k]];
	      int sink = available_sources[circle[(k+1)%circlesize]];
	      migrationRateF[source][sink]=1;
	      migrationRateF[sink][source]=1;
	      sources[sink].push_back(source);
	      sinks[source].push_back(sink);
	      if((sources[source].size()>=nbSources)&&(sinks[source].size()>=nbSources))
		toerase.push_back(source);
	    }
	  for(unsigned k(0);k<toerase.size();++k)
	    { // trouble do not erase with current order
	      for(int t(0);t<available_sources.size();++t)
		if(available_sources[t]==toerase[k])
		  {
		    available_sources.erase(available_sources.begin()+t); // swap with the last one
		    t=1000; // break loop
		  }
	    }
	  nbEdges += circlesize;
	}      
      if(((nbEdges>=(nbPop-1)*nbSources+1)||(iteration>=300))||available_sources.size()<3) //(later = cant make  a circle anymore
	test=false; 
    }
  //then check if one population is missing one source
  unsigned avail_size = available_sources.size();
  vector<int> toconnect;
  for(unsigned i(0);i<avail_size;++i)
    {
      int source = available_sources[i];
      if(sources[source].size()==0)
	toconnect.push_back(source);
    }
  if(toconnect.size()>0)
    {
      for(unsigned i(0);i<toconnect.size();++i)
	{
	  int tobc = toconnect[i];
	  int maillon = -1;
	  for(int j(0);j<nbPop;++j)
	    {
	      if((sources[j].size()>0)&&(sinks[j].size()>0))
		{
		  maillon = j;
		  j=nbPop; // break the loop
		}
	    }
	  if(maillon>=0)
	    {
	      sources[tobc].push_back(maillon);
	      sinks[tobc].push_back(sinks[maillon][0]);
	      sources[sinks[maillon][0]].push_back(tobc);
	      sinks[maillon].push_back(tobc);
	    }
	  else
	    {
	      cerr << "ERROR in the network creation" << endl;
	      exit(42);
	    }
	}
    }
  vector<int> tmp(nbPop);
  int destination,start,end;
  int currentG = generation&1;
  for(int n(0);n<nbPop;++n)
    {
      int nsources = sinks[n].size();
      int fsize = poptab[n]->nbFemale;
      int msize = poptab[n]->popsize - poptab[n]->nbFemale;      
      for(int j(0);j<nsources-1;++j)
	{
	  start=fsize/nsources*j;	  
	  end=fsize/nsources*(j+1);
	  destination = sinks[n][j];
	  for(int i(start);i<end;++i)
	    poptab[n]->female[currentG][i]->destinationPop = destination;
	  start=msize/nsources*j;	  
	  end=msize/nsources*(j+1);
	  destination = sinks[n][j];
	  for(int i(start);i<end;++i)
	    poptab[n]->male[currentG][i]->destinationPop = destination;
	}
      int j=nsources-1;
      start=fsize/nsources*j;	  
      end=fsize;
      destination = sinks[n][j];
      for(int i(start);i<end;++i)
	poptab[n]->female[currentG][i]->destinationPop = destination;
      start=msize/nsources*j;	  
      end=msize;
      destination = sinks[n][j];
      for(int i(start);i<end;++i)
	poptab[n]->male[currentG][i]->destinationPop = destination;    
      poptab[n]->test(n);
    }
}


void Metapop::create_SPA()
{
  assert(nbSources<nbPop);
  bool test = true;
  vector< vector<int> > sources;
  vector< vector<int> > sinks;
  // first circle = simple
  int nbEdges= 0;
  int iteration =0;
  migrationRateF.clear();
  std::vector<int> available_sources(nbPop);
  for(int n(0);n<nbPop;++n)
    {
      //      migrationRateF[n].clear();
      poptab[n]->migrationType = migrationType;
      vector< int > tmp2;
      vector< double > tmp3;
      sources.push_back(tmp2);
      sinks.push_back(tmp2);
      for(int k(0);k<nbPop;++k)
	tmp3.push_back(0);
      migrationRateF.push_back(tmp3);
      available_sources[n]=n;
    }
  // gives a npop * npop matrix
  while(test)
    {
      int avail_size = available_sources.size();
      int circlesize = myrand(avail_size-3)+3; // circle of minimum size 3
      std::vector<int> circle = sample_without_replacement(avail_size,circlesize);
      bool ok = true;
      ++iteration;
      for(int k(0);k<circlesize;++k)
	{
	  int source = available_sources[circle[k]];
	  int sink = available_sources[circle[(k+1)%circlesize]];
	  if((migrationRateF[source][sink]!=0))
	    {
	      ok = false;
	    }
	  if((sources[sink].size()>=nbSources)||(sinks[source].size()>=nbSources))
	    {
	      ok = false;
	    }
	}
      if(ok)
	{
	  vector<int> toerase;
	  for(int k(0);k<circlesize;++k)
	    {
	      int source = available_sources[circle[k]];
	      int sink = available_sources[circle[(k+1)%circlesize]];
	      migrationRateF[source][sink]=1;
	      migrationRateF[sink][source]=1;
	      sources[sink].push_back(source);
	      sources[source].push_back(sink);
	      sinks[source].push_back(sink);
	      sinks[sink].push_back(source);
	      if((sources[source].size()>=nbSources)&&(sinks[source].size()>=nbSources))
		toerase.push_back(source);
	    }
	  for(int k(0);k<toerase.size();++k)
	    { // trouble do not erase with current order
	      for(int t(0);t<available_sources.size();++t)
		if(available_sources[t]==toerase[k])
		  {
		    available_sources.erase(available_sources.begin()+t); // swap with the last one
		    t=1000; // break loop
		  }
	    }
	  nbEdges += circlesize*2;
	}      
      if(((nbEdges>=(nbPop-1)*nbSources+1)||(iteration>=300))||available_sources.size()<3) //(later = cant make  a circle anymore
	test=false; 
    }
  //then check if one population is missing one source
  unsigned avail_size = available_sources.size();
  vector<int> toconnect;
  for(unsigned i(0);i<avail_size;++i)
    {
      int source = available_sources[i];
      if(sources[source].size()==0)
	toconnect.push_back(source);
    }
  if(toconnect.size()>0)
    {
      for(unsigned i(0);i<toconnect.size();++i)
	{
	  int tobc = toconnect[i];
	  int maillon = -1;
	  for(int j(0);j<nbPop;++j)
	    {
	      if((sources[j].size()>0)&&(sinks[j].size()>0))
		{
		  maillon = j;
		  j=nbPop; // break the loop
		}
	    }
	  if(maillon>=0)
	    {
	      sources[tobc].push_back(maillon);
	      sinks[tobc].push_back(maillon);
	      sources[maillon].push_back(tobc);
	      sinks[maillon].push_back(tobc);
	    }
	  else
	    {
	      cerr << "ERROR in the network creation" << endl;
	      exit(42);
	    }
	}
    }
  for(int i(0); i<nbPop;++i)
    cout << "source" << i << "  -- " <<  sinks[i] <<  " destination " << sources[i] << endl;  
  vector<int> tmp(nbPop);
  int destination,start,end;
  int currentG = generation&1;
  for(int n(0);n<nbPop;++n)
    {
      int nsources = sinks[n].size();
      int fsize = poptab[n]->nbFemale;
      int msize = poptab[n]->popsize - poptab[n]->nbFemale;      
      for(int j(0);j<nsources-1;++j)
	{
	  start=fsize/nsources*j;	  
	  end=fsize/nsources*(j+1);
	  destination = sinks[n][j];
	  for(int i(start);i<end;++i)
	    poptab[n]->female[currentG][i]->destinationPop = destination;
	  start=msize/nsources*j;	  
	  end=msize/nsources*(j+1);
	  destination = sinks[n][j];
	  for(int i(start);i<end;++i)
	    poptab[n]->male[currentG][i]->destinationPop = destination;
	}
      int j=nsources-1;
      start=fsize/nsources*j;	  
      end=fsize;
      destination = sinks[n][j];
      for(int i(start);i<end;++i)
	poptab[n]->female[currentG][i]->destinationPop = destination;
      start=msize/nsources*j;	  
      end=msize;
      destination = sinks[n][j];
      for(int i(start);i<end;++i)
	poptab[n]->male[currentG][i]->destinationPop = destination;    
      poptab[n]->test(n);
    }
}


void Metapop::set_mu(double a, double b, double c,double d,double e, double f, double g,double h)
{
  RATES.muMt1 = a;
  RATES.muMt2 = b;
  RATES.muX1 = c;
  RATES.muX2 = d;
  RATES.muY1 = e;
  RATES.muY2 = f;
  RATES.muA1 = g;
  RATES.muA2 = h;
}

void Metapop::set_pmate(double a)
{
  pmate = a;
}

void Metapop::set_pmig(double a)
{
  pmig = a;
}


vector<Human> Metapop::sample_local(int nindiv, int idpop) const
{ // local sample - > sample 1 population only with idpop
  if(popsize>0)
    {
      vector<double> prop;
      vector<double> propf;
      vector<int> pops;
      prop.push_back(1.0);
      propf.push_back(0);
      pops.push_back(idpop);
      return sample(nindiv, pops, prop,propf);     
    }
  else
    return vector<Human>();
}

vector<Human> Metapop::sample_local_random(int nindiv) const
{ // local sample - > sample 1 population only with idpop
  if(popsize>0)
    {
      int idpop = myrand(nbPop);
      vector<double> prop;
      vector<double> propf;
      vector<int> pops;
      prop.push_back(1.0);
      propf.push_back(0);
      pops.push_back(idpop);
      return sample(nindiv, pops, prop,propf);     
    }
  else
    return vector<Human>();
}

vector<Human> Metapop::sample_scattered(int nindiv) const
{ // scatterd sample - > sample nindiv indiv for each population 
  if(popsize>0)
    {
      vector<Human> A(nindiv*nbPop);
      int counter = 0;
      //      vector<int> samplelist;
      int pp,r;
      for(int pp(0);pp<nbPop;++pp)
	{
	  for(int i(0);i<nindiv;++i)
	    {
	      A[counter] =*(poptab[pp]->male[generation&1][i]);
	      ++counter;
	    }      
	}
      assert(counter==nindiv*nbPop);
      return A;
    }
  else
    return vector<Human>();
}

vector<Human> Metapop::sample_pooled(int nindiv,int npop) const
{ // pooled sample - > sample nindiv indiv for npop population 
  if(popsize>0)
    {
      return sample(nindiv,npop);
    }
  else
    return vector<Human>();
}

vector<Human> Metapop::sample(int nindiv, int npop) const
{
  if(popsize>0)
    {
      vector<double> prop(npop);
      for(int i(0); i<npop; ++i)
	prop[i] = 1.0 / (double) npop;
      return sample(nindiv, npop, prop);     
    }
  else
    return vector<Human>();
}

vector<Human> Metapop::sample(int nindiv, int npop,int popID) const
{
  if(popsize>0)
    {
      vector<double> prop;
      vector<int> props;
      prop.push_back(1);
      props.push_back(popID);
      return sample(nindiv, props, prop);
    }
  else
    return vector<Human>();
}

vector<Human> Metapop::sample(int nindiv, int npop,std::vector<double> prop) const // nb of sample + nb of pop + proportion in each pop
{
  if(popsize>0)
    {
      vector<int> sampledpop;
      sampledpop = sample_without_replacement(nbPop,npop);
      return sample(nindiv,sampledpop,prop);
    }
  else
    return vector<Human>();
}

vector<Human> Metapop::sample(int nindiv, std::vector<int> pops,std::vector<double> prop) const // nb of sample + IDs of pop + proportion in each pop
{
  if(popsize>0)
    {
      int npop = pops.size();
      vector<double> propf(npop);
      for(int i(0);i<npop;++i)
	{
	  propf[i] = (double) poptab[pops[i]]->nbFemale / (double) poptab[pops[i]]->popsize; // as many chance as women prportion -> makes it 0 if no women and 1 if no men
	}  
      return sample(nindiv, pops, prop, propf);
    }
  else
    return vector<Human>();
}

vector<Human> Metapop::sample(int nindiv, std::vector<int> pops,std::vector<double> prop,std::vector<double> propf) const// nb of sample + IDs of pop + proportion in each pop + prop of females
{
  if(popsize>0)
    {
      cout << nindiv << endl;
      vector<Human> A(nindiv);
      int counter = 0;
      //      vector<int> samplelist;
      int npop = pops.size();
      int r;
      int p,pp;
      if(npop != (int) prop.size())
	{
	  cerr << " Error: the sampling scheme is malfunctionning " << endl;
	}
      if(npop != (int) propf.size())
	{
	  cerr << " Error: the sampling scheme is malfunctionning " << endl;
	}  
      for(int i(0);i<nindiv;++i)
	{
	  p = myrand(npop);
	  pp = pops[p]; // sample only male !!
	  r = myrand(poptab[pp]->popsize- poptab[pp]->nbFemale);
	  A[counter] =*(poptab[pp]->male[generation&1][r]);
	  ++counter;
	}
      cout << counter << " " << nindiv << endl;
      assert(counter==nindiv);
      return A;
    }
  else
    return vector<Human>();
}


void Metapop::output_pop_size() const
{
  std::cout << " Total population size : " << popsize << " - female " << nbFemale << std::endl;
  for(int i(0);i<nbPop;++i)
    std::cout << " Population " << i << " has size : " << poptab[i]->popsize << " - female " << poptab[i]->nbFemale << std::endl;
}

void Metapop::write_smartFormat(std::string const& filename) const
{
  
}

ostream & operator<<(ostream & os, Metapop const & a)
{
  a.print(os);
  return os;
}

ostream & Metapop::print(ostream & os) const
{
  output_pop_size();
  for (int i(0); i<popsize;++i)
    {
      os << "f " << i << " " << *female[generation&1][i] ;
      //      os << "m " << i << " " << *male[generation&1][i] ;
    }
  return os;
}

void Metapop::move_indiv(int pop1, int id1, int pop2,bool sex)
{
  Human * tmp;
  int currentG = generation &1;
  if(sex)
    {  
      int indexlastpop2=poptab[pop2]->nbFemale-1;
      int indexlastpop1=poptab[pop1]->nbFemale-1;
      tmp = poptab[pop2]->female[currentG][indexlastpop2+1];
      poptab[pop2]->female[currentG][indexlastpop2+1] = poptab[pop1]->female[currentG][id1];
      poptab[pop1]->female[currentG][id1] = poptab[pop1]->female[currentG][indexlastpop1];
      poptab[pop1]->female[currentG][indexlastpop1] = tmp;
      poptab[pop2]->nbFemale+=1;
      poptab[pop1]->nbFemale-=1;
      poptab[pop2]->popsize+=1;
      poptab[pop1]->popsize-=1;
    }
  else
    {  
      int indexlastpop2=poptab[pop2]->popsize-poptab[pop2]->nbFemale-1;
      int indexlastpop1=poptab[pop1]->popsize-poptab[pop1]->nbFemale-1;
      tmp = poptab[pop2]->male[currentG][indexlastpop2+1];
      poptab[pop2]->male[currentG][indexlastpop2+1] = poptab[pop1]->male[currentG][id1];
      poptab[pop1]->male[currentG][id1] = poptab[pop1]->male[currentG][indexlastpop1];
      poptab[pop1]->male[currentG][indexlastpop1] = tmp;
      poptab[pop2]->popsize+=1;
      poptab[pop1]->popsize-=1;
    }
}


