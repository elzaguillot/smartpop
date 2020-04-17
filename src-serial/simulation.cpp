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

#include "simulation.h"

using namespace std;

Simulation::Simulation()
{
  counter = 0;
  maxpop = 0;
  popValues.push_back(defaultPopVal);
}

Simulation::Simulation(unsigned randomSeed)
{
  counter = 0;
  simuP.seed=randomSeed;
  maxpop = metapopValues.popsize;
  popValues.push_back(defaultPopVal);
  simuP.fileoutput = "smartpop_"+to_string(simuP.seed);
}

Simulation::Simulation(Simulation const & a)
{
  counter = 0;
  maxpop = a.maxpop;
  popValues = a.popValues;
  metapopValues = a.metapopValues;
  simuP = a.simuP;
}

Simulation::~Simulation()
{
}

Simulation& Simulation::operator = (Simulation const & a)
{
  if(this != &a)
    {
      counter = 0;
      maxpop = a.maxpop;
      popValues = a.popValues;
      metapopValues = a.metapopValues;
      simuP = a.simuP;
    }
  return *this;
}

void Simulation::load_file(string filename)
{
  try
    {
      simulation_settings options;
      options.load(filename);//"debug_settings.xml");
      options.save("SMARTPOP_parameters.xml");
      std::cout << "Success\n";
      metapopValues = options.metapop;
      simuP = options.simuP;
      popValues = options.populations;
      set_mu(options.mu.muMt1,options.mu.muMt2,options.mu.muX1, options.mu.muX2, options.mu.muY1, options.mu.muY2, options.mu.muA1, options.mu.muA2);
      DNAS = options.pdna;
    }
  catch (std::exception &e)
    {
      std::cout << "Error: " << e.what() << "\n";
    }

}

void Simulation::load_directory(string dirName)
{
  simuP.load = true;
  simuP.directoryLoad = dirName;
}

void Simulation::export_file(string filename)
{
  simulation_settings options;
  options.load(metapopValues,popValues,simuP,RATES,DNAS);
  options.save(filename);
}

void Simulation::check_parameters()
{
  metapopValues.popsize = 0;
  for(int n(0);n<metapopValues.nbPop;++n)
    metapopValues.popsize += popValues[n].popsize;
  if(simuP.diversity.size()>4)
    {
      cerr << " ERROR : something is wrong in your choice of diversity to output " << endl;
      exit(2);
    }
  if(DNAS.nbChromosomesA == 0)
    DNAS.sizePerA = 0;
  if(DNAS.nbChromosomesX == 0)
    DNAS.sizePerX = 0;
  if(DNAS.sizePerX == 0)
    DNAS.nbChromosomesX = 0;
  if(DNAS.sizePerA == 0)
    DNAS.nbChromosomesA = 0;
  DNAS.nbCrumbNuc = DNAS.nbCrumbA + DNAS.nbCrumbX + std::max(DNAS.nbCrumbX,DNAS.nbCrumbY);
  if(simuP.diversity[0]==0)
    {
      if((DNAS.sizePerA == 0)||(DNAS.sizePerX == 0)||(DNAS.sizeY == 0)||(DNAS.sizeMt == 0))
	{
	  simuP.diversity.clear();
	  if(DNAS.sizeMt>0)
	    simuP.diversity.push_back(1);
	  if(DNAS.sizePerX>0)
	    simuP.diversity.push_back(2);
	  if(DNAS.sizeY>0)
	    simuP.diversity.push_back(3);
	  if(DNAS.sizePerA>0)
	    simuP.diversity.push_back(4);
	}
    }
  int index=-1;
  if(DNAS.sizeMt==0)
    for(unsigned int i(0);i<simuP.diversity.size();++i)
      if(simuP.diversity[i]==1)
	index = i;
  if(index>-1)
    simuP.diversity.erase(simuP.diversity.begin()+index);
  index=-1;
  if(DNAS.sizePerX==0)
    for(unsigned int i(0);i<simuP.diversity.size();++i)
      if(simuP.diversity[i]==2)
	index = i;
  if(index>-1)
    simuP.diversity.erase(simuP.diversity.begin()+index);
  index=-1;
  if(DNAS.sizeY==0)
    for(unsigned int i(0);i<simuP.diversity.size();++i)
      if(simuP.diversity[i]==3)
	index = i;
  if(index>-1)
    simuP.diversity.erase(simuP.diversity.begin()+index);
  index=-1;
  if(DNAS.sizePerA==0)
    for(unsigned int i(0);i<simuP.diversity.size();++i)
      if(simuP.diversity[i]==4)
	index = i;
  if(index>-1)
    simuP.diversity.erase(simuP.diversity.begin()+index);
  for(int n(0);n<metapopValues.nbPop;++n)
    {
      if(popValues[n].matingSystem.size()>1)
	{
	  std::vector<int> newm;
	  for(unsigned i(0);i<popValues[n].matingSystem.size();++i)
	    if(popValues[n].matingSystem[i]>1)
	      newm.push_back(popValues[n].matingSystem[i]);
	  popValues[n].matingSystem = newm;
	  if(popValues[n].matingSystem.size()==0)
	    popValues[n].matingSystem.push_back(1);
	}
    }
  if((metapopValues.nbPop>1)&&(metapopValues.migrationType==1))
    {
      if((int) metapopValues.migrationRatesF.size()<metapopValues.nbPop)
	set_migrationRate(metapopValues.nbPop-1,0,1,0);
      else
	if((int) metapopValues.migrationRatesF.size()>metapopValues.nbPop)
	  {
	    std::cerr << "ERROR: migration rates set for non existing population " << std::endl;
	    exit(2);
	  }
      if((int) metapopValues.migrationRatesM.size()<metapopValues.nbPop)
	set_migrationRate(metapopValues.nbPop-1,0,0,0);
      else
	if((int) metapopValues.migrationRatesM.size()>metapopValues.nbPop)
	  {
	    std::cerr << "ERROR: migration rates set for non existing population " << std::endl;
	    exit(2);
	  }
      if((int) metapopValues.migrationRatesF.size()<metapopValues.nbPop)
	{
	  std::cerr << "ERROR: migration rates poorly set - contact the developper " << std::endl;
	  exit(2);
	}
      if((int) metapopValues.migrationRatesM.size()<metapopValues.nbPop)
	{
	  std::cerr << "ERROR: migration rates poorly set - contact the developer" << std::endl;
	  exit(2);
	}
      for(int n(0);n<metapopValues.nbPop;++n)
	{
	  if((int) metapopValues.migrationRatesF[n].size()!=metapopValues.nbPop)
	    {
	      std::cerr << "ERROR: migration rates poorly set - contact the developper " << std::endl;
	      exit(2);
	    }
	  if((int) metapopValues.migrationRatesM[n].size()!=metapopValues.nbPop)
	    {
	      std::cerr << "ERROR: migration rates poorly set - contact the developer" << std::endl;
	      exit(2);
	    }
	}
    }
  if(metapopValues.migrationType>=3)
    if(metapopValues.nbSources==0)
      metapopValues.nbSources =3;
  if(metapopValues.nbPop<1)
    cerr << "ERROR : you must simulate at least 1 population" << endl;
}

void Simulation::set_seed(unsigned int a)
{
  simuP.seed = a;
}

void Simulation::set_popsize(int a)
{
  metapopValues.popsize = 0;
  for(int n(0);n<metapopValues.nbPop;++n)
    {
      popValues[n].popsize = a;
      metapopValues.popsize += a;
    }
}

void Simulation::set_popsize(int a,int popID)
{
  if((popID<0) || (popID>= metapopValues.nbPop))
    {
      cerr << "Error: population " << popID << " does not exist " << endl;
      exit(2);
    }
  popValues[popID].popsize = a;
}

void Simulation::set_nbSources(int a)
{
  metapopValues.nbSources = a;
}

void Simulation::set_nbPop(int a)
{
  if(a>metapopValues.nbPop)
    for(int i(0);i<a-metapopValues.nbPop;++i)
      {
	popVal tmp = popValues[0];
	popValues.push_back(tmp);
      }
  else
    if(a<metapopValues.nbPop)
      for(int i(0);i<metapopValues.nbPop - a;++i)
	{
	  popValues.pop_back();
	}
  metapopValues.nbPop = a;
}

void Simulation::set_sample(int a)
{
  metapopValues.sample = a;
}

void Simulation::set_nbsamplepop(int a)
{
  metapopValues.nbsamplepop = a;
}

void Simulation::set_sampleType(int a)
{
  metapopValues.sampleType = a;
}

void Simulation::set_sampleBetween_on()
{
  metapopValues.sampleBetween = true;
}

void Simulation::set_step(int a)
{
  simuP.step = a;
}

void Simulation::set_nstep(int a)
{
  simuP.nstep = a;
}

void Simulation::set_matsys(int a,int popID)
{
  if(popID<metapopValues.nbPop)
    {
      popValues[popID].matingSystem.push_back(a);
    }
  else
    {
      cerr << " ERROR : attempt to set a mating system to a non existing population " << endl;
      exit(2);
    }
}

void Simulation::set_matsys(int a)
{
  for(int n(0); n<metapopValues.nbPop;++n)
    {
      popValues[n].matingSystem.push_back(a);
    }
}

void Simulation::set_inbreeding(int a)
{
  if(a>2)
    {
      cerr << "Error: inbreeding must take a value 0, 1 or 2 " << endl;
      exit(1);
    }
  for(int n(0); n<metapopValues.nbPop;++n)
    {
      popValues[n].inbreeding = a;
    }
}

void Simulation::set_inbreeding(int a,int popID)
{
  if(popID<metapopValues.nbPop)
    {
      popValues[popID].inbreeding=a;
    }
  else
    {
      cerr << " ERROR : attempt to set a mating system inbreeding factor to a non existing population " << endl;
      exit(2);
    }
}

void Simulation::set_varChild(double a)
{
  for(int n(0); n<metapopValues.nbPop;++n)
    {
      popValues[n].varChild = a;
    }
}

void Simulation::set_varChild(double a,int popID)
{
  if(popID<metapopValues.nbPop)
    {
      popValues[popID].varChild=a;
    }
  else
    {
      cerr << " ERROR : attempt to set a mating system inbreeding factor to a non existing population " << endl;
      exit(2);
    }
}

void Simulation::set_polygamy(int a)
{
  if(a<1)
    {
      cerr << "Error : the polygamy factor must be at least 1 (= monogamy)" << endl;
      exit(2);
    }
  for(int n(0); n<metapopValues.nbPop;++n)
    {
      popValues[n].polygamy = a;
    }
}

void Simulation::set_polygamyall(int a)
{
  if(a<1)
    {
      cerr << "Error : the polygamy factor must be at least 1 (= monogamy)" << endl;
      exit(2);
    }
  for(int n(0); n<metapopValues.nbPop;++n)
    {
      popValues[n].polygamy = a;
      popValues[n].polygamyf = a;
    }
}

void Simulation::set_polygamyf(int a)
{
  if(a<1)
    {
      cerr << "Error : the polygamy factor must be at least 1 (= monogamy)" << endl;
      exit(2);
    }
  for(int n(0); n<metapopValues.nbPop;++n)
    {
      popValues[n].polygamyf = a;
    }
}

void Simulation::set_migrationType(int a)
{
  metapopValues.migrationType = a;
  if((a==3)||(a==4))
    if(metapopValues.nbSources==0) // put default valye
      metapopValues.nbSources=3;
}

void Simulation::set_pmate(double a)
{
  if((a<0)||(a>1))
    {
      cerr << "Error: pmate must be between 0 and 1" << endl;
      exit(2);
    }
  metapopValues.pmate = a;
}

void Simulation::set_pmig(double a)
{
  if((a<0)||(a>1))
    {
      cerr << "Error: pmig must be between 0 and 1" << endl;
      exit(2);
    }
  metapopValues.pmig = a;
}

void Simulation::set_polygamy(int a,int popID)
{
  if(popID<metapopValues.nbPop)
    {
      popValues[popID].polygamy=a;
    }
  else
    {
      cerr << " ERROR : attempt to set a mating system inbreeding factor to a non existing population " << endl;
      exit(2);
    }
}

void Simulation::set_polygamyf(int a,int popID)
{
  if(popID<metapopValues.nbPop)
    {
      popValues[popID].polygamyf=a;
    }
  else
    {
      cerr << " ERROR : attempt to set a mating system inbreeding factor to a non existing population " << endl;
      exit(2);
    }
}

void Simulation::set_nSimu(int a)
{
  simuP.nsimu = a;
}

void Simulation::set_sizeMt(int a)
{
  if(a%CRUMB>0)
    {
      DNAS.nbCrumbMt = a/CRUMB + 1;
      DNAS.sizeMt = DNAS.nbCrumbMt * CRUMB;
    }
  else
    {
      DNAS.sizeMt = a;
      DNAS.nbCrumbMt = a/CRUMB;
    }
}

void Simulation::set_sizePerX(int a)
{
  if(a%CRUMB>0) // if the size is not proportionnal to the crumb, then increase the size to the higher crumb
    {
      DNAS.nbCellPerX = a/CRUMB + 1;
      DNAS.sizePerX = DNAS.nbCellPerX * CRUMB;
      DNAS.nbCrumbX = DNAS.nbCellPerX*DNAS.nbChromosomesX;
    }
  else
    {
      DNAS.sizePerX = a;
      DNAS.nbCellPerX = a/CRUMB;
      DNAS.nbCrumbX = DNAS.nbCellPerX*DNAS.nbChromosomesX;
    }
}

void Simulation::set_sizeY(int a)
{
  if(a%CRUMB>0)
    {
      DNAS.nbCrumbY = a/CRUMB + 1;
      DNAS.sizeY = DNAS.nbCrumbY * CRUMB;
    }
  else
    {
      DNAS.sizeY = a;
      DNAS.nbCrumbY = a/CRUMB;
    }
}

void Simulation::set_sizePerA(int a)
{
  if(a%CRUMB>0) // if the size is not proportionnal to the crumb, then increase the size to the higher crumb
    {
      DNAS.nbCellPerA = a/CRUMB + 1;
      DNAS.sizePerA = DNAS.nbCellPerA * CRUMB;
      DNAS.nbCrumbA = DNAS.nbCellPerA*DNAS.nbChromosomesA;
    }
  else
    {
      DNAS.sizePerA = a;
      DNAS.nbCellPerA = a/CRUMB;
      DNAS.nbCrumbA = DNAS.nbCellPerA*DNAS.nbChromosomesA;
    }
}

void Simulation::set_nbLociA(int a)
{
  DNAS.nbChromosomesA = a;
  DNAS.nbCrumbA = DNAS.nbCellPerA*DNAS.nbChromosomesA*2;
}

void Simulation::set_nbLociX(int a)
{
  DNAS.nbChromosomesX = a;
  DNAS.nbCrumbX = DNAS.nbCellPerX*DNAS.nbChromosomesX;
}

void Simulation::set_mu(double a,double b, double c, double d)
{
  RATES.muMt1 = a;
  RATES.muMt2 = a;
  RATES.muX1 = b;
  RATES.muX2 = b;
  RATES.muY1 = c;
  RATES.muY2 = c;
  RATES.muA1 = d;
  RATES.muA2 = d;
}

void Simulation::set_mu(double a,double b, double c, double d,double e,double f, double g, double h)
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

void Simulation::set_growthFactor(double a,double b,double c)
{
  for(int n(0);n<metapopValues.nbPop;++n)
    {
      popValues[n].growthFactor[0] = a;
      popValues[n].growthFactor[1] = b;
      popValues[n].growthFactor[2] = c;
    }
}

void Simulation::set_growthFactor(double a,double b,double c,int popID)
{
  if(popID<metapopValues.nbPop)
    {
      popValues[popID].growthFactor[0] = a;
      popValues[popID].growthFactor[1] = b;
      popValues[popID].growthFactor[2] = c;
    }
  else
    {
      cerr << "ERROR: trying to set a population demography for non existing population" << endl;
      exit(2);
    }
}

void Simulation::set_burninT(double a)
{
  simuP.burninT = a;
}

void Simulation::set_burninDna(int a)
{
  if((a>0)&&(a<5))
    simuP.burninDna = a;
  else
    {
      cerr << "ERROR: The dna type must be 1 (mtDNA) 2(X) 3(Y) or 4 (autosomes)" << endl;
      exit(1);
    }
}

void Simulation::set_diversity(int a)
{
  if(simuP.diversity[0] <= 0)
    simuP.diversity[0] = a;
  else
    simuP.diversity.push_back(a);
  if(a==-1)
    {
      simuP.diversity.resize(1);
      simuP.diversity[0] = a;
    }
}

void Simulation::set_fileoutput(string a)
{
  simuP.fileoutput = a;
}

void Simulation::set_verbose_on()
{
  simuP.verbose = true;
}

void Simulation::set_fasta_on()
{
  simuP.fasta = true;
}

void Simulation::set_save_on(string dirName)
{
  simuP.smartsave = true;
  simuP.directorySave = dirName;
}



void Simulation::set_arlequin_on()
{
  simuP.arl = true;
}

void Simulation::set_haplotype_on()
{
  simuP.hap = true;
}

void Simulation::set_sfs_on()
{
  simuP.sfs = true;
}

void Simulation::set_APA_mig()
{
  metapopValues.migrationType=4;
}

void Simulation::set_header_on()
{
  simuP.header = true;
}

void Simulation::set_migrationRate(int pop1,int pop2, bool sex, double rate)
{
  if((rate>=0)&&(rate<=1))
    {
      if(rate<1)
	metapopValues.migrationType = 1;
      else
	metapopValues.migrationType = 2;
      int maxpop=std::max(pop1,pop2) +1 ;
      int dim1 = metapopValues.migrationRatesF.size(); // must always squared !
      if(dim1<maxpop)
	{
	  for(int i(0);i<(maxpop-dim1);++i)
	    metapopValues.migrationRatesF.push_back(std::vector<double>(maxpop,0));
	  for(int i(0);i<(dim1);++i)
	    for(int j(0);j<maxpop-dim1;++j)
	      metapopValues.migrationRatesF[i].push_back(0);
	}
      assert( metapopValues.migrationRatesF.size() == metapopValues.migrationRatesF[0].size());
      assert( (int) metapopValues.migrationRatesF.size() >= maxpop);
      dim1 =   metapopValues.migrationRatesM.size(); // must always be  squared !
      if(dim1<maxpop)
	{
	  for(int i(0);i<(maxpop-dim1);++i)
	    metapopValues.migrationRatesM.push_back(std::vector<double>(maxpop,0));
	  for(int i(0);i<(dim1);++i)
	    for(int j(0);j<maxpop-dim1;++j)
	      metapopValues.migrationRatesM[i].push_back(0);
	}
      assert( metapopValues.migrationRatesM.size() == metapopValues.migrationRatesM[0].size());
      assert((int) metapopValues.migrationRatesM.size() >= maxpop);
      if(sex==true)
	{
	  metapopValues.migrationRatesF[pop1][pop2] = rate;
	}
      else
	{
	  metapopValues.migrationRatesM[pop1][pop2] = rate;
	}
    }
  else
    {
      cerr << "ERROR : Migration rate must be between 0 and 1 "<< endl;
      exit(2);
    }
}

int Simulation::get_nSimu()
{
  return simuP.nsimu;
}

bool Simulation::get_verbose()
{
  return simuP.verbose;
}

unsigned int Simulation::get_seed()
{
  return simuP.seed;
}

int Simulation::get_popsize()
{
  return metapopValues.popsize;
}

int Simulation::get_nbPop()
{
  return metapopValues.nbPop;
}

int Simulation::get_sample()
{
  return metapopValues.sample;
}

int Simulation::get_step()
{
  return simuP.step;
}

int Simulation::get_nstep()
{
  return simuP.nstep;
}

std::vector<int> Simulation::get_matsys(int popID)
{
  if(popID<metapopValues.nbPop)
    return popValues[popID].matingSystem;
  else
    return std::vector<int>();
}

int Simulation::get_inbreeding(int popID)
{
  if(popID<metapopValues.nbPop)
    return popValues[popID].inbreeding;
  else
    return -1;
}

double Simulation::get_varChild(int popID)
{
  if(popID<metapopValues.nbPop)
    return popValues[popID].varChild;
  else
    return -1;
}

int Simulation::get_polygamy(int popID)
{
  if(popID<metapopValues.nbPop)
    return popValues[popID].polygamy;
  else
    return -1;
}

int Simulation::get_polygamyf(int popID)
{
  if(popID<metapopValues.nbPop)
    return popValues[popID].polygamyf;
  else
    return -1;
}

double* Simulation::get_growthFactor(int popID)
{
  if(popID<metapopValues.nbPop)
    return popValues[popID].growthFactor;
  else
    return NULL;
}

double Simulation::get_burninT()
{
  return simuP.burninT;
}

int Simulation::get_burninDna()
{
  return simuP.burninDna;
}


std::vector<int> Simulation::get_diversity()
{
  return simuP.diversity;
}

std::string Simulation::get_fileoutput()
{
  return simuP.fileoutput;
}

bool Simulation::get_fasta()
{
  return simuP.fasta;
}

bool Simulation::get_save()
{
  return simuP.smartsave;
}

bool Simulation::get_arlequin()
{
  return simuP.arl;
}

bool Simulation::get_header()
{
  return simuP.header;
}


void Simulation::header_write()
{
  std::string out = "popsize samplesize time pmate pmig ";
  if(DNAS.sizeMt>0)
    out += "mt.size mt.S mt.thetapi mt.nbhap mt.heterozygosityHap mt.heterozygosityNei mt.thetawatterson mt.tajimaD mt.eta  mt.meanSFS mt.varSFS ";
  if(DNAS.sizePerX>0)
    out += "x.size x.S x.thetapi x.nbhap x.heterozygosityHap x.heterozygosityNei x.thetawatterson x.tajimaD x.eta  x.meanSFS x.varSFS ";
  if(DNAS.sizeY>0)
    out += "y.size y.S y.thetapi y.nbhap y.heterozygosityHap y.heterozygosityNei y.thetawatterson y.tajimaD y.eta  y.meanSFS y.varSFS ";
  if(DNAS.sizePerA>0)
    out += "a.size a.S a.thetapi a.nbhap a.heterozygosityHap a.heterozygosityNei a.thetawatterson a.tajimaD a.eta  a.meanSFS a.varSFS ";
  out += "\n";
  ofstream fs;
  fs.open((simuP.fileoutput+"all.txt").c_str(),ios::app);
  fs << out ;
  fs.close();
  if(metapopValues.sampleBetween)
    {
      std::string out2 = "popsize samplesize time pmate pmig pop1 pop2 ";
      if(DNAS.sizeMt>0)    
	out2 += "mt.mean_pw_diff mt.Fst mt.neiD mt.chordD mt.rwcD ";
      if(DNAS.sizePerX>0)
	out2 += "x.mean_pw_diff x.Fst x.neiD x.chordD x.rwcD ";
      if(DNAS.sizeY>0)
	out2 += "y.mean_pw_diff y.Fst y.neiD y.chordD y.rwcD ";
      if(DNAS.sizePerA>0)
	out2 += "a.mean_pw_diff a.Fst a.neiD a.chordD a.rwcD ";
      out2 += "\n";
      ofstream fs2;
      fs2.open((simuP.fileoutput+"_Between.txt").c_str(),ios::app);
      fs2 << out2;
      fs2.close();
    }
}

int Simulation::classic()
{
  std::vector<int> popsizes;
  metapopValues.popsize=0;
  for(int n(0); n<metapopValues.nbPop; ++n)
    {
      popsizes.push_back(popValues[n].popsize);
      metapopValues.popsize+=popValues[n].popsize;
    }
  Metapop epsilon(metapopValues.nbPop,metapopValues.popsize,popsizes);
  epsilon.migrationRateF = metapopValues.migrationRatesF;
  epsilon.migrationRateM = metapopValues.migrationRatesM;
  epsilon.migrationType = metapopValues.migrationType;
  epsilon.set_nbSources(metapopValues.nbSources);
  epsilon.pmig =  metapopValues.pmig;
  epsilon.pmate =  metapopValues.pmate;
  for( int n(0);n<metapopValues.nbPop;++n)
    epsilon.poptab[n]->set_growthFactor(popValues[n].growthFactor);
  epsilon.Init(simuP.nstep*simuP.step);
  for(int i(0);i<simuP.nsimu;i++)
    {
      if(simuP.verbose)
	cout << " Simulation " << i << "/" << simuP.nsimu << endl;
      if(i>0)
	epsilon.reset(metapopValues.popsize,popsizes);
      evolve(epsilon);
      ++counter;
    }
  return 0;
}

int Simulation::successive_classic()
{

  std::vector<int> popsizes;
  for(int n(0); n<metapopValues.nbPop; ++n)
    popsizes.push_back(popValues[n].popsize);
  string filename;
  for(int i(0);i<simuP.nsimu;i++)
    {
      filename = simuP.directoryLoad + "_" + to_string(i) + ".sim";
      if(simuP.verbose)
	cout << " Simulation " << i << "/" << simuP.nsimu << " from file " << filename << endl;
      Metapop epsilon(metapopValues.nbPop,metapopValues.popsize,popsizes);
      cout << filename << endl;
      evolve(epsilon);
      ++counter;
    }
  return 0;
}

int Simulation::print_help()
{
  string outstr="";
  outstr += "SMARTPOP 2.0  by  E.G.Guillot\n This program comes with ABSOLUTELY NO WARRANTY; for details see smartpop.sourceforge.net\n This is free software, and you are welcome to redistribute it  under certain conditions;\n For ANY usage please cite Guillot and Cox, BMC Bioinformatics 2014;\n for details see smartpop.sourceforge.net\n \n";
  outstr += "usage: ./smartpop -p popsize\n";
  outstr += "\t other options (default values)\n";
  outstr += "\t\t -h --help           input file with parameters\n";
  outstr += "\t\t -i --input          input file with parameters\n";
  outstr += "\t\t -v --verbose        verbose\n";
  outstr += "\t\t -o --output XX      name for output files\n";
  outstr += "\t\t --header            write the header in the diversity files\n";
  outstr += "\t\t --fasta             output a fasta file of the sample for each simulation\n";
  outstr += "\t\t --arl               output an arlequin file of the sample for each simulation\n";
  outstr += "\t\t --child2            forces all couple to have 2 children (i.e. no variuance in reproductive success\n"; 
  outstr += "\t\t --pooled            pooled sample (X indiv from Y demes)\n";
  outstr += "\t\t --scattered         scattered sample (X indiv from all demes)\n";
  outstr += "\t\t --local             local sample (X indiv from 1 demes)\n";
  outstr += "\t\t --blind             random sample without knowing the structure \n";
  outstr += "\t\t --sample X          number of individuals to sample per deme (according to the scheme) \n";
  outstr += "\t\t --nbsamplepop X     number of population to sample for pooled scheme and for inter-population diversity comparison \n";
  outstr += "\t\t --seed X            random seed\n";
  outstr += "\t\t --step X            nb of generation X to evolve in one step\n";
  outstr += "\t\t --nstep X           nb of steps X to evolve: 1 simulation = 'nbStep' times checking diversity every 't' generations\n";
  outstr += "\t\t --nsimu X           number of simulations X\n";
  outstr += "\t\t --sample X          sample size X\n";
  outstr += "\t\t --nbLociX X         number of loci on X chromosome\n";
  outstr += "\t\t --sizeX X           size (number of sites) per locus on X chromosome\n";
  outstr += "\t\t --nbLociA X         number of loci on autosome\n";
  outstr += "\t\t --sizeA X           size (number of sites) per locus on autosome\n";
  outstr += "\t\t --sizeMt X          size mtDNA sequence \n";
  outstr += "\t\t --sizeY X           number of sites X per locus on autosome\n";
  outstr += "\t\t -p --popsize X      set the population size to X for deme\n";
  outstr += "\t\t --popsizeN X Y      set the population size to X for deme Y\n";
 
  outstr += "\t\t --mat X             mating system (X=:1 random, 2 MBD/FZS, 3 MZD/MZS, 4 FZD/MBS, FBD/FBS)\n";
  outstr += "\t\t --polygamy X        set the maximum number of mates X per indiv (if X=1 it becomes monagamy) \n";
  outstr += "\t\t --polyandry X       set the maximum number of husbands X per woman (if X=1 it becomes monagamy) \n";
  outstr += "\t\t --polygyny X        set the maximum number of wives X per man (if X=1 it becomes monagamy) \n";
  outstr += "\t\t --polygamyN X Y     set the maximum number of mates X per indiv (if 1 it becomes monagamy) for population Y \n";
  outstr += "\t\t --polygynyN X Y     set the maximum number of wives X per man (if 1 it becomes monagamy) for population Y \n";
  outstr += "\t\t --inbreeding X      X=0 nothing; X=1 mating between half siblings is forbidden; X=2 mathing between full sibling is forbidden\n";
  outstr += "\t\t --mu X X X X (X X X X) if 4 arguments: rate_mitochondrial_transition rate_X_chr_transition rate_Y_chr_transition rate_autosome_transition \n\t\t\t WARNING : mutation rate is per site per generation \n";
  outstr += "\t\t                        if 8 arguments: rate_mitochondrial_transition rate_mitochondrial_transversion rate_X_chr_transition rate_X_chr_transversion rate_Y_chr_transition rate_Y_chr_transversion rate_autosome_transition rate_autosome_transversion \n\t\t\t WARNING : mutation rate per site per generation \n";
  outstr += "\t\t --burnin X  number of generations in the burnin phasewith normal mutation rate\n";
  outstr += "\t\t --burninHigh X  number of generations in the burnin phase with high mutation rate\n";
  outstr += "\t\t --demog a b c       with a,b,c decimal numbers which describes the change in population size through time like this (applied to all demes) :\n\t\t\t                    popsize(t+1) = a + b*popsize(t) + c*popsize(t)*popsize(t)\n";
  outstr += "\t\t --demogn a b c X    as previously but applied to deme X ";

  outstr += "See smartpop.sourceforge.net for more details\n";
  cout << outstr;
  exit(0);
  return 0;
}

int Simulation::verbose_start()
{
  cout << "SMARTPOP 2.0 by  E.G.Guillot\n This program comes with ABSOLUTELY NO WARRANTY; for details see smartpop.sourceforge.net\n This is free software, and you are welcome to redistribute it  under certain conditions;\nFor ANY usage please cite Guillot and Cox, BMC Bioinformatics 2014;\n for details see smartpop.sourceforge.net\n \n";
  cout << " Random seed : " << simuP.seed << endl;
  cout << "You have chosen the following options : " << endl;
  cout << "\t output file : " << simuP.fileoutput;
  cout << endl;
  cout << "\t t : " << simuP.step;
  cout << "\t -nbStep : " << simuP.nstep;
  cout << "\t -burnin : " << simuP.burninT << endl;
  cout << "\t -burninH : " << simuP.burninDna << endl;
  if(simuP.burninT == 0)
    cout << "\t WARNING : you have not input any burnin, you simulation will start with a diversity 0 withon your population, be aware it has poor biological meaning" << endl;
  cout << "\t -nsimu : " << simuP.nsimu << endl;
  //  cout << "\t -mat : " << matsys;
  for(int i(0);i<metapopValues.nbPop;++i)
    {
      if(popValues[i].inbreeding == 0)
	cout << " Siblings alliances are fully allowed " << endl;
      else
	cout << " Siblings alliances are forbidden (between those who share at least one parent) " << endl;
    }
  cout << "\t Complete DNA simulations" << endl;
  cout << "\t Mitochondiral sequences have " << DNAS.sizeMt << " base pair of nucleotides" << endl;
  cout << "\t Y chromosome sequences have " << DNAS.sizeY << " base pair of nucleotides" << endl;
  cout << "\t " << DNAS.nbChromosomesX << " X loci sequences have " << DNAS.sizePerX << " base pair of nucleotides" << endl;
  cout << "\t " << DNAS.nbChromosomesA << " Autosome loci sequences have " << DNAS.sizePerA << " base pair of nucleotides" << endl;
  cout << "\t with Autosome transition rate " << RATES.muA1 << "  - transversion rate " << RATES.muA2 <<endl;
  cout << "\t with Autosome transition rate " << RATES.muX1 << "  - transversion rate " << RATES.muX2 <<endl;
  cout << "\t with Autosome transition rate " << RATES.muY1 << "  - transversion rate " << RATES.muY2 <<endl;
  cout << "\t with mtDNA transition rate " << RATES.muMt1 << "  - transversion rate " << RATES.muMt2 <<endl;
  for(unsigned k(0); k < simuP.diversity.size(); ++k)
    {
      switch(simuP.diversity[k]){
      case(-1):
	cout << "outputting no DNA diversity";
	break;
      case(0):
	cout << "outputting complete DNA diversity";
	break;
      case(1):
	cout << "outputting mtdna diversity only\n WARNING : if you only whish to compute mtdna diversity use the flag -mtdna for faster simulations";
	break;
      case(2):
	cout << "outputting Y chromosome diversity ";
	break;
      case(3):
	cout << "outputting X chromosome diversity ";
	break;
      case(4):
	cout << "outputting Autosome diversity ";
	break;
      default:
	cerr << "ERROR in the diversity profile"<<endl;
	exit(1);
      }
      if(simuP.fasta)
	cout << " and FASTA file";
      if(simuP.arl)
	cout << " and ARLEQUIN file";
      if(simuP.hap)
	cout << " and haplotypes";
      if(simuP.sfs)
	cout << " and Site Frequency Spectrum";
      if(simuP.smartsave)
	cout << " and SMARTPOP population file";
      cout << endl;
    }
  cout << " Starting now " << simuP.nsimu<< " simulations "<<endl;
  return 0;
}

int Simulation::verbose_end()
{
  return 0;
}

int Simulation::output_main( Metapop & epsilon) // X is diversity to output
{
  std::vector<double> out;
  if(epsilon.popsize>0)
    {
      ofstream fs;
      int p = epsilon.popsize;
      int timeG = epsilon.generation;
      int nF = 0;
      double pmate = epsilon.pmate;
      double pmig = epsilon.pmig;
      int samplesize = metapopValues.sample;
      int sampleType = metapopValues.sampleType;
      if(samplesize>0)
	{
	  vector<Human> epsilonSub;
	  std::string outputTowrite;
	  nF = 0;
	  if(sampleType == 0) // local
	    epsilonSub = epsilon.sample_local_random(samplesize);
	  else{ if(sampleType == 1) // pooled
	      epsilonSub = epsilon.sample_pooled(samplesize,metapopValues.nbsamplepop);
	    else{ 
	      if(sampleType == 2) // scattered
		{
		  epsilonSub = epsilon.sample_scattered(samplesize);  
		  samplesize *= metapopValues.nbPop;
		}
	      else //(sampleType == 3) // scatterd
		epsilonSub = epsilon.sample(samplesize,metapopValues.nbPop);
	    }
	  }
	  // mtdna on males only
	  out.push_back(p);
	  out.push_back(samplesize);
	  out.push_back(timeG);
	  out.push_back(pmate);
	  out.push_back(pmig);  
	  if(DNAS.sizeMt>0)
	    {
	      DNA sampleall(epsilonSub,samplesize,nF,3);
	      out += sampleall.output_diversity_with_time(simuP.fileoutput+"_all_div.txt",p,timeG);
	      if(simuP.fasta)
		sampleall.output_fasta(simuP.fileoutput+"_mt_"+to_string(counter));
	      if(simuP.arl)
		sampleall.output_arlequin(simuP.fileoutput+"_mt_"+to_string(counter));
	      if(simuP.sfs)
		sampleall.SFS_folded(simuP.fileoutput+"_mt_"+to_string(timeG)+".sfs");
	    }
	  if(DNAS.sizePerX>0)
	    {
	      DNA sampleall(epsilonSub,samplesize,nF,0);
	      out += sampleall.output_diversity_with_time(simuP.fileoutput+"_x_div.txt",p,timeG);
	      if(simuP.fasta)
		sampleall.output_fasta(simuP.fileoutput+"_x_"+to_string(counter));
	      if(simuP.arl)
		sampleall.output_arlequin(simuP.fileoutput+"_x_"+to_string(counter));
	      if(simuP.sfs)
		sampleall.SFS_folded(simuP.fileoutput+"_x_"+to_string(timeG)+".sfs");
	    }
	  if(DNAS.sizeY>0)
	    {
	      DNA sampleall(epsilonSub,samplesize,nF,1);
	      out += sampleall.output_diversity_with_time(simuP.fileoutput+"_y_div.txt",p,timeG);
	      if(simuP.fasta)
		sampleall.output_fasta(simuP.fileoutput+"_y_"+to_string(counter));
	      if(simuP.arl)
		sampleall.output_arlequin(simuP.fileoutput+"_y_"+to_string(counter));
	      if(simuP.sfs)
		sampleall.SFS_folded(simuP.fileoutput+"_y_"+to_string(timeG)+".sfs");
	    }
	  if(DNAS.sizePerA>0)
	    {
	      DNA sampleall(epsilonSub,samplesize,nF,2);
	      out += sampleall.output_diversity_with_time(simuP.fileoutput+"_a_div.txt",p,timeG);
	      if(simuP.fasta)
		sampleall.output_fasta(simuP.fileoutput+"_a_"+to_string(counter));
	      if(simuP.arl)
		sampleall.output_arlequin(simuP.fileoutput+"_a_"+to_string(counter));
	      if(simuP.sfs)
		sampleall.SFS_folded(simuP.fileoutput+"_a_"+to_string(timeG)+".sfs");
	    }
	  fs.open((simuP.fileoutput+"all.txt").c_str(),ios::app);
	  fs << outputTowrite ;
	  fs << out ;
	  fs.close();
	  if((metapopValues.nbPop>1) && (metapopValues.sampleBetween))
	    {
	      vector<DNA> dnaMT;
	      vector<DNA> dnaX;
	      vector<DNA> dnaY;
	      vector<DNA> dnaA;
	      for(int popid(0);popid<metapopValues.nbsamplepop;++popid) //only sample half, too slow otherwise
		{
		  epsilonSub = epsilon.sample(metapopValues.sample,metapopValues.nbPop,popid);
		  if(DNAS.sizeMt>0) // mtdna or complete
		    {
		      DNA tmp(epsilonSub, metapopValues.sample, nF, 3);
		      dnaMT.push_back(tmp);
		    }
		  if(DNAS.sizePerX>0)
		    {
		      dnaX.push_back(DNA(epsilonSub, metapopValues.sample, nF, 0));
		    }
		  if(DNAS.sizeY>0)
		    {
		      dnaY.push_back(DNA(epsilonSub, metapopValues.sample, nF, 1));
		    }
		  if(DNAS.sizePerA>0)
		    {
		      dnaA.push_back(DNA(epsilonSub, metapopValues.sample, nF, 2));
		    }
		}
	      for(int popid(0);popid<metapopValues.nbsamplepop;++popid)
		{
		  for(int popid2(popid+1);popid2<metapopValues.nbsamplepop;++popid2)
		    {
		      std::vector<double> out2;
		      out2.push_back(p);
		      out2.push_back(samplesize);
		      out2.push_back(timeG);
		      out2.push_back(pmate);
		      out2.push_back(pmig);
		      out2.push_back(popid);
		      out2.push_back(popid2);
		      if(DNAS.sizeMt>0)
			{
			  out2 += dnaMT[popid].output_diversity(dnaMT[popid2],simuP.fileoutput+"_mt_between_div.txt",p,timeG);
			}
		      if(DNAS.sizePerX>0)
			{
			  out2 += dnaX[popid].output_diversity(dnaX[popid2],simuP.fileoutput+"_x_between_div.txt",p,timeG);
			}
		      if(DNAS.sizeY>0)
			{
			  out2 += dnaY[popid].output_diversity(dnaY[popid2],simuP.fileoutput+"_y_between_div.txt",p,timeG);
			}
		      if(DNAS.sizePerA>0)
			{
			  out2 += dnaA[popid].output_diversity(dnaA[popid2],simuP.fileoutput+"_a_between_div.txt",p,timeG);
			}
		      ofstream fs2;
		      fs2.open((simuP.fileoutput+"_Between.txt").c_str(),ios::app);
		      fs2 << out2 ;
		      fs2.close();
		    }
		}
	    }
	}
    }
  else
    cerr << "Sample size is null " << endl;
  return 0;
}


int Simulation::run()
{
  if(simuP.verbose)
    verbose_start();
  if(simuP.header)
    header_write();
  if(simuP.load == false)
    {
      classic();
    }
  else
    {
      successive_classic();
    }
  if(simuP.verbose)
    verbose_end();
  return 0;
}

void Simulation::output_sfs(vector<Human> & epsilonSub,int counter,int generation)
{
  //  vector<Human> epsilonSub = epsilon.sample(metapopValues.sample,metapopValues.nbPop);
  std::string filenamesfsmt = simuP.fileoutput+"_"+to_string(counter)+".sfsmt";
  std::string filenamesfsx = simuP.fileoutput+"_"+to_string(counter)+".sfsx";
  std::string filenamesfsa = simuP.fileoutput+"_"+to_string(counter)+".sfsa";
  std::string filenamesfsy = simuP.fileoutput+"_"+to_string(counter)+".sfsy";
  if(DNAS.sizeMt>0)
    {
      DNA happy(epsilonSub,metapopValues.sample,(int) 0.5*metapopValues.sample,3);
      happy.SFS_DNA(filenamesfsmt,generation);
      happy.SFS_folded();
    }
  if(DNAS.sizePerX>0)
    {
      DNA happyX(epsilonSub,metapopValues.sample,(int) 0.5*metapopValues.sample,0);
      happyX.SFS_DNA(filenamesfsx,generation);
    }
  if(DNAS.sizePerA>0)
    {
      DNA happyA(epsilonSub,metapopValues.sample,(int) 0.5*metapopValues.sample,2);
      happyA.SFS_DNA(filenamesfsa,generation);
    }
  if(DNAS.sizeY)
    {
      DNA happyY(epsilonSub,metapopValues.sample,(int) 0.5*metapopValues.sample,1);
      happyY.SFS_DNA(filenamesfsy,generation);
    }
}

void Simulation::output_afs(Metapop & epsilon,int counter)
{
  vector<Human> epsilonSub = epsilon.sample(metapopValues.sample,metapopValues.nbPop);
  std::string filenameafsmt = simuP.fileoutput+"_"+to_string(counter)+".afsmt";
  std::string filenameafsx = simuP.fileoutput+"_"+to_string(counter)+".afsx";
  std::string filenameafsa = simuP.fileoutput+"_"+to_string(counter)+".afsa";
  std::string filenameafsy = simuP.fileoutput+"_"+to_string(counter)+".afsy";
  if(DNAS.sizeMt>0)
    {
      DNA happy(epsilonSub,metapopValues.sample,(int) 0.5*metapopValues.sample,3);
      happy.AFS_DNA(filenameafsmt,epsilon.generation);
    }
  if(DNAS.sizePerX>0)
    {
      DNA happyX(epsilonSub,metapopValues.sample,(int) 0.5*metapopValues.sample,0);
      happyX.AFS_DNA(filenameafsx,epsilon.generation);
      for(int j(0);j<DNAS.nbChromosomesA;++j)
	{
	  std::string filenameafsxj = simuP.fileoutput+"_"+to_string(counter)+".afsx"+to_string(j);
	  DNA happy(epsilonSub,metapopValues.sample,(int) 0.5*metapopValues.sample,0,j);
	  happy.AFS_DNA(filenameafsxj,epsilon.generation);
	}
    }
  if(DNAS.sizePerA>0)
    {
      DNA happyA(epsilonSub,metapopValues.sample,(int) 0.5*metapopValues.sample,2);
      happyA.AFS_DNA(filenameafsa,epsilon.generation);
      for(int j(0);j<DNAS.nbChromosomesA;++j)
	{
	  std::string filenameafsaj = simuP.fileoutput+"_"+to_string(counter)+".afsa"+to_string(j);
	  DNA happy(epsilonSub,metapopValues.sample,(int) 0.5*metapopValues.sample,2,j);
	  happy.AFS_DNA(filenameafsaj,epsilon.generation);
	}
    }
  if(DNAS.sizeY)
    {
      DNA happyY(epsilonSub,metapopValues.sample,(int) 0.5*metapopValues.sample,1);
      happyY.AFS_DNA(filenameafsy,epsilon.generation);
    }
}

void Simulation::output_haplotype(Metapop & epsilon, int counter)
{
  std::string filenamehapmt = simuP.fileoutput+"_"+to_string(counter)+"mt.hap";
  std::string filenamehapy = simuP.fileoutput+"_"+to_string(counter)+"y.hap";

  DNA happyMt(epsilon,3);
  happyMt.output_haplotypes(filenamehapmt,epsilon.generation);
  DNA happyY(epsilon,1);
  happyY.output_haplotypes(filenamehapy,epsilon.generation);
}

void Simulation::save(Metapop & epsilon)
{
  std::string filename = simuP.directorySave + "_" + to_string(counter) + ".sim";
  cout << filename << endl;
  std::ofstream ofs(filename.c_str());
}

void Simulation::evolve(Metapop & epsilon)
{
  for(int n(0);n<metapopValues.nbPop;++n)
    {
      epsilon.poptab[n]->set_matingSystem(popValues[n].matingSystem);
      epsilon.poptab[n]->set_inbreeding(popValues[n].inbreeding);
      epsilon.poptab[n]->set_varChild(popValues[n].varChild);
      epsilon.poptab[n]->set_polygamy(popValues[n].polygamy);
      epsilon.poptab[n]->set_polygamyf(popValues[n].polygamyf);
      epsilon.poptab[n]->set_growthFactor(popValues[n].growthFactor);
    }
  if(simuP.burninT>0)
    {
      epsilon.burnin(simuP.burninT,simuP.burninDna,10,simuP.burninDna);
    }
  for(int j(0);j<simuP.nstep;++j)
    {
      epsilon.evolve(simuP.step,simuP.verbose);
      output_main(epsilon);
      if(simuP.hap)
	{
	  output_haplotype(epsilon,counter);
	}
      if(simuP.smartsave)
	{
	  save(epsilon);
	}
    }
}
