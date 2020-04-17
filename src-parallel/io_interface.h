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

#ifndef SIMULATION_SETTINGS_H
#define SIMULATION_SETTINGS_H


#include <set>
#include "utils.h"
#include "msnetwork.h"
#include "dnastructure.h"


struct simuParam
{
  std::vector<int> diversity;
  bool header;
  bool fasta;
  bool nexus;
  bool arl;
  bool hap;
  bool sfs;
  bool smartsave;
  bool load;
  bool verbose;
  //simulation
  unsigned int seed;
  int step;
  int nstep;
  int nsimu;
  double burninT; /** Diversity threshold for the burnin phase */
  int burninDna; /** 1:mtdna ; 2:x ; 3:y ; 4:a */
  std::string fileoutput;
  std::string directoryLoad; /** Name of the directory from where the population have tobe loaded at the bening */ // TODO still useful ?
  std::string directorySave; /** Name of the directory where the population will be saved at the end */
  simuParam()
  {
    diversity.push_back(-1);
    header = false;
    fasta = false;
    nexus = false;
    arl = false;
    hap =false;
    sfs =false;
    load= false;
    smartsave = false;
    verbose = false;
  //simulation
    seed = time(0);
    step = 1000;
    nstep = 1;
    nsimu = 100;
    fileoutput = "smartpop_"+to_string(seed);
    burninT = 0;
    burninDna = 1;
    directoryLoad = "";
    directorySave = "";
  }
  simuParam(simuParam const & a)
  {
    diversity = a.diversity;
    header = a.header;
    fasta = a.fasta;
    nexus = a.nexus;
    arl = a.arl;
    hap = a.hap;
    sfs = a.sfs;
    load = a.load;
    smartsave = a.smartsave;
    verbose = a.verbose;
  //simulation
    seed = a.seed;
    step = a.step;
    nstep = a.nstep;
    nsimu = a.nsimu;
    fileoutput = a.fileoutput;;
    burninT = a.burninT;
    burninDna = a.burninDna;
    directoryLoad = a.directoryLoad;
    directorySave = a.directorySave;
  }
  simuParam & operator=(simuParam const & a)
  {
    if(&a!=this)
      {
	diversity = a.diversity;
	header = a.header;
	fasta = a.fasta;
	nexus= a.nexus;
	arl = a.arl;
	hap = a.hap;
	sfs = a.sfs;
	load = a.load;
	smartsave = a.smartsave;
	verbose = a.verbose;
	//simulation
	seed = a.seed;
	step = a.step;
	nstep = a.nstep;
	nsimu = a.nsimu;
	fileoutput = a.fileoutput;;
	burninT = a.burninT;
	burninDna = a.burninDna;
	directoryLoad = a.directoryLoad;
	directorySave = a.directorySave;
      }
    return *this;
  }  
};

struct popVal
{
  int popsize;
  int inbreeding;
  double varChild;
  int polygamy;
  int polygamyf;
  double growthFactor[3];
  std::vector<int> matingSystem;
  popVal()
  {
    popsize = 200;
    inbreeding = 0;
    varChild = 1;
    polygamy = 1;
    polygamyf = 1;
    growthFactor[0]=0;
    growthFactor[1]=1;
    growthFactor[2]=0;
    matingSystem.push_back(1);
  }
  popVal(popVal const & a)
  {
    popsize = a.popsize;
    inbreeding = a.inbreeding;
    varChild = a.varChild;
    polygamy = a.polygamy;
    polygamyf = a.polygamyf;
    growthFactor[0]=a.growthFactor[0];
    growthFactor[1]=a.growthFactor[1];
    growthFactor[2]=a.growthFactor[2];
    matingSystem = a.matingSystem;      
  }
  popVal & operator=(popVal const & a)
  {
    if(&a != this)
      {
	popsize = a.popsize;
	inbreeding = a.inbreeding;
	varChild = a.varChild;
	polygamy = a.polygamy;
	polygamyf = a.polygamyf;
	growthFactor[0]=a.growthFactor[0];
	growthFactor[1]=a.growthFactor[1];
	growthFactor[2]=a.growthFactor[2];
	matingSystem = a.matingSystem; 
      }     
    return *this;
  }
};

const struct popVal defaultPopVal;
const struct popVal zeroPopVal;

struct metapopVal
{
  int popsize;
  int sample;
  int sampleType;
  bool sampleBetween;
  int nbsamplepop;
  double migrationType;
  int nbPop;
  int nbSources;
  double pmate;
  double pmig;
  std::vector< std::vector<double> > migrationRatesF;
  std::vector< std::vector<double> > migrationRatesM;
  metapopVal()
  {
    popsize = 200;
    sample =50;
    migrationType = 0;
    nbPop = 1;
    pmate = 0;
    pmig = 0;
    nbSources = 0;
    sampleType = 2;
    nbsamplepop = 1;
    sampleBetween = false;
  }
  metapopVal(metapopVal const & a )
  {
    popsize = a.popsize;
    sample = a.sample;
    migrationType = a.migrationType;
    nbPop = a.nbPop;
    pmate = a.pmate;
    pmig = a.pmig;
    migrationRatesF = a.migrationRatesF;
    migrationRatesM = a.migrationRatesM;
    nbSources = a.nbSources;
    nbSources = a.nbSources;
    sampleType = a.sampleType;
    nbsamplepop = a.nbsamplepop;
    sampleBetween = a.sampleBetween;
  }
  metapopVal& operator=(metapopVal const & a )
  {
    if(&a != this)
      {
	popsize = a.popsize;
	sample = a.sample;
	migrationType = a.migrationType;
	nbPop = a.nbPop;
	pmate = a.pmate;
	pmig = a.pmig;
	migrationRatesF = a.migrationRatesF;
	migrationRatesM = a.migrationRatesM;
	nbSources = a.nbSources;
	sampleType = a.sampleType;
	nbsamplepop = a.nbsamplepop;
	sampleBetween = a.sampleBetween;
      }
    return *this;
  }
};

const struct metapopVal defaultMetapopVal;
//const struct metapopVal zeroMetapopVal = {.popsize = 0, .sample=50,.migrationType = 0, .nbPop = 0};
const struct metapopVal zeroMetapopVal;


struct simulation_settings
{
  metapopVal metapop;
  simuParam simuP;
  std::vector<popVal> populations;
  MutationRates mu;
  DnaStructure pdna;
  simulation_settings()
  {   
  }
  simulation_settings(simulation_settings const & a)
  {   
    mu = a.mu;
    pdna = a.pdna;
    simuP = a.simuP;
    metapop = a.metapop;
    populations = a.populations;
  }
  void load(const std::string &filename);
  void save(const std::string &filename);
  void load(metapopVal mtpv,std::vector<popVal> pv,simuParam psimu, MutationRates pmu, DnaStructure ppdna)
  {
    metapop = mtpv;
    populations = pv;
    mu = pmu;
    pdna = ppdna;
    simuP = psimu;
  }
};

#endif
