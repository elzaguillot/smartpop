
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

#include "dna.h"
using namespace std;


DNA::DNA()
{
  N = 0;
  nbXcell = 0;
  nbS = 0;
  nbSreal = 0;
  hh = -1;
  pw = -1;
  tw = -1;
  ps = -1;
  hs = -1;
  nbh = -1;
}

DNA::DNA(Population & pop)
{
  N = pop.get_popsize();
  nbXcell = DNAS.nbCrumbMt;
  nbS = DNAS.sizeMt;
  nbSreal = DNAS.sizeMt;
  assert(nbXcell*32 == nbS);
  if(N==1)
    {
      tab = new typeDNA*[1];
      tab[1] = new typeDNA[nbXcell];
    }
  else
    {
      tab = new typeDNA*[N];
      for(int i(0);i<N;++i)
	tab[i] = new typeDNA[nbXcell];
    }
  int nbFemale = pop.get_nbFemale();
  int currentG = pop.get_generation();
  currentG &= 1;
  for(int i(0); i<N; ++i)
    {
      for(int j(0) ;j<nbXcell ;++j)
	{
	  if(i < nbFemale)
	    tab[i][j] = pop.female[currentG][i]->mtDna[j];
	  else
	    tab[i][j] = pop.male[currentG][i-nbFemale]->mtDna[j];
	}
    }
  hh = -1;
  pw = -1;
  tw = -1;
  ps = -1;
  hs = -1;
  nbh = -1;
}

DNA::DNA(Population  & pop,int tdna)
{
  int nbFemale = pop.get_nbFemale();
  int Npop = pop.get_popsize();
  int currentG = pop.get_generation();
  currentG &= 1;
  switch(tdna)
    {
    case 0: // X chromosome
      {
	N = Npop + nbFemale;
	nbXcell = DNAS.nbCrumbX;
	nbS = nbXcell * CRUMB;
	tab = new typeDNA*[N];
	for(int i(0);i<N;++i)
	  tab[i] = new typeDNA[nbXcell];
 	int index = DNAS.nbCrumbA;
	for(int i(0); i<Npop; ++i)
	  {	  
	    for(int j(0); j<nbXcell; ++j)
	      {
		if(i < nbFemale)
		  {
		    tab[i][j] = pop.female[currentG][i]->dna[index+j]; 
		  }
		else
		  {
		    tab[i][j] = pop.male[currentG][i-nbFemale]->dna[index+j];
		  }
	      }}	      
	for(int i(0); i<nbFemale; ++i)
	  {
	    for(int j(0); j<nbXcell; ++j)
	      {
		tab[i+Npop][j] = pop.female[currentG][i]->dna[index+nbXcell+j];
	      }	      
	  }
	break;
      }
    case 1: // Y chromosome
      {
	N = Npop - nbFemale;
	nbXcell = DNAS.nbCrumbY;
	nbS = nbXcell * CRUMB;
	tab = new typeDNA*[N];
	for(int i(0);i<N;++i)
	  tab[i] = new typeDNA[nbXcell];
	int index = DNAS.nbCrumbA + DNAS.nbCrumbX;
	for(int i(0); i<N; ++i)
	  {
	    for(int j(0); j<nbXcell; ++j)
	      {
		tab[i][j] = pop.male[currentG][i]->dna[index+j];
	      }
	  }
	break;
      }
    case 2: // autosomes
      {
	N = 2*Npop;
	nbXcell = DNAS.nbCrumbA/2; // only half because it is diploid
	nbS = nbXcell * CRUMB;
	tab = new typeDNA*[N];
	for(int i(0);i<N;++i)
	  tab[i] = new typeDNA[nbXcell];
	int index1 =0;
	int index2 =0;
	int cellPerX = DNAS.sizePerA/32;
	int nbChromosomes = DNAS.nbChromosomesA;
	for(int i(0); i<Npop; ++i)
	  {	  
	    for(int j(0); j<nbChromosomes; j+=2)
	      {
		for(int k(0); k<cellPerX; ++k)
		  {
		    index1 = j*cellPerX;
		    index2 = (j/2)*cellPerX;
		    if(i < nbFemale)
		      {
			tab[i][index2+k] = pop.female[currentG][i]->dna[index1+k];
			tab[i+Npop][index2+k] = pop.female[currentG][i]->dna[index1+cellPerX+k];
		      }
		    else
		      {
			tab[i][index2+k] = pop.male[currentG][i-nbFemale]->dna[index1+k];
			tab[i+Npop][index2+k] = pop.male[currentG][i-nbFemale]->dna[index1+cellPerX+k];
		      }
		  }
	      }
	  } 
	break;
      }
    case 3:
      {
	N = Npop;
	nbXcell = DNAS.nbCrumbMt;
	nbS = DNAS.sizeMt;
	assert(nbXcell*32 == nbS);
	tab = new typeDNA*[N];
	for(int i(0);i<N;++i)
	  tab[i] = new typeDNA[nbXcell];
	int nbFemale = pop.get_nbFemale();
	int currentG = pop.get_generation();
	currentG &= 1;
	for(int i(0); i<N; ++i)
	  {
	    for(int j(0) ;j<nbXcell ;++j)
	      {
		if(i < nbFemale)
		  tab[i][j] = pop.female[currentG][i]->mtDna[j];
		else
		  tab[i][j] = pop.male[currentG][i-nbFemale]->mtDna[j];
	      }
	  }
      }
    default:
      {
	cerr << " Call DNA () for non existing chromosome:" << tdna << endl;
	exit(3);
	break;
      }
    }  
  nbSreal = nbS;
  hh = -1;
  pw = -1;
  tw = -1;
  ps = -1;
  hs = -1;
  nbh = -1;
}


DNA::DNA(Metapop  & metapop,int tdna)
{
  int nbFemale = metapop.nbFemale;
  int Npop = metapop.popsize;
  int currentG = metapop.generation;
  int offset = 0;
  currentG &= 1;
  switch(tdna)
    {
    case 0: // X chromosome
      {
	N = Npop + nbFemale;
	nbXcell = DNAS.nbCrumbX;
	nbS = nbXcell * CRUMB;
	tab = new typeDNA*[N];
	for(int i(0);i<N;++i)
	  tab[i] = new typeDNA[nbXcell];
	int index = DNAS.nbCrumbA;	
	for(int n(0);n<metapop.nbPop;++n)
	  {
	    int nf = metapop.poptab[n]->nbFemale;
	    int ntot = metapop.poptab[n]->popsize;
	    for(int i(0); i<metapop.poptab[n]->popsize; ++i)
	      {	  
		for(int j(0); j<nbXcell; ++j)
		  {
		    if(i < metapop.poptab[n]->nbFemale)
		      {
			tab[i+offset][j] = metapop.poptab[n]->female[currentG][i]->dna[index+j];
		      }
		    else
		      {
			tab[i+offset][j] = metapop.poptab[n]->male[currentG][i-nf]->dna[index+j];
		      }
		  }}	      
	    for(int i(0); i<metapop.poptab[n]->nbFemale; ++i)
	      {
		for(int j(0); j<nbXcell; ++j)
		  {
		    tab[i+ntot+offset][j] = metapop.poptab[n]->female[currentG][i]->dna[index+nbXcell+j];
		  }
	      }
	    offset += ntot + metapop.poptab[n]->nbFemale;
	  }
	break;
      }
    case 1: // Y chromosome
      {
	N = Npop - nbFemale;
	nbXcell = DNAS.nbCrumbY;
	nbS = nbXcell * CRUMB;
	tab = new typeDNA*[N];
	for(int i(0);i<N;++i)
	  tab[i] = new typeDNA[nbXcell];
	int index = DNAS.nbCrumbA + DNAS.nbCrumbX;
	for(int n(0);n<metapop.nbPop;++n)
	  {
	    for(int i(0); i<(metapop.poptab[n]->popsize-metapop.poptab[n]->nbFemale); ++i)
	      {
		for(int j(0); j<nbXcell; ++j)
		  {
		    tab[i+offset][j] = metapop.poptab[n]->male[currentG][i]->dna[index+j];
		  }
	      }
	    offset += metapop.poptab[n]->popsize-metapop.poptab[n]->nbFemale; // += nbmale
	  }
	break;
      }
    case 2: // autosomes
      {
	N = 2*Npop;
	nbXcell = DNAS.nbCrumbA/2; // only half because it is diploid
	nbS = nbXcell * CRUMB;
	tab = new typeDNA*[N];
	for(int i(0);i<N;++i)
	  tab[i] = new typeDNA[nbXcell];
	int index1 =0;
	int index2 =0;
	int cellPerX = DNAS.sizePerA/32;
	int nbChromosomes = DNAS.nbChromosomesA;
	for(int n(0);n<metapop.nbPop;++n)
	  {
	    int nf = metapop.poptab[n]->nbFemale;
	    int ntot = metapop.poptab[n]->popsize;
	    for(int i(0); i<metapop.poptab[n]->popsize; ++i)
	      {	
		for(int j(0); j<nbChromosomes; j+=2)
		  {
		    for(int k(0); k<cellPerX; ++k)
		      {
			index1 = j*cellPerX;
			index2 = (j/2)*cellPerX;
			if(i < metapop.poptab[n]->nbFemale)
			  {
			    tab[i+offset][index2+k] = metapop.poptab[n]->female[currentG][i]->dna[index1+k];
			    tab[i+ntot+offset][index2+k] = metapop.poptab[n]->female[currentG][i]->dna[index1+cellPerX+k];
			  }
			else
			  {
			    tab[i+offset][index2+k] = metapop.poptab[n]->male[currentG][i-nf]->dna[index1+k];
			    tab[i+ntot+offset][index2+k] = metapop.poptab[n]->male[currentG][i-nf]->dna[index1+cellPerX+k];
			  }
		      }
		  }
	      }
	    offset += ntot * 2;
	  } 
	break;       
      }
    case 3: // mtdna
      {
	N = Npop;
	nbXcell = DNAS.nbCrumbMt; // only half because it is diploid
	nbS = nbXcell * CRUMB;
	tab = new typeDNA*[N];
	for(int i(0);i<N;++i)
	  tab[i] = new typeDNA[nbXcell];
	for(int n(0);n<metapop.nbPop;++n)
	  {
	    int nf = metapop.poptab[n]->nbFemale;
	    for(int i(0); i<metapop.poptab[n]->popsize; ++i)
	      {
		if(i < metapop.poptab[n]->nbFemale)
		  {
		    for(int x(0);x<nbXcell;++x)
		      tab[i+offset][x]=metapop.poptab[n]->female[currentG][i]->mtDna[x];
		  }
		else
		  {
		    for(int x(0);x<nbXcell;++x)
		      tab[i+offset][x]=metapop.poptab[n]->male[currentG][i-nf]->mtDna[x];
		  }
	      }
	    offset += metapop.poptab[n]->popsize;
	  }
	break;       
      }
    default:
      {
	cerr << " Call DNA () for non existing chromosome:" << tdna << endl;
	exit(3);
	break;
      }
    }  
  nbSreal = nbS;
  hh = -1;
  pw = -1;
  tw = -1;
  ps = -1;
  hs = -1;
  nbh = -1;
}

DNA::DNA(Population  & pop,int tdna,int loci)
{ // extract a loci from X chromosome or autosome
  int nbFemale = pop.get_nbFemale();
  int Npop = pop.get_popsize();
  int currentG = pop.get_generation();
  currentG &= 1;
  switch(tdna)
    {
    case 0: // X chromosome
      {
	N = Npop + nbFemale;
	nbXcell = DNAS.sizePerX/CRUMB;
	nbS = DNAS.sizePerX;
	int offset = loci * nbXcell + DNAS.nbCrumbA;
	tab = new typeDNA*[N];
	for(int i(0);i<N;++i)
	  tab[i] = new typeDNA[nbXcell];
	for(int i(0); i<Npop; ++i)
	  {
	    for(int j(0); j<nbXcell; ++j)
	      {
		if(i < nbFemale)
		  {
		    tab[i][j] = pop.female[currentG][i]->dna[offset+j];
		  }
		else
		  {
		    tab[i][j] = pop.male[currentG][i-nbFemale]->dna[offset+j];
		  }
	      }
	  }
	for(int i(0); i<nbFemale; ++i)
	  {
	    for(int j(0); j<nbXcell; ++j)
	      {
		tab[i+Npop][j] = pop.female[currentG][i]->dna[offset+nbXcell+j];	      
	      }
	  }
	break;
      }
    case 1: // Y chromosome
      {
	cerr << "ERROR : there is only one loci on Y chromosome " << endl;
	exit(2);
	break;
      }
    case 2: // autosomes
      {
	N = 2 * Npop;
	nbS = DNAS.sizePerA;
	nbXcell = nbS/CRUMB; // only half because it is diploid
	tab = new typeDNA*[N];
	for(int i(0);i<N;++i)
	  tab[i] = new typeDNA[nbXcell];
	int offset = loci * nbXcell;
	for(int i(0); i<Npop; ++i)
	  {
	    for(int k(0); k<nbXcell; ++k)
	      {
		if(i < nbFemale)
		  {
		    tab[i][k] = pop.female[currentG][i]->dna[offset+k];
		    tab[i+Npop][k] = pop.female[currentG][i]->dna[offset+nbXcell+k];
		  }
		else
		  {
		    tab[i][k] = pop.male[currentG][i-nbFemale]->dna[offset+k];
		    tab[i+Npop][k] = pop.male[currentG][i-nbFemale]->dna[offset+nbXcell+k];
		  }
	      }
	  } 
	break;       
      }
    case 3:
      cerr << "ERROR there is only one loci on the mtdna " << endl;
    default:
      {
	cerr << " Call DNA () for non existing chromosome:" << tdna << endl;
	exit(3);
	break;
      }
    }
  nbSreal = nbS;
  hh = -1;
  pw = -1;
  tw = -1;
  ps = -1;
  hs = -1;
  nbh = -1;
}

DNA::DNA(vector< vector<typeDNA> > input)
{
  N = input.size();
  nbXcell = input[0].size();
  if((N >1)&&(nbXcell>1))
    {
      tab = new typeDNA*[N];
      for(int i(0);i<N;++i)
	tab[i] = new typeDNA[nbXcell];
      nbS = nbXcell * CRUMB;
      for(int i(0); i<N; ++i)
	{
	  for(int j(0); j<nbXcell; ++j)
	    {
	      tab[i][j] = input[i][j];
	    }
	}
    }
  else
    {
      cerr << "Error in the input of genetics - N : " << N << " and nbXcell " << nbXcell << endl;
      exit(3);
    }
  nbSreal = nbS;
  hh = -1;
  pw = -1;
  tw = -1;
  ps = -1;
  hs = -1;
  nbh = -1;
}

DNA::DNA(std::vector<Human> & indivList,int nHuman, int nFemale, int tdna)
{
  assert(nHuman == (int) indivList.size());
  int nbFemale = nFemale;
  int Npop = nHuman;
  switch(tdna)
    {
    case 0: // X chromosome
      {
	N = Npop + nbFemale;
	nbXcell = DNAS.nbCrumbX;
	nbS = nbXcell * CRUMB;
	tab = new typeDNA*[N];
	for(int i(0);i<N;++i)
	  tab[i] = new typeDNA[nbXcell];
	int index = DNAS.nbCrumbA;
	for(int i(0); i<Npop; ++i)
	  {	  
	    for(int j(0); j<nbXcell; ++j)
	      {
		tab[i][j] = indivList[i].dna[index+j]; 
	      }
	  }
	for(int i(0); i<nbFemale; ++i)
	  {
	    for(int j(0); j<nbXcell; ++j)
	      {
		tab[i+Npop][j] = indivList[i].dna[index+nbXcell+j];
	      }
	  }
	break;
      }
    case 1: // Y chromosome
      {
	N = Npop - nbFemale;
	nbXcell = DNAS.nbCrumbY;
	nbS = nbXcell * CRUMB;
	tab = new typeDNA*[N];
	for(int i(0);i<N;++i)
	  tab[i] = new typeDNA[nbXcell];
	int index = DNAS.nbCrumbA + DNAS.nbCrumbX;
	for(int i(0); i<N; ++i)
	  {
	    for(int j(0); j<nbXcell; ++j)
	      {
		tab[i][j] = indivList[i+nFemale].dna[index+j];
	      }
	  }
	break;
      }
    case 2: // autosomes
      {
	N = 2*Npop;
	nbXcell = DNAS.nbCrumbA/2; // only half because it is diploid
	nbS = nbXcell * CRUMB;
	tab = new typeDNA*[N];
	for(int i(0);i<N;++i)
	  tab[i] = new typeDNA[nbXcell];
	int index1 =0;
	int index2 =0;
	int cellPerX = DNAS.nbCellPerA;
	int nbChromosomes = DNAS.nbChromosomesA;
	for(int i(0); i<Npop; ++i)
	  {	  
	    for(int j(0); j<nbChromosomes; ++j)
	      {
		index1 = j*2*cellPerX;
		index2 = j*cellPerX;
		for(int k(0); k<cellPerX; ++k)
		  {
		    tab[i][index2+k] = indivList[i].dna[index1+k];
		    tab[i+Npop][index2+k] = indivList[i].dna[index1+cellPerX+k];
		  }
	      }
	  } 
	break;       
      }
    case 3: // mtdna
      {
	N = nHuman;
	nbXcell = DNAS.nbCrumbMt;
	nbS = DNAS.sizeMt;
	assert(nbXcell*32 == nbS);
	tab = new typeDNA*[N];
	for(int i(0);i<N;++i)
	  tab[i] = new typeDNA[nbXcell];
	for(int i(0); i<N; ++i)
	  {
	    for(int j(0) ;j<nbXcell ;++j)
	      {
		tab[i][j] = indivList[i].mtDna[j];
	      }
	  } 
	break;
      }
    default:
      {
	cerr << " Call DNA () for non existing chromosome :" << tdna  << endl;
	exit(3);
	break;
      }
    }
  nbSreal = nbS;
  hh = -1;
  pw = -1;
  tw = -1;
  ps = -1;
  hs = -1;
  nbh = -1;
}

DNA::DNA(std::vector<Human> & indivList,int nHuman, int nFemale, int tdna,int loci)
{
  assert(nHuman == (int) indivList.size());
  int nbFemale = nFemale;
  int Npop = nHuman;
  switch(tdna)
    {
    case 0: // X chromosome
      {
	N = Npop + nbFemale;
	nbXcell = DNAS.sizePerX/CRUMB;
	nbS = DNAS.sizePerX;
	int offset = loci * nbXcell + DNAS.nbCrumbA;
	tab = new typeDNA*[N];
	for(int i(0);i<N;++i)
	  tab[i] = new typeDNA[nbXcell];
	for(int i(0); i<Npop; ++i)
	  {
	    for(int j(0); j<nbXcell; ++j)
	      {
		if(i < nbFemale)
		  {
		    tab[i][j] = indivList[i].dna[offset+j];
		  }
		else
		  {
		    tab[i][j] = indivList[i].dna[offset+j];
		  }
	      }
	  }
	for(int i(0); i<nbFemale; ++i)
	  {
	    for(int j(0); j<nbXcell; ++j)
	      {
		tab[i+Npop][j] = indivList[i].dna[offset+nbXcell+j];	      
	      }
	  }
	break;
      }
    case 1: // Y chromosome
      {
	cerr << "ERROR : there is only one loci on Y chromosome " << endl;
	exit(2);
	break;
      }
    case 2: // autosomes
      {
	N = 2 * Npop;
	nbS = DNAS.sizePerA;
	nbXcell = nbS/CRUMB; // only half because it is diploid
	tab = new typeDNA*[N];
	for(int i(0);i<N;++i)
	  tab[i] = new typeDNA[nbXcell];
	int offset = loci * nbXcell;
	//	int index1 = offset;
	//	int cellPerX = pop->get_sizePerA()/32;
	//	int nbChromosomes = pop->get_nbChromosomesA();
	for(int i(0); i<Npop; ++i)
	  {
	    for(int k(0); k<nbXcell; ++k)
	      {
		//		index1 = offset;
		if(i < nbFemale)
		  {
		    tab[i][k] = indivList[i].dna[offset+k];
		    tab[i+Npop][k] = indivList[i].dna[offset+nbXcell+k];
		  }
		else
		  {
		    tab[i][k] = indivList[i].dna[offset+k];
		    tab[i+Npop][k] = indivList[i].dna[offset+nbXcell+k];
		  }
	      }
	  } 
	break;       
      }
    case 3: // mtdna
      {
	cerr << " ERROR: there is only one loci on the mtdna " << endl;
      }
    default:
      {
	cerr << " Call DNA () for non existing chromosome:" << tdna << endl;
	exit(3);
	break;      
      }
    }
  nbSreal = nbS;
  hh = -1;
  pw = -1;
  tw = -1;
  ps = -1;
  hs = -1;
  nbh = -1;
}

DNA::DNA(DNA const & a)
{
  N = a.N;
  nbXcell = a.nbXcell;
  nbS = a.nbS;
  if((a.N>0)&&(a.nbXcell>0))
    {
      tab = new typeDNA*[N];
      for(int i(0);i<N;++i)
	{
	  tab[i] = new typeDNA[nbXcell];
	  memcpy(tab[i],a.tab[i],nbXcell*sizeof(typeDNA));
	}
    }
  else
    {
      cerr << "ERROR attempt to create an empty DNA() " << endl;
      cerr <<  " N " << a.N << " - nbXcell " << a.nbXcell << endl;
    }
  nbSreal = a.nbSreal;
  hh = a.hh;
  pw = a.pw;
  tw = a.tw;
  ps = a.ps;
  hs = a.hs;
  nbh = a.nbh;
}

DNA::DNA(array_type const  & input,int sizeN, int sizeX)//:tab(boost::extents[2][NB_MTDNA_CELL])
{
  nbS = sizeX * CRUMB;
  N = sizeN;
  nbXcell = sizeX;  
  tab = new typeDNA*[N];
  for(int i(0);i<N;++i)
    {
      tab[i] = new typeDNA[nbXcell];
      memcpy(&tab[i],&input[i],sizeof(typeDNA));
    }
  nbSreal = nbS;
  hh = -1;
  pw = -1;
  tw = -1;
  ps = -1;
  hs = -1;
  nbh = -1;
}


DNA::~DNA()
{
  if(N>0)
    {
      if(nbXcell>0){
	for(int i(0);i<N;++i)
	  delete [] tab[i];
      }
      delete [] tab;
    }
}

ostream &operator<<(ostream & os, DNA const& a)
{
  a.print(os);
  return os;
}

ostream& DNA::print(ostream & os) const
{
  for(int i(0);i<N;++i) 
    {
      for(int j(0);j<nbXcell;++j)
	{
	  os << convertIntToQuaternaryString(tab[i][j],CRUMB);
	}
      os << "\n";
    }
  return os;
}

/*DNA operator+(DNA const & a, DNA const & b)
  {
  DNA out(a,b);
  return out;
  }*/

DNA& DNA::operator=(DNA const & a)
{
  if(this != &a)
    {
      N = a.N;
      nbXcell = a.nbXcell;
      nbS = a.nbS;
      tab = new typeDNA*[N];
      for(int i(0);i<N;++i)
	{
	  tab[i] = new typeDNA[nbXcell];
	  memcpy(tab[i],a.tab[i],nbXcell*sizeof(typeDNA));
	}
      nbSreal = nbS;
      hh = a.hh;
      pw = a.pw;
      tw = a.tw;
      ps = a.ps;
      hs = a.hs;
      nbh = a.nbh;
    }
  return *this;
}

int DNA::get_N() const
{
  return N;
}

int DNA::get_nbXcell() const
{
  return nbXcell;
}

int DNA::get_nbS() const
{
  return nbS;
}

int DNA::get_nbSreal() const
{
  return nbS;
}

void DNA::set_nbSreal(int a)
{
  if(a>nbS)
    cerr << "ERROR in the inuput, sequences length do not match";
  nbSreal = a;
}

long DNA::poly_sites() const
{
  typeDNA tot(0),c(0);
  typeDNA tmp(0),tmp2(0),tmp3(0); // storing the number is faster
  for(int k(0); k<nbXcell; ++k) 
    {
      c = 0;
      tmp = tab[0][k];
      for(int j(1); j<N; ++j)
	{
	  tmp2 = tmp ^ tab[j][k];
	  c |= tmp2;
	}
      tmp3 = dna_distance(c);
      tot += tmp3;
    }
  return tot;  
}

long DNA::poly_sites2(std::vector<int> tsfs) const
{
  int tot=0;
  for(unsigned int i(1);i<tsfs.size()-1;++i)
    tot+= tsfs[i];
  return tot;  
}

double DNA::prop_poly_sites() const
{
  double out;
  out = (double) poly_sites()/(double) nbSreal;
  return out;
}

double DNA::mean_pw_diff() const
{
  typeDNA x = 0;
  double tot(0);
  double factoriel = 0.0;
  for(int i(0); i<N-1; ++i)
    for(int j(i+1); j<N; ++j)
      {
	x += dist_two(i,j);
	++factoriel;
      }
  tot += (double) x /(factoriel);  
  return tot ;// (double) nbS;
}

double DNA::mean_pw_diff2(std::vector<int> SFScounts) const
{
  double tot = 0;
  int n = SFScounts.size();
  for (int i(1);i<N;++i)
    {
      tot +=  i*(N-i)*SFScounts[i];
    }
  tot *= 2.0 / (double) (N*(N-1));
  return tot ;
}


double DNA::mean_pw_diff_between(DNA & dna2) const
{
  if(nbXcell == dna2.nbXcell)
    {
      typeDNA x = 0;
      double tot = 0;
      int factoriel = 0;
      for(int i(0); i<N; ++i)
	for(int j(0); j<dna2.N; ++j)
	  {
	    x += dist_two(i,j,dna2);
	    ++factoriel;
	  }
      tot += (double) x /((double) factoriel);
      return tot ;// (double) nbS;
    }
  else
    cerr << " ERROR : attempt to compare diversity of different DNA " << endl;
  return 0;
}

double DNA::Nei_D(DNA & dna2) const
{
  if(nbXcell == dna2.nbXcell)
    {
      double tot = 0;
      double D1 =0;
      double D2 = 0;
      double D3 = 0;
      std::vector<double> a,b;
      for(int i(0); i<nbXcell; ++i)
	{
	  a=site_frequency_DNA(i);
	  b=dna2.site_frequency_DNA(i);
	  for(int j(0);j<4;++j)
	    {
	      D1+=a[j]*b[j];
	      D2+=pow(a[j],2);
	      D3+=pow(b[j],2);
	    }
	}
      tot=-log(D1/(sqrt(D2)*sqrt(D3)));
      return tot;
    }
  else
    cerr << " ERROR : attempt to compare diversity of different DNA " << endl;
  return 0;
}

double DNA::privateSNP(DNA & dna2) const
{
  if(nbXcell == dna2.nbXcell)
    {
      double tot = 0;
      double nbPriv =0;
      std::vector<double> a,b;
      for(int i(0); i<nbXcell; ++i)
	{
	  a=site_frequency_DNA(i);
	  b=dna2.site_frequency_DNA(i);
	  for(int j(0);j<4;++j)
	    {
	      if(a[j]>0)
		if(b[j]==0)
		  ++nbPriv;
	    }
	}
      return nbPriv;
    }
  else
    cerr << " ERROR : attempt to compare diversity of different DNA " << endl;
  return 0;
}

double DNA::Chord_D(DNA & dna2) const
{
  if(nbXcell == dna2.nbXcell)
    {
      double tot = 0;
      double D1 =0;
      std::vector<double> a,b;
      for(int i(0); i<nbXcell; ++i)
	{
	  a=site_frequency_DNA(i);
	  b=dna2.site_frequency_DNA(i);
	  for(int j(0);j<4;++j)
	    {
	      D1+=sqrt(a[j]*b[j]);
	    }
	}      
      if(D1/nbXcell>1) // problem of round 1 becomes 1.000000000000000001
	tot=0;
      else
	tot=(2.0/3.141593)*sqrt(2-2*D1/nbXcell);
      return tot;
    }
  else
    cerr << " ERROR : attempt to compare diversity of different DNA " << endl;
  return 0;
}


double DNA::RWC_D(DNA & dna2) const
{
  if(nbXcell == dna2.nbXcell)
    {
      double tot = 0;
      double D1 =0;
      double D2 = 0;
      std::vector<double> a,b;
      for(int i(0); i<nbXcell; ++i)
	{
	  a=site_frequency_DNA(i);
	  b=dna2.site_frequency_DNA(i);
	  for(int j(0);j<4;++j)
	    {
	      D1+=pow((a[j]-b[j]),2);
	      D2+=a[j]*b[j];
	    }
	}
      tot=sqrt(D1/(2.0*nbXcell-2.0*D2));
      return tot;
    }
  else
    cerr << " ERROR : attempt to compare diversity of different DNA " << endl;
  return 0;
}

int DNA::nb_of_haplotypes() const
{
  vector<double> hapFreq = hap_frequency();
  return hapFreq.size();
}

vector<double> DNA::hap_frequency() const
{
  if(N>2 && nbXcell>0)
    {
      int nbHap = 1;
      int k = 0; // loop on haplotypes
      int j = 0; // loop on chromosome cell
      int matchHap;
      array_type haplotypes;
      vector<double> pi; // frequency of the haplotypes
      haplotypes.push_back(tab[0]);
      //     haplotypes.push_back(new typeDNA[nbXcell]);
      //      memcpy(&haplotypes[0],&tab[0],sizeof(typeDNA)*nbXcell);
      pi.push_back(1.0 / (double) N);
      if(N>1) 
	{
	  for(int i(1); i<N; ++i) // for each indiv
	    {
	      j = 0;
	      k = 0;
	      do{ // for each haplotype
		matchHap = 1;
		j = 0;
		do{ // for each chromosome cell
		  if(tab[i][j] ^ haplotypes[k][j])
		    {
		      matchHap = 0;	      
		      break;	      	      
		    }
		  ++ j;
		}
		while((j < nbXcell));
		if (matchHap == 1)
		  {// went to the end and found no difference, if diff in the last call j = nbXcell -1 
		    pi[k] += 1.0 / (double) N;
		  }
		else
		  {
		    ++ k;	
		  }
	      }
	      while((k < nbHap) && (matchHap == 0));
	      if (matchHap == 0)
		{
		  haplotypes.push_back(tab[i]);
		  //		  haplotypes.push_back(new typeDNA[nbXcell]);
		  //		  memcpy(&haplotypes[nbHap],&tab[i],sizeof(typeDNA)*nbXcell);
		  //		  haplotypes[nbHap] = tab[i];
		  pi.push_back(1.0 / (double) N);
		  ++ nbHap;
		}
	    }    
	}
      return pi;
    }
  else
    {
      vector<double> out(1,0);
      return out;
    }
}

DNA DNA::haplotypes() const
{
  if(N>1 && nbXcell>0)
    {
      int nbHap = 1;
      int k = 0; // loop on haplotypes
      int j = 0; // loop on chromosome cell
      int matchHap;
      array_type haplotypes;
      haplotypes.push_back(tab[0]);
      //     haplotypes.push_back(new typeDNA[nbXcell]);
      //      memcpy(&haplotypes[0],&tab[0],sizeof(typeDNA)*nbXcell);
      vector<double> pi; // frequency of the haplotypes
      for(int i(1); i<N; ++i) // for each indiv
	{
	  j = 0;
	  k = 0;
	  do{ // for each haplotype
	    matchHap = 1;
	    j = 0;
	    do{ // for each chromosome cell
	      if(tab[i][j] ^ haplotypes[k][j])
		{
		  matchHap = 0;	      
		  break;	      	      
		}
	      ++ j;
	    }
	    while((j < nbXcell));
	    if(matchHap == 0)
	      {
		++ k;	
	      }
	  }
	  while((k < nbHap) && (matchHap == 0));
	  if (matchHap == 0)
	    {
	      haplotypes.push_back(tab[i]);
	      // haplotypes.push_back(new typeDNA[nbXcell]);
	      // memcpy(&haplotypes[nbHap],&tab[i],sizeof(typeDNA)*nbXcell);
	      ++ nbHap;
	    }
	}
      DNA out(haplotypes,nbHap,nbXcell);
      return out;
    }
  else
    {
      DNA out;
      return out;
    }
}

vector<double> DNA::site_frequency(int site) const
{
  int nbHap = 1;
  int k = 0; // loop on haplotypes
  int matchHap;
  int haplotypes[N];
  vector<double> pi; // frequency of the haplotypes
  if( nbXcell>0){
    int cell =  site / (CRUMB);
    int position = 2 *(site - CRUMB * cell);
    typeDNA screen;
    typeDNA tmp;
    screen = 3;
    screen <<= position;
    tmp = (tab[0][cell] & screen);
    haplotypes[0] = (int) (tmp >> position);
    pi.push_back(1.0 / (double) N);
    for(int i(1); i<N; ++i) // for each indiv
      {
	matchHap = 0;
	k = 0;
	do{ // for each haplotype
	  tmp = ((tab[i][cell] & screen));
	  tmp >>= position;
	  if(tmp>3)
	    {
	      cerr<<" ERROR "<<endl;
	      exit(3);
	    }
	  if( (int) tmp == haplotypes[k])
	    {  // matching haplotypes
	      matchHap = 1;
	      pi[k] += 1.0 / (double) N;
	    }
	  ++ k;
	}
	while((k < nbHap ) && (matchHap == 0));
	if (matchHap == 0)
	  {
	    tmp = (tab[i][cell] & screen);
	    haplotypes[nbHap] = (int) (tmp >> position);
	    pi.push_back(1.0 / (double) N);
	    ++ nbHap;
	  }
      }    
  }
  return pi;  
}

vector<int> DNA::site_frequency_int(int site) const
{
  int nbHap = 1;
  int k = 0; // loop on haplotypes
  int matchHap;
  int haplotypes[N];
  vector<int> pi; // frequency of the haplotypes
  if( nbXcell>0){
    int cell =  site / (CRUMB);
    int position = 2 *(site - CRUMB * cell);
    typeDNA screen;
    typeDNA tmp;
    screen = 3;
    screen <<= position;
    tmp = (tab[0][cell] & screen);
    haplotypes[0] = (int) (tmp >> position);
    pi.push_back(1);
    for(int i(1); i<N; ++i) // for each indiv
      {
	matchHap = 0;
	k = 0;
	do{ // for each haplotype
	  tmp = ((tab[i][cell] & screen));
	  tmp >>= position;
	  if(tmp>3)
	    {
	      cerr<<" ERROR "<<endl;
	      exit(3);
	    }
	  if( (int) tmp == haplotypes[k])
	    {  // matching haplotypes
	      matchHap = 1;
	      ++pi[k] ;
	    }
	  ++ k;
	}
	while((k < nbHap ) && (matchHap == 0));
	if (matchHap == 0)
	  {
	    tmp = (tab[i][cell] & screen);
	    haplotypes[nbHap] = (int) (tmp >> position);
	    pi.push_back(1 );
	    ++ nbHap;
	  }
      }    
  }
  return pi;  
}

double DNA::heterozigosity_hap() const
{
  double H = 0;
  vector<double> pi;
  pi = hap_frequency();
  int nb_hap = pi.size();
  if(nb_hap > 1)
    {
      for(int i(0); i <nb_hap; ++i)
	{
	  H += pi[i] * pi[i];
	}
      H = (double) N / (N - 1.0) * (1 - H);
    }
  else
    H = 0;
  return H;
}

double DNA::heterozigosity_sites() const
{
  double H = 0;
  double tmp;
  for(int k(0); k<nbSreal; k++)
    {
      tmp = 0;
      vector<double> pi;
      pi = site_frequency_DNA(k);
      int nb_hap = pi.size();
      if(nb_hap>1) // if 1 hap H=0
	{
	  for(int i(0); i <nb_hap; ++i)
	    {
	      tmp += pi[i] * pi[i];
	    }
	  H +=  (1 - tmp);
	}
    }
  if (H > 0)
    {
      H *= (double) N / (N - 1.0);
      H /= nbSreal;
    }
  return H;
}

double DNA::heterozigosity_sites2( std::vector<int> tsfs) const
{
  double H = 0;
  double tmp;
  for(int k(1); k<N; k++)
    {
      tmp = 0;
      double freq = tsfs[k];
      tmp += k*k + (N-k)*(N-k);
      H -= freq*tmp*1.0/(double) (N*N);
      H+= freq;
    }
  H *= (double) N / (N - 1.0);
  H /= nbSreal;  
  return (1-H);
}

double DNA::shannon_index() const
{
  double H = 0;
  double tmp;
  for(int k(0); k<nbSreal; k++)
    {
      tmp = 0;
      vector<double> pi;
      pi = site_frequency_DNA(k);
      int nb_hap = pi.size();
      if(nb_hap>1) // if 1 hap H=0
	{
	  for(int i(0); i <nb_hap; ++i)
	    {
	      tmp += pi[i]*log( pi[i]);
	    }
	  H += tmp;
	}
    }
  if (H > 0)
    {
      H *= (double) N / (N - 1.0);
      H /= nbSreal;
    }
  return H;
}

double DNA::Neff(double theta, double mu) const
{
  return theta / (2.0 * mu);
}

double DNA::theta_waterson() const
{
  double theta = harmonic(N-1);
  theta = (double) poly_sites() / theta ;
  return theta;
}

double DNA::theta_waterson(double polysites) const
{
  double theta = harmonic(N-1);
  theta = polysites / theta ;
  return theta;
}

double DNA::theta_hom() const
{
  double H = heterozigosity_hap();
  double theta = 1.0 / (double) (1.0 - H) -1;
  return theta;
}

double DNA::theta_pi() const
{
  double H = 0;
  vector<double> pi;
  DNA haplo;  
  haplo = haplotypes();
  pi = hap_frequency();
  int nb_hap = pi.size();
  for(int i(0); i <nb_hap; ++i)
    {
      for(int j(0); j <nb_hap; ++j)
	{
	  H += pi[i] * pi[j] * haplo.dist_two(i,j);
	}
    }  
  H = (double) N * H / ((double) (N - 1.0));
  return H;
}

double DNA::tajima_D(long ps, long pw, long tw) const // polymorphic sites, pairwise difference, theta_watterson
{
  long S = ps;
  double a1 = tajima_a1(N-1);
  double a2 = tajima_a2(N-1);
  double c1 = tajima_c1(N);
  double c2 = tajima_c2(N);
  return (pw - tw) / sqrt( c1 /a1 *S + c2 / (a1*a1 + a2) * S * (S -1)) ;
}

typeDNA DNA::dist_two(int a, int b) const
{
  typeDNA out = 0;
  typeDNA tmp,c;
  for(int i(0); i<nbXcell;  ++i)
    {
      c = tab[a][i]^tab[b][i];
      tmp = dna_distance(c);
      out = tmp + out;
    }
  return out;
}

typeDNA DNA::dist_two(int a, int b,DNA & dna2) const
{
  typeDNA out = 0;
  typeDNA tmp,c;
  for(int i(0); i<nbXcell;  ++i)
    {
      c = tab[a][i]^dna2.tab[b][i];
      tmp = dna_distance(c);
      out = tmp + out;
    }
  return out;
}

vector< vector<typeDNA> > input_fasta(string filename)
{ // need to change the size of chromosome if not CRUMB divided
  ifstream fs;
  string inputstring;
  int s;
  string name;
  string dna("");
  vector< vector<typeDNA> > dnatab;
  vector<typeDNA> dnacrumby;
  fs.open(filename.c_str());
  if(fs.is_open())
    {
      while(fs>>inputstring)
	{
	  if(strchr(inputstring.c_str(),'>')) // new indiv
	    {
	      name = inputstring;
	      if(dna.size()>0)
		{
		  // transform string to int
		  dnacrumby = dna_string_to_int(dna);
		  dnatab.push_back(dnacrumby);
		  dna="";
		}
	    }
	  else // dna ?
	    {
	      s = inputstring.size();
	      if (s>0) // it must be some dna bits
		{
		  dna += inputstring;
		}
	    }
	}
      dnatab.push_back(dnacrumby);
      return dnatab;
    }
  else
    {
      cerr << " ERROR : problem was encountered reading the input file"<<endl;
      exit(3);
    }
}

vector< vector<typeDNA> > input_ped(string filename)
{ // need to change the size of chromosome if not CRUMB divided
  ifstream fs;
  string inputstring;
  int s;
  string indiv;
  string dna1("");
  string dna2("");
  vector< vector<typeDNA> > dnatab;
  vector<typeDNA> dnacrumby;
  vector<string> tmp;
  bool test=true;
  fs.open(filename.c_str());
  if (fs.good())
    {
      string str;
      while(getline(fs, str)) 
	{
	  test=true;
	  istringstream ss(str);
	  while(ss >> indiv)
	    {	      
	      if((indiv=="A")||(indiv=="C")||(indiv=="T")||(indiv=="G"))
		{
		  if(test)
		    {
		      dna1+=indiv;
		      test=false;
		    }
		  else
		    {
		      dna2+=indiv;
		      test=true;
		    }
		}
	    }
	  dnacrumby = dna_string_to_int(dna1);
	  dnatab.push_back(dnacrumby);
	  dnacrumby = dna_string_to_int(dna2);
	  dnatab.push_back(dnacrumby);
	  dna1="";
	  dna2="";
	}
      fs.close();
      return dnatab;
    }
  else
    {
      cerr << " ERROR : problem was encountered reading the input file"<<endl;
      exit(3);
    }
}


vector< vector<typeDNA> > input_pedY(string filename)
{ // need to change the size of chromosome if not CRUMB divided
  ifstream fs;
  string inputstring;
  int s;
  string indiv;
  string dna1("");
  string dna2("");
  vector< vector<typeDNA> > dnatab;
  vector<typeDNA> dnacrumby;
  vector<string> tmp;
  bool test=true;
  fs.open(filename.c_str());
  if (fs.good())
    {
      string str;
      while(getline(fs, str)) 
	{
	  test=true;
	  istringstream ss(str);
	  while(ss >> indiv)
	    {	      
	      if((indiv=="A")||(indiv=="C")||(indiv=="T")||(indiv=="G"))
		{

		  if(test)
		    {
		      dna1+=indiv;
		      test=false;
		    }
		  else
		    { // ghost signal of ped
		      test=true;
		    }
		}
	    }
	  dnacrumby = dna_string_to_int(dna1);
	  dnatab.push_back(dnacrumby);
	  dna1="";
	}
      fs.close();
      return dnatab;
    }
  else
    {
      cerr << " ERROR : problem was encountered reading the input file"<<endl;
      exit(3);
    }
}


int input_ped_nbS(string filename)
{ // need to change the size of chromosome if not CRUMB divided
  ifstream fs;
  string inputstring;
  int s;
  string indiv;
  int size =0;
  string dna2("");
  vector< vector<typeDNA> > dnatab;
  vector<typeDNA> dnacrumby;
  vector<string> tmp;
  unsigned int test=0;
  fs.open(filename.c_str());
  if (fs.good())
    {
      string str;
      getline(fs, str);
      {
	istringstream ss(str);	 
	while(ss >> indiv)
	  {	      
	    if((indiv=="A")||(indiv=="C")||(indiv=="T")||(indiv=="G"))
	      {
		++size;
	      }
	  }
      }
    }
  return size/2;
}


int input_ped_nbSY(string filename)
{ // need to change the size of chromosome if not CRUMB divided
  ifstream fs;
  string inputstring;
  int s;
  string indiv;
  int size =0;
  string dna2("");
  vector< vector<typeDNA> > dnatab;
  vector<typeDNA> dnacrumby;
  vector<string> tmp;
  unsigned int test=0;
  fs.open(filename.c_str());
  if (fs.good())
    {
      string str;
      getline(fs, str);
      {
	istringstream ss(str);	 
	while(ss >> indiv)
	  {	      
	    if((indiv=="A")||(indiv=="C")||(indiv=="T")||(indiv=="G"))
	      {
		++size;
	      }
	  }
      }
    }
  return size/2;
}

std::vector<double> DNA::output_diversity_with_time(string filename, int Ntheo, int timeG)
{
  nbh = nb_of_haplotypes();
  hh = heterozigosity_hap();
  vector<int> tsfs = SFS_folded();
  pw = mean_pw_diff2(tsfs);
  ps = poly_sites2(tsfs);
  tw = theta_waterson(ps);
  hs = heterozigosity_sites2(tsfs);
  vector<double> mfs = MFS(tsfs);
  double eta = (double) tsfs[1] / (double) nbSreal; // frequency of single polymorphic site = nb of snp / nbsites
  double eta2 = (double) tsfs[2] / (double) nbSreal; // frequency of single polymorphic site = nb of snp / nbsites
  double eta3 = (double) tsfs[3] / (double) nbSreal; // frequency of single polymorphic site = nb of snp / nbsites
  //  double th = theta_hom();
  double td = tajima_D(ps,pw,tw);
  std::vector<double> out;
  out.push_back(nbSreal);
  out.push_back(ps);
  out.push_back(pw);
  out.push_back(nbh);
  out.push_back(hh);
  out.push_back(hs);
  out.push_back(tw);
  out.push_back(td);
  out.push_back(eta);
  //  out.push_back(eta2);
  //  out.push_back(eta3);
  out.push_back(mfs[0]);
  out.push_back(mfs[1]);
  return out;
}

void DNA::output_diversity_with_time(string filename, int Ntheo, int timeG,double pmate, double pmig) 
{
  //  int ps = poly_sites();
  //  double pps = poly_sites()/(double) nbS;
  //  double pw = mean_pw_diff();
   nbh = nb_of_haplotypes();
   hh = heterozigosity_hap();
  vector<int> tsfs = SFS_folded();
  double pw2 = mean_pw_diff2(tsfs);
  double ps2 = poly_sites2(tsfs);
   tw = theta_waterson(ps2);
  double hs2 = heterozigosity_sites2(tsfs);
  vector<double> mfs = MFS(tsfs);
  double eta = (double) tsfs[1] / (double) nbSreal; // frequency of single polymorphic site = nb of snp / nbsites
  double eta2 = (double) tsfs[2] / (double) nbSreal; // frequency of single polymorphic site = nb of snp / nbsites
  double eta3 = (double) tsfs[3] / (double) nbSreal; // frequency of single polymorphic site = nb of snp / nbsites
  string out;
  out += to_string(Ntheo) + "\t" +  to_string(nbSreal) + "\t" +  to_string(ps2) + "\t" +  to_string(pw2) + "\t" +  to_string(nbh) +  "\t" +  to_string(hh) + "\t"+   to_string(hs2) +  "\t"+  to_string(tw) + "\t" + to_string(eta) + "\t"+ to_string(eta2) + "\t" + to_string(eta3) + "\t"  + to_string(mfs[0]) + "\t"  + to_string(mfs[1]) + "\t"  +  to_string(timeG) +  "\t"  +  to_string(pmate) + "\t"  +  to_string(pmig) +"\n";
  ofstream fs;
  fs.open(filename.c_str(),ios::app);
  fs << out;
  //  fs << tsfs;
  fs.close();
}

void DNA::output_diversity_with_time_nosimu(string filename, int Ntheo) 
{
  nbh = nb_of_haplotypes();
  hh = heterozigosity_hap();
  vector<int> tsfs = SFS_folded();
  pw = mean_pw_diff2(tsfs);
  ps = poly_sites2(tsfs);
  tw = theta_waterson(ps);
  hs = heterozigosity_sites2(tsfs);
  vector<double> mfs = MFS(tsfs);
  double eta = (double) tsfs[1] / (double) nbSreal; // frequency of single polymorphic site = nb of snp / nbsites
  double eta2 = (double) tsfs[2] / (double) nbSreal; // frequency of single polymorphic site = nb of snp / nbsites
  double eta3 = (double) tsfs[3] / (double) nbSreal; // frequency of single polymorphic site = nb of snp / nbsites
  //  double th = theta_hom();
  //  double tp = mean_pw_diff();
  //  double td = tajima_D(ps,pw,tw);
  string out;
  //  out += to_string(Ntheo) + "\t" +  to_string(nbS) + "\t" +  to_string(ps) + "\t" +  to_string(pps) + "\t" +  to_string(pw) + "\t" +  to_string(nbh) +  "\t" +  to_string(hh) + "\t" +  to_string(hs) + "\t" +  to_string(tw) +  "\t" +  to_string(th) + "\t" +  to_string(tp) + "\t" +  to_string(td) + "\t" +  to_string(timeG) + "\n";
  out += to_string(Ntheo) + "\t" +  to_string(nbSreal) + "\t" +  to_string(ps) + "\t" +  to_string(pw) + "\t" +  to_string(nbh) +  "\t" +  to_string(hh) + "\t"+   to_string(hs) +  "\t"+  to_string(tw) + "\t" + to_string(eta) + "\t"+ to_string(eta2) + "\t" + to_string(eta3) + "\t"  + to_string(mfs[0]) + "\t"  + to_string(mfs[1])  + "\n";
  ofstream fs;
  fs.open(filename.c_str(),ios::app);
  fs << out;
  fs.close();
}


std::vector<double> DNA::output_diversity(DNA & dna2,string filename, int Ntheo, int timeG) 
{
  double pwb = mean_pw_diff_between(dna2);
  pw = mean_pw_diff();
  double Fst = (pwb-pw)/pwb;
  double neiD = Nei_D(dna2);
  double chordD = Chord_D(dna2);
  double rwcD = RWC_D(dna2);
  std::vector<double> out;
  out.push_back(pwb);
  out.push_back(Fst);
  out.push_back(neiD);
  out.push_back(chordD);
  out.push_back(rwcD);
  return out;
}

void DNA::output_diversity(string filename, int Ntheo) 
{
  ps = poly_sites();
  double pps = poly_sites()/(double) nbSreal;
  pw = mean_pw_diff();
  nbh = nb_of_haplotypes();
  hh = heterozigosity_hap();
  hs = heterozigosity_sites();
  tw = theta_waterson();
  double th = theta_hom();
  double tp = theta_pi();
  double td = tajima_D(ps,pw,tw);
  string out;
  out += to_string(Ntheo) + "\t" +  to_string(nbSreal) + "\t" +  to_string(ps) + "\t" +  to_string(pps) + "\t" +  to_string(pw) + "\t" +  to_string(nbh) +  "\t" +  to_string(hh) + "\t" +  to_string(hs) + "\t" +  to_string(tw) +  "\t" +  to_string(th) + "\t" +  to_string(tp) + "\t" +  to_string(td) +  "\n";
  ofstream fs;
  fs.open(filename.c_str(),ios::app);
  fs << out;
  fs.close();
}

void output_diversity_headline(string  filename)
{
  string out =  "N \t nbS \t poly_sites \t prop_poly_sites \t mean_pw_diff \t nb_of_haplotypes \t heterozygosity_hap \t heterozygosity_sites \t theta_waterson \t theta_hom \t theta_pi \t tajima_D \n";
  ofstream fs;
  fs.open(filename.c_str(),ios::app);
  fs << out;
  fs.close();
}

void output_diversity_headline_nosimu(string  filename)
{
  string out = "N \t nbS \t S \t mean_pw_diff \t nb_of_haplotypes \t h_hap \t h_sites \t theta_waterson \t eta \t eta2 \t eta3 \t msfs \t sdmsfs \n";
  ofstream fs;
  fs.open(filename.c_str(),ios::app);
  fs << out;
  fs.close();
}

void output_diversity_headline_with_time(string filename)
{
  string out = "N \t nbS \t S \t mean_pw_diff \t nb_of_haplotypes \t h_hap \t h_sites \t theta_waterson \t eta \t eta2 \t eta3 \t msfs \t sdmsfs \t time_in_generation \t pmate \t pmig  \n";
  ofstream fs;
  fs.open(filename.c_str(),ios::app);
  fs << out;
  fs.close();
}

void output_diversity_headline_with_time(string filename,string a,string b)
{

  string out = "N \t nbS \t S \t mean_pw_diff \t nb_of_haplotypes \t h_hap \t h_sites \t theta_waterson \t eta \t eta2 \t eta3 \t msfs \t sdmsfs \t time_in_generation \t"+a+"\t "+b+"\n";
  ofstream fs;
  fs.open(filename.c_str(),ios::app);
  fs << out;
  fs.close();
}

void output_diversity_headline_with_time2(string filename)
{
  string out = "  pw \t neiD \t chorD \t rwcD \t time  \t pmate \t pmig \n";
  ofstream fs;
  fs.open(filename.c_str(),ios::app);
  fs << out;
  fs.close();
}

void DNA::output_arlequin(string filename) const
{
  ofstream fs;
  fs.open((filename+".arp").c_str());
  fs << " # Arlequin file outputted by the program SMARTPOP \n";
  fs << "[Profile] \n";
  fs << "\t Title=\"Simulation done on " << time(0) << "\"\n";
  fs << "\t NbSamples=1 \n";
  fs << "\t DataType=DNA \n";
  fs << "\t GenotypicData=0 \n";
  fs << "\t LocusSeparator=NONE \n";
  fs << "\t MissingData=\"?\" \n \n";
  fs << "[Data]\n";
  fs << "\t [[Samples]]\n \n";
  fs << "\t \t SampleName=\"Sample1\" \n";
  fs << "\t \t SampleSize=" << N << "\n";
  fs << "\t \t SampleData ={\n";
  for(int i(0); i<N; ++i)
    {
      fs << i+1 << " 1 ";
      for(int j(0); j<nbXcell; ++j)
	fs << convertIntToQuaternaryString(tab[i][j]);
      fs << "\n";
    }
  fs << "\t \t }\n \n";
  fs << "[[Structure]] \n \n";
  fs << "\t StructureName=\"Simulated data\" \n";
  fs << "\t NbGroups=1 \n";
  fs << "\t Group={ \n";
  fs << "\t \"Sample1\" \n";
  fs << "\t } \n" ;
  fs.close();
}

void DNA::output_fasta(string filename) const
{
  ofstream fs;
  fs.open((filename+".fasta").c_str());
  for(int i(0); i<N; ++i)
    {
      fs << ">Indiv" << i+1 << " \n";
      for(int j(0); j<nbXcell; ++j)
	fs << convertIntToQuaternaryString(tab[i][j]) << "\n";
    }
  fs << "\n" ;
  fs.close();
}

void DNA::output_nexus(string filename) const
{
  ofstream fs;
  fs.open((filename+".nex").c_str());
  std::string out = "#NEXUS\n [Example of a Nexus file used in PopABC1.0]\n [Individuals defined here]\n BEGIN DATA;\n DIMENSIONS NTAX=" +to_string(N)+" NCHAR="+to_string(nbS)+";\n FORMAT DATATYPE=DNA;\n MATRIX\n";
  fs<<out;
  for(int i(0); i<N; ++i) 
    {
      fs << "Indiv" << i+1 << " ";
      for(int j(0); j<nbXcell; ++j)
	fs << convertIntToQuaternaryString(tab[i][j]) << "\n";
    }
  fs << "\n" ;
  out=";\n END;\n ";
  out+="[Populations defined here]\n BEGIN SETS;\n  TAXSET 'population1' = 1-"+to_string(N)+";\n END;\n";
  fs<<out;
  fs.close();
}

void DNA::output_haplotypes(string filename, int timeG)
{
  DNA H = haplotypes();
  vector<double> freq = hap_frequency();
  ofstream fs;
  fs.open(filename.c_str(),ios::app);
  fs<<"Generation "<<timeG<<" " << H.N<< "\n";
  for(int k(0);k<H.N;++k)
    {
      for(int i(0);i<H.nbXcell;++i)
	fs<<H.tab[k][i];
      fs << " " << freq[k] <<" \n";
    }
  fs << "\n";
  fs.close();
}



void DNA::SFS(string filename) const
{ 
  vector<int> tsfs = SFS();
  ofstream fs;
  fs.open(filename.c_str(),ios::app);
  for(unsigned i(0);i<tsfs.size();++i)
    fs << tsfs[i] << "\t";
  fs<<endl;
  fs.close();
}


void DNA::SFS_folded(string filename) const
{ 
  vector<int> tsfs = SFS_folded();
  ofstream fs;
  fs.open(filename.c_str(),ios::app);
  for(unsigned i(0);i<tsfs.size();++i)
    fs << tsfs[i] << "\t";
  fs<<endl;
  fs.close();
}


std::vector<int> DNA::SFS() const
{
  std::vector<int> sfs(N+1,0); // a vector of N cells filled with 0
  vector<int> pi;
  for(int i(0); i<nbSreal; ++i)
    {
      int freq = 0;
      pi = site_frequency_DNA_int(i);
      freq = N - pi[0];
      sfs[freq]+=1;
    }
  return sfs;
}

std::vector<int> DNA::SFS_folded() const
{ 
  std::vector<int> sfs(N+1,0); // a vector of N cells filled with 0
  vector<int> pi;
  for(int i(0); i<nbSreal; ++i)
    {
      int freq = N;
      pi = site_frequency_DNA_int(i);
      if(pi[0]>0)
	freq=pi[0];
      if(pi[1]>0)
	if(pi[1]<freq)
	  freq=pi[1];
      if(pi[2]>0)
	if(pi[2]<freq)
	  freq=pi[2];
      if(pi[3]>0)
	if(pi[3]<freq)
	  freq=pi[3];
      if((freq<0)||(freq>=N))
	freq=0;
      sfs[freq]+=1;
    }
  return sfs;
}


std::vector<double> DNA::SFS_folded_dbl(std::vector<int> tsfs) const
{
  std::vector<double> out(tsfs.size(),0);
  for(unsigned int i(0);i<tsfs.size();++i)
    {
      out[i] = (double) tsfs[i]/ (double) ( nbS);
    }
  return out;
}

double ssq_sfs(std::vector<double> tsfs1, std::vector<double> tsfs2)
{
  assert(tsfs1.size()==tsfs2.size());
  double out = 0;
  for(unsigned int i(0);i<tsfs1.size();++i)
    {
      if((tsfs1[i]>0)||( tsfs2[i]>0))
	{
	  double tmp;
	  tmp = tsfs1[i] - tsfs2[i];	  
	  out += tmp * tmp/(tsfs1[i] + tsfs2[i]);
	}
    }
  return out;
}

double sabs_sfs(std::vector<double> tsfs1,std::vector<double> tsfs2) 
{
  assert(tsfs1.size()==tsfs2.size());
  double out =0;
  for(unsigned int i(0);i<tsfs1.size();++i)
    {
      out +=  abs(tsfs1[i]- tsfs2[i]);
    }
  return out;
}

std::vector<double> resize_sfs(std::vector<double> tsfs)
{
  int nsize = std::floor(tsfs.size()/2)+1;
  std::vector<double> out(nsize,0);
  out[0]=tsfs[0]+tsfs[1]*0.5;
  int index = 0;
  for(unsigned int i(1);i<nsize-1;++i)
    {
      index += 2;
      out[i] = (tsfs[index-1]*0.5+ tsfs[index]+ tsfs[index+1]*0.5);
    }
  out[nsize-1] = tsfs[tsfs.size()-1]+tsfs[tsfs.size()-2]*0.5;
  return out;
}

std::vector<double> resize_simple_sfs(std::vector<double> tsfs)
{
  int nsize = std::floor(tsfs.size()/2)+1;
  std::vector<double> out(nsize,0);
  out[0]=tsfs[0];
  int index = 0;
  for(unsigned int i(1);i<nsize-1;++i)
    {
      index += 2;
      out[i] =  tsfs[index];
    }
  out[nsize-1] = tsfs[tsfs.size()-1];
  return out;
}
std::vector<double> DNA::MFS(std::vector<int> msfs) const // return mfs and sdmfs
{
  int n = msfs.size();
  double meansfs=0;
  double meansfs2 = 0;
  for(int i(1);i<n;++i)
    {
      double a = msfs[i];
      meansfs += a;
      meansfs2 += a*a; 
    }
  meansfs /= (double) n;
  meansfs2 /= (double) n;
  meansfs2 = sqrt(meansfs2 -meansfs*meansfs ); // square or mean - mean of square = std
  std::vector<double> out;
  out.push_back(meansfs);
  out.push_back(meansfs2);
  return out;
}


vector<double> DNA::site_frequency_DNA(int site) const 
{ // take into account DNA sequence - output frequency for base ACTG -> pi[0,1,2,3]
  //  int k = 0; // loop on haplotypes
  //  int matchHap;
  vector<double> pi(4,0); // frequency of the haplotypes
  int cell =  site / (CRUMB);
  int position = 2 *(site - CRUMB * cell);
  typeDNA screen;
  typeDNA tmp;
  screen = 3;
  screen <<= position;
  for(int i(0); i<N; ++i) // for each indiv
    {
      //      matchHap = 0;
      //      k = 0;
      //      do{ // for each haplotype
      tmp = ((tab[i][cell] & screen));
      tmp >>= position;
      //	if( (int) tmp == k)
      //	  {  // matching haplotypes
      //	    matchHap = 1;
      pi[tmp] += 1.0 / (double) N;
      //	  }
      //	++ k;
      //      }
      //      while((k < 4 ) && (matchHap == 0));
    }    
  return pi;
}

vector<double> DNA::site_frequency_DNA_crumb(int site) const // optimization ?
{ // take into account DNA sequence - output frequency for base ACTG -> pi[0,1,2,3]
  int k = 0; // loop on haplotypes
  int matchHap;
  vector<double> pi(4,0); // frequency of the haplotypes
  int cell =  site / (CRUMB);
  int position = 2 *(site - CRUMB * cell);
  typeDNA screen;
  typeDNA tmp;
  screen = 3;
  screen <<= position;
  //  typeDNA hap1 = tab[1][cell];
  for(int i(0); i<N; ++i) // for each indiv
    {
      matchHap = 0;
      k = 0;
      //      if((hap1)^tab[i][cell]!=0)
      do{ // for each haplotype
	tmp = ((tab[i][cell] & screen));
	tmp >>= position;
	if( (int) tmp == k)
	  {  // matching haplotypes
	    matchHap = 1;
	    pi[k] += 1.0 / (double) N;
	  }
	++ k;
      }
      while((k < 4 ) && (matchHap == 0));
    }    
  return pi;
}

vector<int> DNA::site_frequency_DNA_int(int site) const 
{ // take into account DNA sequence - output frequency for base ACTG -> pi[0,1,2,3]
  //  int k = 0; // loop on haplotypes
  //  int matchHap;
  vector<int> pi(4,0); // frequency of the haplotypes
  typeDNA cell =  site / (CRUMB);
  typeDNA position = 2 *(site - CRUMB * cell);
  typeDNA screen;
  typeDNA tmp;
  screen = 3;
  screen <<= position;
  for(int i(0); i<N; ++i) // for each indiv
    {
      tmp = ((tab[i][cell] & screen));
      tmp >>= position;
      pi[tmp] += 1.0 ;
    }
  return pi;
}

void DNA::SFS_DNA(string filename, int generation) const
{ // use site_frequency_DNA
  ofstream fs;
  fs.open(filename.c_str(),ios::app);
  vector<double> pi;
  for(int i(0); i<nbSreal; ++i)
    {
      pi = site_frequency_DNA(i);
      fs<<generation<< " " << i << " " ;
      for(int j(0); j<4; ++j)
	if (pi[j]<=1) // checl that is less than 1, sometime sums up to a little bit more, because of rounding
	  fs<< pi[j] << " " ;
	else
	  fs<< 1 << " " ;
      fs<<endl;
    }
  fs.close();
}

void DNA::AFS_DNA(string filename, int generation) const
{ // use site_frequency_DNA
  ofstream fs;
  fs.open(filename.c_str(),ios::app);
  vector<double> pi;
  pi = hap_frequency();
  for(int j(0); j<4; ++j)
    if(pi[j]<=1)
      fs<<generation<< " " << pi[j] << endl;
    else
      fs<<generation<< " " << 1 << endl;
  fs.close();
}
