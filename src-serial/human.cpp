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

#include "human.h"

using namespace std;



Human::Human()
{
  nbChild = 0;
  nbMate = 0;
  IDNetwork = -1;
  destinationPop = -1;
  if(DNAS.nbCrumbMt>0)
    {
      mtDna = new typeDNA[DNAS.nbCrumbMt]();
    }
  if(DNAS.nbCrumbNuc>0)
    {
      dna = new typeDNA[DNAS.nbCrumbNuc]();
    }
}

Human::Human(typeDNA *tabX,typeDNA *pMtDna)
{
  nbChild = 0;
  nbMate = 0;
  destinationPop = -1;
  IDNetwork = -1;
  if(DNAS.nbCrumbMt>0)
    {
      mtDna = new typeDNA[DNAS.nbCrumbMt]();
      memcpy(mtDna,pMtDna,DNAS.nbCrumbMt * sizeof(typeDNA));
    }
  if(DNAS.nbCrumbNuc>0)
    {
      dna = new typeDNA[DNAS.nbCrumbNuc]();
      memcpy(dna,tabX,DNAS.nbCrumbNuc * sizeof(typeDNA));
    }
 }

Human::Human(Human const& a)
{
  if(DNAS.nbCrumbNuc>0)
    {
      dna = new typeDNA[DNAS.nbCrumbNuc]();
      memcpy(dna,(a.dna),DNAS.nbCrumbNuc * sizeof(typeDNA));
    }
  if(DNAS.nbCrumbMt>0)
    {
      mtDna = new typeDNA[DNAS.nbCrumbMt]();
      memcpy(mtDna,(a.mtDna),DNAS.nbCrumbMt * sizeof(typeDNA));
    }
  nbChild = a.nbChild;
  nbMate = a.nbMate;
  IDNetwork = a.IDNetwork;
  destinationPop = a.destinationPop;
}

Human::~Human()
{  
   if(DNAS.nbCrumbNuc >0)
    {
      delete [] dna;
    }
  if(DNAS.nbCrumbMt >0)
    {
      delete [] mtDna;
    }
}

ostream & operator<<(ostream& os, Human const& a)
{
  a.print(os);
  return os;
}

Human& Human::operator = (Human const& a)
{
  if(this != &a)
    {
      //      delete from the old Human
      if(DNAS.nbCrumbNuc>0)
	{
	  if(dna)
	    delete [] dna;
	}
      if(DNAS.nbCrumbMt>0)
	if(mtDna)
	  delete [] mtDna;
      //we create an all brand new Human
      if(DNAS.nbCrumbNuc>0)
	{
	  //	  if(dna==NULL)
	  //	    {
	  typeDNA * tmp = new typeDNA[DNAS.nbCrumbNuc]();
	  dna = tmp;
	  //	    }	  
	  memcpy(dna,(a.dna),DNAS.nbCrumbNuc * sizeof(typeDNA));
	}
      if(DNAS.nbCrumbMt>0)
	{
	  //	  if(mtDna==NULL)
	  //	    {	  
	  typeDNA *tmp2 = new typeDNA[DNAS.nbCrumbMt]();
	  mtDna = tmp2;
	  //	    }
	  memcpy((mtDna),(a.mtDna),DNAS.nbCrumbMt * sizeof(typeDNA));
	}
      nbChild = a.nbChild;
      nbMate = a.nbMate;
      IDNetwork = a.IDNetwork;
      destinationPop =  a.destinationPop;
    }
  return *this;
}

ostream& Human::print(ostream& os) const
{
  string tmp2;
  typeDNA bob;
  //  const int crumb = CRUMB;
  os << " " <<  IDNetwork << " " << nbChild << " " << nbMate << " " ;
  for(int j(0);j<DNAS.nbCrumbNuc;++j)
    {
      bob = dna[j];	  
      if(bob>0)
	os << j << " " << bob << endl;
      //      tmp2 = convertIntToQuaternaryString(bob,CRUMB);
      //      os<<" "<<tmp2<<" *";
    }
  os<<endl;
  return os;
}


typeDNA Human::distance_dna(Human const& a)
{
  typeDNA distance(0);
  for(int i(0); i < DNAS.nbCrumbNuc; ++i)
    {
      distance += dist_dna(dna[i],a.dna[i]);
    }
  return distance;
}

void Human::reborn(typeDNA *newDna,typeDNA* newMtDna)
{ // faster to modify an existing Human than t create a new one
  if(TRACK_MATE_AND_CHILDREN)
    {
      nbChild = 0;
      nbMate = 0;
    }
  if(DNAS.nbCrumbMt>0)
    memcpy(mtDna,newMtDna,DNAS.nbCrumbMt * sizeof(typeDNA));
  if(DNAS.nbCrumbNuc>0)
    memcpy(dna,newDna,DNAS.nbCrumbNuc * sizeof(typeDNA));
  IDNetwork = -1;
  destinationPop = -1;
}


void Human::reborn(Human * & a)
{
  if(TRACK_MATE_AND_CHILDREN)
    {
      nbMate = 0;
      nbChild = 0;
    }
  if(DNAS.nbCrumbMt>0)
    memcpy(mtDna,(a->mtDna), DNAS.nbCrumbMt * sizeof(typeDNA));
  if(DNAS.nbCrumbNuc>0)
    memcpy((dna),(a->dna),DNAS.nbCrumbNuc * sizeof(typeDNA));
  IDNetwork = a->IDNetwork;
  destinationPop = a->destinationPop;
}

void Human::reset()
{
  if(TRACK_MATE_AND_CHILDREN)
    {      
      nbMate = 0;
      nbChild = 0;
    }
  if(DNAS.nbCrumbMt>0)
    std::fill(mtDna,mtDna + DNAS.nbCrumbMt, 0);
  if(DNAS.nbCrumbNuc>0)
    std::fill(dna,dna + DNAS.nbCrumbNuc, 0);
  IDNetwork = -1; // so that it bugs if IDNetwork is not well set
  destinationPop = -1;
}


int Human::getNbMate() const
{
  return nbMate;
}

int Human::getNbChild() const
{
  return nbChild;
}

void Human::hadChild(int n)
{
  nbChild += n;
}

void Human::newMate()
{
  ++nbMate;
}
