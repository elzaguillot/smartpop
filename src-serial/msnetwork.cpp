
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
#include "msnetwork.h"

using namespace boost;

#define NSTG 4

// +2 => +2*nbPop
// fatherNode / meotherNode changes
// delete generation ?

MSNetwork::MSNetwork():myg(2)
{ // dull constructor
  nbNodes = 2;
  storedGeneration = 0;
  nnpg=0;
  nbPop = 0;
}

MSNetwork::MSNetwork(int a):myg(NSTG*a+2)
{ // constructor from population size
  nbNodes = NSTG*a+2;
  storedGeneration = NSTG;
  nnpg = a;
  nbPop = 1;
  popsizes.push_back(a);
}

MSNetwork::MSNetwork(int a, int b):myg(b*a+2)
{ // constructor from population size and the number of stored generations
  nbNodes = b*a+2;
  storedGeneration = b;
  nnpg = a;
  nbPop = 1;
  popsizes.push_back(a);
}

MSNetwork::MSNetwork( MSNetwork const & a):myg(a.myg)
{ // constructor of copy
  nbNodes = a.nbNodes;
  storedGeneration = a.storedGeneration;
  nnpg = a.nnpg;
  popsizes = a.popsizes;
  nbPop = a.nbPop;
}

MSNetwork::~MSNetwork()
{
}

MSNetwork& MSNetwork::operator = (MSNetwork const& a)
{
  if(this != &a)
    {
      myg = a.myg;
      nbNodes = a.nbNodes;
      storedGeneration = a.storedGeneration;
      nnpg = a.nnpg;
      popsizes = a.popsizes;
      nbPop = a.nbPop;
    }
  return *this;
}

void MSNetwork::Init(std::vector<int> nbFemales, std::vector<int> popsizes,int generationZero) // return the index of females
{  // Initiate the Class depending on the number of female
  int nbAssignedIndiv = 0;
  int offset = generationZero*nnpg;
  assert(nbNodes==num_vertices(myg));
  assert(nbNodes==nnpg*storedGeneration+nbPop*2);
  for(int i(0); i<offset; ++i)
    {
      //      myg[vertex(i,myg)].ID = -1;
      myg[i].ID = -1;
      //      myg[vertex(i,myg)].generation = 0;
      //      myg[vertex(i,myg)].sex = false;
      //      myg[vertex(i,myg)].population = -1;
      myg[i].population = -1;
      //      std::cout << i << " " ;
    }
  //  std::cout << std::endl;
  for(unsigned i(offset+nnpg); i<nbNodes; ++i)
    {
      myg[i].ID = -1;
      myg[i].population = -1;
      //      myg[vertex(i,myg)].ID = -1;
      //      myg[vertex(i,myg)].generation = 0;
      //      myg[vertex(i,myg)].sex = false;
      //      myg[vertex(i,myg)].population = -1;
      //      std::cout << i << ". " ;
    }
  //  std::cout <<std::endl;
  nbAssignedIndiv=offset;
  for(int k(0);k<nbPop;++k)
    {
      int nbFemale = nbFemales[k];
      for(int i(0); i<nbFemale; ++i)
	{
	  //	  std:: cout << nbAssignedIndiv << "-" << k << " ";
	  /*	  myg[vertex(nbAssignedIndiv,myg)].sex = true;
	  myg[vertex(nbAssignedIndiv,myg)].ID = i;
	  myg[vertex(nbAssignedIndiv,myg)].generation = generationZero;
	  myg[vertex(nbAssignedIndiv,myg)].population = k;*/
	  myg[nbAssignedIndiv].sex = true;
	  myg[nbAssignedIndiv].ID = i;
	  myg[nbAssignedIndiv].generation = generationZero;
	  myg[nbAssignedIndiv].population = k;
	  ++nbAssignedIndiv;
	}
      //      std::cout << std::endl;
      int popsize = popsizes[k];
      for(int i(nbFemale); i<popsize; ++i)
	{
	  //	  std:: cout << nbAssignedIndiv << "-" << k << " ";
	  myg[vertex(nbAssignedIndiv,myg)].sex = false;
	  myg[vertex(nbAssignedIndiv,myg)].ID = i-nbFemale;
	  myg[vertex(nbAssignedIndiv,myg)].generation = generationZero;
	  myg[vertex(nbAssignedIndiv,myg)].population = k;
	  ++nbAssignedIndiv;
	}
      //      std::cout << std::endl;
      //      std::cout <<  nbNodes-1-(2*k) << std::endl;
      myg[vertex(nbNodes-1-(2*k),myg)].ID=-2;
      myg[vertex(nbNodes-1-(2*k),myg)].population = k;
      myg[vertex(nbNodes-2-(2*k),myg)].ID=-2;
      myg[vertex(nbNodes-2-(2*k),myg)].population = k;
      //      assert(nbAssignedIndiv == nnpg); // watchout
    }
  for(int j(nbAssignedIndiv);j<(nnpg+offset);++j)
    {
      //      std:: cout << j << "-";
      myg[vertex(j,myg)].sex = false;
      myg[vertex(j,myg)].ID = -1;
      myg[vertex(j,myg)].generation = -1;
      myg[vertex(j,myg)].population = -1;
    }
  //  std::cout << std::endl;
}

void MSNetwork::Init2() // return the index of females
{  // Initiate the Class depending on the number of female
  assert(nbNodes==num_vertices(myg));
  assert(nbNodes==nnpg*storedGeneration+nbPop*2);
  for(int i(0); i<nbNodes-2*nbPop; ++i)
    {
      /*	  myg[vertex(i,myg)].sex = false;
	  myg[vertex(i,myg)].ID = -1;
	  myg[vertex(i,myg)].generation = -1;
	  myg[vertex(i,myg)].population = -1;*/
      myg[i].sex = false;
      myg[i].ID = -1;
      myg[i].generation = -1;
      myg[i].population = -1;
    }
  for(int k(0);k<nbPop;++k)
    {
      /*      myg[vertex(nbNodes-1-(2*k),myg)].ID=-2;
      myg[vertex(nbNodes-1-(2*k),myg)].population = k;
      myg[vertex(nbNodes-2-(2*k),myg)].ID=-2;
      myg[vertex(nbNodes-2-(2*k),myg)].population = k;*/
      myg[nbNodes-1-(2*k)].ID=-2;
      myg[nbNodes-1-(2*k)].population = k;
      myg[nbNodes-2-(2*k)].ID=-2;
      myg[nbNodes-2-(2*k)].population = k;
      //      assert(nbAssignedIndiv == nnpg); // watchout
    }
}

/*void MSNetwork::resize(int population, int newsize)
{ // resize the graph as a function of the new population size
  int nbNewNodes = newsize - popsizes[population]*nnpg;
  if(nbNewNodes>0)
    {
      for(int i(0);i<nbNewNodes * storedGeneration;++i)
	{
	  add_vertex(myg);
	}
      nnpg += nbNewNodes;
      nbNodes = nnpg * storedGeneration + 2*nbPop;
      Init(0.5,0);
    }
  else
    {
      std::cerr << " ERROR : problem in the resizing of the myg" << std::endl;
      exit(2);
    }
    }*/

void MSNetwork::reset(std::vector<int> nbFemales, std::vector<int> popsizes)
{ // reset the whole classe instead of creating a new one if same number of nodes
  // remove all the edges
  GraphTraits::vertex_iterator v,vend;
  /*  for (tie(v,vend) = vertices(myg);v!=vend;++v)
    {
      clear_vertex(*v,myg);
      }*/
  for (int i(0);i<nbNodes;++i)
    {
      clear_vertex(vertex(i,myg),myg);
    }

  Init(nbFemales,popsizes,0);
}

void MSNetwork::clear_premating(int gparent)
{
  int indexc = gparent * nnpg;
  std::vector<GraphTraits::vertex_descriptor> fatherNode;// = vertex(nbNodes-1-2*l,myg);
  std::vector<GraphTraits::vertex_descriptor> motherNode;// = vertex(nbNodes-2-2*l,myg);
  GraphTraits::vertex_descriptor v;
  for(int l(0);l<nbPop;++l)
    {
      fatherNode.push_back(vertex(nbNodes-1-2*l,myg));
      motherNode.push_back(vertex(nbNodes-2-2*l,myg));
    }
  for(int k(0); k<nnpg; ++k)
    { // for all nodes in this generation, clear the out edges and connect the node to the fatherNOde and motherNode
      v = vertex(indexc,myg);
      //      if(
      clear_out_edges(v,myg);
      if(myg[v].population>-1)
	{
	  if(myg[v].sex == true)
	    {
	      add_edge(v,motherNode[myg[v].population],myg);
	      //	      std::cout << myg[v].population << " - " << myg[v].ID << " " << myg[v].sex << " " << k << std::endl;
	    }
	  else
	    add_edge(v,fatherNode[myg[v].population],myg);
	}
      ++indexc;
    }
  int nexG = (gparent+1)%storedGeneration;
  indexc = nexG * nnpg;
  for(int k(0); k<nnpg; ++k)
    { // for all nodes in this generation, clear the out edges and connect the node to the fatherNOde and motherNode
      v = vertex(indexc,myg);
      myg[v].population = -1;
      myg[v].ID = -1;
      ++indexc;
    }
}

void MSNetwork::clear_postmating()
{
  GraphTraits::vertex_descriptor fatherNode;
  GraphTraits::vertex_descriptor motherNode;
  for(int l(0);l<nbPop;++l)
    {
      fatherNode = vertex(nbNodes-1-2*l,myg);
      motherNode = vertex(nbNodes-2-2*l,myg);
      clear_vertex(fatherNode,myg);
      clear_vertex(motherNode,myg);
    }
}

void MSNetwork::set_nbPop(int a,std::vector<int> maxpopsize,std::vector<int> popsizes, std::vector<int> nbFemales)
{
  int diffNbPop = a-nbPop;
  int newnnpg = 0;
  for(int i(0);i<a;++i)
    newnnpg += maxpopsize[i];
  int diffSize = (newnnpg-nnpg)*storedGeneration;
  for(int i(0);i<diffSize;++i)
    {
      //	  nbNodes +=2;
      Node hNode;
      //	  hNode.ID = -1;
      //	  hNode.population = -1;
      //	  hNode.generation = -1;
      add_vertex(hNode,myg);
    }
  if(newnnpg>nnpg)
    {
      nnpg = newnnpg;
      nbNodes += diffSize;
    }
  if(diffNbPop >0)
    {      
      nbPop = a;
      popsizes = maxpopsize;
      nbNodes = nnpg*storedGeneration+nbPop*2;
      for(int i(0); i< diffNbPop;++i)
	{
	  Node newFatherNode;
	  //newFatherNode.ID = -nbPop+(diffNbPop-i);
	  /*	  newFatherNode.ID = nbPop-(diffNbPop-i);
	  newFatherNode.sex = false;
	  newFatherNode.generation = -2;*/
	  add_vertex(newFatherNode,myg);
	  Node newMotherNode;
	  //newMotherNode.ID = -nbPop+(diffNbPop-i);
	  /*newMotherNode.ID = nbPop-(diffNbPop-i);
	  newMotherNode.sex = true;
	  newMotherNode.generation = -2;*/
	  add_vertex(newMotherNode,myg);
	}
      Init(nbFemales,popsizes,0);
    }
  assert(num_vertices(myg)==nbNodes);
  // if the number of population is less, do nothing, at least for now
}

void MSNetwork::mateOne(int generation,int propsize,int popDemo, int nbFemale )
{}

void MSNetwork::writeGraphviz(std::string filename,int generation)
{ // write the graph in a .dot file
  // color represent the sex
  // the label represents the generation
  //  clear_premating(storedGeneration-1);
  clear_postmating();
  GraphTraits::vertex_iterator v,vend;
  int k = 0;
  for( tie(v,vend) = vertices(myg);v!=vend;++v)
    {
      myg[*v].ID= ++k;
      if(myg[*v].sex == 1)
	myg[*v].color = "#ff0000"; // red
      else
	myg[*v].color = "#0000ff"; // blue
      int a = myg[*v].generation;
      int b = myg[*v].population;
      int c = myg[*v].ID;
      //      name = to_string(a)+"-"+to_string(k)+" "+to_string(b)+" "+to_string(c);
      std::stringstream ss;
      ss << a<<"-"<<k<<"-"<< b<<"-"<< c;
      myg[*v].name= ss.str();
    }
  boost::dynamic_properties dp;
  dp.property("color", get(&Node::color, myg));
  dp.property("label", get(&Node::name, myg));
  dp.property("node_id", get(&Node::ID, myg)); 
  std::ofstream fout(filename.c_str(),std::ios::trunc);
  // Write out the graph
  //  boost::write_graphviz_dp(std::cout, myg, dp);
  boost::write_graphviz_dp(fout, myg, dp);
  fout.close();
}

int MSNetwork::get_indivID(int population, int generation, int index)
{
  int gfactor = generation%storedGeneration;
  int out =gfactor*nnpg;
  for(int i(0);i<nbPop;++i)
    out += popsizes[i];
  out += index;
  return out;
}

bool MSNetwork::share_one_parent(int ID1, int ID2)
{
  std::vector<GraphTraits::vertex_descriptor> p1 =  get_parents(ID1);
  std::vector<GraphTraits::vertex_descriptor> p2 = get_parents(ID2);
  for(unsigned i(0); i <p1.size();++i)
    for(unsigned j(0); j <p2.size();++j)
      if (p1[i]==p2[j])
	return true;
  return false;	  
}

bool MSNetwork::share_two_parents(int ID1, int ID2)
{
  std::vector<GraphTraits::vertex_descriptor> p1 = get_parents(ID1);
  std::vector<GraphTraits::vertex_descriptor> p2 = get_parents(ID2);
  //  for(int i(0); i <p1.size();++i)
  if(p1.size()==2)
    {
      for(unsigned j(0); j <p2.size();++j)
	if (p1[0]==p2[j])
	  for(unsigned j(0); j <p2.size();++j)
	    if (p1[1]==p2[j])
	      return true;
    }
  return false;
}

GraphTraits::vertex_descriptor MSNetwork::get_indivNode(int population, int generation, int index)
{
  int out = get_indivID(population,generation,index);
  return vertex(out,myg);
}

std::vector<int> MSNetwork::get_femalesID(int generation, int popID)
{
  int offset = (generation%storedGeneration) * nnpg;
  GraphTraits::vertex_descriptor v;
  std::vector<int> femalesID;
  for(int i(0); i< nnpg; ++i)
    {
      v = vertex(i+offset,myg);
      if(myg[v].population == popID)
	if(myg[v].sex == true)
	  femalesID.push_back(i+offset);
    }
  return femalesID;
}

std::vector<int> MSNetwork::get_femalesID(int generation)
{
  int offset = (generation%storedGeneration) * nnpg;
  GraphTraits::vertex_descriptor v;
  std::vector<int> femalesID;
  for(int i(0); i< nnpg; ++i)
    {
      v = vertex(i+offset,myg);
      if(myg[v].sex == true)
	femalesID.push_back(i+offset);
    }
  return femalesID;
}

std::vector<int> MSNetwork::get_malesID(int generation, int popID)
{
  int offset = (generation%storedGeneration) * nnpg;
  GraphTraits::vertex_descriptor v;
  std::vector<int> malesID;
  for(int i(0); i< nnpg; ++i)
    {
      v = vertex(i+offset,myg);
      if(myg[v].population == popID)
	if(myg[v].sex == false)
	  malesID.push_back(i+offset);
    }
  return malesID;
}

std::vector<int> MSNetwork::get_malesID(int generation)
{
  int offset = (generation%storedGeneration) * nnpg;
  GraphTraits::vertex_descriptor v;
  std::vector<int> malesID;
  for(int i(0); i< nnpg; ++i)
    {
      v = vertex(i+offset,myg);
      if(myg[v].sex == false)
	malesID.push_back(i+offset);
    }
  return malesID;
}

GraphTraits::vertex_descriptor MSNetwork::get_father(int index)
{
  Tinv_adjacency_iterator vup;
  Tinv_adjacency_iterator vup2;
  for ( tie(vup,vup2) = inv_adjacent_vertices(vertex(index,myg),myg); vup != vup2; ++vup)
    {
      if (myg[*vup].sex == 0)
	return *vup;
    }
  std::cerr << " ERROR : access to a non existing generation" << std::endl;
  return NULL;// seg fault
}

GraphTraits::vertex_descriptor MSNetwork::get_mother(int index)
{
  return get_mother(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor > MSNetwork::get_parents(int index)
{
  return get_parents(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor > MSNetwork::get_children(int index)
{
  return get_children(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor > MSNetwork::get_sons(int index)
{
  return get_sons(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_daughters(int index)
{
  return get_daughters(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_brothers(int index)
{
  return get_brothers(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_half_brothers(int index)
{
  return get_half_brothers(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_mother_sons(int index)
{
  return get_mother_sons(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_father_sons(int index)
{
  return get_father_sons(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_sisters(int index)
{
  return get_sisters(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_mother_daughters(int index)
{
  return get_mother_daughters(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_father_daughters(int index)
{
  return get_father_daughters(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_grandparents(int index)
{
  return get_grandparents(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_mother_parents(int index)
{
  return get_mother_parents(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_father_parents(int index)
{
  return get_father_parents(vertex(index,myg));
}

GraphTraits::vertex_descriptor MSNetwork::get_mother_mother(int index)
{
  return get_mother_mother(vertex(index,myg));
}

GraphTraits::vertex_descriptor MSNetwork::get_mother_father(int index)
{
  return get_mother_father(vertex(index,myg));
}

GraphTraits::vertex_descriptor MSNetwork::get_father_mother(int index)
{
  return get_father_mother(vertex(index,myg));
}

GraphTraits::vertex_descriptor MSNetwork::get_father_father(int index)
{
  return get_father_father(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_mother_brothers(int index)
{
  return get_mother_brothers(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_father_sisters(int index)
{
  return get_father_sisters(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_grandchildren(int index)
{
  return get_grandchildren(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_grandsons(int index)
{
  return get_grandsons(vertex(index,myg));
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_granddaughters(int index)
{
  return get_granddaughters(vertex(index,myg));
}

GraphTraits::vertex_descriptor MSNetwork::get_father(GraphTraits::vertex_descriptor vdes)
{
  Tinv_adjacency_iterator vup;
  Tinv_adjacency_iterator vup2;
  for( tie(vup,vup2) = inv_adjacent_vertices(vdes,myg); vup != vup2; ++vup)
    {
      if (myg[*vup].sex == 0)
	return *vup;
    }
  std::cerr << "ERROR: can not find the father" << std::endl;
  return GraphTraits::vertex_descriptor(); // default constructor?
}

GraphTraits::vertex_descriptor MSNetwork::get_mother(GraphTraits::vertex_descriptor vdes)
{
  Tinv_adjacency_iterator vup;
  Tinv_adjacency_iterator vup2;
  for  ( tie(vup,vup2) = inv_adjacent_vertices(vdes,myg); vup != vup2; ++vup)
    {
      if (myg[*vup].sex == 1)
	return *vup;
    }
  std::cerr << "ERROR: can not find the mother" << std::endl;
  return GraphTraits::vertex_descriptor(); // default constructor?
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_parents(GraphTraits::vertex_descriptor vdes)
{
  Tinv_adjacency_iterator vup;
  Tinv_adjacency_iterator vup2;
  std::vector<GraphTraits::vertex_descriptor> out;
  for(tie(vup,vup2) = inv_adjacent_vertices(vdes,myg); vup!=vup2; ++vup)
    {
      out.push_back(*vup);
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_children(GraphTraits::vertex_descriptor vdes)
{
  GraphTraits::adjacency_iterator vdown;
  GraphTraits::adjacency_iterator vdown2;
  std::vector<GraphTraits::vertex_descriptor> out;
  for(tie(vdown,vdown2) = adjacent_vertices(vdes,myg); vdown!=vdown2; ++vdown)
    {
      if (myg[*vdown].ID >=0) // take of mother / fatherNode
	out.push_back(*vdown);
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_sons(GraphTraits::vertex_descriptor vdes)
{
  GraphTraits::adjacency_iterator vdown;
  GraphTraits::adjacency_iterator vdown2;
  std::vector<GraphTraits::vertex_descriptor> out;
  for(tie(vdown,vdown2) = adjacent_vertices(vdes,myg); vdown!=vdown2; ++vdown)
    {
      if((myg[*vdown].sex==0)&&(myg[*vdown].ID>=0))
	out.push_back(*vdown);
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_daughters(GraphTraits::vertex_descriptor vdes)
{
  GraphTraits::adjacency_iterator vdown;
  GraphTraits::adjacency_iterator vdown2;
  std::vector<GraphTraits::vertex_descriptor> out;
  for(tie(vdown,vdown2) = adjacent_vertices(vdes,myg); vdown!=vdown2; ++vdown)
    {
      if((myg[*vdown].sex==1)&&(myg[*vdown].ID>=0))
	out.push_back(*vdown);
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_half_brothers(GraphTraits::vertex_descriptor vdes)
{
  Tinv_adjacency_iterator vup;
  Tinv_adjacency_iterator vup2;
  std::vector<GraphTraits::vertex_descriptor> out;
  std::vector<GraphTraits::vertex_descriptor> sons;
  for( tie(vup,vup2) = inv_adjacent_vertices(vdes,myg); vup != vup2; ++vup)
    {
      sons = get_sons((*vup));
      out = merge(out,sons);
    }
  return out;
}


std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_brothers(GraphTraits::vertex_descriptor vdes)
{
  std::vector<GraphTraits::vertex_descriptor> out;
  std::vector<GraphTraits::vertex_descriptor> v1;
  std::vector<GraphTraits::vertex_descriptor> v2;
  std::vector<GraphTraits::vertex_descriptor> matribro = get_matribrothers(vdes);
  std::vector<GraphTraits::vertex_descriptor> patribro = get_patribrothers(vdes);
  for(std::vector<GraphTraits::vertex_descriptor>::iterator itm = matribro.begin(); itm != matribro.end(); ++itm)
    {
      for(std::vector<GraphTraits::vertex_descriptor>::iterator itp = matribro.begin(); itp != matribro.end(); ++itp)
	{
	  if(*itm==*itp)
	    out.push_back(*itm);
	}
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_sisters(GraphTraits::vertex_descriptor vdes)
{
  std::vector<GraphTraits::vertex_descriptor> out;
  std::vector<GraphTraits::vertex_descriptor> v1;
  std::vector<GraphTraits::vertex_descriptor> v2;
  std::vector<GraphTraits::vertex_descriptor> matrisis = get_matrisisters(vdes);
  std::vector<GraphTraits::vertex_descriptor> patrisis = get_patrisisters(vdes);
  for(std::vector<GraphTraits::vertex_descriptor>::iterator itm = matrisis.begin(); itm != matrisis.end(); ++itm)
    {
      for(std::vector<GraphTraits::vertex_descriptor>::iterator itp = matrisis.begin(); itp != matrisis.end(); ++itp)
	{
	  if(*itm==*itp)
	    out.push_back(*itm);
	}
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_matribrothers(GraphTraits::vertex_descriptor vdes)
{ // get brother with which vdes shares a mother
  std::vector<GraphTraits::vertex_descriptor> out;
  if(in_degree(vdes,myg)>0)
    out = get_sons(get_mother(vdes));
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_patribrothers(GraphTraits::vertex_descriptor vdes)
{ // get brother with which vdes share a father
  std::vector<GraphTraits::vertex_descriptor> out;
  if(in_degree(vdes,myg))
    {
      out = get_sons(get_father(vdes));
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_mother_sons(GraphTraits::vertex_descriptor vdes)
{ // same as matriborthers
  std::vector<GraphTraits::vertex_descriptor> out;
  GraphTraits::vertex_descriptor mum;
  if(in_degree(vdes,myg)>0)
    {
      mum = get_mother(vdes);
      out = get_sons(mum);
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_father_sons(GraphTraits::vertex_descriptor vdes)
{ // same as patribrothers
  std::vector<GraphTraits::vertex_descriptor> out;
  GraphTraits::vertex_descriptor dad;
  if(in_degree(vdes,myg)>0)
    {
      dad = get_father(vdes);
      out = get_sons(dad);
    }
  return out;
}


std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_half_sisters(GraphTraits::vertex_descriptor vdes)
{
  Tinv_adjacency_iterator vup;
  Tinv_adjacency_iterator vup2;
  std::vector<GraphTraits::vertex_descriptor> out;
  std::vector<GraphTraits::vertex_descriptor> daughters;
  for( tie(vup,vup2) = inv_adjacent_vertices(vdes,myg); vup != vup2; ++vup)
    { // from both parents
      daughters = get_daughters((*vup));
      out = merge(out,daughters);
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_matrisisters(GraphTraits::vertex_descriptor vdes)
{
  if(in_degree(vdes,myg)>0)
    return get_daughters(get_mother(vdes));
  std::vector<GraphTraits::vertex_descriptor> out;
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_patrisisters(GraphTraits::vertex_descriptor vdes)
{
  if(in_degree(vdes,myg)>0)
    return get_daughters(get_father(vdes));
  std::vector<GraphTraits::vertex_descriptor> out;
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_mother_daughters(GraphTraits::vertex_descriptor vdes)
{
  std::vector<GraphTraits::vertex_descriptor> out;
  GraphTraits::vertex_descriptor mum;
  if(in_degree(vdes,myg)>0)
    {  
      mum = get_mother(vdes);
      out = get_daughters(mum);
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_father_daughters(GraphTraits::vertex_descriptor vdes)
{
  std::vector<GraphTraits::vertex_descriptor> out;
  GraphTraits::vertex_descriptor dad;
  if(in_degree(vdes,myg)>0)
    {
      dad = get_father(vdes);
      out = get_daughters(dad);
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_mother_parents(GraphTraits::vertex_descriptor vdes)
{
  std::vector<GraphTraits::vertex_descriptor> out;
  GraphTraits::vertex_descriptor mum;
  if(in_degree(vdes,myg)>0)
    {
      mum = get_mother(vdes);
      if(in_degree(mum,myg)>0)
	out = get_parents(mum);
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_father_parents(GraphTraits::vertex_descriptor vdes)
{
  std::vector<GraphTraits::vertex_descriptor> out;
  GraphTraits::vertex_descriptor dad;
  if(in_degree(vdes,myg)>0)
    {
      dad = get_father(vdes);
      if(in_degree(dad,myg)>0)
	out = get_parents(dad);
    }
  return out;
}

GraphTraits::vertex_descriptor MSNetwork::get_mother_mother(GraphTraits::vertex_descriptor vdes)
{
  GraphTraits::vertex_descriptor out;
  GraphTraits::vertex_descriptor mum;
  if(in_degree(vdes,myg)>0)
    {
      mum = get_mother(vdes);
      if(in_degree(mum,myg)>0)
	out = get_mother(mum);
    }
  return out;
}

GraphTraits::vertex_descriptor MSNetwork::get_mother_father(GraphTraits::vertex_descriptor vdes)
{
  GraphTraits::vertex_descriptor out;
  GraphTraits::vertex_descriptor mum;
  if(in_degree(vdes,myg)>0)
    {
      mum = get_mother(vdes);
      if(in_degree(mum,myg)>0)
	{
	  out = get_father(mum);
	}
    }
  return out;
}

  
GraphTraits::vertex_descriptor MSNetwork::get_father_mother(GraphTraits::vertex_descriptor vdes)
{
  GraphTraits::vertex_descriptor out;
  GraphTraits::vertex_descriptor dad;
  if(in_degree(vdes,myg)>0)
    {
      dad = get_father(vdes);
      if(in_degree(dad,myg)>0)
	out = get_mother(dad);
    }
  return out;
}

GraphTraits::vertex_descriptor MSNetwork::get_father_father(GraphTraits::vertex_descriptor vdes)
{
  GraphTraits::vertex_descriptor out;
  GraphTraits::vertex_descriptor dad;
  if(in_degree(vdes,myg)>0)
    {
      dad = get_father(vdes);
      if(in_degree(dad,myg)>0)
	out = get_father(dad);
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_mother_brothers(GraphTraits::vertex_descriptor vdes)
{
  std::vector<GraphTraits::vertex_descriptor> out;
  GraphTraits::vertex_descriptor mum;
  if(in_degree(vdes,myg)>0)
    {
      mum = get_mother(vdes);
      out = get_brothers(mum);
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_father_sisters(GraphTraits::vertex_descriptor vdes)
{
  std::vector<GraphTraits::vertex_descriptor> out;
  GraphTraits::vertex_descriptor dad;
  if(in_degree(vdes,myg)>1)
    {
      dad = get_father(vdes);
      out = get_sisters(dad);
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_father_sisters_sons(GraphTraits::vertex_descriptor vdes)
{
  std::vector<GraphTraits::vertex_descriptor> sons;
  std::vector<GraphTraits::vertex_descriptor> sisters;
  std::vector<GraphTraits::vertex_descriptor> out;
  GraphTraits::vertex_descriptor dad;
  if(in_degree(vdes,myg)>1)
    {
      dad = get_father(vdes);
      sisters = get_half_sisters(dad);
      for(unsigned int i(0);i<sisters.size();++i)
	{
	  sons = get_sons(sisters[i]);
	  out += sons;
	}
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_father_brothers_sons(GraphTraits::vertex_descriptor vdes)
{
  std::vector<GraphTraits::vertex_descriptor> sons;
  std::vector<GraphTraits::vertex_descriptor> brothers;
  std::vector<GraphTraits::vertex_descriptor> out;
  GraphTraits::vertex_descriptor dad;
  if(in_degree(vdes,myg)>0)
    {
      dad = get_father(vdes);
      brothers = get_half_brothers(dad);
      for(unsigned int i(0);i<brothers.size();++i)
	{
	  sons = get_sons(brothers[i]);
	  out += sons;
	}
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_father_sisters_sons_nob_nofbs_nomzs(GraphTraits::vertex_descriptor vdes)
{
  std::vector<GraphTraits::vertex_descriptor> sons;
  std::vector<GraphTraits::vertex_descriptor> sisters;
  std::vector<GraphTraits::vertex_descriptor> out;
  std::vector<GraphTraits::vertex_descriptor> forbidden;
  GraphTraits::vertex_descriptor dad;
  int inDg = in_degree(vdes,myg);
  if(inDg>0)
    {
      dad = get_father(vdes);
      sisters = get_matrisisters(dad);
      for(unsigned int i(0);i<sisters.size();++i)
	{
	  sons = get_sons(sisters[i]);
	  out += sons;
	}
      forbidden = get_brothers(vdes);
      out = take_off(out,forbidden);
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_mother_brothers_daughters(GraphTraits::vertex_descriptor vdes)
{
  std::vector<GraphTraits::vertex_descriptor> daughters;
  std::vector<GraphTraits::vertex_descriptor> brothers;
  std::vector<GraphTraits::vertex_descriptor> out;
  GraphTraits::vertex_descriptor mum;
  mum = get_mother(vdes);
  brothers = get_brothers(mum);
  for(unsigned int i(0);i<brothers.size();++i)
    {
      daughters = get_daughters(brothers[i]);
      out += daughters; // can not be double, a girl only has one father
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_mother_brothers_sons_nob(GraphTraits::vertex_descriptor vdes)
{
  std::vector<GraphTraits::vertex_descriptor> sons;
  std::vector<GraphTraits::vertex_descriptor> brothers;
  std::vector<GraphTraits::vertex_descriptor> out;
  std::vector<GraphTraits::vertex_descriptor> forbidden;
  GraphTraits::vertex_descriptor mum;
  int inDg = in_degree(vdes,myg);
  if(inDg>0)
    {
      mum = get_mother(vdes);
      brothers = get_matribrothers(mum);
      for(unsigned int i(0);i<brothers.size();++i)
	{
	  sons = get_sons(brothers[i]);
	  out += sons;
	}
      forbidden = get_brothers(vdes);
      out = take_off(out,forbidden);
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_mother_brothers_sons(GraphTraits::vertex_descriptor vdes)
{
  std::vector<GraphTraits::vertex_descriptor> sons;
  std::vector<GraphTraits::vertex_descriptor> brothers;
  std::vector<GraphTraits::vertex_descriptor> out;
  GraphTraits::vertex_descriptor mum;
  int inDg = in_degree(vdes,myg);
  if(inDg>0)
    {
      mum = get_mother(vdes);
      brothers = get_half_brothers(mum);
      for(unsigned int i(0);i<brothers.size();++i)
	{
	  sons = get_sons(brothers[i]);
	  out += sons;
	}
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_mother_sisters_sons(GraphTraits::vertex_descriptor vdes)
{
  std::vector<GraphTraits::vertex_descriptor> sons;
  std::vector<GraphTraits::vertex_descriptor> sisters;
  std::vector<GraphTraits::vertex_descriptor> out;
  GraphTraits::vertex_descriptor mum;
  int inDg = in_degree(vdes,myg);
  if(inDg>0)
    {
      mum = get_mother(vdes);
      sisters = get_half_sisters(mum);
      for(unsigned int i(0);i<sisters.size();++i)
	{
	  sons = get_sons(sisters[i]);
	  out += sons;
	}
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_cousins_nob(GraphTraits::vertex_descriptor vdes)
{
  std::vector<GraphTraits::vertex_descriptor> cousins;
  std::vector<GraphTraits::vertex_descriptor> grandparents;
  std::vector<GraphTraits::vertex_descriptor> out;
  std::vector<GraphTraits::vertex_descriptor> forbidden;
  GraphTraits::vertex_descriptor dad;
  grandparents=get_grandparents(vdes);
  for(unsigned int i(0);i<grandparents.size();++i)
    {
      cousins = get_grandchildren(grandparents[i]);
      out += cousins;
    }
  forbidden = get_brothers(vdes);
  out = take_off(out,forbidden);
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_cousins_male_nob(GraphTraits::vertex_descriptor vdes)
{
  std::vector<GraphTraits::vertex_descriptor> cousins;
  std::vector<GraphTraits::vertex_descriptor> grandparents;
  std::vector<GraphTraits::vertex_descriptor> out;
  std::vector<GraphTraits::vertex_descriptor> forbidden;
  GraphTraits::vertex_descriptor dad;
  grandparents=get_grandparents(vdes);
  for(unsigned int i(0);i<grandparents.size();++i)
    {
      cousins = get_grandsons(grandparents[i]);
      out += cousins;
    }
  forbidden = get_brothers(vdes);
  out = take_off(out,forbidden);
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_grandparents(GraphTraits::vertex_descriptor vdes)
{
  Tinv_adjacency_iterator vup;
  Tinv_adjacency_iterator vup2;
  std::vector<GraphTraits::vertex_descriptor> out;
  std::vector <GraphTraits::vertex_descriptor> parents = get_parents(vdes);
  std::vector <GraphTraits::vertex_descriptor> parentsup;
  for (std::vector<GraphTraits::vertex_descriptor>::iterator it = parents.begin(); it != parents.end();++it)
    {
      parentsup = get_parents(*it);
      out = merge(out,parentsup);
    }
  return out;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_grandchildren(GraphTraits::vertex_descriptor vdes)
{
  Tinv_adjacency_iterator vup;
  Tinv_adjacency_iterator vup2;
  std::vector<GraphTraits::vertex_descriptor> out;
  std::vector<GraphTraits::vertex_descriptor> out2;
  std::vector<GraphTraits::vertex_descriptor> children;
  std::vector<GraphTraits::vertex_descriptor> grandchildren;
  for( tie(vup,vup2) = inv_adjacent_vertices(vdes,myg); vup != vup2; ++vup)
    {
      children = get_children((*vup));
      out += children;
    }
  for(std::vector<GraphTraits::vertex_descriptor>::iterator it = out.begin(); it != out.end(); ++it)
    for( tie(vup,vup2) = inv_adjacent_vertices(*it,myg); vup != vup2; ++vup)
      {
	grandchildren = get_children((*vup));
	out2 = merge(out2,grandchildren);
      }
  return out2;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_granddaughters(GraphTraits::vertex_descriptor vdes)
{
  Tinv_adjacency_iterator vup;
  Tinv_adjacency_iterator vup2;
  std::vector<GraphTraits::vertex_descriptor> out;
  std::vector<GraphTraits::vertex_descriptor> out2;
  std::vector<GraphTraits::vertex_descriptor> children;
  std::vector<GraphTraits::vertex_descriptor> grandchildren;
  for( tie(vup,vup2) = inv_adjacent_vertices(vdes,myg); vup != vup2; ++vup)
    {
      children = get_children((*vup));
      out.insert(out.begin(),children.begin(),children.end());
    }
  for(std::vector<GraphTraits::vertex_descriptor>::iterator it = out.begin(); it != out.end(); ++it)
    for( tie(vup,vup2) = inv_adjacent_vertices(*it,myg); vup != vup2; ++vup)
      {
	grandchildren = get_daughters((*vup));
	out2 = merge(out2,grandchildren);
      }
  return out2;
}

std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_grandsons(GraphTraits::vertex_descriptor vdes)
{
  Tinv_adjacency_iterator vup;
  Tinv_adjacency_iterator vup2;
  std::vector<GraphTraits::vertex_descriptor> out;
  std::vector<GraphTraits::vertex_descriptor> out2;
  std::vector<GraphTraits::vertex_descriptor> children;
  std::vector<GraphTraits::vertex_descriptor> grandchildren;
  for( tie(vup,vup2) = inv_adjacent_vertices(vdes,myg); vup != vup2; ++vup)
    {
      children = get_children((*vup));
      out.insert(out.begin(),children.begin(),children.end());
    }
  for(std::vector<GraphTraits::vertex_descriptor>::iterator it = out.begin(); it != out.end(); ++it)
    for( tie(vup,vup2) = inv_adjacent_vertices(*it,myg); vup != vup2; ++vup)
      {
	grandchildren = get_sons((*vup));
	out2 = merge(out2,grandchildren);
      }
  return out2;
}

GraphTraits::vertex_descriptor MSNetwork::get_random_male()
{
  //  Tinv_adjacency_iterator vup;
  //  Tinv_adjacency_iterator vup2;
  GraphTraits::in_edge_iterator in_start,in_end;
  int degIn = in_degree(vertex(nbNodes-1,myg),myg);
  if(degIn>0)
    {
      int chosen = myrand(degIn);
      tie(in_start,in_end) = in_edges(vertex(nbNodes-1,myg),myg);
      in_start = move_iterator_on(in_start,in_end,chosen);
      return source(*in_start,myg);
    }
  std::cerr << " ERROR: No male in this generation" << std::endl;
  return GraphTraits::vertex_descriptor();
}

GraphTraits::vertex_descriptor MSNetwork::get_random_female()
{
  Tinv_adjacency_iterator vup;
  Tinv_adjacency_iterator vup2;
  int degIn = in_degree(vertex(nbNodes-2,myg),myg);
  if(degIn>0)
    {
      int chosen = myrand(degIn);
      tie(vup,vup2) = inv_adjacent_vertices(vertex(nbNodes-2,myg),myg);
      vup = move_iterator_on(vup,vup2,chosen);
      return *vup;
    }
  std::cerr << " ERROR: No female in this generation" << std::endl;
  return NULL;
}

GraphTraits::vertex_descriptor MSNetwork::get_random_male(int popID)
{
  GraphTraits::in_edge_iterator in_start,in_end;
  GraphTraits::vertex_descriptor fatherNode = vertex(nbNodes-1-2*popID,myg);
  int degIn = in_degree(fatherNode,myg);
  if(degIn>0)
    {
      int chosen = myrand(degIn);
      tie(in_start,in_end) = in_edges(fatherNode,myg);
      in_start = move_iterator_on(in_start,in_end,chosen);
      return source(*in_start,myg);
    }
  std::cerr << " ERROR: No male in this generation" << std::endl;
  return GraphTraits::vertex_descriptor();
}

GraphTraits::vertex_descriptor MSNetwork::get_random_female(int popID)
{
  //  Tinv_adjacency_iterator vup;
  //  Tinv_adjacency_iterator vup2;
  GraphTraits::in_edge_iterator in_start,in_end;
  GraphTraits::vertex_descriptor motherNode = vertex(nbNodes-2-2*popID,myg);
  int degIn = in_degree(motherNode,myg);
  if(degIn>0)
    {
      int chosen = myrand(degIn);
      tie(in_start,in_end) = in_edges(motherNode,myg);
      in_start = move_iterator_on(in_start,in_end,chosen);
      return source(*in_start,myg);
    }
  std::cerr << " ERROR: No male in this generation" << std::endl;
  return GraphTraits::vertex_descriptor();
}

void MSNetwork::test(int generation, int nbhuman)
{ // test that each of the current generation has two diff parents
  int deg;
  for(int i(0); i<nnpg;++i)
    {
      deg = in_degree(vertex(i+(generation%storedGeneration)*nnpg,myg),myg);
      if(deg>2)
	std::cout << " index " << i+(generation%storedGeneration)*nnpg << " i " << i << " degree " << deg << " generation " << generation << " sex " << myg[vertex(i,myg)].sex << " generation%storedGeneration " << generation%storedGeneration << " nnpg " << nnpg << std::endl;
    }
  assert(num_vertices(myg)==nbNodes);
}


std::vector<GraphTraits::vertex_descriptor> MSNetwork::population_filter(std::vector<GraphTraits::vertex_descriptor> & vdesList, int populationID) // TODO
{
  std::vector<GraphTraits::vertex_descriptor> out;
  for (std::vector<GraphTraits::vertex_descriptor>::iterator it = vdesList.begin(); it != vdesList.end();++it)
    {
      if(myg[*it].population == populationID)
	{
	  out.push_back(*it);
	}
    }
  return out;
}

bool MSNetwork::check_population(GraphTraits::vertex_descriptor vdes,int popID)
{
  return (myg[vdes].population==popID);
}


std::vector<GraphTraits::vertex_descriptor> MSNetwork::get_male(GraphTraits::vertex_descriptor vdes,bool inbreeding, std::vector<int> mating,int popID)
{
  std::vector<GraphTraits::vertex_descriptor> out;
  std::vector<GraphTraits::vertex_descriptor> forbidden;
  int mat;
  for(std::vector<int>::iterator itmat = mating.begin(); itmat != mating.end(); ++itmat)
    {
      mat =*itmat;
      switch(mat)
	{	  
	case(1):
	  {
	    out.push_back(get_random_male(popID));
	    break;
	  }
	case(2):
	  {
	    out += get_father_sisters_sons(vdes);
	    break;
	  }
	case(3):
	  {
	    out += get_mother_sisters_sons(vdes);
	    break;
	  }
	case(4):
	  {
	    out += get_mother_brothers_sons(vdes);
	    break;
	  }
	case(5):
	  {
	    out += get_father_brothers_sons(vdes);
	    break;
	  }
	default:
	  {
	    out.push_back(get_random_male(popID));
	    break;
	  }
	}
    }
  if(inbreeding>0)
    {  
      if(inbreeding==1)
	forbidden = get_half_brothers(vdes);
      else if(inbreeding==2)
	forbidden = get_brothers(vdes);
      out = take_off(out,forbidden);  
    }
  //  out = unique(out);
  return out;
}



/*void clear_vertex(GraphTraits::vertex_descriptor u,Graph& g)
{
  Graph::edge_parallel_category Cat = disallow_parallel_edge_tag;
  GraphTraits::vertex_iterator vi, viend;
  for (boost::tie(vi, viend) = vertices(g); vi != viend; ++vi)
    detail::erase_from_incidence_list(g.out_edge_list(*vi), u,Cat);
  g.out_edge_list(u).clear();
  // clear() should be a req of Sequence and AssociativeContainer,
  // or maybe just Container
  }*/
