
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
#ifndef GENETWORK_H
#define  GENETWORK_H
#include <boost/version.hpp>
#include <iostream>
#include <iomanip>
#include <utility> // std::pair
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/directed_graph.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/tuple/tuple.hpp>


#include "utils.h"

struct Node
{
  int population ;
  int generation ;
  int ID ;
  bool sex; // 0 male - 1 female
  std::string color;
  std::string name;
};

typedef boost::adjacency_list<boost::vecS,boost::vecS, boost::bidirectionalS, Node, boost::no_property > Graph;
//typedef boost::adjacency_list<boost::vecS,boost::listS, boost::bidirectionalS, Node, boost::no_property > Graph;
//typedef boost::adjacency_list<boost::listS,boost::listS, boost::bidirectionalS, Node, boost::no_property > Graph;
typedef boost::graph_traits<Graph> GraphTraits;
typedef boost::inv_adjacency_iterator_generator<Graph,GraphTraits::vertex_descriptor, GraphTraits::in_edge_iterator >::type Tinv_adjacency_iterator; 

class MSNetwork
{
 public:
  MSNetwork();
  MSNetwork(int a);
  MSNetwork(int a, int b);
  MSNetwork(MSNetwork const & a);
  ~MSNetwork();
  MSNetwork& operator = (MSNetwork const &);
  // protected:
  unsigned nbNodes;
  int storedGeneration;
  int nnpg; // number of node per generation
  int nbPop;
  std::vector<int> popsizes;
  Graph myg;
  void test(int, int);

  void Init(std::vector<int> nbFemales,std::vector<int> popsizes,int generation); // return the indexes of females
  void Init2();// called by burnin
  void mateOne(int,int,int,int);
  void writeGraphviz(std::string filename,int generation);
  void set_nbPop(int a,std::vector<int> maxpopsize,std::vector<int> popsizes,std::vector<int> nbFemales);
  //  void resize(int population,int newsize);
  void reset(std::vector<int> nbFemales,std::vector<int> popsizes);
  void clear_premating(int gparent);
  void clear_postmating();
  int get_indivID(int population,int generation,int index);
  bool share_one_parent(int ID1, int ID2);
  bool share_two_parents(int ID1, int ID2);
  GraphTraits::vertex_descriptor get_indivNode(int population,int generation,int index);
  std::vector<int> get_femalesID(int generation, int popID);
  std::vector<int> get_femalesID(int generation);
  std::vector<int> get_malesID(int generation, int popID);
  std::vector<int> get_malesID(int generation);
  GraphTraits::vertex_descriptor get_father(int index);
  GraphTraits::vertex_descriptor get_mother(int index);
  std::vector<GraphTraits::vertex_descriptor> get_parents(int index);
  std::vector<GraphTraits::vertex_descriptor> get_children(int index);
  std::vector<GraphTraits::vertex_descriptor> get_sons(int index);
  std::vector<GraphTraits::vertex_descriptor> get_daughters(int index);
  std::vector<GraphTraits::vertex_descriptor> get_brothers(int index);
  std::vector<GraphTraits::vertex_descriptor> get_half_brothers(int index);
  std::vector<GraphTraits::vertex_descriptor> get_mother_sons(int index);
  std::vector<GraphTraits::vertex_descriptor> get_father_sons(int index);
  std::vector<GraphTraits::vertex_descriptor> get_sisters(int index);
  std::vector<GraphTraits::vertex_descriptor> get_mother_daughters(int index);
  std::vector<GraphTraits::vertex_descriptor> get_father_daughters(int index);
  std::vector<GraphTraits::vertex_descriptor> get_grandparents(int index);
  std::vector<GraphTraits::vertex_descriptor> get_mother_parents(int index);
  std::vector<GraphTraits::vertex_descriptor> get_father_parents(int index);
  GraphTraits::vertex_descriptor get_mother_mother(int index);
  GraphTraits::vertex_descriptor get_mother_father(int index);
  GraphTraits::vertex_descriptor get_father_mother(int index);
  GraphTraits::vertex_descriptor get_father_father(int index);
  std::vector<GraphTraits::vertex_descriptor> get_mother_brothers(int index);
  std::vector<GraphTraits::vertex_descriptor> get_father_sisters(int index);
  std::vector<GraphTraits::vertex_descriptor> get_grandchildren(int index);
  std::vector<GraphTraits::vertex_descriptor> get_grandsons(int index);
  std::vector<GraphTraits::vertex_descriptor> get_granddaughters(int index);

  GraphTraits::vertex_descriptor get_father(GraphTraits::vertex_descriptor);
  GraphTraits::vertex_descriptor get_mother(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_parents(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_children(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_sons(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_daughters(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_brothers(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_half_brothers(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_matribrothers(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_patribrothers(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_mother_sons(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_father_sons(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_sisters(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_half_sisters(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_matrisisters(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_patrisisters(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_mother_daughters(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_father_daughters(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_mother_parents(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_father_parents(GraphTraits::vertex_descriptor);
  GraphTraits::vertex_descriptor get_mother_mother(GraphTraits::vertex_descriptor);
  GraphTraits::vertex_descriptor get_mother_father(GraphTraits::vertex_descriptor);
  GraphTraits::vertex_descriptor get_father_mother(GraphTraits::vertex_descriptor);
  GraphTraits::vertex_descriptor get_father_father(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_mother_brothers(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_father_sisters(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_grandparents(GraphTraits::vertex_descriptor vdes);
  std::vector<GraphTraits::vertex_descriptor> get_grandchildren(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_grandsons(GraphTraits::vertex_descriptor);
  std::vector<GraphTraits::vertex_descriptor> get_granddaughters(GraphTraits::vertex_descriptor);

  std::vector<GraphTraits::vertex_descriptor> get_father_sisters_sons_nob_nofbs_nomzs(GraphTraits::vertex_descriptor vdes);
  std::vector<GraphTraits::vertex_descriptor> get_father_sisters_sons(GraphTraits::vertex_descriptor vdes);
  std::vector<GraphTraits::vertex_descriptor> get_father_brothers_sons(GraphTraits::vertex_descriptor vdes);
  std::vector<GraphTraits::vertex_descriptor> get_mother_brothers_sons(GraphTraits::vertex_descriptor vdes);
  std::vector<GraphTraits::vertex_descriptor> get_mother_sisters_sons(GraphTraits::vertex_descriptor vdes);
  std::vector<GraphTraits::vertex_descriptor> get_mother_brothers_sons_nob(GraphTraits::vertex_descriptor vdes);
  std::vector<GraphTraits::vertex_descriptor> get_mother_brothers_daughters(GraphTraits::vertex_descriptor vdes);
  std::vector<GraphTraits::vertex_descriptor> get_cousins_nob(GraphTraits::vertex_descriptor vdes);
  std::vector<GraphTraits::vertex_descriptor> get_cousins_male_nob(GraphTraits::vertex_descriptor vdes);

  GraphTraits::vertex_descriptor get_random_male();
  GraphTraits::vertex_descriptor get_random_female();
  GraphTraits::vertex_descriptor get_random_male(int popID);
  GraphTraits::vertex_descriptor get_random_female(int popID);

  std::vector<GraphTraits::vertex_descriptor> population_filter(std::vector<GraphTraits::vertex_descriptor> & vdesList, int populationID); // TODO
  bool check_population(GraphTraits::vertex_descriptor vdes,int popID);
  std::vector<GraphTraits::vertex_descriptor> get_male(GraphTraits::vertex_descriptor vdes,bool inbreeding, std::vector<int> mating,int popID);
};

template<class T>
T move_iterator_on(T myIterator, T myIteratorend,int nstep)
{
  for (int i(0); i<nstep; ++i)
    {
      ++myIterator; // careful this is not robust to seg fault itertaing over the edge, meh ...
      /*      if(myIterator == myIteratorend)
	{
	  break;
	  }*/
    }
  return myIterator;
}

#endif
