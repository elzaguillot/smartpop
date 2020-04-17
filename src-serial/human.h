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

#ifndef DEF_HUMAN2
#define DEF_HUMAN2
#include <string>
#include <string.h> // need this one for memcpy
#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <stdio.h>
#include "utils.h"
#include "dnastructure.h"
#include <fstream>
#include "dnastructure.h"

#define TRACK_MATE_AND_CHILDREN 1

class Human
{
 public:
  Human(); //
  Human(typeDNA *tabX,typeDNA *pMtDna); /** Constructor from DNA sequences  */
  Human(Human const &);// copy constructor
  ~Human();
  friend std::ostream &operator<<(std::ostream&, Human const&); /** default output */
  Human& operator = (Human const&);
  void reborn(typeDNA *newDna,typeDNA* newMtDna); /** Change an old instance instead of creating a new one each generation */
  void reborn(Human  * &); /** Clone into an existing instance */
  void reset();
  int getNbMate() const; /** Return Number of Mates for this Human ID */
  int getNbChild() const; /** Return Number of children for this Human ID */
  void hadChild(int n); /** Increase the nb of children of this human by n */
  void newMate(); /** Increase the nb of mates by 1 */
  typeDNA distance_dna(Human const& a);
  int IDNetwork;
  typeDNA* dna; /** DNA sequences, all stored in an array of 64 bits . each chromosome has a fixed number of cell */
  typeDNA* mtDna;
  // protected:
  int nbChild; /** number of children */
  int nbMate; /** number of mates */
  int destinationPop;
 private:
  std::ostream& print(std::ostream&) const;

};



#endif
