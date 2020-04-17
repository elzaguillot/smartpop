
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
#ifndef DEF_DNA_STRUCT_H
#define DEF_DNA_STRUCT_H

#include "utils.h"

struct DnaStructure
{
  int nbCrumbMt; /** Number of cell taken by mt */
  int nbCrumbA; /** Number of cell taken by autosome */
  int nbCrumbY; /** Number of cell taken by Y chromosome */
  int nbCrumbX; /** Number of cell taken by X chromosome */
  int nbCrumbNuc; /** Number of cell taken A +X +Y */
  int nbCellPerA; /** Number of cell for each autosomal loci */
  int nbCellPerX; /** Number of cell for each X loci */
  int nbChromosomesA; /** Number of autosomal loci */
  int nbChromosomesX; /** Number of X loci */
  int sizeMt; /** Size in bp of mito */
  int sizeY; /** Size in bp of Y */
  int sizePerX; /** Size in bp of X */
  int sizePerA; /** Size in bp of A */
  DnaStructure()
  {
    nbCrumbMt = 100; /** Number of cell taken by mt */
    nbCrumbA = 2000; /** Number of cell taken by autosome */
    nbCrumbY = 100; /** Number of cell taken by Y chromosome */
    nbCrumbX = 100; /** Number of cell taken by X chromosome */
    nbCrumbNuc = 2200; /** Number of cell taken A +X +Y */
    nbCellPerA = 100; /** Number of cell for each autosomal loci */
    nbCellPerX = 100; /** Number of cell for each X loci */
    nbChromosomesA = 10; /** Number of autosomal loci */
    nbChromosomesX = 1; /** Number of X loci */
    sizeMt = 3200; /** Size in bp of mito */
    sizeY = 3200; /** Size in bp of Y */
    sizePerX = 3200; /** Size in bp of X */
    sizePerA = 3200; /** Size in bp of A */
  }
  DnaStructure(DnaStructure const &  a)
  {
    nbCrumbMt = a.nbCrumbMt; /** Number of cell taken by mt */
    nbCrumbA = a.nbCrumbA; /** Number of cell taken by autosome */
    nbCrumbY = a.nbCrumbY; /** Number of cell taken by Y chromosome */
    nbCrumbX = a.nbCrumbX; /** Number of cell taken by X chromosome */
    nbCrumbNuc = a.nbCrumbNuc; /** Number of cell taken A +X +Y */
    nbCellPerA = a.nbCellPerA; /** Number of cell for each autosomal loci */
    nbCellPerX = a.nbCellPerX; /** Number of cell for each X loci */
    nbChromosomesA = a.nbChromosomesA ; /** Number of autosomal loci */
    nbChromosomesX = a.nbChromosomesX; /** Number of X loci */
    sizeMt = a.sizeMt; /** Size in bp of mito */
    sizeY = a.sizeY; /** Size in bp of Y */
    sizePerX = a.sizePerX; /** Size in bp of X */
    sizePerA = a.sizePerA; /** Size in bp of A */
  }
  DnaStructure & operator=(DnaStructure const &  a)
  {
    if(&a!=this)
      {
	nbCrumbMt = a.nbCrumbMt; /** Number of cell taken by mt */
	nbCrumbA = a.nbCrumbA; /** Number of cell taken by autosome */
	nbCrumbY = a.nbCrumbY; /** Number of cell taken by Y chromosome */
	nbCrumbX = a.nbCrumbX; /** Number of cell taken by X chromosome */
	nbCrumbNuc = a.nbCrumbNuc; /** Number of cell taken A +X +Y */
	nbCellPerA = a.nbCellPerA; /** Number of cell for each autosomal loci */
	nbCellPerX = a.nbCellPerX; /** Number of cell for each X loci */
	nbChromosomesA = a.nbChromosomesA ; /** Number of autosomal loci */
	nbChromosomesX = a.nbChromosomesX; /** Number of X loci */
	sizeMt = a.sizeMt; /** Size in bp of mito */
	sizeY = a.sizeY; /** Size in bp of Y */
	sizePerX = a.sizePerX; /** Size in bp of X */
	sizePerA = a.sizePerA; /** Size in bp of A */
      }
    return *this;
  }
};

struct MutationRates
{
  double muMt1;
  /** Transition rate for mtDNA */
  double muMt2;
  /** Transversion rate for mtDNA */
  double muX1;
  /** Transition rate for X chromosome */
  double muX2;
  /** Transversion rate for X chromosome */
  double muY1;
  /** Transition rate for Y chromosome */
  double muY2;
  /** Transversion rate for Y chromosome */
  double muA1; // Transition A<->G C<->T i.e. change on 1 bit
  /** Transition rate for autosomes */
  double muA2; // Transversion A<->C A<->T G<->C G<->T i.e. change on 2 bits
  /** Transition rate for autosomes */
  MutationRates()
  {
    /*    muMt1 = 1e-6;
    muMt2  = 1e-7;
    muX1 = 1e-8;
    muX2 = 1e-9;
    muY1 = 1e-8;
    muY2 = 1e-9;
    muA1 = 1e-8;
    muA2 = 1e-9;*/
    muMt1 = 3e-6;
    muMt2  = 0;
    muX1 = 1e-7;
    muX2 = 0;
    muY1 = 1e-7;
    muY2 = 0;
    muA1 = 1e-7;
    muA2 = 0;
  }
  MutationRates(MutationRates const & a)
  {
    muMt1 = a.muMt1;
    muMt2 = a.muMt2;
    muX1 = a.muX1;
    muX2 = a.muX2;
    muY1 = a.muY1;
    muY2 = a.muY2;
    muA1 = a .muA1;
    muA2 = a.muA2;
  }
  MutationRates& operator=(MutationRates const & a)
  {
    if(&a != this)
      {
	muMt1 = a.muMt1;
	muMt2 = a.muMt2;
	muX1 = a.muX1;
	muX2 = a.muX2;
	muY1 = a.muY1;
	muY2 = a.muY2;
	muA1 = a.muA1;
	muA2 = a.muA2;
      }
    return *this;
  }
};

extern struct DnaStructure DNAS;
extern struct MutationRates RATES;

#endif
