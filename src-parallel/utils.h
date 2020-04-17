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

#ifndef DEF_UTIL
#define DEF_UTIL
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <functional>
#include <sstream>
#include <climits>
#include "boost/random/binomial_distribution.hpp"
#include "boost/random/poisson_distribution.hpp"
#include "boost/random/variate_generator.hpp"
#include "boost/random/mersenne_twister.hpp"
#include <boost/foreach.hpp>

// Try to sniff out whether or not a long is 64 bits, or 32.  This is needed to decide whether to apply a "ULL" or "UL" suffix to 64-bit numeric literals.
//#ifdef __SIZEOF_LONG__
//// __SIZEOF_LONG__ is what g++ uses, according to http://gcc.gnu.org/onlinedocs/cpp/Common-Predefined-Macros.html.
//#if __SIZEOF_LONG__ == 8
//#define LONG_IS_64_BITS
//#endif
//#endif
#if ULONG_MAX > 0xFFFFFFFFUL
#define LONG_IS_64_BITS
#endif

// If compiling with g++ on a 32-bit machine, g++ demands that 64-bit unsigned constants have a special "ULL" ("unsigned long long") suffix.
// (Hopefully this nonstandard suffix works on other compilers.)  But if we're compiling on a 64-bit machine, these constants need only a "UL" suffix.
// This macro attempts to figure out the right suffix, and glue it on.  LONG_IS_64_BITS needs to be
// either set by the user or guessed from other predefined constants.
#ifdef LONG_IS_64_BITS
#define U64(x) x##UL
#else
#define U64(x) x##ULL
#endif

#if defined(_MSC_VER) || defined(__BORLANDC__) // windobe version 64bits int
typedef unsigned __int64 ulong64;
typedef signed __int64 long64;
#define CRUMB 32 // number of base that can be coded in a number, 64 for uint128, 32 for ulong64
#define typeDNA  __ulong64
#else //linux version 128 bits int
typedef unsigned long long ulong64;
typedef signed long long long64;
#define CRUMB 32 // number of base that can be coded in a number, 64 for uint128, 32 for ulong64
#define typeDNA ulong64

#endif

#define rand2 (int)((rand()/((double)RAND_MAX+1)*2))

double factorial(int);
int rBinom(int, double );
int rPoisson(double);
int * poisson_array(int N,double lambda);
int myrand(int const);
int logrand(int);
ptrdiff_t myrand(ptrdiff_t);
int dist_dna(int, int);
unsigned long long int dist_dna(unsigned long long int, unsigned long long int);
unsigned long long convertDnaToInt(std::string dna);
std::vector<unsigned long long> convertDnaToIntSplit(std::string dna);
#ifdef HAS_UINT128
__uint128_t dist_dna(__uint128_t a, __uint128_t b);// NOT WORKING
#endif	// HAS_UINT128
double harmonic(int n);
int * poisson_array_smart(int Nparents, int Nchildren);
char convertToChar(unsigned number);
std::string convertIntToQuaternaryString(typeDNA dna,int nbcrumb = CRUMB);
double tajima_a1(int n); // warning in tajima D computaion call with n-1
double tajima_a2(int n); // warning in tajima D computaion call with n-1
double tajima_b1(int n);
double tajima_b2(int n);
double tajima_c1(int n);
double tajima_c2(int n);
unsigned long long dna_distance(unsigned long long c);
std::vector<std::string> split_string(std::string tosplit,std::string sep);
std::vector<typeDNA> dna_string_to_int(std::string);
double harm_square(int i);
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);
bool isFloat(std::string myString);
bool isInt(std::string myString);
std::vector<int> sample_without_replacement(int maxnb,int n);


template <class T>
std::ostream & operator<<(std::ostream & fs, const std::vector<T>  & a)
{
  for(typename std::vector<T>::const_iterator it = a.begin(); it!= a.end(); ++it)
    fs << *it << " " ;
  fs << std::endl;
  return fs;
}


template <class T>
inline std::string to_string (const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

template<class T >
std::string convertToStr(T number)
{
  switch(number){
  case(0):    
    return "A";
  case(1):
    return "G";
  case(2):
    return "C";
  case(3):
    return "T";
  default:
    return "-";
  }  
}

template<class T >
int convertToInt(T letter)
{
  switch(letter){
  case('A'):
    return 0;
  case('G'):
    return 1;
  case('C'):
    return 2;
  case('T'):
    return 3;
  default:
    return 0;
  }
}


template <typename T>
std::vector<T>& operator+=(std::vector<T>& a, const std::vector<T>& b)
{
    a.insert(a.end(), b.begin(), b.end());
    return a;
}

template <typename T>
std::vector<T>& operator+=(std::vector<T>& aVector, const T& aObject)
{
  aVector.push_back(aObject);
  return aVector;
}

template <typename T>
std::vector<T> unique(std::vector<T>& a)
{
  std::vector<T> out;
  int test = 0;
  for(unsigned int i(0); i<a.size();++i)
    {
      test = 0;
      for(unsigned int j(0); j<out.size();++j)
	  if(a[i]==out[j])
	    test=1;
      if (test==0)
	out.push_back(a[i]);
    }
  return out;
}

template <typename T>
std::vector<T> take_off(std::vector<T>& a,std::vector<T>& b) // return a without the elements of b
{
  std::vector<T> out;
  int test = 0;
  for(unsigned int i(0); i<a.size();++i)
    {
      test = 0;
      for(unsigned int j(0); j<b.size(); ++j)
	if(a[i] == b[j])
	  test = 1;
      if (test == 0)
	out.push_back(a[i]);
    }
  return out;
}

template <typename T>
std::vector<T> merge(std::vector<T>& a,std::vector<T>& b) // return a without the elements of b
{
  std::vector<T> out;
  int test = 0;
  if(a.size() == 0)
    out = b;
  else if (b.size()==0)
    out = a;
  else if( a.size() > b.size() )
    {
      out = a;
      for(unsigned int i(0); i<b.size();++i)
	{
	  test = 0;
	  for(unsigned int j(0); j<a.size(); ++j)
	    if(b[i] == a[j])
	      test = 1;
	  if (test == 0)
	    out.push_back(b[i]);
	}
    }
  else
    {
      out = b;
      for(unsigned int i(0); i<a.size();++i)
	{
	  test = 0;
	  for(unsigned int j(0); j<b.size(); ++j)
	    if(a[i] == b[j])
	      test = 1;
	  if (test == 0)
	    out.push_back(a[i]);
	}
    }
  return out;
}

#endif
