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

#include "utils.h"
boost::mt19937 rng;

double factorial(int n) 
{
  if ((n == 1)||(n == 0))
    return 1;
  else
    {
      double k =  factorial(n - 1);
      std::cout << n << " " << k << std::endl;
      return n * k;//factorial(n - 1);
    }
}


int rBinom(int n, double p) { 
	boost::binomial_distribution<> dist(n, p);
	return dist(rng);
}

int rPoisson(double lambda)
{
boost::variate_generator<boost::mt19937&, boost::poisson_distribution<> > 
dist(rng,boost::poisson_distribution<>(lambda));
  return dist();	
}

int * poisson_array(int N,double lambda)
{
  int* poissonArray=new int[N];  
  for(int i(0);i<N;i++)
    {
      poissonArray[i]=rPoisson(lambda);
    }  
  return poissonArray;
}

int * poisson_array_smart(int Nparents, int Nchildren)
{
  int* poissonArray=new int[Nparents];
  int children(Nchildren),parents(Nparents+1);
  double raison;
  for(int i(0);i<Nparents-1;i++)
    {
      parents--;
      raison=1.0/parents;
      poissonArray[i]=rBinom(children,raison);
      children-=poissonArray[i];
    }
  poissonArray[Nparents-1]=children;
  return poissonArray;
}

int myrand(int const m)
{
  return (int)(rand()/((double)RAND_MAX+1)*m);
}

int logrand(int m)
{
  return pow(10,myrand(m))*(myrand(10)+1);
}


int dist_dna(int a, int b)
{  
  // WTJW's code to calculate distance btween 2seq coded in quaternary
  int c(0),d(0),e(0),f(0),g(0);
  c=a^b; //a crumb will be 00 iff it is equal in a and b
  if(c) // speed up if 32 sites are often found similar (which should be !)
    {
      d=c&0x55555555; //keep only first of the 2 bits // 5=0101
      e=c|(d<<1); // Now the upper bit in each crumb is on if either the upper or lower bit in c is
      f = (e<<1) & 0x55555555; // move the upper bit to th elower bit and zero the upper bit
      g=f+(f>>2); // Crumb i contains sum of original crumb i and original crumb i+1
      g&=0x33333333; // Zero out odd numbered crumbs to give us space to accumulate larger sums without overflow
      g=g+(g>>4); // Crumb i contains sum of originak crumbs i, i+1, i+2, i+3 (for even i)
      g=g+(g>>8); // Crumb i contains some of original crumbs i..i+7
      g&=0x000F0000F; // Zero out some more unimportant bits to give us more space again
      g=g+(g>>16); // Crumb i contains sum of original crumb i..i+7 for i==0 or i==8
      g&=0x0000001F;// Get rid of higher bits, leaving the total count (which must be <=16) in g
      return g;
    }
  else
    {
      return 0;
    }
}

unsigned long long int dist_dna(unsigned long long int a, unsigned long long int b)
{
  // WTJW's code to calculate distance btween 2seq coded in quaternary
  unsigned long long int c(0),d(0),e(0),f(0),g(0);
  c=a^b; //a crumb will be 00 iff it is equal in a and b
  if(c)
    {
      d=c&U64(0x5555555555555555); //keep only first of the 2 bits // 5=0101
      e=c|(d<<1); // Now the upper bit in each crumb is on if either the upper or lower bit in c is
      f=(e<<1)&U64(0x5555555555555555); // move the upper bit to th elower bit and zero the upper bit
      g=f+(f>>2); // Crumb i contains sum of original crumb i and original crumb i+1
      g&=U64(0x3333333333333333); // Zero out odd numbered crumbs to give us space to accumulate larger sums without overflow //3=0011
      g=g+(g>>4); // Crumb i contains sum of originak crumbs i, i+1, i+2, i+3 (for even i)
      g=g+(g>>8); // Crumb i contains some of original crumbs i..i+7
      g&=U64(0x000F000F000F000F); // Zero out some more unimportant bits to give us more space again
      g=g+(g>>16); // Crumb i contains sum of original crumb i..i+7 for i==0 or i==8
      g&=U64(0x000000FF000000FF);// Get rid of higher bits, leaving the total count (which must be <=16) in g
      g=g+(g>>32);
      g&=U64(0x00000000000000FF);// Get rid of higher bits, leaving the total count (which must be <=32) in g
      return g;
    }
  else 
    {
      return 0;
    }
}

#ifdef HAS_UINT128
__uint128_t dist_dna(__uint128_t a, __uint128_t b) // NOT WORKING ><
{
  // WTJW's code to calculate distance btween 2seq coded in quaternary
  __uint128_t c(0),d(0),e(0),f(0),g(0);
  c=a^b;
  d=c&0x55555555555555555555555555555555;
  e=c|(d<<1); 
  f = (e<<1) & 0x55555555555555555555555555555555;
  g=f+(f>>2); 
  g&=0x33333333333333333333333333333333;
  g=g+(g>>4); 
  g=g+(g>>8); 
  g&=0x000F000F000F000F000F000F000F000F;
  g=g+(g>>16); 
  g&=0x000000FF000000FF000000FF000000FF;
  g=g+(g>>32);
  g=g+(g>>64);
  g&=0x0000000000000000000000000001FFFF;
  return g;
}
#endif	// HAS_UINT128

unsigned long long convertDnaToInt(std::string dna)
{
  unsigned long long d,base(0);
  unsigned long long dnaInt(0);
  d=dna.std::string::size();
  if(d>0)
    {
      while(d>0)
	{
	  d--;
	  base=convertToInt(dna[d]);
	  dnaInt=(base<<2*d)+dnaInt;
	}
    }  
  return dnaInt;  
}


double harmonic(int n)
{
  if(n==0)
    return 0;
  if(n==1)
    {
      return 1;
    }
  else
    {
      return((1/double(n)) +harmonic(n-1));
    }
}

char convertToChar(unsigned number)
{
  switch(number){
  case(0):    
    return 'A';
  case(1):
    return 'G';
  case(2):
    return 'C';
  case(3):
    return 'T';
  default:
    return '-';
  }  
}

std::string convertIntToQuaternaryString(typeDNA dna,int nbcrumb)
{
  std::string strDNA;
  strDNA.reserve(nbcrumb);
  // bit manipulation to to translate the in to the sequence in base 4
  // a crumb = 2 bits
  while (nbcrumb >0) {
    strDNA+=convertToChar((unsigned) (dna&3));
    dna >>= 2;
    nbcrumb--;
  }
  return strDNA;  
}


double tajima_a1(int n)
{
  if(n == 1)
    return 1;
  else
    return (1.0 / (double) n) + tajima_a1(n-1);
}

double tajima_a2(int n)
{
  if(n == 1)
    return 1;
  else
    return (1.0 / (double) (n * n)) + tajima_a2(n-1);
}

double tajima_b1(int n)
{
  return (double) (n + 1) / (double) (3 * (n - 1));
}

double tajima_b2(int n)
{
  return (double) (2 * (n*n + n + 3)) / (double) (9 * n * (n-1));
}

double tajima_c1(int n)
{
  return tajima_b1(n) - 1.0 / tajima_a1(n-1);
}

double tajima_c2(int n)
{
  return tajima_b2(n) - (double) (n +1) / (tajima_a1(n-1) * n) + tajima_a2(n-1) / (tajima_a1(n-1) * tajima_a1(n-1));
}

unsigned long long dna_distance(unsigned long long c)
{
  unsigned long long g = 0;
  unsigned long long f = 0;
  unsigned long long d = 0;
  unsigned long long e = 0;
  if(c)
    {
      d = c & U64(0x5555555555555555); //keep only first of the 2 bits // 5 = 0101
      e = c | (d << 1); // Now the upper bit in each crumb is on if either the upper or lower bit in c is
      f = (e >> 1 ) & U64(0x5555555555555555); // move the upper bit to th elower bit and zero the upper bit
      e = f >> 2;
      g = f + e; // Crumb i contains sum of original crumb i and original crumb i+1      
      g &= U64(0x3333333333333333); // Zero out odd numbered crumbs to give us space to accumulate larger sums without overflow //3 = 0011
      g = g + (g >> 4); // Crumb i contains sum of originak crumbs i, i+1, i+2, i+3 (for even i)
      g = g + (g >> 8); // Crumb i contains some of original crumbs i..i+7
      g &= U64(0x000F000F000F000F); // Zero out some more unimportant bits to give us more space again
      g = g + (g >> 16); // Crumb i contains sum of original crumb i..i+7 for i == 0 or i == 8
      g &= U64(0x000000FF000000FF);// Get rid of higher bits, leaving the total count (which must be < = 16) in g
      g = g + (g >> 32);
      g &= U64(0x00000000000000FF);// Get rid of higher bits, leaving the total count (which must be < = 32) in g
      return g;
    }
  else
    {
      return 0;
    }
}

std::vector<std::string> split_string(std::string tosplit,std::string sep)
{
  unsigned found;
  int sizesep=sep.size();
  std::vector<std::string> out;
  found=tosplit.find(sep.c_str());
  if(found!=std::string::npos)
    {
      std::string tmp=tosplit.substr(0,found);
      if(found<tosplit.size()-sizesep)
	{
	  out=split_string(tosplit.substr(found+sizesep,tosplit.size()),sep);
	}
      out.push_back(tmp);
    }
  else
    {
      if ( !((tosplit[0]<33)))
	out.push_back(tosplit);
    }
  return out;
}

std::vector<unsigned long long> convertDnaToIntSplit(std::string dna)
{ // a dna a string, output a vector of numbers corresponding to this dna  
  std::vector<unsigned long long> out;
  std::string sub;
  int sz = dna.size();
  int nbcell = sz / CRUMB;
  for(int i=0; i<nbcell-1 ; ++i)
    {
      sub = dna.substr(i*CRUMB,CRUMB);// substr (position,length)
      out.push_back(convertDnaToInt(sub));
    }
  // TODO, handle DNA not exactly CRUMB size sort out the end
  return out;
}

std::vector<typeDNA> dna_string_to_int(std::string stringput)
{
  int s=stringput.size();
  int l = floor((double) s/CRUMB);
  std::vector<typeDNA> out(l);
  typeDNA newCell;
  for(int i(0);i<l;++i)
    {
      newCell = convertDnaToInt(stringput.substr(i*CRUMB,CRUMB));
      out[i] = newCell;
    }
  int incomplete=s-l*CRUMB;
  if(incomplete>0)
    {
      newCell = convertDnaToInt(stringput.substr(l*CRUMB,incomplete));
      out.push_back(newCell);
    }
  return out;
}

double harm_square(int i)
{
  if(i<2)
    return 0;
  else 
    return 1.0 / (double) (i*i) + harm_square(i-1);
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
{
  std::stringstream ss(s);
  std::string item;
  while(std::getline(ss, item, delim)) 
    {
      elems.push_back(item);
    }
  return elems;
}

std::vector<std::string> split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

bool isFloat( std::string myString ) 
{
  std::istringstream iss(myString);
  float f;
  iss >> std::noskipws >> f; // noskipws considers leading whitespace invalid
  // Check the entire string was consumed and if either failbit or badbit is set
  return iss.eof() && !iss.fail();
  // stackoverflow.com/questions/447206/c-isfloat-function
}

bool isInt( std::string myString ) 
{
  std::istringstream iss(myString);
  int f;
  iss >> std::noskipws >> f; // noskipws considers leading whitespace invalid
  // Check the entire string was consumed and if either failbit or badbit is set
  return iss.eof() && !iss.fail();
  // stackoverflow.com/questions/447206/c-isfloat-function
} 

std::vector<int> sample_without_replacement(int maxnb,int n)
{
  if(maxnb<n)
    {
      std::cerr << "ERROR : sample without replacement of more elements than existing" << std::endl;
      return std::vector<int>();
    }
  std::vector<int> out;
  std::vector<int> tosample(maxnb);
  for(int i(0);i<maxnb;++i)
    tosample[i]=i;
  if(maxnb==n)
    return tosample;
  --maxnb;
  for(int i(0);i<n;++i)
    {
      int r = myrand(maxnb);
      int a = tosample[r];
      out.push_back(a);
      tosample[r] = tosample[maxnb];
      tosample[maxnb] = a;
      --maxnb;
    }
  return out;
  // algorithm inspired by Knuth Fisher found on :
  // http://stackoverflow.com/questions/196017/unique-random-numbers-in-o1#196065
}
