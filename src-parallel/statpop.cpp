
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
#include "myprogram_optionsstat.h"
using namespace std;

int main(int argc, char* argv[])
{
  std::string filename =argv[1];
  simuParam para;
  DNA lambda = parse_mycommand_line(argc,argv,para);
  if(para.header)
    {
      if(para.verbose)
	cout << "Creating file with header " << filename+"_div.txt" << endl;
      output_diversity_headline_nosimu(filename+"_div.txt");
    }
  if(para.smartsave)
    {
      if(para.verbose)
	cout << "Outputting diversity into " << filename+"_div.txt" << endl;
      lambda.output_diversity_with_time_nosimu(filename+"_div.txt",lambda.get_N());
    }
  if(para.sfs)
    {
      if(para.verbose)
	cout << "Creating file with folded site frequency spectrum " << filename+"smart.sfs" << endl;
      lambda.SFS_folded(filename+"smart.sfs");
    }
  if(para.fasta)
    {
      if(para.verbose)
	cout << "Creating fasta file " << filename+".fasta" << endl;
      lambda.output_fasta(filename);
    }
  if(para.nexus)
    {
      if(para.verbose)
	cout << "Creating nexus file " << filename+".nex" << endl;
      lambda.output_nexus(filename);
    }
  if(para.arl)
    {
      if(para.verbose)
	cout << "Creating arlequin file " << filename+".arp" << endl;
      lambda.output_arlequin(filename);
    }
  return 0;
}
