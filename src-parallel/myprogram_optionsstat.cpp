
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
#include "myprogram_optionsstat.h"

DNA parse_mycommand_line(int argc,char * argv[],simuParam & param)
{
  std::string filename =argv[1];
  std::string expFile="SMARTPOP_parameters.xml";
  po::options_description desc;
  desc.add_options()
    ("help,h", "produce help")
    ("verbose,v", "verbose")
    ("div", "ouput diversity")
    ("header", "ouput headers in files")
    ("fasta", "ouput fasta files")
    ("infasta", "ouput fasta files")
    ("inped", "ouput fasta files")
    ("nexus", "ouput nexus files")
    ("arl", "ouput arlequin files")
    ("sfs", "ouput Site Frequency Spectrum files")
    ("haplo", "haploid marker e.g. X or Y on male")
    ("out",po::value<std::string >(), "ouput Site Frequency Spectrum files")
    ;
  po::variables_map vm;
  try
    {
      po::store(po::parse_command_line(argc, argv, desc), vm);    
      po::notify(vm);
    }  
  catch(boost::program_options::error& e)
    { 
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
      exit(2) ; 
    } 
  if (vm.count("verbose")) 
    param.verbose = true;
  if (vm.count("header")) 
    param.header = true;
  if (vm.count("fasta")) 
    param.fasta = true;
  if (vm.count("nexus")) 
    param.nexus = true;
  if (vm.count("arl")) 
    param.arl = true;
  if (vm.count("sfs")) 
    param.sfs = true;
  if (vm.count("div")) 
    param.smartsave = true;
  if (vm.count("out"))
    param.fileoutput = vm["out"].as< std::string >();
  if(vm.count("infasta"))
    {
      std::cout << "reading input " << filename << std::endl;
      std::vector< std::vector<typeDNA> > mypeddna=  input_fasta(filename);
      DNA lambda(mypeddna);
      return lambda;
    }
  else
    {
      if(vm.count("haplo"))
	{
	  std::cout << "reading input (haploid) " << filename << std::endl;
	  std::vector< std::vector<typeDNA> > mypeddna=  input_pedY(filename);
	  int mysize =  input_ped_nbSY(filename);
	  DNA lambda(mypeddna);
	  lambda.set_nbSreal(mysize);
	  return lambda;
	}
      else
	{
	  std::cout << "reading input " << filename << std::endl;
	  std::vector< std::vector<typeDNA> > mypeddna=  input_ped(filename);
	  int mysize =  input_ped_nbS(filename);
	  DNA lambda(mypeddna);
	  lambda.set_nbSreal(mysize);
	  return lambda;
	}
    }
  std::cout << "WTF" << std::endl;
}
