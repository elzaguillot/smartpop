
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
#include "myprogram_options.h"

void parse_mycommand_line(int argc,char * argv[],Simulation & sim,int nbCores,int id)
{
  std::string expFile="SMARTPOP_parameters.xml";
  po::options_description desc;
  desc.add_options()
    ("help,h", "Produce help")
    ("verbose,v", "Verbose")
    ("header", "Ouput headers in files")
    ("fasta", "Ouput fasta files")   
    ("arl", "Ouput arlequin files")
    ("hap", "Ouput haplotype files")
    ("sfs", "Ouput Site Frequency Spectrum files")
    ("child2", "Perfect theoretical model")
    ("pooled", "Pooled sample")
    ("scattered", "Scattered sample")
    ("local", "Local sample")
    ("blind", "Sample without knowing the structure")
    ("between", "Compute the diversity between populations")
    ("inbreeding", po::value<int>()->notifier(boost::bind(&Simulation::set_inbreeding, &sim,_1)), "Prevent inbreeding")
    ("seed,s", po::value<unsigned int>()->notifier(boost::bind(&Simulation::set_seed, &sim,_1)), "Random Seed")
    ("step,t", po::value<int>()->notifier(boost::bind(&Simulation::set_step, &sim,_1)), "Time per step")
    ("nstep", po::value<int>()->notifier(boost::bind(&Simulation::set_nstep, &sim,_1)), "Number of step")
    ("popsize,p", po::value<int>()->notifier(boost::bind(&Simulation::set_popsize, &sim,_1)), "Population size")
    ("npop", po::value<int>()->notifier(boost::bind(&Simulation::set_nbPop, &sim,_1)), "Number of populations")
    ("nbSources", po::value<int>()->notifier(boost::bind(&Simulation::set_nbSources, &sim,_1)), "Number of populations to which the population sends women")
    ("mat,m", po::value<  std::vector<int> >()->multitoken(), "Mating System Preferred") // boost::bind(&Simulation::set_matsys, &sim,_1)
    ("pmate", po::value<double>()->notifier(boost::bind(&Simulation::set_pmate, &sim,_1)), "Mating System Preferred")
    ("pmig", po::value<double>()->notifier(boost::bind(&Simulation::set_pmig, &sim,_1)), "Mating System Preferred")
    ("polygamy", po::value<int>()->notifier(boost::bind(&Simulation::set_polygamyall, &sim,_1)), "Maximum polygamy authorized")
    //    ("polygamy", po::value<int>()->notifier(boost::bind(&Simulation::set_polygamyf, &sim,_1)), "Maximum polygamy authorized")
    ("polyandry", po::value<int>()->notifier(boost::bind(&Simulation::set_polygamyf, &sim,_1)), "Maximum polygamy authorized")
    ("polygyny", po::value<int>()->notifier(boost::bind(&Simulation::set_polygamy, &sim,_1)), "Maximum polygamy authorized")
    ("popsizeN", po::value< std::vector<int> >()->multitoken(), "Population size of population n") //TOLOOK
    ("matingN", po::value< std::vector<int> >()->multitoken(), "mating system of population n") //TOLOOK
    ("polygamyN,", po::value< std::vector<int> >()->multitoken(), "Polygamy of population n") //TOLOOK
    ("polygynyN,", po::value< std::vector<int> >()->multitoken(), "Polygamy of population n") //TOLOOK
    //    ("polygamyN,polyn", po::value<int>()->multitoken(), "Polygamy of population n") //TOLOOK
    ("mig", po::value< std::vector<double> >()->multitoken(), "Migration rate from a to b with rate r ") //TOLOOK
    ("migf", po::value< std::vector<double> >()->multitoken(), "Female Migration rate from a to b with rate r ") //TOLOOK
    ("migm", po::value< std::vector<double> >()->multitoken(), "Male Migration rate from a to b with rate r ") //TOLOOK
    ("mu", po::value< std::vector<double> >()->multitoken(), "mutation rates") //TOLOOK
    ("sizeMt", po::value<int>()->notifier(boost::bind(&Simulation::set_sizeMt, &sim,_1)), "Size of mtdna sequence in number of bp")
    ("sizeA", po::value<int>()->notifier(boost::bind(&Simulation::set_sizePerA, &sim,_1)), "Size of autosome sequence per loci in number of bp")
    ("sizeX", po::value<int>()->notifier(boost::bind(&Simulation::set_sizePerX, &sim,_1)), "Size of X sequence per loci in number of bp")
    ("sizeY", po::value<int>()->notifier(boost::bind(&Simulation::set_sizeY, &sim,_1)), "Size of Y sequence in number of bp")
    ("nbLociA", po::value<int>()->notifier(boost::bind(&Simulation::set_nbLociA, &sim,_1)), "Number of autosomal loci")
    ("nbLociX", po::value<int>()->notifier(boost::bind(&Simulation::set_nbLociX, &sim,_1)), "Number of X chromosome loci")
    ("output,o", po::value< std::string >()->notifier(boost::bind(&Simulation::set_fileoutput, &sim,_1)), "Output file name")
    //    ("save", po::value< std::string >()->notifier(boost::bind(&Simulation::set_save_on, &sim,_1)), "Save simulation to file ")
    //    ("load", po::value< std::string >()->notifier(boost::bind(&Simulation::load_directory, &sim,_1)), "Load simulation from file ")
    ("input,i", po::value< std::string >()->notifier(boost::bind(&Simulation::load_file, &sim,_1)), "Input parameter file ")
    ("nsimu", po::value< int >()->notifier(boost::bind(&Simulation::set_nSimu, &sim,_1)), "Number of simulations ")
    ("sample", po::value< int >()->notifier(boost::bind(&Simulation::set_sample, &sim,_1)), "Sample size ")
    ("nbsamplepop", po::value< int >()->notifier(boost::bind(&Simulation::set_nbsamplepop, &sim,_1)), "Sample size ")
    ("burnin", po::value< double >()->notifier(boost::bind(&Simulation::set_burninT, &sim,_1)), "Burnin until reaching a threshold value ")
    ("burninH", po::value<int >()->notifier(boost::bind(&Simulation::set_burninDna, &sim,_1)), "Set the DNA maker on which must be tested the burnin ")
    ("exp", po::value<std::string >(&expFile), "Export file name")
    ("demog", po::value< std::vector<double> >()->multitoken(), "Demography function") //TOLOOK
    ("demogn", po::value< std::vector<double> >()->multitoken(), "Demography function for population n") //TOLOOK
    ("migrationType", po::value< int >()->notifier(boost::bind(&Simulation::set_migrationType, &sim,_1)), "Export file name")
    ;
  po::variables_map vm;
  try
    {
      po::store(po::parse_command_line(argc, argv, desc), vm);    
      if( vm.count("help")) 
	{ 
	  sim.print_help();
	  return; 
	} 
      po::notify(vm);
    }  
  catch(boost::program_options::error& e)
    { 
      sim.print_help();
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
      exit(2) ; 
    } 
  if (vm.count("verbose")) 
    if(id==0)
      sim.set_verbose_on();
  if (vm.count("header")) 
    if(id==0)
      sim.set_header_on();
  if (vm.count("fasta")) 
    sim.set_fasta_on();
  if (vm.count("arl")) 
    sim.set_arlequin_on();
  if (vm.count("hap")) 
    sim.set_haplotype_on();
  if (vm.count("sfs")) 
    sim.set_sfs_on();
  if (vm.count("child2")) 
    sim.set_varChild(0);
  if (vm.count("pooled")) 
    sim.set_sampleType(1);
  if (vm.count("local")) 
    sim.set_sampleType(0);
  if (vm.count("scattered")) 
    sim.set_sampleType(2);
  if (vm.count("blind")) 
    sim.set_sampleType(3);
  if (vm.count("between")) 
    sim.set_sampleBetween_on();
  if (vm.count("popsizeN")) {
    std::vector<int> popn=vm["popsizeN"].as< std::vector<int> >();
    if(popn.size()==2)
      sim.set_popsize(popn[0],popn[1]);}
  if (vm.count("matingN")) {
    std::vector<int> matn=vm["matingN"].as< std::vector<int> >();
    if(matn.size()==2)
      sim.set_matsys(matn[0],matn[1]);}
  if (vm.count("mat")) {
    std::vector<int> matn=vm["mat"].as< std::vector<int> >();
    for(uint i(0);i<matn.size();++i)
      sim.set_matsys(matn[i]);}
  if (vm.count("polygamyN")) {
    std::vector<int> polyn=vm["polygamyN"].as< std::vector<int> >();
    if(polyn.size()==2){
      sim.set_polygamy(polyn[0],polyn[1]);
      sim.set_polygamyf(polyn[0],polyn[1]);}
    else
      {
	std::cerr << "Error: --polygamyN takes 3 arguments" << std::endl;
	exit(0);
      }}
  if (vm.count("polygynyN")) {
    std::vector<int> polyn=vm["polygamyN"].as< std::vector<int> >();
    if(polyn.size()==2)
      sim.set_polygamy(polyn[0],polyn[1]);
    else
      {
	std::cerr << "Error: --polygamyN takes 3 arguments" << std::endl;
	exit(0);
      }}
  if (vm.count("mig")) {
    std::vector<double> mig=vm["mig"].as< std::vector<double> >();
    if(mig.size()==3)
      {
	sim.set_migrationRate((int) mig[0],(int) mig[1],1,mig[2]);     
	sim.set_migrationRate((int) mig[0],(int) mig[1],0,mig[2]);
      }
    else
      {
	if(mig.size()==1)
	  {
	    int np = sim.get_nbPop();
	    double newmig = (double) mig[0] /(-1+ np);
	    for(int i(0);i<np;++i)
	      for(int j(0);j<np;++j)
		if(j!=i)
		  {
		    sim.set_migrationRate(i,j,1,newmig);
		    sim.set_migrationRate(i,j,0,newmig);
		  }
	  }
	else
	  {
	    std::cerr << "Error: --mig takes 3 arguments" << std::endl;
	    exit(0);
	  }
      }
  }
  if (vm.count("migf")) {
    std::vector<double> migf=vm["migf"].as< std::vector<double> >();
    if(migf.size()==3)
      sim.set_migrationRate((int) migf[0],(int) migf[1],1,migf[2]);
    else
      {
	std::cerr << "Error: --migf takes 3 arguments" << std::endl;
	exit(0);
      }
  }
  if (vm.count("migm")) {
    std::vector<double> migm=vm["migm"].as< std::vector<double> >();
    if(migm.size()==3)
      sim.set_migrationRate((int) migm[0],(int) migm[1],0,migm[2]);
    else
      {
	std::cerr << "Error: --migm takes 3 arguments" << std::endl;
	exit(0);
      }
  }
  if (vm.count("mu")) {
    std::vector<double> mu=vm["mu"].as< std::vector<double> >();
    if(mu.size()==4)
      sim.set_mu(mu[0],0,mu[1],0,mu[2],0,mu[3],0);
    if(mu.size()==8)
      sim.set_mu(mu[0],mu[1],mu[2],mu[3],mu[4],mu[5],mu[6],mu[7]);
    if((mu.size()!=8)&&(mu.size()!=4))
      {
	std::cerr << "Error: --mu takes 4 or 8 arguments as mutations rates" << std::endl;
	exit(2);
      }
  }
  if (vm.count("demogn")) {
    std::vector<double> demogn=vm["demogn"].as< std::vector<double> >();
    if(demogn.size()==4)
      sim.set_growthFactor(demogn[0],demogn[1],demogn[2],(int) demogn[3]);
    else
      {
	std::cerr << "Error: --demogn takes 4 arguments" << std::endl;
	exit(0);
      }
  }
  if (vm.count("demog")) {
    std::vector<double> demog=vm["demog"].as< std::vector<double> >();
    if(demog.size()==3)
      sim.set_growthFactor(demog[0],demog[1],demog[2]);
    else
      {
	std::cerr << "Error: --demog takes 3 arguments" << std::endl;
	exit(0);
      }
  }
  if(nbCores>1)
    {
      if(vm.count("save"))
	{
	  std::string saveon = vm["save"].as< std::string >();
	  sim.set_save_on(saveon+to_string(id));
	}
      if(vm.count("load"))
	{
	  std::string loadon = vm["load"].as< std::string >();
	  sim.load_directory(loadon+to_string(id));
	}
    }
  srand(sim.get_seed()); 
  sim.check_parameters();
  if(id==0)
    sim.export_file(expFile);
}     
