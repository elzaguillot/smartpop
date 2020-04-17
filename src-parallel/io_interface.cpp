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

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "io_interface.h"

// Loads debug_settings structure from the specified XML file
void simulation_settings::load(const std::string &filename)
{
    // Create an empty property tree object
    using boost::property_tree::ptree;
    ptree pt;
    read_xml(filename, pt,boost::property_tree::xml_parser::trim_whitespace );
    //output
    simuP.header = pt.get("output.header",false);
    simuP.fasta = pt.get("output.fasta",false);
    simuP.nexus = pt.get("output.nexus",false);
    simuP.arl = pt.get("output.arlequin",false);
    simuP.hap = pt.get("output.haplotype",false);
    simuP.verbose = pt.get("output.verbose",false);
    simuP.fileoutput = pt.get("output.filename","simu"+to_string(simuP.seed));
    //simulation
    simuP.seed = pt.get("simulation.seed",time(0));
    simuP.step = pt.get("simulation.step",1000);
    simuP.nstep = pt.get("simulation.nstep",1);
    simuP.nsimu = pt.get("simulation.nsimu",500);
    simuP.burninT = pt.get("simulation.burnin",0.0);
    simuP.burninDna = pt.get("simulation.burninH",1);
    //model
    metapop.nbPop = pt.get("model.nbPop",1);
    metapop.sample = pt.get("simulation.sample",50);
    metapop.sampleBetween = pt.get("simulation.betweenDiv",false);
    std::string stype = pt.get("simulation.sampleType","blind");
    if(stype == "local")
      metapop.sampleType = 0;
    if(stype == "pooled")
      metapop.sampleType = 1;
    if(stype == "scattered")
      metapop.sampleType = 2;
    if(stype == "blind")
      metapop.sampleType = 3;
    metapop.nbsamplepop=pt.get("simulation.nbpopsample",1);
    metapop.nbSources = pt.get("model.nbSources",0);
    metapop.popsize = pt.get("model.popsize",100);
    metapop.migrationType = pt.get("model.migrationType",0);
    metapop.pmig = pt.get("model.pmig",0.0);
    metapop.pmate = pt.get("model.pmate",0.0);
    for(int i(0);i<metapop.nbPop;++i)
      for(int j(0);j<metapop.nbPop;++j)
	{
	  if( pt.count("model.migration.F."+to_string(i)+"."+to_string(j)) > 0 )
	    metapop.migrationRatesF[i][j] = pt.get("model.migration.F."+to_string(i)+"."+to_string(j),0.0);	  
	}
    for(int i(0);i<metapop.nbPop;++i)
      for(int j(0);j<metapop.nbPop;++j)
	{
	  if( pt.count("model.migration.M."+to_string(i)+"."+to_string(j)) > 0 )
	    metapop.migrationRatesM[i][j] = pt.get("model.migration.M."+to_string(i)+"."+to_string(j),0.0);	  
	}
    mu.muMt1 = pt.get("model.mu.mt.transition",1e-7);
    mu.muMt2 = pt.get("model.mu.mt.transversion",1e-8);
    mu.muX1 = pt.get("model.mu.x.transition",1e-8);
    mu.muX2 = pt.get("model.mu.x.transversion",1e-9);
    mu.muY1 = pt.get("model.mu.y.transition",1e-8);
    mu.muY2 = pt.get("model.mu.y.transversion",1e-9);
    mu.muA1 = pt.get("model.mu.a.transition",1e-8);
    mu.muA2 = pt.get("model.mu.a.transversion",1e-8);
    pdna.sizePerX = pt.get("model.dna.sizeX",3200);
    pdna.sizePerA = pt.get("model.dna.sizeA",3200);
    pdna.sizeY = pt.get("model.dna.sizeY",3200);
    pdna.sizeMt = pt.get("model.dna.sizeMT",3200);
    pdna.nbChromosomesA = pt.get("model.dna.nbLociA",1);
    pdna.nbChromosomesX = pt.get("model.dna.nbLociX",1);
    BOOST_FOREACH( ptree::value_type const& v, pt.get_child("model") ) {
      if( v.first == "population" ) {
	popVal f;
	f.popsize = v.second.get("popsize",100);
	f.inbreeding = v.second.get("inbreeding",0);
	f.varChild = v.second.get("varianceNbChildren",1);
	f.polygamyf = v.second.get("polyandry",1);
	f.polygamy = v.second.get("polygyny",1);
	f.polygamyf = v.second.get("polyandry",1);
	f.growthFactor[0] = v.second.get("demography.a",0.0);
	f.growthFactor[1] = v.second.get("demography.b",1.0);
	f.growthFactor[2] = v.second.get("demography.c",0.0);
	BOOST_FOREACH(ptree::value_type const& val, v.second)
	  {	
	    if (val.first== "MS")
	      {
		f.matingSystem.push_back(val.second.get("",1));
	      }
	  }
	populations.push_back(f);
      }}
}


// Saves the debug_settings structure to the specified XML file
void simulation_settings::save(const std::string &filename)
{
   // Create an empty property tree object
   using boost::property_tree::ptree;
   ptree pt;
   ptree ptchild;
   pt.put("output.header",simuP.header);
   pt.put("output.fasta",simuP.fasta);
   pt.put("output.nexus",simuP.nexus);
   pt.put("output.arlequin",simuP.arl);
   pt.put("output.haplotype",simuP.hap);
   pt.put("output.verbose",simuP.verbose);
   pt.put("output.filename",simuP.fileoutput);
   pt.put("simulation.seed",simuP.seed);
   pt.put("simulation.step",simuP.step);
   pt.put("simulation.nstep",simuP.nstep);
   pt.put("simulation.nsimu",simuP.nsimu);
   pt.put("simulation.burnin",simuP.burninT);
   pt.put("simulation.burninH",simuP.burninDna);

    //model
   pt.put("model.popsize",metapop.popsize);
   pt.put("simulation.sample",metapop.sample);
   pt.put("simulation.nbpopsample",metapop.nbsamplepop);
   pt.put("simulation.betweenDiv",metapop.sampleBetween);
   if(metapop.sampleType == 0)
     pt.put("simulation.sample","local");
   if(metapop.sampleType == 1)
     pt.put("simulation.sample","pooled");
   if(metapop.sampleType == 2)
     pt.put("simulation.sample","scattered");
   if(metapop.sampleType == 3)
     pt.put("simulation.sample","blind");
   pt.put("model.migrationType",metapop.migrationType);
   pt.put("model.nbPop",metapop.nbPop);
   pt.put("model.nbSources",metapop.nbSources);
   pt.put("model.pmate",metapop.pmate);
   pt.put("model.pmig",metapop.pmig);
   if(metapop.migrationRatesF.size()>0)
     for(int i(0);i<metapop.migrationRatesF.size();++i)
       for(int j(0);j<metapop.migrationRatesF[i].size();++j)
	 pt.put("model.migration.F."+to_string(i)+"."+to_string(j),metapop.migrationRatesF[i][j]);
   if(metapop.migrationRatesM.size()>0)
     for(int i(0);i<metapop.migrationRatesM.size();++i)
       for(int j(0);j<metapop.migrationRatesM[i].size();++j)
	 pt.put("model.migration.M."+to_string(i)+"."+to_string(j),metapop.migrationRatesM[i][j]);
   pt.put("model.mu.mt.transition",mu.muMt1);
   pt.put("model.mu.mt.transversion",mu.muMt2);
   pt.put("model.mu.x.transition",mu.muX1);
   pt.put("model.mu.x.transversion",mu.muX2);
   pt.put("model.mu.y.transition",mu.muY1);
   pt.put("model.mu.y.transversion",mu.muY2);
   pt.put("model.mu.a.transition",mu.muA1);
   pt.put("model.mu.a.transversion",mu.muA2);

   pt.put("model.dna.sizeX",pdna.sizePerX);
   pt.put("model.dna.sizeA",pdna.sizePerA);
   pt.put("model.dna.sizeY",pdna.sizeY);
   pt.put("model.dna.sizeMT",pdna.sizeMt);
   pt.put("model.dna.nbLociA",pdna.nbChromosomesA);
   pt.put("model.dna.nbLociX",pdna.nbChromosomesX);
      
   BOOST_FOREACH( popVal  pop,populations) {
     ptree & node = pt.add("model.population", "");
     node.put("popsize", pop.popsize);
     node.put("inbreeding", pop.inbreeding);
     node.put("varianceNbChildren", pop.varChild);
     node.put("polygyny", pop.polygamy);
     node.put("polyandry", pop.polygamyf);
     node.put("demography.a", pop.growthFactor[0]);
     node.put("demography.b", pop.growthFactor[1]);
     node.put("demography.c", pop.growthFactor[2]);
     for(unsigned j(0);j<pop.matingSystem.size();++j)
       {
	 node.add("MS",pop.matingSystem[j]);
       }
   }   
   // Write the property tree to the XML file.
   boost::property_tree::xml_writer_settings<char> settings('\t', 1);
   write_xml(filename, pt,std::locale(),settings);
}
