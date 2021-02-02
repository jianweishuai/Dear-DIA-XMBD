#ifndef MULTITHREADS_H
#define MULTITHREADS_H

#include <thread>

#include "core.h"
#include "utility.h"
#include "configure.h"
#include "namespace.h"
#include "parser_fasta.h"
#include "mzxml_reader.h"

using namespace deardia;

void deardia::MultiThreadsRun(int argc, char* argv[])
{
	//=========================================================================
	// Load configure file
	//=========================================================================
	deardia::ReadConfigureFile(argc, argv);
	deardia::LoadWindowFile();

	std::vector<std::string> nameSuffix = deardia::GetInputFileNameAndSuffix(deardia::inputFile);
	deardia::mgfFile = deardia::outputDir + nameSuffix[0] + ".mgf";

	//=========================================================================
	// allocate threads
	//=========================================================================
	int nThreads = 0;

	if (deardia::totalThreads == 0)
	{
		nThreads = std::thread::hardware_concurrency() - 1;
	}
	else
	{
		nThreads = deardia::totalThreads;
	}

	deardia::parserFastaThreads = nThreads;
	deardia::analyzeSliderThreads = (int)(nThreads / 3);
	deardia::clusterThreads = nThreads - deardia::analyzeSliderThreads;

	//=========================================================================
	// Initialze deep neural network
	//=========================================================================
	std::thread deepLearning = std::thread([&]()
	{
		deardia::DeepLearning();
	});

	while (1)
	{
		if (init_network_done)
		{
			break;
		}
	}

	//=========================================================================
	// Load fasta file
	//=========================================================================
	std::thread laodFastaFile = std::thread([&]()
	{
		deardia::LoadFastaFile();
	});
	
	std::vector< std::thread > parserFasta(deardia::parserFastaThreads);

	for (int i = 0; i < deardia::parserFastaThreads; ++i)
	{
		parserFasta[i] = std::thread([&]()
		{
			deardia::ParserFasta();
		});
	}

	std::thread creatInvertedIndex = std::thread([&]()
	{
		deardia::CreatInvertedIndex();
	});
		
	laodFastaFile.join();
	for (int i = 0; i < deardia::parserFastaThreads; ++i)
	{
		parserFasta[i].join();
	}
	creatInvertedIndex.join();

	delete deardia::proteinQueue;
	deardia::proteinQueue = nullptr;

	delete deardia::peptideQueue;
	deardia::peptideQueue = nullptr;

	//=========================================================================
	// Parser input file
	//=========================================================================
	std::thread parserInputFile;

	if (nameSuffix[1] == "mzXML")
	{
		deardia::vendorFile = 0;

		parserInputFile = std::thread([&]()
		{
			deardia::ParserMzXMLFile();
		});
	}
	else
	{
		std::cout << "Dear-DIA-XMBD only support .mzXML format.\n" 
			<<"Dear-DIA-XMBD does not support this file: " << nameSuffix[0] + nameSuffix[1] << '\n'
			<< "Please cheak your file." << std::endl;

		exit(0);
	}

	//=========================================================================
	// Process sliders
	//=========================================================================
	std::vector<std::thread> processSliders(deardia::analyzeSliderThreads);

	for (int i = 0; i < deardia::analyzeSliderThreads; ++i)
	{
		processSliders[i] = std::thread([&]()
		{
			deardia::ProcessSliders();
		});
	}

	std::vector<std::thread> clusteringAndSearchingDatabase(deardia::clusterThreads);

	for (int i = 0; i < deardia::clusterThreads; ++i)
	{
		clusteringAndSearchingDatabase[i] = std::thread([&]()
		{
			deardia::ClusteringAndSearchingDatabase();
		});
	}

	std::thread writeMGF = std::thread([&]()
	{
		deardia::WriteMGF();
	});

	writeMGF.join();
	parserInputFile.join();

	for (int i = 0; i < deardia::analyzeSliderThreads; ++i)
	{
		processSliders[i].join();
	}

	for (int i = 0; i < deardia::clusterThreads; ++i)
	{
		clusteringAndSearchingDatabase[i].join();
	}

	deepLearning.join();
}


#endif