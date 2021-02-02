#ifndef CONFIGURE_H
#define CONFIGURE_H

#include <string>
#include <fstream>
#include <string.h>
#include "namespace.h"

using namespace deardia;

std::string deardia::ConfigureValue(std::string value)
{
	std::string v = "";
	for (int i = 0; i < value.length(); ++i)
	{
		if (value[i] == '\t' || value[i] == ' ' || value[i] == '#') break;
		v += value[i];
	}
	return v;
}

void deardia::ReadConfigureFile(int argc, char* argv[])
{
	
	if (argc == 1)
	{
		printf(deardia::SOFTWARE_INFO.c_str());
		exit(0);
	}

	if (argc == 2)
	{
		if (!strcmp(argv[1], "--help") ||
			!strcmp(argv[1], "--h") ||
			!strcmp(argv[1], "-help") ||
			!strcmp(argv[1], "-h"))
		{
			printf(deardia::SOFTWARE_INFO.c_str());
			exit(0);
		}
		else if (!strcmp(argv[1], "--p"))
		{
			deardia::GernerateDefaultConfigureFile();
			printf("\nCreated:  deardia.configure.new\n");
			exit(0);
		}
		else
		{
			printf("\nPlease use ./deardia --help.\n");
			exit(0);
		}
	}
	
	int cond = -1;
	deardia::configureFile = "";
	deardia::outputDir = "";
	
	if (argc > 2)
	{
		
		for (int i = 1; i < argc; ++i)
		{
			if (!strcmp(argv[i], "--cf"))
			{
				cond = 1;
				continue;
			}
			if (!strcmp(argv[i], "--in"))
			{
				cond = 2;
				continue;
			}
			if (!strcmp(argv[i], "--out_dir"))
			{
				cond = 3;
				continue;
			}

			switch (cond)
			{

			case 1:

				deardia::configureFile = argv[i];
				break;
			case 2:
				deardia::inputFile = argv[i];
				break;
			case 3:
				deardia::outputDir = argv[i];

				if (deardia::outputDir.back() != '/')
				{
					deardia::outputDir += '/';
				}
				
				break;
			default:
				break;
			}
		}

		if (deardia::configureFile.length() == 0)
		{
			printf("\nThere is no input parameter of configure file: --cf\n");
			printf(deardia::SOFTWARE_INFO.c_str());
			exit(0);
		}
		if (deardia::outputDir.length() == 0)
		{
			printf("\nThere is no input parameter of output directory: --out_dir\n");
			printf(deardia::SOFTWARE_INFO.c_str());
			exit(0);
		}
		if (deardia::inputFile.length() == 0)
		{
			printf("\nThere is no input parameter of input file: --in\n");
			printf(deardia::SOFTWARE_INFO.c_str());
			exit(0);
		}
		
		std::ifstream ConfigureFile(deardia::configureFile);

		if (!ConfigureFile)
		{
			printf("Cannot open %s \nNo such file or directory.\n", deardia::configureFile.c_str());
			exit(0);
		}

		std::string line;

		while (std::getline(ConfigureFile, line))
		{
			if (line.empty()) continue;
			
			if (line[0] == '#') continue;

			while (line.back() == '\r') line.pop_back(); // avoid windows file

			int pos1 = line.find(" =");

			if (pos1 == line.npos) continue;

			std::string key, value;

			key = line.substr(0, pos1);
			value = line.substr(pos1 + 3);

			if (key == "database_name")
			{
				deardia::fastaFile = value;
			}
			if (key == "window_file_name")
			{
				deardia::windowFile = value;
			}
			if (key == "mxnet_json_file_name")
			{
				deardia::mxnetJsonFile = value;
			}
			if (key == "mxnet_params_file_name")
			{
				deardia::mxnetParamFile = value;
			}

			if (key == "GPU")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%d", &deardia::useGPU);
			}

			if (key == "num_threads")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%d", &deardia::totalThreads);
			}

			if (key == "variable_mod")
			{
				int pos2 = value.find("@");
				float mod_mass;
				sscanf(value.substr(0, pos2).c_str(), "%f", &mod_mass);
				std::string mod_AA = value.substr(pos2 + 1);

				for (int i = 0; i < mod_AA.length(); ++i)
				{
					deardia::modInfo[mod_AA[i]] = mod_mass;
				}
			}

			if (key == "buffer_size")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%d", &deardia::mzxmlBufferSize);
			}
			if (key == "slider_strides")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%d", &deardia::sliderStrides);
			}
			if (key == "batch_size")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%d", &deardia::batchSize);
			}
			if (key == "n_clusters")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%d", &deardia::kClusters);
			}
			if (key == "min_hit_num")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%d", &deardia::minHitNumber);
			}
			if (key == "precursor_min_charge")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%d", &deardia::ms1MinCharge);
			}
			if (key == "precursor_max_charge")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%d", &deardia::ms1MaxCharge);
			}
			if (key == "fragment_min_charge")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%d", &deardia::ms2MinCharge);
			}
			if (key == "fragment_max_charge")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%d", &deardia::ms2MaxCharge);
			}
			if (key == "fragment_min_mz")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%f", &deardia::ms2MinMz); 
			}
			if (key == "fragment_max_mz")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%f", &deardia::ms2MaxMz);
			}
			if (key == "pep_min_length")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%d", &deardia::pepMinLength);
			}
			if (key == "pep_max_length")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%d", &deardia::pepMaxLength);
			}
			if (key == "miss_cleavage")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%d", &deardia::missCleavage);
			}
			if (key == "ms1_tol_ppm")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%f", &deardia::ms1TolPPM);
				
				deardia::ms1TolPPM *= 1e-6;
			}
			if (key == "ms2_tol_ppm")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%f", &deardia::ms2TolPPM);
				
				deardia::ms2TolPPM *= 1e-6;
			}
			if (key == "slider_queue_size")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%d", &deardia::sliderQueueSize);
			}
			if (key == "cluster_queue_size")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%d", &deardia::clusterQueueSize);
			}
			if (key == "deeplearning_queue_size")
			{
				sscanf(deardia::ConfigureValue(value).c_str(), "%d", &deardia::deepLearningQueueSize);
			}
			if (key == "decoy_prefix")
			{
				deardia::decoyPrefix = deardia::ConfigureValue(value);
			}
		}
	}
}

void deardia::GernerateDefaultConfigureFile()
{
	std::ofstream ConfigureFile("deardia.configure.new");

	if (!ConfigureFile)
	{
		printf("Cannot write default configure file. Please cheak your current dirctory.\n");
		exit(0);
	}

	ConfigureFile << "# deardia version v" + deardia::SOFTWARE_VERSION << '\n'
		<< "# deardia parameters file. " << '\n'
		<< "# Everything following the '#' symbol means comment." << '\n'
		<< '\n'
		<< "database_name = /path/to/db.fasta" << '\n'
		<< "window_file_name = /path/to/win.100.tsv" << '\n'
		<< "mxnet_json_file_name = /path/to/deep_vae_mxnet-symbol.json" << '\n'
		<< "mxnet_params_file_name = /path/to/deep_vae_mxnet.params" << '\n'
		<< '\n'
		<< "GPU = 0 \t#1: use GPU acceleration; 0: not use GPU device." << '\n'
		<< "num_threads = 0 \t#0=use all CPU threads; else specify num threads directly." << '\n'
		<< "# Up to 3 variable modifications per peptide are supported. " << '\n'
		<< "# format <mass>@<residues>" << '\n'
		<< "# Note: 42.0106@n specify represent Acetylation" << '\n'
		<< "variable_mod = 42.0106@n" << '\n'
		<< "variable_mod = 15.9949@M" << '\n'
		<< "#variable_mod = 79.96633@STY" << '\n'
		<< '\n'
		<< "buffer_size = 4096\t#buffer space size when loading mzXML file." << '\n'
		<< "slider_strides = 1\t#overlap of neighboring sliders." << '\n'
		<< "batch_size = 512\t#batch size of deep auto-encoder inputs." << '\n'
		<< "n_clusters = 20\t#number of clusters in k-means." << '\n'
		<< "min_hit_num = 6\t#hit b ions + hit y ions >= 6." << '\n'
		<< "precursor_min_charge = 2\t#precursor min charge." << '\n'
		<< "precursor_max_charge = 4\t#precursor max charge." << '\n'
		<< "fragment_min_charge = 1\t#minimum charge of fragment." << '\n'
		<< "fragment_max_charge = 2\t#maximum charge of fragment." << '\n'
		<< "fragment_min_mz = 100\t#fragment min m/z" << '\n'
		<< "fragment_max_mz = 1800\t#fragment max m/z" << '\n'
		<< "pep_min_length = 7\t#minimum length of peptides" << '\n'
		<< "pep_max_length = 40\t#maximum length of peptides" << '\n'
		<< "miss_cleavage = 2\t#miss cleavage of enzyme digestion " << '\n'
		<< "ms1_tol_ppm = 50\t#tolerance of ms1 m/z (unit: ppm)" << '\n'
		<< "ms2_tol_ppm = 100\t#tolerance of ms2 m/z (unit: ppm)" << '\n'
		<< "decoy_prefix = DECOY_\t#prefix of decoy proteins in fasta database" << '\n'
		<< '\n'
		<< "#Advance options" << '\n'
		<< "#Warning: larger queue size will cost more memory." << '\n'
		<< "slider_queue_size = 16\t#the size of slider concurrency queue." << '\n'
		<< "cluster_queue_size = 512\t#the size of clustering concurrency queue" << '\n'
		<< "deeplearning_queue_size = 256\t#the size of deeplearning concurrency queue" << '\n';

	ConfigureFile.close();
}

#endif