#ifndef NAMESPACE_H
#define NAMESPACE_H

#include <atomic>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>

#include "../lib/zlib.h"
#include "../lib/Base64.h"
#include "../lib/blockingconcurrentqueue.h"

//=========================================================================
// Globle variables
//=========================================================================
std::atomic<int> cluster_done(0);
std::atomic<int> deeplearning_done(0);
std::atomic<int> process_slider_done(0);
std::atomic<bool> init_network_done(false);

namespace deardia
{
	//=========================================================================
	// Information of software
	//=========================================================================
	std::string SOFTWARE_VERSION = "1.0.0";
	std::string SOFTWARE_INFO =
		"\n Welcome to use deardia v" + SOFTWARE_VERSION + " software.\n" +
		"\n Usage: deardia.exe --cf deardia.configure --out_dir /path/to/dir/ --in xx.mzXML\n" +
		"\nIntroduction of args:\n" + 
		"\n--cf:\tconfigure file of deardia.\n" +
		"\n--out_dir:\toutput directory or path.\n" +
		"\n--in:\tinput file, deardia supports *.mzXML, *.wiff, *.raw, *.d/ .\n" +
		"\n If there is no deardia.configure, please run the following command:\n" +
		"\n\t ./deardia --p\n";

	//=========================================================================
	// Default parameters
	//=========================================================================
	int vendorFile = 0;
	
	int useGPU = 0;
	int totalThreads = 0;

	int clusterThreads;
	int parserFastaThreads;
	int analyzeSliderThreads;

	int sliderWidth = 20;
	int latentFeaturesDim = 16;

	float ms1BinSize = 30.0f;
	float ms2BinSize = 30.0f;

	int sliderQueueSize;
	int clusterQueueSize;
	int deepLearningQueueSize;

	int sliderStrides;
	int mzxmlBufferSize;

	std::vector< std::vector<float> > winRange;

	int batchSize;
	int kClusters;

	// Parameters of PIndex
	float ms1TolPPM;
	float ms2TolPPM;
	
	int numFragmentBin;
	int numPrecursorBin;

	int minHitNumber;

	int ms1MinCharge;
	int ms1MaxCharge;

	int ms2MinCharge;
	int ms2MaxCharge;

	float ms2MinMz;
	float ms2MaxMz;

	int pepMinLength;
	int pepMaxLength;

	int missCleavage;

	std::string decoyPrefix;

	std::unordered_map<char, float> modInfo;

	// setting of directories
	std::string inputFile;
	std::string outputDir;

	std::string mgfFile;
	std::string fastaFile;
	std::string windowFile;

	std::string configureFile;
	std::string mxnetJsonFile;
	std::string mxnetParamFile;

	//=========================================================================
	// Data structures
	//=========================================================================
	struct DigestionInfo
	{
		bool stop = false;

		std::vector<int> byIons;
		std::vector< std::vector<int> > ms1Info;
	};

	struct SpectrumData
	{
		int winId;

		int precision;
		int peaksCount;

		std::string zlib;
		std::string peakString;
		std::vector< std::vector<float> > spectrum;

		void DecodeMzXMLData();
	};

	struct SpectraData
	{
		bool stop = false;

		float retentionTime;

		std::vector< SpectrumData > ms1Spectra;
		std::vector< std::vector<SpectrumData> > ms2Spectra;
	};

	struct SliderData
	{
		bool stop = false;

		int winId;
		float retentionTime;
		
		std::vector<float> kurtosis;

		std::vector<float> ms2Mz;
		std::vector<float> ms2Intensity;

		std::vector< std::vector<float> > ms1XICs;
		std::vector< std::vector<float> > ms2XICs;

		std::vector<float> ms1Mz;
		std::vector<int> ms1Charge;

		std::vector< std::vector<float> > ms2XICFeatures;
	};

	struct MgfData
	{
		bool stop = false;

		float ms1Mz;
		int ms1Charge;
		float retentionTime;

		std::vector< std::vector<float> > ms2Spectrum;
	};

	//=========================================================================
	// Global variables
	//=========================================================================
	std::vector< std::vector<int> > ppidMapToFragment;
	std::vector< std::vector< std::vector<int> > > ms1MapToPpid;

	moodycamel::BlockingConcurrentQueue<std::string>* proteinQueue
		= new moodycamel::BlockingConcurrentQueue<std::string>();

	moodycamel::BlockingConcurrentQueue<DigestionInfo>* peptideQueue
		= new moodycamel::BlockingConcurrentQueue<DigestionInfo>();

	moodycamel::BlockingConcurrentQueue<SpectraData> sliderQueue;

	moodycamel::BlockingConcurrentQueue<SliderData> deepLearningQueue;

	moodycamel::BlockingConcurrentQueue<SliderData> clusterQueue;

	moodycamel::BlockingConcurrentQueue<MgfData> mgfQueue;

	//=========================================================================
	// Utilities
	//=========================================================================
	void LoadWindowFile();

	std::string FindKeyValue(const std::string& in_str, const std::string& key);
	std::vector<std::string> GetInputFileNameAndSuffix(std::string inputFileName);

	float Skewness(std::vector<float>& x);
	float Kurtosis(std::vector<float>& x);
	float KL_divergence(std::vector<float>& x, std::vector<float>& y);
	float PearsonCorrcoef(std::vector<float>& X, std::vector<float>& Y);
	float AutoCorrelation(std::vector<float>& x);

	class MXNetBufferFile
	{
	public:

		std::string file_path_;
		int length_;

		char* buffer_;

		explicit MXNetBufferFile(std::string file_path)
			:file_path_(file_path)
		{
			std::ifstream ifs(file_path.c_str(), std::ios::in | std::ios::binary);
			if (!ifs)
			{
				std::cerr << "Can't open the parameters file of deep neural network. Please cheak" << file_path << ". \n";
				length_ = 0;
				buffer_ = NULL;
				return;
			}

			ifs.seekg(0, std::ios::end);
			length_ = ifs.tellg();
			ifs.seekg(0, std::ios::beg);

			buffer_ = new char[sizeof(char) * length_];
			ifs.read(buffer_, length_);
			ifs.close();
		}

		int GetLength()
		{
			return length_;
		}

		char* GetBuffer()
		{
			return buffer_;
		}

		~MXNetBufferFile()
		{
			if (buffer_)
			{
				delete[] buffer_;
				buffer_ = NULL;
			}
		}
	};

	//=========================================================================
	// Major functions
	//=========================================================================
	void MultiThreadsRun(int argc, char* argv[]);

	void GernerateDefaultConfigureFile();
	std::string ConfigureValue(std::string value);
	void ReadConfigureFile(int argc, char* argv[]);

	void ParserFasta();
	void LoadFastaFile();
	void CreatInvertedIndex();

	void ParserMzXMLFile();
	
	void WriteMGF();
	void DeepLearning();
	void ProcessSliders();
	void ClusteringAndSearchingDatabase();

	void Deisotopes(std::vector<SpectrumData>& ms1Spectra,
		std::vector< std::vector<float> >& ms1Mz,
		std::vector< std::vector<int> >& ms1Charge,
		std::vector< std::vector< std::vector<float> > >& ms1XICs);
}
#endif