#ifndef FASTAPARSER_H
#define FASTAPARSER_H

#include <string>
#include <vector>
#include <numeric>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include "namespace.h"
using namespace deardia;

void deardia::LoadFastaFile()
{
	printf("Start loadding db.fasta file ...\n");
	
	std::ifstream DBFile(deardia::fastaFile.c_str());
	
	if (!DBFile)
	{
		printf("Cannot open db.fasta file. Please cheak %s\n", deardia::fastaFile.c_str());
		exit(0);
	}

	std::string line;
	std::string proteinStr = "";
	
	int numProtein = 0;
	int numDecoy = 0;

	while (std::getline(DBFile, line))
	{
		if (line.find(deardia::decoyPrefix) != line.npos) ++numDecoy;

		if (line.find('>') != line.npos)
		{
			if (proteinStr.length() > 0)
			{
				++numProtein;
				deardia::proteinQueue->enqueue(proteinStr);
			}
			proteinStr.clear();
		}
		else if (line.find('>') == line.npos)
		{
			while (line.back() == '\n' || line.back() == '\r')
			{
				line.pop_back();
			}

			proteinStr += line;
		}
	}
	++numProtein;
	deardia::proteinQueue->enqueue(proteinStr);

	for (int t = 0; t < deardia::parserFastaThreads; ++t)
	{
		deardia::proteinQueue->enqueue("stop");
	}

	DBFile.close();

	std::cout << "Load fasta file done." << '\n'
		<< "This database includes " << numProtein - numDecoy
		<< " proteins and " << numDecoy << " decoy proteins." << '\n'
		<< "Start creating inverted index of MS1 and MS2 ..." << '\n';
}

void deardia::ParserFasta()
{
	std::vector<float> AminoAcidsMAss = {	89.047679f, // A
											0.0f, 		 // B
											121.019751f + 57.021464f, // C
											133.037509f, // D
											147.053159f, // E
											165.078979f, // F
											75.032029f,  // G
											155.069477f, // H
											131.094629f, // I
											0.0f, 		 // J
											146.105528f, // K
											131.094629f, // L
											149.051051f, // M
											132.053493f, // N
											0.0f, 		 // O
											115.063329f, // P
											146.069143f, // Q
											174.111676f, // R
											105.042594f, // S
											119.058244f, // T
											0.0f,		 // U
											117.078979f, // V
											204.089878f, // W
											0.0f, 		 // X
											181.073894f, // Y
											0.0f 		 // Z
										};

	float H2O = 18.010565f;
	float PROTON = 1.007825f;
	std::unordered_map<char, float> modInfo = deardia::modInfo;

	std::string proteinStr;

	while (1)
	{
		deardia::proteinQueue->wait_dequeue(proteinStr);
		
		if (proteinStr == "stop") break;
		
		if (proteinStr[0] == 'M') proteinStr = proteinStr.substr(1);
		
		std::vector<int> splitPos;
		splitPos.push_back(-1);

		for (int i = 0; i < proteinStr.length(); ++i)
		{
			if (proteinStr[i] == 'K' || proteinStr[i] == 'R')
			{
				splitPos.push_back(i);
			}
			if (proteinStr[i] == 'P')
			{
				if (!splitPos.empty())
				{
					if (splitPos.back() + 1 == i)
					{
						splitPos.pop_back();
					}
				}
			}
		}

		splitPos.push_back(proteinStr.length() - 1);

		std::vector< std::vector<int> > cutPos;

		for (int i = 0; i < splitPos.size() - 1; ++i)
		{
			for (int j = i + 1; j <= i + 1 + deardia::missCleavage; ++j)
			{
				if (j < splitPos.size())
				{
					std::vector<int> tmp(2, 0);
					tmp[0] = splitPos[i] + 1;
					tmp[1] = splitPos[j] - tmp[0] + 1;

					if (tmp[1] < deardia::pepMinLength ||
						tmp[1] > deardia::pepMaxLength) continue;

					cutPos.push_back(tmp);
				}
			}
		}

		std::string peptideStr = "";
		
		for (int i = 0; i < cutPos.size(); ++i)
		{
			peptideStr = proteinStr.substr(cutPos[i][0], cutPos[i][1]);
			
			std::vector<int> modPos;
			std::vector<float> origAcidsMass;

			for (int j = 0; j < peptideStr.length(); ++j)
			{
				if (modInfo.find(peptideStr[j]) != modInfo.end())
				{
					modPos.push_back(j);
				}

				origAcidsMass.push_back(AminoAcidsMAss[peptideStr[j] - 'A']);
			}

			std::vector< std::vector<float> > modAcidsMass;

			modAcidsMass.push_back(origAcidsMass); // 0 modification

			if (modPos.size() >= 1)// 1 modification
			{
				for (int j = 0; j < modPos.size(); ++j)
				{
					std::vector<float> modAcids = origAcidsMass;
					modAcids[modPos[j]] += modInfo[peptideStr[modPos[j]]];
					modAcidsMass.push_back(modAcids);
				}
			}

			if (modPos.size() >= 2)// 2 modifications
			{
				for (int j = 0; j < modPos.size() - 1; ++j)
				{
					for (int k = j + 1; k < modPos.size(); ++k)
					{
						std::vector<float> modAcids = origAcidsMass;
						modAcids[modPos[j]] += modInfo[peptideStr[modPos[j]]];
						modAcids[modPos[k]] += modInfo[peptideStr[modPos[k]]];

						modAcidsMass.push_back(modAcids);
					}
				}
			}

			if (modPos.size() >= 3)// 2 modifications
			{
				for (int j = 0; j < modPos.size() - 2; ++j)
				{
					for (int k = j + 1; k < modPos.size() - 1; ++k)
					{
						for (int m = k + 1; m < modPos.size(); ++m)
						{
							std::vector<float> modAcids = origAcidsMass;
							modAcids[modPos[j]] += modInfo[peptideStr[modPos[j]]];
							modAcids[modPos[k]] += modInfo[peptideStr[modPos[k]]];
							modAcids[modPos[m]] += modInfo[peptideStr[modPos[m]]];

							modAcidsMass.push_back(modAcids);
						}
					}
				}
			}

			float origPepMass = -((float)peptideStr.length() - 1.0f) * H2O;
			
			for (int n = 0; n < modAcidsMass.size(); ++n)
			{
				float pepMass = origPepMass + std::accumulate(modAcidsMass[n].begin(), modAcidsMass[n].end(), 0.0f);

				deardia::DigestionInfo digestionInfo;

				digestionInfo.ms1Info.resize(2);

				for (int z = deardia::ms1MinCharge; z <= deardia::ms1MaxCharge; ++z)
				{
					float ms1Mz = pepMass / (float)z + PROTON;

					for (int k = 0; k < deardia::winRange.size(); ++k)
					{
						if (ms1Mz >= deardia::winRange[k][0] && ms1Mz <= deardia::winRange[k][1])
						{
							int ms1BinMz = (int)(ms1Mz * deardia::ms1BinSize);

							digestionInfo.ms1Info[0].push_back(ms1BinMz);
							digestionInfo.ms1Info[1].push_back(z - deardia::ms1MinCharge);

							break;
						}
					}
				}

				if (digestionInfo.ms1Info[0].size() == 0) continue;

				std::unordered_set<int> skipIons; // skip indistinguishable ions

				for (int z = deardia::ms2MinCharge; z <= deardia::ms2MaxCharge; ++z)
				{
					float pepMass = origPepMass + std::accumulate(modAcidsMass[n].begin(), modAcidsMass[n].end(), 0.0f);
					float sumMass = modAcidsMass[n][0];

					float yFirstIonMz = (pepMass + z * PROTON) / z;

					if (yFirstIonMz >= deardia::ms2MinMz &&
						yFirstIonMz <= deardia::ms2MaxMz)
					{
						int iyFirstIonMz = (int)(yFirstIonMz * deardia::ms2BinSize);

						if (skipIons.find(iyFirstIonMz) == skipIons.end())
						{
							skipIons.insert(iyFirstIonMz);
							digestionInfo.byIons.push_back(iyFirstIonMz);
						}
					}

					for (int j = 1; j < peptideStr.length(); ++j)
					{
						float bIonMass = sumMass - H2O + PROTON;
						float yIonMass = pepMass - bIonMass + 2.0f * PROTON;

						float bIonMz = (bIonMass + PROTON * (z - 1)) / z;
						float yIonMz = (yIonMass + PROTON * (z - 1)) / z;

						if (bIonMz >= deardia::ms2MinMz && 
							bIonMz <= deardia::ms2MaxMz)
						{
							int ibIonMz = (int)(bIonMz * deardia::ms2BinSize);

							if (skipIons.find(ibIonMz) == skipIons.end())
							{
								skipIons.insert(ibIonMz);
								digestionInfo.byIons.push_back(-ibIonMz);
							}
						}

						if (yIonMz >= deardia::ms2MinMz && 
							yIonMz <= deardia::ms2MaxMz)
						{
							int iyIonMz = (int)(yIonMz * deardia::ms2BinSize);

							if (skipIons.find(iyIonMz) == skipIons.end())
							{
								skipIons.insert(iyIonMz);
								digestionInfo.byIons.push_back(iyIonMz);
							}
						}

						sumMass += modAcidsMass[n][j] - H2O;
					}

					float bEndIonMz = (sumMass - H2O + z * PROTON) / z;

					if (bEndIonMz >= deardia::ms2MinMz && 
						bEndIonMz <= deardia::ms2MaxMz)
					{
						int ibEndIonMz = (int)(bEndIonMz * deardia::ms2BinSize);

						if (skipIons.find(ibEndIonMz) == skipIons.end())
						{
							skipIons.insert(ibEndIonMz);
							digestionInfo.byIons.push_back(-ibEndIonMz);
						}
					}
				}

				deardia::peptideQueue->enqueue(digestionInfo);
			}
			
		}
	}

	deardia::DigestionInfo digestionInfo;
	digestionInfo.stop = true;

	deardia::peptideQueue->enqueue(digestionInfo);
}

void deardia::CreatInvertedIndex()
{
	
	deardia::DigestionInfo digestionInfo;

	int ppid = 0;
	int stop_count = 0;

	deardia::numFragmentBin = (int)(deardia::ms2MaxMz * deardia::ms2BinSize) + 8;
	deardia::numPrecursorBin = (int)(deardia::winRange.back()[1] * deardia::ms1BinSize) + 8;

	int delta_z = deardia::ms1MaxCharge - deardia::ms1MinCharge + 1;
	
	deardia::ms1MapToPpid.resize(deardia::numPrecursorBin, std::vector< std::vector<int> >(delta_z));
	
	while (1)
	{
		deardia::peptideQueue->wait_dequeue(digestionInfo);

		if (digestionInfo.stop)
		{
			++stop_count;
			if (stop_count == deardia::parserFastaThreads) break;
		}
		else
		{
			for (int i = 0; i < digestionInfo.ms1Info[0].size(); ++i)
			{
				int mz = digestionInfo.ms1Info[0][i];
				int z  = digestionInfo.ms1Info[1][i];

				deardia::ms1MapToPpid[mz][z].push_back(ppid);
			}

			deardia::ppidMapToFragment.push_back(digestionInfo.byIons);

			++ppid;
		}
	}
	
	printf("Creating inverted index done. This database includes %d peptides.\n", ppid);
	
}

#endif