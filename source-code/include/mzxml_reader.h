#ifndef VENDORREADER_H
#define VENDORREADER_H

#include <string>
#include <vector>

#include "../lib/zlib.h"
#include "../lib/Base64.h"

#if defined _WIN32
#pragma comment(lib, "./lib/zlibstat.lib")
#endif

#include "namespace.h"

//=========================================================================
// Parser mzXML file (open source format)
//=========================================================================
void deardia::ParserMzXMLFile()
{
	FILE* xml_file;
	xml_file = fopen(deardia::inputFile.c_str(), "r");

	char* buffer_char = (char*)malloc((deardia::mzxmlBufferSize + 1) * sizeof(char));
	buffer_char[deardia::mzxmlBufferSize] = 0;
	
	int pos;
	bool flag = false;
	int slider_num;
	int slider_id = 0;
	int cycle_id = -1;
	int middle_id = -1;
	int middle_rt = (int)(deardia::sliderWidth / 2);

	int strides = deardia::sliderStrides;
	int slider_width = deardia::sliderWidth;

	std::string name = "";
	std::string contents = "";

	float rt;
	int scanCount;
	int precision;
	int peaksCount;
	std::string zlib = "";
	std::string mslevel = "";
	std::string peak_string = "";
	std::string precursor_mz = "";

	std::vector<deardia::SpectrumData> ms1_spectra;
	std::vector< std::vector<deardia::SpectrumData> > ms2_spectra(slider_width);

	deardia::SpectraData spectra_data;

	while (1)
	{
		if (feof(xml_file)) break;
		size_t len = fread(buffer_char, 1, deardia::mzxmlBufferSize, xml_file);

		for (int i = 0; i < deardia::mzxmlBufferSize; ++i)
		{
			if (buffer_char[i] < ' ') buffer_char[i] = ' ';

			if (!flag && buffer_char[i] == '<')
			{
				flag = true;
			}
			else if (flag && buffer_char[i] != '<' && buffer_char[i] != '>')
			{
				contents += buffer_char[i];
			}
			else if (!flag && buffer_char[i] != '<' && buffer_char[i] != '>')
			{
				if (name == "precursorMz")
				{
					precursor_mz += buffer_char[i];
				}
				else if (name == "peaks")
				{
					peak_string += buffer_char[i];
				}
			}
			else if (flag && buffer_char[i] == '>')
			{
				flag = false;
				contents += buffer_char[i];

				pos = contents.find(' ');
				name = contents.substr(0, pos);

				if (name == "msRun")
				{
					sscanf(deardia::FindKeyValue(contents, "scanCount").c_str(), "%d", &scanCount);
					std::cout << "scanCount: " << scanCount << std::endl;
				}

				if (name == "scan")
				{
					mslevel = deardia::FindKeyValue(contents, "msLevel");
					sscanf(deardia::FindKeyValue(contents, "peaksCount").c_str(), "%d", &peaksCount);

					if (mslevel == "1")
					{
						++cycle_id;
						++middle_id;

						std::string rt_str = FindKeyValue(contents, "retentionTime");
						rt = atof(rt_str.substr(2, rt_str.length() - 3).c_str());
					}
				}
				else if (name == "peaks")
				{
					zlib = FindKeyValue(contents, "compressionType");
					sscanf(FindKeyValue(contents, "precision").c_str(), "%d", &precision);
				}

				if (name == "/peaks>")
				{
					if (cycle_id == -1) continue;

					if (mslevel == "1")
					{
						deardia::SpectrumData spectrumdata;

						spectrumdata.zlib = zlib;
						spectrumdata.peaksCount = peaksCount;
						spectrumdata.precision = precision / 8;
						spectrumdata.peakString = peak_string;

						ms1_spectra.push_back(spectrumdata);
					}
					if (mslevel == "1" && middle_id == middle_rt)
					{
						spectra_data.retentionTime = rt;
					}
					else if (mslevel == "2")
					{
						deardia::SpectrumData spectrumdata;

						spectrumdata.zlib = zlib;
						spectrumdata.peaksCount = peaksCount;
						spectrumdata.precision = precision / 8;
						spectrumdata.peakString = peak_string;

						float precursormz;

						sscanf(precursor_mz.c_str(), "%f", &precursormz);

						for (int n = 0; n < deardia::winRange.size(); ++n)
						{
							if (precursormz >= deardia::winRange[n][0] &&
								precursormz <= deardia::winRange[n][1])
							{
								spectrumdata.winId = n;
								break;
							}
						}

						ms2_spectra[cycle_id].push_back(spectrumdata);
					}

					peak_string.clear();
					precursor_mz.clear();
				}

				if (name == "/peaks>" && cycle_id == slider_width)
				{
					cycle_id -= strides;
					middle_id = middle_rt - strides;

					spectra_data.ms1Spectra.assign(ms1_spectra.begin(), ms1_spectra.begin() + slider_width);
					spectra_data.ms2Spectra.assign(ms2_spectra.begin(), ms2_spectra.begin() + slider_width);

					++slider_id;

					std::cout << slider_id << std::endl;

					if (slider_id < deardia::sliderQueueSize)
					{
						deardia::sliderQueue.enqueue(spectra_data);
					}

					while (slider_id >= deardia::sliderQueueSize)
					{
						if (process_slider_done.load(std::memory_order_acquire) > 0)
						{
							deardia::sliderQueue.enqueue(spectra_data);
							process_slider_done.fetch_add(-1, std::memory_order_release);
							break;
						}
						else if (process_slider_done.load(std::memory_order_acquire) == 0)
						{
							continue;
						}
					}

					ms1_spectra.erase(ms1_spectra.begin(), ms1_spectra.begin() + strides);
					ms2_spectra.erase(ms2_spectra.begin(), ms2_spectra.begin() + strides);
					ms2_spectra.resize(slider_width);
				}
				contents.clear();
			}
		}
	}

	fclose(xml_file);
	free(buffer_char);

	for (int i = 0; i < deardia::analyzeSliderThreads; ++i)
	{
		deardia::SpectraData spectra_data;
		spectra_data.stop = true;
		deardia::sliderQueue.enqueue(spectra_data);
	}
}

void deardia::SpectrumData::DecodeMzXMLData()
{
	spectrum.resize(2);
	spectrum[0].resize(peaksCount);
	spectrum[1].resize(peaksCount);

	unsigned char* uncomp_result;
	size_t size = textToBinarySize(peakString.length());

	if (zlib == "zlib")
	{
		unsigned char* result = (unsigned char*)malloc(size * sizeof(unsigned char));
		textToBinary(peakString.c_str(), peakString.length(), result);

		unsigned long retLen = 2 * peaksCount * precision;
		uncomp_result = (unsigned char*)malloc(retLen * sizeof(unsigned char));

		int OK = uncompress(uncomp_result, &retLen, result, size);

		free(result);

		if (OK != 0)
		{
			std::cerr << "Something was wrong in decode peak data.\n\n";
			free(uncomp_result);
			exit(3);
		}
	}
	else
	{
		uncomp_result = (unsigned char*)malloc(size * sizeof(unsigned char));
		textToBinary(peakString.c_str(), peakString.length(), uncomp_result);
	}

	int p = 0;
	for (int i = 0; i < peaksCount; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			int beg = p;
			int end = p + precision - 1;
			unsigned char tmp = 0;

			while (beg < end)
			{
				tmp = uncomp_result[beg];
				uncomp_result[beg] = uncomp_result[end];
				uncomp_result[end] = tmp;

				++beg;
				--end;
			}
			p += precision;
		}
	}

	if (precision == 4)
	{
		float* arr = (float*)uncomp_result;
		for (int i = 0; i < peaksCount; ++i)
		{
			spectrum[0][i] = arr[2 * i];
			spectrum[1][i] = arr[2 * i + 1];
		}
	}
	else
	{
		double* arr = (double*)uncomp_result;
		for (int i = 0; i < peaksCount; ++i)
		{
			spectrum[0][i] = arr[2 * i];
			spectrum[1][i] = arr[2 * i + 1];
		}
	}
	free(uncomp_result);
}

#endif