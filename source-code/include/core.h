#ifndef CORE_H
#define CORE_H

#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include "kmeanspp.h"
#include "findpeaksx.h"

#include "namespace.h"
#include "../lib/c_predict_api.h"

void deardia::Deisotopes(	std::vector<SpectrumData>& ms1Spectra,
							std::vector< std::vector<float> >& ms1Mz,
							std::vector< std::vector<int> >& ms1Charge,
							std::vector< std::vector< std::vector<float> > >& ms1XICs)
{
	int halfWidth = deardia::sliderWidth / 2;
	int numWindows = deardia::winRange.size();

	ms1Mz.resize(numWindows);
	ms1XICs.resize(numWindows);
	ms1Charge.resize(numWindows);

	if (deardia::vendorFile == 0)
	{
		for (int t = 0; t < deardia::sliderWidth; ++t)
		{
			ms1Spectra[t].DecodeMzXMLData();
		}
	}
	
	std::unordered_set<int> skipIons;

	for (int t = halfWidth - 3; t < halfWidth + 3; ++t)
	{
		std::vector<float> mzPeaks = FindPeaksSX::findpeaks(ms1Spectra[t].spectrum[0], ms1Spectra[t].spectrum[1]);

		if (mzPeaks.size() == 0) continue;

		for (int i = 0; i < mzPeaks.size() - 1; ++i)
		{
			if (mzPeaks[i] < deardia::winRange[0][0] ||
				mzPeaks[i] > deardia::winRange.back()[1])
			{
				continue;
			}

			int winId = -1;
			for (int n = 0; n < numWindows; ++n)
			{
				if (mzPeaks[i] > deardia::winRange[n][0] &&
					mzPeaks[i] < deardia::winRange[n][1])
				{
					winId = n;
					break;
				}
			}

			if (winId == -1) continue;

			bool isotope = false;

			for (int z = deardia::ms1MinCharge; z <= deardia::ms1MaxCharge; ++z)
			{
				if (isotope) break;

				float diff = 1.00287f / (float)z;

				for (int j = i + 1; j < mzPeaks.size(); ++j)
				{
					float mzError = mzPeaks[j] - mzPeaks[i] - diff;

					if (mzError >= -0.1f && mzError <= 0.1f)
					{
						int mzCharge = (int)(mzPeaks[i] * deardia::ms1BinSize) * 10 + z;

						if (skipIons.find(mzCharge) == skipIons.end())
						{
							skipIons.insert(mzCharge);

							float mzTolerance = mzPeaks[i] * deardia::ms1TolPPM;

							float leftMz  = mzPeaks[i] - mzTolerance;
							float rightMz = mzPeaks[i] + mzTolerance;

							if (rightMz > deardia::winRange.back()[1])
							{
								rightMz = deardia::winRange.back()[1];
							}

							std::vector<float> curve(deardia::sliderWidth, 0.0f);
							std::vector<float> mzTrace(deardia::sliderWidth, 0.0f);

							int nonZeroCount = 0;
							for (int tt = 0; tt < deardia::sliderWidth; ++tt)
							{
								int leftIndex = std::lower_bound(ms1Spectra[tt].spectrum[0].begin(),
								ms1Spectra[tt].spectrum[0].end(), leftMz) - ms1Spectra[tt].spectrum[0].begin();

								if (leftIndex == ms1Spectra[tt].spectrum[0].size()) continue;
								
								float maxMz = ms1Spectra[tt].spectrum[0][leftIndex];
								float maxIntensity = ms1Spectra[tt].spectrum[1][leftIndex];

								if (maxMz > rightMz) continue;

								for (int k = leftIndex; k < ms1Spectra[tt].spectrum[0].size(); ++k)
								{
									if (ms1Spectra[tt].spectrum[0][k] <= rightMz)
									{
										curve[tt] += ms1Spectra[tt].spectrum[1][k];

										if (ms1Spectra[tt].spectrum[1][k] > maxIntensity)
										{
											maxMz = ms1Spectra[tt].spectrum[0][k];
											maxIntensity = ms1Spectra[tt].spectrum[1][k];
										}
									}
									else
									{
										break;
									}
								}

								if (curve[tt] > 1.0f)
								{
									++nonZeroCount;
									mzTrace[tt] = maxMz;
								}
							}

							if (nonZeroCount < 5) continue;

							float detectedMz = 0.0f;
							float sumIntensity = 0.0f;

							for (int k = 0; k < deardia::sliderWidth; ++k)
							{
								sumIntensity += curve[k];
								detectedMz += mzTrace[k] * curve[k];
							}
							detectedMz = detectedMz / sumIntensity;

							int binMz = (int)(detectedMz * deardia::ms1BinSize);
							if (binMz >= deardia::numPrecursorBin) continue;

							ms1Mz[winId].push_back(detectedMz);
							ms1Charge[winId].push_back(z - deardia::ms1MinCharge);

							std::vector<float> scaleCurve(deardia::sliderWidth);
							auto m = std::minmax_element(curve.begin(), curve.end());

							for (int k = 0; k < deardia::sliderWidth; ++k)
							{
								scaleCurve[k] = (curve[k] - *m.first) / (*m.second - *m.first);
							}

							ms1XICs[winId].push_back(scaleCurve);
						}

						isotope = true;
						break;
					}
					else if (mzError > 0.1f)
					{
						break;
					}
				}
			}
		}
	}
}

#if 0
void deardia::ProcessSliders()
{
	int doneCount = 0;
	int numWindows = deardia::winRange.size();

	deardia::SpectraData spectraData;

	while (1)
	{
		deardia::sliderQueue.wait_dequeue(spectraData);

		if (spectraData.stop) break;

		std::vector< std::vector< std::vector<int> > > winDict(numWindows);

		for (int i = 0; i < spectraData.ms2Spectra.size(); ++i)
		{
			for (int j = 0; j < spectraData.ms2Spectra[i].size(); ++j)
			{
				int winId = spectraData.ms2Spectra[i][j].winId;

				std::vector<int> Idx = { i, j };
				winDict[winId].push_back(Idx);
			}
		}

		std::vector< std::vector<float> > ms1Mz;
		std::vector< std::vector<int> > ms1Charge;
		std::vector< std::vector< std::vector<float> > > ms1XICs;
		Deisotopes(spectraData.ms1Spectra, ms1Mz, ms1Charge, ms1XICs);

		for (int n = 0; n < numWindows; ++n)
		{
			if (winDict[n].size() < 5 || ms1Mz[n].size() < 1) continue;

			std::vector< std::vector<float> > ms2TempSpec(deardia::sliderWidth);

			for (int i = 0; i < winDict[n].size(); ++i)
			{
				int cycleId = winDict[n][i][0];
				int windowId = winDict[n][i][1];

				spectraData.ms2Spectra[cycleId][windowId].DecodeMzXMLData();

				std::vector<float> tempSpec;
				tempSpec.resize((int)(spectraData.ms2Spectra[cycleId][windowId].spectrum[0].back() * 32.0f), 0.0f);

				float element[2] = { 0.0f, 0.0f };
				for (int j = 0; j < spectraData.ms2Spectra[cycleId][windowId].spectrum[0].size(); ++j)
				{
					element[0] = spectraData.ms2Spectra[cycleId][windowId].spectrum[0][j];
					element[1] = spectraData.ms2Spectra[cycleId][windowId].spectrum[1][j];

					if (element[1] < 0.1f ||
						element[0] < deardia::ms2MinMz ||
						element[0] > deardia::ms2MaxMz)
					{
						continue;
					}

					int pos = (int)(element[0] * 30.0f);

					tempSpec[pos] += element[1];
					tempSpec[pos - 1] += element[1];
					tempSpec[pos + 1] += element[1];
				}
				ms2TempSpec[cycleId] = tempSpec;
			}

			int maxSize = 0;
			for (int i = 0; i < ms2TempSpec.size(); ++i)
			{
				if (ms2TempSpec[i].size() > maxSize)
				{
					maxSize = ms2TempSpec[i].size();
				}
			}

			deardia::SliderData sliderData;

			for (int i = 1; i < maxSize - 1; ++i)
			{
				if (i > (int)(deardia::winRange[n][0] * 30.0f) &&
					i < (int)(deardia::winRange[n][1] * 30.0f))
				{
					continue;
				}

				std::vector<float> curve(deardia::sliderWidth, 0.0f);

				int t = 0;
				int nonZeroCount = 0;

				float maxIntensity = 0.0f;
				float nonZeroMinIntensity = 0.0f;

				for (auto it = ms2TempSpec.begin(); it != ms2TempSpec.end(); ++it)
				{
					if (it->size() == 0)
					{
						++t;
						continue;
					}

					float temp[3] = { 0.0f, 0.0f, 0.0f };
					int t_size = it->size();

					if (i < t_size) temp[1] = (*it)[i];
					if (i - 1 < t_size) temp[0] = (*it)[i - 1];
					if (i + 1 < t_size) temp[2] = (*it)[i + 1];

					if (temp[0] > 1e-3 && temp[1] > 1e-3 && temp[2] > 1e-3)
					{
						curve[t] = temp[0] + temp[1] + temp[2];

						if (nonZeroCount == 0)
						{
							maxIntensity = curve[t];
							nonZeroMinIntensity = curve[t];
						}
						else
						{
							maxIntensity = std::max(maxIntensity, curve[t]);
							nonZeroMinIntensity = std::min(nonZeroMinIntensity, curve[t]);
						}

						++nonZeroCount;
					}

					++t;
				}

				if (nonZeroCount < 5) continue;
				if (maxIntensity / nonZeroMinIntensity < 4.0f) continue;

				float kurtosis = deardia::Kurtosis(curve);

				std::vector<float> scaleCurve(deardia::sliderWidth);
				auto m = std::minmax_element(curve.begin(), curve.end());
				float scale = *m.second - *m.first;

				for (int j = 0; j < deardia::sliderWidth; ++j)
				{
					scaleCurve[j] = (curve[j] - *m.first) / scale;
				}

				sliderData.ms2Mz.push_back(i);
				sliderData.kurtosis.push_back(kurtosis);
				sliderData.ms2XICs.push_back(scaleCurve);
				sliderData.ms2Intensity.push_back(*m.second);
			}

			if (sliderData.ms2Mz.size() < deardia::minHitNumber * deardia::kClusters)
			{
				continue;
			}

			sliderData.winId = n;
			sliderData.retentionTime = spectraData.retentionTime;

			sliderData.ms1Mz = ms1Mz[n];
			sliderData.ms1XICs = ms1XICs[n];
			sliderData.ms1Charge = ms1Charge[n];

			++doneCount;

			if (doneCount < deardia::deepLearningQueueSize)
			{
				deardia::deepLearningQueue.enqueue(sliderData);
			}

			while (doneCount >= deardia::deepLearningQueueSize)
			{
				if (deeplearning_done.load(std::memory_order_acquire) > 0)
				{
					deardia::deepLearningQueue.enqueue(sliderData);
					deeplearning_done.fetch_add(-1, std::memory_order_release);
					break;
				}
				else if (deeplearning_done.load(std::memory_order_acquire) == 0)
				{
					continue;
				}
			}
		}

		process_slider_done.fetch_add(1, std::memory_order_release);
	}

	deardia::SliderData sliderData;
	sliderData.stop = true;
	deardia::deepLearningQueue.enqueue(sliderData);
}
#endif

void deardia::ProcessSliders()
{
	int doneCount = 0;
	int numWindows = deardia::winRange.size();

	deardia::SpectraData spectraData;

	while (1)
	{
		deardia::sliderQueue.wait_dequeue(spectraData);

		if (spectraData.stop) break;

		std::vector< std::vector< std::vector<int> > > winDict(numWindows);

		for (int i = 0; i < spectraData.ms2Spectra.size(); ++i)
		{
			for (int j = 0; j < spectraData.ms2Spectra[i].size(); ++j)
			{
				int winId = spectraData.ms2Spectra[i][j].winId;

				std::vector<int> Idx = { i, j };
				winDict[winId].push_back(Idx);
			}
		}

		std::vector< std::vector<float> > ms1Mz;
		std::vector< std::vector<int> > ms1Charge;
		std::vector< std::vector< std::vector<float> > > ms1XICs;
		Deisotopes(spectraData.ms1Spectra, ms1Mz, ms1Charge, ms1XICs);

		for (int n = 0; n < numWindows; ++n)
		{
			if (winDict[n].size() < 5 || ms1Mz[n].size() < 1) continue;

			std::unordered_set<int> candidatePpId;
			for (int i = 0; i < ms1Mz[n].size(); ++i)
			{
				int z = ms1Charge[n][i];
				int binMz = (int)(ms1Mz[n][i] * deardia::ms1BinSize);
				
				if (deardia::ms1MapToPpid[binMz][z].size() > 0)
				{
					candidatePpId.insert(deardia::ms1MapToPpid[binMz][z].begin(),
						deardia::ms1MapToPpid[binMz][z].end());
				}
			}
			if (candidatePpId.size() == 0) continue;

			std::unordered_set<int> candidateMS2;
			for (auto it = candidatePpId.begin(); it != candidatePpId.end(); ++it)
			{
				for (int i = 0; i < deardia::ppidMapToFragment[*it].size(); ++i)
				{
					candidateMS2.insert(std::abs(deardia::ppidMapToFragment[*it][i]));
				}
			}

			std::vector< std::vector< std::vector<float> > > ms2TempSpec(deardia::sliderWidth);

			for (int i = 0; i < winDict[n].size(); ++i)
			{
				int cycleId = winDict[n][i][0];
				int windowId = winDict[n][i][1];

				spectraData.ms2Spectra[cycleId][windowId].DecodeMzXMLData();

				ms2TempSpec[cycleId] = spectraData.ms2Spectra[cycleId][windowId].spectrum;
			}

			deardia::SliderData sliderData;

			std::unordered_set<int> skipMz;

			for (auto mz = candidateMS2.begin(); mz != candidateMS2.end(); ++mz)
			{
				float mzFloat = (float)(*mz) / deardia::ms2BinSize;

				if (mzFloat > deardia::winRange[n][0] &&
					mzFloat < deardia::winRange[n][1])
				{
					continue;
				}

				int nonZeroCount = 0;
				std::vector<float> curve(deardia::sliderWidth, 0.0f);
				std::vector<float> mzTrace(deardia::sliderWidth, 0.0f);

				float leftMz  = mzFloat * (1.0f - deardia::ms2TolPPM);
				float rightMz = mzFloat * (1.0f + deardia::ms2TolPPM);

				for (int t = 0; t < deardia::sliderWidth; ++t)
				{
					if (ms2TempSpec[t].size() == 0) continue;

					int leftId = std::lower_bound(ms2TempSpec[t][0].begin(), 
						ms2TempSpec[t][0].end(), leftMz) - ms2TempSpec[t][0].begin();

					if (leftId == ms2TempSpec[t][0].size()) continue;

					float weightMz = 0.0f;
					for (int i = leftId; i < ms2TempSpec[t][0].size(); ++i)
					{
						if (ms2TempSpec[t][0][i] <= rightMz)
						{
							curve[t] += ms2TempSpec[t][1][i];
							weightMz += ms2TempSpec[t][0][i] * ms2TempSpec[t][1][i];
						}
						else
						{
							break;
						}
					}

					if (curve[t] > 1e-3f)
					{
						++nonZeroCount;
						mzTrace[t] = weightMz / curve[t];
					}
				}

				if (nonZeroCount < 5) continue;

				float sumIntensity = std::accumulate(curve.begin(), curve.end(), 0.0f);
				float detectedMz = std::inner_product(mzTrace.begin(), mzTrace.end(), curve.begin(), 0.0f);
				detectedMz /= sumIntensity;

				/*
				float detectedMz = 0.0f;
				float sumIntensity = 0.0f;
				for (int j = 0; j < deardia::sliderWidth; ++j)
				{
					sumIntensity += curve[j];
					detectedMz += mzTrace[j] * curve[j];
				}
				detectedMz /= sumIntensity;
				*/

				int binMz = (int)(detectedMz * deardia::ms2BinSize);

				if (binMz >= deardia::numFragmentBin) continue;

				if (skipMz.find(binMz) == skipMz.end())
				{
					skipMz.insert(binMz);

					float kurtosis = deardia::Kurtosis(curve);

					std::vector<float> scaleCurve(deardia::sliderWidth);
					auto m = std::minmax_element(curve.begin(), curve.end());
					float scale = *m.second - *m.first;

					for (int j = 0; j < deardia::sliderWidth; ++j)
					{
						scaleCurve[j] = (curve[j] - *m.first) / scale;
					}

					sliderData.ms2Mz.push_back(detectedMz);
					sliderData.kurtosis.push_back(kurtosis);
					sliderData.ms2XICs.push_back(scaleCurve);
					sliderData.ms2Intensity.push_back(*m.second);
				}
			}

			if (sliderData.ms2Mz.size() < deardia::minHitNumber * deardia::kClusters)
			{
				continue;
			}

			sliderData.winId = n;
			sliderData.retentionTime = spectraData.retentionTime;

			sliderData.ms1Mz = ms1Mz[n];
			sliderData.ms1XICs = ms1XICs[n];
			sliderData.ms1Charge = ms1Charge[n];

			++doneCount;

			if (doneCount < deardia::deepLearningQueueSize)
			{
				deardia::deepLearningQueue.enqueue(sliderData);
			}

			while (doneCount >= deardia::deepLearningQueueSize)
			{
				if (deeplearning_done.load(std::memory_order_acquire) > 0)
				{
					deardia::deepLearningQueue.enqueue(sliderData);
					deeplearning_done.fetch_add(-1, std::memory_order_release);
					break;
				}
				else if (deeplearning_done.load(std::memory_order_acquire) == 0)
				{
					continue;
				}
			}
		}

		process_slider_done.fetch_add(1, std::memory_order_release);
	}

	deardia::SliderData sliderData;
	sliderData.stop = true;
	deardia::deepLearningQueue.enqueue(sliderData);
}


void deardia::DeepLearning()
{
	std::cout << "\n Start initilazing deep auto-encoder ......\n";

	mx_uint xDim = (mx_uint)deardia::sliderWidth;
	mx_uint batchSize = (mx_uint)deardia::batchSize;
	mx_uint zDim = (mx_uint)deardia::latentFeaturesDim;

	deardia::MXNetBufferFile jsonData(deardia::mxnetJsonFile);
	deardia::MXNetBufferFile paramData(deardia::mxnetParamFile);

	int deviceId = 0;

	int deviceType;
	if (deardia::useGPU == 1)
	{
		deviceType = 2;
	}
	else
	{
		deviceType = 1;
	}

	mx_uint numInputNodes = 1;

	const char* inputKey[1] = { "data" };
	const char** inputKeys = inputKey;

	const mx_uint inputShapeIndPtr[2] = { 0, 2 };
	const mx_uint inputShapeData[2] = { batchSize, xDim };

	//PredictorHandle predHnd = 0;
	PredictorHandle predHnd = 0;

	if (jsonData.GetLength() == 0 || paramData.GetLength() == 0)
	{
		std::cout << '\n'
			<< "Something was wrong in loadding mxnet json and params files."
			<< "Please check if the files or directory are correct."
			<< '\n';

		exit(0);
	}

	//Create Predictor
	int status = MXPredCreate((const char*)jsonData.GetBuffer(),
		(const char*)paramData.GetBuffer(),
		static_cast<size_t> (paramData.GetLength()),
		deviceType,
		deviceId,
		numInputNodes,
		inputKeys,
		inputShapeIndPtr,
		inputShapeData,
		&predHnd);

	if (status != 0)
	{
		std::cout << "Initialize auto-encoder failed." << '\n'
			<< "Note: in your configure file, the input dimension of deep auto-encoder equals "
			<< deardia::sliderWidth
			<< " and the latent dimension of deep auto-encoder equals "
			<< deardia::latentFeaturesDim
			<< ". Plearse cheak if the mxnet .json file and .params file are corresponded."
			<< std::endl;
	}

	assert(status == 0);
	assert(predHnd);

	std::cout << "Initialize deep variational auto-encoder done.\n" << std::endl;

	init_network_done = true;

	// Apply deep auto-encoder
	
	int stopCount = 0;
	int doneCount = 0;
	
	deardia::SliderData sliderData;
	int inputSize = batchSize * xDim;

	while (1)
	{
		deardia::deepLearningQueue.wait_dequeue(sliderData);

		if (sliderData.stop)
		{
			++stopCount;
			if (stopCount == deardia::analyzeSliderThreads)
			{
				break;
			}
		}
		else
		{
			int nSample = sliderData.ms2XICs.size();
			int batchNumber = (int)(nSample / deardia::batchSize);

			if (nSample % deardia::batchSize != 0) ++batchNumber;

			sliderData.ms2XICFeatures.resize(nSample, std::vector<float>(zDim));

			int n = 0;

			for (int i = 0; i < batchNumber; ++i)
			{
				std::vector<mx_float> inputData(inputSize, 0.0f);

				int m = 0;

				for (int j = i * batchSize; j < (i + 1) * batchSize; ++j)
				{
					if (j >= nSample) break;

					for (int k = 0; k < xDim; ++k)
					{
						inputData[m] = sliderData.ms2XICs[j][k];
						++m;
					}
				}
				
				MXPredSetInput(predHnd, "data", inputData.data(), inputSize);

				// Do predict forward
				MXPredForward(predHnd);

				mx_uint outputIndex = 0;
				mx_uint* shape = 0;
				mx_uint shapeLength;

				//Get output results
				MXPredGetOutputShape(predHnd, outputIndex, &shape, &shapeLength);

				size_t t_size = 1;
				for (mx_uint k = 0; k < shapeLength; ++k)
				{
					t_size *= shape[k];
				}

				std::vector<float> feat(t_size);

				MXPredGetOutput(predHnd, outputIndex, &(feat[0]), t_size);

				for (int j = 0; j < batchSize; ++j)
				{
					if (n >= nSample) break;

					for (int k = 0; k < zDim; ++k)
					{
						sliderData.ms2XICFeatures[n][k] = feat[j * zDim + k];
					}
					++n;
				}
			}

			++doneCount;

			if (doneCount < deardia::clusterQueueSize)
			{
				deardia::clusterQueue.enqueue(sliderData);
			}

			while (doneCount >= deardia::clusterQueueSize)
			{
				if (cluster_done.load(std::memory_order_acquire) > 0)
				{
					deardia::clusterQueue.enqueue(sliderData);
					cluster_done.fetch_add(-1, std::memory_order_release);
					break;
				}
				else if (cluster_done.load(std::memory_order_acquire) == 0)\
				{
					continue;
				}
			}
		}

		deeplearning_done.fetch_add(1, std::memory_order_release);
	}

	//Release Predictor
	MXPredFree(predHnd);

	for (int i = 0; i < deardia::clusterThreads; ++i)
	{
		deardia::SliderData sliderData;

		sliderData.stop = true;
		deardia::clusterQueue.enqueue(sliderData);
	}
}

void deardia::ClusteringAndSearchingDatabase()
{
	double LogFactorial[64] = {		0.000000, 0.000000, 0.693147, 1.791759,
									3.178054, 4.787492, 6.579251, 8.525161,
									10.604603, 12.801827, 15.104413, 17.502308,
									19.987214, 22.552164, 25.191221, 27.899271,
									30.671860, 33.505073, 36.395445, 39.339884,
									42.335616, 45.380139, 48.471181, 51.606676,
									54.784729, 58.003605, 61.261702, 64.557539,
									67.889743, 71.257039, 74.658236, 78.092224,
									81.557959, 85.054467, 88.580828, 92.136176,
									95.719695, 99.330612, 102.968199,106.631760,
									110.320640, 114.034212, 117.771881, 121.533082,
									125.317271, 129.123934, 132.952575, 136.802723,
									140.673924, 144.565744, 148.477767, 152.409593,
									156.360836, 160.331128, 164.320112, 168.327445,
									172.352797, 176.395848, 180.456291, 184.533829,
									188.628173, 192.739047, 196.866182, 201.009316
								}; // be used in hyper-score

	deardia::SliderData sliderData;

	while (1)
	{
		deardia::clusterQueue.wait_dequeue(sliderData);

		if (sliderData.stop) break;

		cluster_done.fetch_add(1, std::memory_order_release);

		int winId = sliderData.winId;
		int nSamples = sliderData.ms2Mz.size();

		std::unordered_set<int> ms1HitPepId;
		std::unordered_map<int, int> ppidMapToMS1;

		for (int i = 0; i < sliderData.ms1Mz.size(); ++i)
		{
			int binMz = (int)(sliderData.ms1Mz[i] * deardia::ms1BinSize);
			int charge = sliderData.ms1Charge[i];

			if (deardia::ms1MapToPpid[binMz][charge].size() > 0)
			{
				for (int j = 0; j < deardia::ms1MapToPpid[binMz][charge].size(); ++j)
				{
					int ppid = deardia::ms1MapToPpid[binMz][charge][j];

					ms1HitPepId.insert(ppid);
					ppidMapToMS1[ppid] = i;
				}
			}
		}

		if (ms1HitPepId.size() == 0) continue;

		std::vector< std::vector<int> > ms2MapToPepId(deardia::numFragmentBin);

		for (auto it = ms1HitPepId.begin(); it != ms1HitPepId.end(); ++it)
		{
			for (int i = 0; i < deardia::ppidMapToFragment[*it].size(); ++i)
			{
				int binMz = deardia::ppidMapToFragment[*it][i];

				if (binMz < 0)
				{
					ms2MapToPepId[-binMz].push_back(-(*it));
				}
				else
				{
					ms2MapToPepId[binMz].push_back(*it);
				}
			}
		}

		for (int i = 0; i < nSamples; ++i)
		{
			float corr = deardia::AutoCorrelation(sliderData.ms2XICs[i]);
			sliderData.ms2XICFeatures[i].push_back(corr);
		}

		std::vector<int> labels = KMEANS::Kmeans(sliderData.ms2XICFeatures, deardia::kClusters);
		int maxLabel = *std::max_element(labels.begin(), labels.end());

		for (int n = 0; n < maxLabel; ++n)
		{
			std::vector<int> clusterIdx;
			for (int i = 0; i < nSamples; ++i)
			{
				if (labels[i] == n) clusterIdx.push_back(i);
			}
			if (clusterIdx.size() < deardia::minHitNumber) continue;

			int clusterSize = clusterIdx.size();

			float maxPearson = 0.0f;
			int maxRow = 0;
			int maxCol = 0;

			for (int i = 0; i < clusterSize - 1; ++i)
			{
				for (int j = i + 1; j < clusterSize; ++j)
				{
					float p = deardia::PearsonCorrcoef(sliderData.ms2XICs[clusterIdx[i]],
														sliderData.ms2XICs[clusterIdx[j]]);

					if (p > maxPearson &&
						std::abs(sliderData.ms2Mz[clusterIdx[i]] - 
								 sliderData.ms2Mz[clusterIdx[j]]) > 100.0f)
					{
						maxPearson = p;
						maxRow = i;
						maxCol = j;
					}
				}
			}

			if (maxPearson < 0.6f) continue;

			std::vector<int> clusterSet;
			clusterSet.push_back(clusterIdx[maxRow]);
			clusterSet.push_back(clusterIdx[maxCol]);

			for (int i = 0; i < clusterSize; ++i)
			{
				if (i == maxRow || i == maxCol) continue;

				float p1 = deardia::PearsonCorrcoef(sliderData.ms2XICs[clusterIdx[i]], sliderData.ms2XICs[clusterIdx[maxRow]]);
				float p2 = deardia::PearsonCorrcoef(sliderData.ms2XICs[clusterIdx[i]], sliderData.ms2XICs[clusterIdx[maxCol]]);

				if (p1 > 0.4f || p2 > 0.4f)
				{
					clusterSet.push_back(clusterIdx[i]);
				}
			}
			if (clusterSet.size() < deardia::minHitNumber) continue;

			if (clusterSet.size() > 250)
			{
				std::vector< std::pair<int, float> > sortKurtosis;
				for (auto it = clusterSet.begin(); it != clusterSet.end(); ++it)
				{
					sortKurtosis.push_back(std::make_pair(*it, sliderData.kurtosis[*it]));
				}

				std::sort(sortKurtosis.begin(), sortKurtosis.end(),
					[](const std::pair<int, float>& x, const std::pair<int, float>& y)
				{
					return x.second > y.second;
				});

				clusterSet.clear();
				for (int i = 0; i < 250; ++i)
				{
					clusterSet.push_back(sortKurtosis[i].first);
				}
			}
			std::sort(clusterSet.begin(), clusterSet.end());
			
			std::unordered_map<int, int> ms2HitFragment;
			std::unordered_map<int, std::vector<int> > clusterHitPepId;

			for (int i = 0; i < clusterSet.size(); ++i)
			{
				int idx = clusterSet[i];

				int mz = (int)(sliderData.ms2Mz[idx] * deardia::ms2BinSize);
				
				if (ms2MapToPepId[mz].size() == 0) continue;

				ms2HitFragment[mz] = idx;

				for (int j = 0; j < ms2MapToPepId[mz].size(); ++j)
				{
					int ppid = std::abs(ms2MapToPepId[mz][j]);

					if (clusterHitPepId.find(ppid) == clusterHitPepId.end())
					{
						std::vector<int> hitInfo(4, 0);

						if (ms2MapToPepId[mz][j] < 0)
						{
							hitInfo[0] = 1;
							hitInfo[1] = sliderData.ms2Intensity[idx];
						}
						else
						{
							hitInfo[2] = 1;
							hitInfo[3] = sliderData.ms2Intensity[idx];
						}
						clusterHitPepId[ppid] = hitInfo;
					}
					else
					{
						if (ms2MapToPepId[mz][j] < 0)
						{
							++clusterHitPepId[ppid][0];
							clusterHitPepId[ppid][1] += sliderData.ms2Intensity[idx];
						}
						else
						{
							++clusterHitPepId[ppid][2];
							clusterHitPepId[ppid][3] += sliderData.ms2Intensity[idx];
						}
					}
				}
				
			}
			if (clusterHitPepId.size() == 0) continue;

			std::vector< std::pair<int, double> > ppidScore;
			for (auto it = clusterHitPepId.begin(); it != clusterHitPepId.end(); ++it)
			{
				if ((it->second[0] + it->second[2]) < deardia::minHitNumber) continue;
				if (it->second[0] == 0 || it->second[2] == 0) continue;

				int Nb = it->second[0];
				int Ny = it->second[2];

				double sumIntensity = it->second[1] + it->second[3];
				double hyperScore = LogFactorial[Nb] + LogFactorial[Ny] + log(sumIntensity);

				ppidScore.push_back(std::make_pair(it->first, hyperScore));
			}
			if (ppidScore.size() == 0) continue;

			if (ppidScore.size() > 1)
			{
				std::sort(ppidScore.begin(), ppidScore.end(),
					[](const std::pair<int, double>& x, const std::pair<int, double>& y)
				{
					return x.second > y.second;
				});
			}

			for (int i = 0; i < ppidScore.size(); ++i)
			{
				if (i > 100) break;

				int ppid = ppidScore[i].first;
				int ms1Idx = ppidMapToMS1[ppid];

				bool accept = false;
				for (int j = 0; j < deardia::ppidMapToFragment[ppid].size(); ++j)
				{
					int fragMz = deardia::ppidMapToFragment[ppid][j];

					if (ms2HitFragment.find(fragMz) != ms2HitFragment.end())
					{
						float pMS1MS2 = deardia::PearsonCorrcoef(sliderData.ms1XICs[ms1Idx], sliderData.ms2XICs[ms2HitFragment[fragMz]]);
						if (pMS1MS2 > 0.4f)
						{
							accept = true;
							break;
						}
					}

				}

				if (accept)
				{
					deardia::MgfData mgfData;

					mgfData.ms1Mz = sliderData.ms1Mz[ms1Idx];
					mgfData.retentionTime = sliderData.retentionTime;
					mgfData.ms1Charge = sliderData.ms1Charge[ms1Idx] + deardia::ms1MinCharge;

					for (int j = 0; j < clusterSet.size(); ++j)
					{
						std::vector<float> temp = { sliderData.ms2Mz[clusterSet[j]], 
													sliderData.ms2Intensity[clusterSet[j]] };

						mgfData.ms2Spectrum.push_back(temp);
					}

					std::sort(mgfData.ms2Spectrum.begin(), mgfData.ms2Spectrum.end(),
						[](std::vector<float>& x, std::vector<float>& y)
						{
							return x[0] < y[0];
						});

					deardia::mgfQueue.enqueue(mgfData);
					break;
				}
			}
		}
	}

	deardia::MgfData mgfData;
	mgfData.stop = true;

	deardia::mgfQueue.enqueue(mgfData);
}

void deardia::WriteMGF()
{
	std::ofstream MgfFile(deardia::mgfFile.c_str());

	if (!MgfFile)
	{
		std::cout << "Cannot open mgf file: " << deardia::mgfFile
			<< "\nPlease cheak your file or directory.\n";
		exit(1);
	}

	int stopCount = 0;
	int spectrumCount = 0;

	deardia::MgfData mgfData;

	while (1)
	{
		deardia::mgfQueue.wait_dequeue(mgfData);

		if (mgfData.stop)
		{
			++stopCount;
			if (stopCount == deardia::clusterThreads) break;
		}
		else
		{
			++spectrumCount;

			MgfFile << "BEGIN IONS\n"
					<< "PEPMASS=" << mgfData.ms1Mz << "\n"
					<< "RTINSECONDS=" << mgfData.retentionTime << "\n"
					<< "CHARGE=" << mgfData.ms1Charge << "\n";

			for (int i = 0; i < mgfData.ms2Spectrum.size(); ++i)
			{
				char buff[32];
				sprintf(buff, "%.2f %.1f\n", mgfData.ms2Spectrum[i][0], mgfData.ms2Spectrum[i][1]);
				MgfFile << buff;
			}

			MgfFile << "END IONS\n" << std::endl;
		}
	}
	MgfFile.close();

	std::cout << "pesudo spectrum number: " << spectrumCount << std::endl;
}

#endif
