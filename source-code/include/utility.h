#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <vector>
#include <fstream>
#include <numeric>
#include <algorithm>

#include "namespace.h"

std::vector<std::string> deardia::GetInputFileNameAndSuffix(std::string inputFileName)
{
	std::string name = "";
	std::string rev_name = "";

	std::string suffix = "";
	std::string rev_suffix = "";

	bool flag = false;
	for (int i = inputFileName.length() - 1; i >= 0; --i)
	{
		if (inputFileName[i] == '/' ||
			inputFileName[i] == '\\')
		{
			break;
		}
		if (inputFileName[i] == '.')
		{
			flag = true;
			continue;
		}
		if (!flag)
		{
			rev_suffix += inputFileName[i];
		}

		if (flag)
		{
			rev_name += inputFileName[i];
		}
	}
	for (int i = rev_suffix.length() - 1; i >= 0; --i)
	{
		suffix += rev_suffix[i];
	}
	for (int i = rev_name.length() - 1; i >= 0; --i)
	{
		name += rev_name[i];
	}

	std::vector<std::string> nameSuffix(2);
	nameSuffix[0] = name;
	nameSuffix[1] = suffix;

	return nameSuffix;
}

void deardia::LoadWindowFile()
{
	std::ifstream WindowFile(deardia::windowFile.c_str());

	if (!WindowFile)
	{
		printf("Cannot open ms1 windows file. Please cheak %s\n", deardia::windowFile.c_str());
		exit(0);
	}

	std::string line;

	int p, pos1, pos2, pos3;

	while (std::getline(WindowFile, line))
	{
		if (line.empty()) continue;

		pos1 = line.find(' ');
		pos2 = line.find('\t');
		pos3 = line.find(',');

		if (pos1 == line.npos &&
			pos2 == line.npos &&
			pos3 == line.npos)
		{
			printf("Something was wrong in loadding window file, please cheak your file format.\n \
				deardia only support the following format:\n \
				' ' or '\t' or ','\n");
			exit(0);
		}

		if (pos1 != line.npos)
		{
			p = pos1;
		}
		else if (pos2 != line.npos)
		{
			p = pos2;
		}
		else if (pos3 != line.npos)
		{
			p = pos3;
		}

		std::vector<float> range(2, 0.0f);
		range[0] = atof(line.substr(0, p).c_str());
		range[1] = atof(line.substr(p + 1).c_str());

		deardia::winRange.push_back(range);
	}
	WindowFile.close();

	if (deardia::winRange.size() == 0)
	{
		printf("Something was wrong in loadding window file, please cheak your file format.\n \
				deardia only support the following format:\n \
				' ' or '\t' or ','\n");
		exit(0);
	}

	printf("\nLoad window file done. This window file includes %d windows.\n", (int)deardia::winRange.size());
}

std::string deardia::FindKeyValue(const std::string& in_str, const std::string& key)
{
	int count = 0;
	int pos;
	int flag = 0;
	std::string value = "";

	pos = in_str.find(key);
	if (pos != in_str.npos)
	{
		for (int i = pos; i < in_str.length(); ++i)
		{
			if (flag == 0 && in_str[i] == '\"')
			{
				flag = 1;
			}
			else if (flag == 1 && in_str[i] != '\"')
			{
				value.push_back(in_str[i]);
			}
			else if (flag == 1 && in_str[i] == '\"')
			{
				return value; //value
			}
		}
	}
	else
	{
		return "";
	}
}

float deardia::PearsonCorrcoef(std::vector<float>& X, std::vector<float>& Y)
{
	float n = (float)X.size();

	float xy_sum = std::inner_product(X.begin(), X.end(), Y.begin(), 0.0f);
	float x2_sum = std::inner_product(X.begin(), X.end(), X.begin(), 0.0f);
	float y2_sum = std::inner_product(Y.begin(), Y.end(), Y.begin(), 0.0f);
	float x_sum = std::accumulate(X.begin(), X.end(), 0.0f);
	float y_sum = std::accumulate(Y.begin(), Y.end(), 0.0f);

	float varX_varY = sqrtf((x2_sum - x_sum * x_sum / n) * (y2_sum - y_sum * y_sum / n));

	return (xy_sum - x_sum * y_sum / n) / varX_varY;

}

float deardia::Skewness(std::vector<float>& x)
{
	float n = (float)x.size();

	float x_mean = std::accumulate(x.begin(), x.end(), 0.0f) / n;

	float sum1 = 0.0f;
	float sum2 = 0.0f;
	for (int i = 0; i < x.size(); ++i)
	{
		sum1 += powf(x[i] - x_mean, 3.0f);
		sum2 += powf(x[i] - x_mean, 2.0f);
	}
	sum1 /= n;
	sum2 /= n;

	float skewness = sqrtf(n * (n - 1.0f)) / (n - 2.0f) * sum1 / powf(sum2, 1.5f);

	return skewness;
}

float deardia::Kurtosis(std::vector<float>& x)
{
	float n = (float)x.size();

	float x_mean = std::accumulate(x.begin(), x.end(), 0.0f) / n;

	float sum1 = 0.0f;
	float sum2 = 0.0f;
	for (int i = 0; i < x.size(); ++i)
	{
		sum1 += powf(x[i] - x_mean, 4.0f);
		sum2 += powf(x[i] - x_mean, 2.0f);
	}

	float kurtosis = n * (n + 1.0f) * (n - 1.0f) / ((n - 2.0f) * (n - 3.0f)) *
		sum1 / powf(sum2, 2.0f) - 3 * powf(n - 1.0f, 2.0f) / ((n - 2.0f) * (n - 3.0f));

	return kurtosis;
}

float deardia::KL_divergence(std::vector<float>& x, std::vector<float>& y)
{
	float value = 0.0f;
	for (int i = 0; i < x.size(); ++i)
	{
		value += x[i] * logf(x[i] / y[i]);
	}

	return value;
}

float deardia::AutoCorrelation(std::vector<float>& x)
{
	float mean = std::accumulate(x.begin(), x.end(), 0.0f);
	float sigma = 0.0f;
	for (int i = 0; i < deardia::sliderWidth; ++i)
	{
		sigma += powf(x[i] - mean, 2.0f);
	}

	sigma /= (float)(deardia::sliderWidth - 1.0f);

	float corr = 0.0f;
	for (int i = 0; i < deardia::sliderWidth - 3; ++i)
	{
		corr += (x[i] - mean) * (x[i + 3] - mean);
	}
	corr /= (float)(deardia::sliderWidth - 3.0f) * sigma;

	return corr;
}

#endif