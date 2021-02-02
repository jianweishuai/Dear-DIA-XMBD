#ifndef FINDPEAKSX_H
#define FINDPEAKSX_H

#include <math.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

namespace FindPeaksSX
{
	int sign(float x);
	std::vector<float> derive(std::vector<float>& x);
	int val2ind(std::vector<float>& x, float val);
	std::vector<float> SA(std::vector<float>& y, int smoothwidth, int ends = 0);
	void fastsmooth(std::vector<float>& y, std::vector<float>& smooth_y, int smoothwidth, int smoothtype, int ends);
	std::vector<float> findpeaks(std::vector<float>& x, std::vector<float>& y);
}

int FindPeaksSX::sign(float x)
{
	if (x < 0)
	{
		return -1;
	}
	else if (x == 0)
	{
		return 0;
	}
	else if (x > 0)
	{
		return 1;
	}
}

std::vector<float> FindPeaksSX::derive(std::vector<float> &x)
{
	unsigned int n = x.size() - 1;
	std::vector<float> d(x.size(), 0.0f);
	
	d[0] = x[1] - x[0];
	d[n] = x[n] - x[n-1];
	
	for (unsigned int i = 1; i < n-1; ++i)
	{
		d[i] = (x[i+1] - x[i-1]) / 2;
	}
	
	return d;
}

int FindPeaksSX::val2ind(std::vector<float> &x, float val)
{
	std::vector<float> dif(x.size(), 0.0f);
	for (int i = 0; i < dif.size(); ++i)
	{
		dif[i] = std::abs(x[i] - val);
	}
	
	auto min_dif = std::min_element(dif.begin(), dif.end());
	
	int index =  std::distance(dif.begin(), min_dif);
	
	return index;
}

std::vector<float> FindPeaksSX::SA(std::vector<float> &y, int smoothwidth, int ends)
{
	float sum_points = std::accumulate(y.begin(), y.begin()+smoothwidth, 0.0f);
	
	int L = y.size();
	int w = int(round(smoothwidth));
	std::vector<float> s(y.size(), 0.0f);
	int half_w = int(round(smoothwidth / 2));
	
	for (int k = 0; k < L-w; ++k)
	{
		s[k+half_w-1] = sum_points;
		sum_points = sum_points - y[k];
		sum_points = sum_points + y[k+w];
	}
	
	s[L-w-1+half_w] = std::accumulate(y.begin()+L-w, y.begin()+L, 0.0f);
	
	std::vector<float> smooth_y(L, 0.0f);
	
	for (int i =0 ; i < s.size(); ++i)
	{
		smooth_y[i] = s[i] / (float)w;
	}
	
	return smooth_y;
}

void FindPeaksSX::fastsmooth(std::vector<float> &y, std::vector<float> &smooth_y, int smoothwidth, int smoothtype, int ends)
{
	if (smoothtype == 1)
	{
		smooth_y = FindPeaksSX::SA(y, smoothwidth, ends);
	}
	else if (smoothtype == 2)
	{
		std::vector<float> t1 = FindPeaksSX::SA(y, smoothwidth, ends);
		smooth_y = FindPeaksSX::SA(t1, smoothwidth, ends);
	}
	else if (smoothtype == 3)
	{
		std::vector<float> t1 = FindPeaksSX::SA(y, smoothwidth, ends);
		std::vector<float> t2 = FindPeaksSX::SA(t1, smoothwidth, ends);
		smooth_y = FindPeaksSX::SA(t2, smoothwidth, ends);
	}
	else if (smoothtype == 4)
	{
		std::vector<float> t1 = FindPeaksSX::SA(y, smoothwidth, ends);
		std::vector<float> t2 = FindPeaksSX::SA(t1, smoothwidth, ends);
		std::vector<float> t3 = FindPeaksSX::SA(t2, smoothwidth, ends);
		smooth_y = FindPeaksSX::SA(t3, smoothwidth, ends);
	}
}

std::vector<float> FindPeaksSX::findpeaks(std::vector<float> &x, std::vector<float> &y)
{
	std::vector<float> peakX;
	
	float WidthPoints = 3.0f;
	float SlopThreashold = 0.7 / (WidthPoints * WidthPoints);
	
	float AmpThreashold = 0.0f;
	
	int smoothwidth = 3;
	int peakgroup = 3;
	int smoothtype = 3;
	
	std::vector<float> d;
	
	if (smoothwidth > 1)
	{
		std::vector<float> derive_y = FindPeaksSX::derive(y);
		FindPeaksSX::fastsmooth(derive_y, d, smoothwidth, smoothtype, 0);
	}
	else
	{
		d = FindPeaksSX::derive(y);
	}

	int n = int(round(peakgroup / 2)) + 1;
	
	int vectorlen = y.size();
	
	int start = 2 * int(round(smoothwidth / 2)) - 2;
	int end = y.size() - smoothwidth;
	
	for (int i = start; i < end; ++i)
	{
		if (FindPeaksSX::sign(d[i]) > FindPeaksSX::sign(d[i+1]))
		{
			if (d[i]-d[i+1] > SlopThreashold)
			{
				if (y[i] > AmpThreashold)
				{
					std::vector<float> xx(peakgroup, 0.0f);
					std::vector<float> yy(peakgroup, 0.0f);
					
					for (int j = 0; j < peakgroup; ++j)
					{
						int groupindex = i + j - n + 1;
						if (groupindex < 0)
						{
							groupindex = 0;
						}
						if (groupindex >= vectorlen)
						{
							groupindex = vectorlen-1;
						}
						xx[j] = x[groupindex];
						yy[j] = y[groupindex];
					}
					
					float peakY = *std::max_element(yy.begin(), yy.end());

					int peak_index = FindPeaksSX::val2ind(yy, peakY);
					peakX.push_back(xx[peak_index]);
				}
			}
		}
	}
	
	return peakX;
}

#endif