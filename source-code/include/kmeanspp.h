#ifndef KMEANSPP_H
#define KMEANSPP_H

#include <vector>
#include <stdio.h>
#include <stdlib.h>

#define MAX_ITER 200
#define FLT_MAX 3.402823466e+38F 

#define N_FEATURES 17

namespace KMEANS
{
	typedef struct
	{
		int group = 0;
		float data[N_FEATURES];
	} Point;

	float randf(float x);

	float dist2(Point& x, Point& y);
	int nearest(Point& pt, std::vector<Point>& centers, int n_centers, float* d2);
	void kpp(std::vector<Point>& pts, int n_samples, std::vector<Point>& centers, int n_clusters);
	std::vector<int> Kmeans(std::vector< std::vector<float> >& features, int n_clusters);
}

float KMEANS::randf(float x)
{
	return x * rand() / (RAND_MAX - 1.0f);
}

float KMEANS::dist2(Point &x, Point &y)
{
	float sum = 0.0f;
	
	for (int i = 0; i < N_FEATURES; ++i)
	{
		sum += (x.data[i] - y.data[i]) * (x.data[i] - y.data[i]);
	}
	
	return sum;
}

int KMEANS::nearest(Point &pt, std::vector<Point> &centers, int n_centers, float *d2)
{
	int min_idx = pt.group;
	float min_dist = FLT_MAX;
	
	float dist = 0.0f;
	
	for (int i = 0; i < n_centers; ++i)
	{
		if (min_dist > (dist = KMEANS::dist2(pt, centers[i])))
		{
			min_dist = dist;
			min_idx = i;
		}
	}
	
	if (d2) *d2 = min_dist;
	
	return min_idx;
}

void KMEANS::kpp(std::vector<Point> &pts, int n_samples, std::vector<Point> &centers, int n_clusters)
{
	centers[0] = pts[ rand() % n_samples ];
	
	float sum;
	float *d = (float *)malloc(n_samples * sizeof(float));
	
	for (int i = 1; i < n_clusters; ++i)
	{
		sum = 0.0f;
		
		for (int j = 0; j < n_samples; ++j)
		{
			KMEANS::nearest(pts[j], centers, i, d+j);
			sum += d[j];
		}
		sum = KMEANS::randf(sum);
		
		for (int j = 0; j < n_samples; ++j)
		{
			if ((sum -= d[j]) > 0) continue;
			
			centers[i] = pts[j];
			
			break;
		}
		
	}
	
	for (int j = 0; j < n_samples; ++j)
	{
		pts[j].group = KMEANS::nearest(pts[j], centers, n_clusters, 0);
	}
	
	free(d);
}

std::vector<int> KMEANS::Kmeans(std::vector< std::vector<float> > &features, int n_clusters)
{
	int n_samples = features.size();
	
	std::vector<Point> pts(n_samples);
	
	for (int i = 0; i < n_samples; ++i)
	{
		for (int j = 0; j < N_FEATURES; ++j)
		{
			pts[i].data[j] = features[i][j];
		}
	}
	
	std::vector<Point> cluster_centers(n_clusters);
	
	KMEANS::kpp(pts, n_samples, cluster_centers, n_clusters);
	
	int changed;
	int stop_msg = n_samples >> 10; //0.1% of n_samples
	
	int iter = 0;
	
	do
	{
		++iter;
		
		for (int i = 0; i < n_clusters; ++i)
		{
			for (int j = 0; j < N_FEATURES; ++j)
			{
				cluster_centers[i].data[j] = 0.0f;
			}
			cluster_centers[i].group = 0;
		}
		
		for (int i = 0; i < n_samples; ++i)
		{
			++cluster_centers[pts[i].group].group;
			
			for (int j = 0; j < N_FEATURES; ++j)
			{
				cluster_centers[pts[i].group].data[j] += pts[i].data[j];
			}
		}
		
		for (int i = 0; i < n_clusters; ++i)
		{
			if (cluster_centers[i].group == 0) continue;
			
			for (int j = 0; j < N_FEATURES; ++j)
			{
				cluster_centers[i].data[j] /= cluster_centers[i].group;
			}
		}
		
		changed = 0;
		
		for (int i = 0; i < n_samples; ++i)
		{
			int min_idx = KMEANS::nearest(pts[i], cluster_centers, n_clusters, 0);
			if (min_idx != pts[i].group)
			{
				++changed;
				pts[i].group = min_idx;
			}
		}
		
	} while (changed > stop_msg && iter <= MAX_ITER);
	
	std::vector<int> labels(n_samples);
	
	for (int i = 0; i < n_samples; ++i)
	{
		labels[i] = pts[i].group;
	}
	
	return labels;
}

#endif