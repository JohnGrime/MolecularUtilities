/*
	Author: John Grime, The University of Chicago.
*/

#if !defined( UTIL_STATS_DEFINED )
#define UTIL_STATS_DEFINED

namespace Util
{
namespace Stats
{

//
// Running statistics, based on algorithms of B. P. Welford
// (via Knuth, "The Art of Computer Programming").
//
// This algorithm provides not only the capability to determine
// variance (and hence, stdev and stderr) from a running input
// with very little storage, it is also robust to catastrophic
// cancellation.
//
struct Stats
{
	size_t N; // number of samples so far
	double S; // N*sigma^2
	double min, mean, max;
	
	Stats()
	{
		Clear();
	}
	void Clear()
	{
		N = 0;
		S = min = mean = max = 0.0;
	}
	void AddSample( double x )
	{
		N++;
		if( N == 1 )
		{
			min = mean = max = x;
			S = 0.0;
			return;
		}
		
		//
		// Update values (new marked with prime):
		//   mean' = mean + (x-mean)/N
		//   S' = S + (x-mean)*(x-mean')
		//

		double delta = (x-mean);
		mean += delta/N;
		S += delta * (x-mean); // <- note: uses updated value of "mean" as wella s old value via "delta".

		if( x < min ) min = x;
		if( x > max ) max = x;
	}
	double Sum() const
	{
		return mean*N;
	}
	double Variance() const
	{
		// SAMPLE variance
		return (N>1) ? (S/(N-1)) : (0.0);
	}
	double StdDev() const
	{
		// SAMPLE standard deviation
		return sqrt( Variance() );
	}
	double StdErr() const
	{
		// Estimated standard error of the sample mean:
		// SE = stdev / sqrt(N) : stdev = sample standard deviation
		// SE = sqrt(variance) / sqrt(N) : as stdev = sqrt(variance)
		// SE = sqrt( variance / N ) : as sqrt() distributive
		return (N>1) ? ( sqrt(Variance()/N) ) : (0.0);
	}
	
	// Combine separate sets of sample stats
	Stats& operator += ( const Stats &rhs )
	{
		// Temporary values, in case &rhs == this
		size_t new_N;
		double new_sum, new_mean, new_S;
		
		// Ignore if no data present in rhs
		if( rhs.N < 1 ) return *this;
				
		new_N = N + rhs.N;
		new_sum = (mean*N) + (rhs.mean*rhs.N);
		new_mean = new_sum / new_N; // safe: rhs.N >= 1, so new_N >= 1.
		
		// This is basically the "parallel algorithm" version of the "online" algorithm
		// for calculating variance in one pass when the sample is partitioned into multiple
		// sets. This is attributed to Chan et al, "Updating Formulae and a Pairwise Algorithm for
		// Computing Sample Variances.", Technical Report STAN-CS-79-773, Stanford CS, (1979).
		//
		// http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Parallel_algorithm
		//
		double delta = (mean - rhs.mean);
		new_S = (S + rhs.S) + (delta*delta)*(N*rhs.N)/new_N;
		
		// If no samples in this Stats structure, or if rhs min/max should replace current:
		if( N < 1 || rhs.min < min ) min = rhs.min;
		if( N < 1 || rhs.max > max ) max = rhs.max;
		
		// Only change N now, as we needed it in the min/max update above.
		N = new_N;
		S = new_S;
		mean = new_mean;
		
		return *this;
	}
};

//
// "Sparse" statistics object for building histograms etc with specified bin width "delta" which acts to
// discretise the sampling domain: this is useful where we don't know the upper and/or lower bounds of the
// sampling domain in advance. This approach also avoids wasted memory where we require a potentially large
// sampling domain but comparatively narrow bins, only some of which may actually be populated.
//
// There is also an optional per-bin coordinate stats object; this allows us to save x columns which better
// reflect the location of the points added to the bin.
//
struct MapStats
{
	std::map<int,Stats> values, coords;
	double delta, min;

	MapStats()                     { Reset( 1.0, 0.0 ); }
	MapStats( double d )           { Reset( d,   0.0 ); }
	MapStats( double d, double m ) { Reset( d,   m   ); }

	void Reset( double delta_, double min_ = 0.0 )
	{
		min = min_;
		Clear( delta_ );
	}
	void Clear( double delta_ = 0.0 )
	{
		if( delta_ != 0.0 )  delta = delta_;
		values.clear();
		coords.clear();
	}

	double CoordToBin( double coord ) const { return (int)floor( (coord-min)/delta ); }
	double BinToCoord( int bin_id )   const { return min+(0.5+bin_id)*delta; }

	//
	// Same params, so need different method names.
	//
	void AddSample( double coord, double val )
	{
		int bin = CoordToBin( coord );
		coords[bin].AddSample( coord );
		values[bin].AddSample( val );
	}

	int Save( FILE* f, bool print_header = true )
	{
		if( f == nullptr ) return -1;
		if( print_header )
		{
			fprintf( f, "# x, y = coordinate and mean value at coordinate\n" );
			fprintf( f, "#\n" );
			fprintf( f, "# N(y)      = number of samples\n" );
			fprintf( f, "# StdDev(y) = standard deviation: sqrt(variance)\n" );
			fprintf( f, "# StdErr(y) = standard error of the estimated mean: sqrt(variance)/N = StdDev(y)/N\n" );
			fprintf( f, "#\n" );
			fprintf( f, "# Both StdDev(y) and StdErr(y) use the SAMPLE variance\n" );
			fprintf( f, "#\n" );
			fprintf( f, "# %12.12s  %12.12s  %12.12s  %12.12s  %12.12s\n", "x", "y", "N(y)", "StdDev(y)", "StdErr(y)" );
			fprintf( f, "#\n" );
		}
		for( const auto& it : values )
		{
			int bin_no = it.first;
			const Stats& v = it.second;
			const Stats& c = coords[bin_no];
			fprintf( f, "  %e  %e  %12d  %e  %e\n", c.mean, v.mean, (int)v.N, v.StdDev(), v.StdErr() );
		}
		return 1;
	}
};

//
// Distribution class: designed to build data in a number of passes. For example, if we were interested in an angle distribution
// measured from simulation:
//
// 1. Load trajectory frame.
// 3. Walk angles in trajectory frame, calling AddSample( theta ) for each angle.
// 4. Call Accumulate()
// 5. Repeat from step 1.
//
// This class basically "adds" each per-pass distribution as a separate sample for the total distribution. The stddev, stderr
// etc are thus calculated for each "bin" using the final per-pass value in the bins when Accumulate() called.
//
struct Distribution
{
	MapStats pass, total;

	Distribution( double delta_ = 1.0 ) { Clear( delta_ );   }

	void Clear( double delta_ = 0.0 )
	{
		if( delta_ == 0.0 ) delta_ = pass.delta;
		pass.Clear( delta_ );
		total.Clear( delta_ );
	}

	double CoordToBin( double coord ) const { return pass.CoordToBin(coord);  }
	double BinToCoord( int bin_id )   const { return pass.BinToCoord(bin_id); }

	void AddSample( double coord ) { pass.AddSample( coord, 1.0 ); }

	// Optional prefactor for weighting accumulation!
	void Accumulate( double prefactor = 1.0 )
	{
		// Accumulate total values (i.e. N*mean) from "pass" into "total".
		for( const auto& it : pass.values )
		{
			const auto& key = it.first;
			const auto& s = it.second;
			if( s.N < 1 ) continue; // ignore unsampled data, rather than accumulating a sample of 0 into "total"!
			total.values[key].AddSample( prefactor*(s.mean*s.N) );
		}
		pass.Clear();
	}

	int Save( FILE* f, bool print_header = true ) const
	{
		if( f == nullptr ) return -1;
		if( print_header )
		{
			fprintf( f, "# x, y = coordinate and mean value at coordinate\n" );
			fprintf( f, "#\n" );
			fprintf( f, "# N(y)      = number of samples\n" );
			fprintf( f, "# StdDev(y) = standard deviation: sqrt(variance)\n" );
			fprintf( f, "# StdErr(y) = standard error of the estimated mean: sqrt(variance)/N = StdDev(y)/N\n" );
			fprintf( f, "# Sum(y) = sum over all y (should be a constant value)\n" );
			fprintf( f, "#\n" );
			fprintf( f, "# Both StdDev(y) and StdErr(y) use the SAMPLE variance\n" );
			fprintf( f, "# P(x) can be plotted as y/Sum(y)\n" );
			fprintf( f, "#\n" );
			fprintf( f, "# %12.12s  %12.12s  %12.12s  %12.12s  %12.12s  %12.12s\n", "x", "y", "N(y)", "StdDev(y)", "StdErr(y)", "Sum(y)" );
			fprintf( f, "#\n" );
		}

		// Accumulate normalizing value, Sum(y)
		double sum_y = 0.0;
		for( const auto& it : total.values ) { sum_y += it.second.mean; }

		// Write data columns. Note Sum(y) as the final column, so y/Sum(x) = P(x): probability of x.
		for( const auto& it : total.values )
		{
			int bin_no = it.first;
			const Stats& s = it.second;
			double coord = begin(total.values)->first + (0.5+bin_no)*total.delta;
			fprintf( f, "  %e  %e  %12d  %e  %e  %e\n", coord, s.mean, (int)s.N, s.StdDev(), s.StdErr(), sum_y );
		}
		return 1;
	}
};

}
}

#endif
