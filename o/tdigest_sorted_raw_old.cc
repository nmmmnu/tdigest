#include <cstdint>
#include <cassert>
#include <cstdio>

#include <algorithm>
#include <iterator>
#include <limits>



namespace t_digest_impl_{
	#if 0
	template<typename T>
	void insertSorted2(T *arr, size_t const n, T &&last){
	    auto it = std::lower_bound(arr, arr + n, last);
	    std::move_backward(it, arr + n, arr + n + 1);
	    *it = std::move(last);
	}
	#endif

	template<typename IT>
	void insertIntoSortedRange(IT first, IT last, typename std::iterator_traits<IT>::value_type &&item){
	    auto it = std::lower_bound(first, last, item);
	    std::move_backward(it, last, std::next(last));

	    *it = std::move(item);
	}
}



struct CentroidData{
	struct Centroid{
		double   mean	= 0;
		uint64_t weight	= 0;

		constexpr Centroid() = default;
		constexpr Centroid(double m, uint64_t w) : mean(m), weight(w){}

		void clear(){
			mean   = 0;
			weight = 0;
		}

		operator bool() const{
			return weight;
		}

		double getWeightedMean() const{
			return mean * static_cast<double>(weight);
		}

		void print() const{
			printf("> Addr %p | mean: %10.4f | weight: %5zu\n", (void *) this, mean, weight);
		}

		friend bool operator<(Centroid const &a, Centroid const &b){
			return a.mean < b.mean;
		}
	};

	uint64_t	size		= 0;
	uint64_t	weight		= 0;

	Centroid	centroids[1];



	constexpr auto begin(){
		return std::begin(centroids);
	}

	constexpr auto end(){
		return std::begin(centroids) + size;
	}

	constexpr auto begin() const{
		return std::begin(centroids);
	}

	constexpr auto end() const{
		return std::begin(centroids) + size;
	}



	void add(double value, uint64_t weight){
		using namespace t_digest_impl_;

		insertIntoSortedRange(begin(), end(), Centroid{ value, weight });

		++size;

		this->weight += weight;
	}
};



template<size_t N>
class RawTDigest{
	static_assert(N >= 2, "if N is less 2, there will be nothing to compress to.");

	double		delta_;

public:
	constexpr RawTDigest(double delta) : delta_(delta){}

	enum class Compression{
		NONE		,
		STANDARD	,
		AGGRESSIVE
	};

	constexpr static size_t capacity(){
		return N;
	}

	static bool check(const CentroidData *cd){
		return cd && cd->size <= N;
	}

	template<Compression C = Compression::AGGRESSIVE>
	void add(CentroidData *cd, double value, uint64_t weight = 1) const{
		assert(weight > 0);

		if (!check(cd))
			return;

		if (cd->size < capacity())
			return cd->add(value, weight);

		if constexpr(C == Compression::NONE)
			return;

		if constexpr(C == Compression::STANDARD)
			compressNormal_(cd);

		if constexpr(C == Compression::AGGRESSIVE)
			compressAggressive_(cd);

		if (cd->size < capacity())
			return cd->add(value, weight);

		// drop the value
		// should be unreachible if Aggressive,
	}

	void compress(CentroidData *cd) const{
		if (!check(cd))
			return;

		return compressNormal_(cd);
	}

	static void print(const CentroidData *cd){
		if (!check(cd))
			return;

		printf("Centroids, size %zu, capacity %zu\n", cd->size, capacity());

		for(auto it = cd->begin(); it != cd->end(); ++it)
			it->print();
	}

	static double percentile_50(CentroidData *cd){
		return percentile(cd, 0.50);
	}

	static double percentile_95(CentroidData *cd){
		return percentile(cd, 0.95);
	}

	static double percentile(CentroidData *cd, double const p){
		assert(p >= 0.00 && p <= 1.00);

		if (!check(cd))
			return 0;

		if (cd->size == 0)
			return 0;

		double const targetRank = p * static_cast<double>(cd->weight);
		double cumulativeWeight = 0;

		for(auto it = cd->begin(); it != cd->end(); ++it){
			cumulativeWeight += static_cast<double>(it->weight);

			if (cumulativeWeight > targetRank)
				return it->mean;
		}

		return std::prev(cd->end())->mean;
	}

private:
	void compressNormal_(CentroidData *cd) const{
		if (cd->size < 2)
			return;

		return compressCentroids_<1>(cd, delta_);
	}

	void compressAggressive_(CentroidData *cd) const{
		if (cd->size < 2)
			return;

		auto const distance = findMinDistance_(cd);

		if (distance > delta_)
			return compressCentroids_<0>(cd, distance);
		else
			return compressCentroids_<1>(cd, delta_);
	}

	static double findMinDistance_(const CentroidData *cd){
		double minDistance = std::numeric_limits<double>::max();

		for(auto it = cd->begin(); it != cd->end(); ++it){
			auto const distance = std::abs(it->mean - std::next(it)->mean);

			if (distance < minDistance)
				minDistance = distance;
		}

		return minDistance;
	}

	template<bool UseWeight>
	void compressCentroids_(CentroidData *cd, double delta) const{
		assert(cd->size > 1);

		if constexpr(!UseWeight){
		//	printf("compressAggressive with delta: %10.5f\n", delta);
		}

		size_t newSize = 0;
		auto   current = cd->centroids[0];

		auto const _ = [](double weight) -> double{
			if constexpr(UseWeight)
				return weight;
			else
				return 1.0;
		};

		for (size_t i = 1; i < cd->size; ++i){
			auto const distance = std::abs(cd->centroids[i].mean - current.mean);
			auto const weight_u = current.weight + cd->centroids[i].weight;
			auto const weight   = static_cast<double>(weight_u);

			if (_(weight) * distance <= delta) {
				current.mean   = (current.getWeightedMean() + cd->centroids[i].getWeightedMean()) / weight;
				current.weight = weight_u;
			}else{
				cd->centroids[newSize++] = current;
				current = cd->centroids[i];
			}
		}

		cd->centroids[newSize++] = current;

		cd->size		= newSize;
	}
};



namespace{
	using MyRawTDigest = RawTDigest<5>;

	constexpr double DELTA = 0.05;
	MyRawTDigest td{ DELTA };

	using C = MyRawTDigest::Compression;

	auto get(CentroidData *cd){
		td.add<C::NONE>(cd, 25.0);

		td.add<C::NONE>(cd, 1.50);
		td.add<C::NONE>(cd, 1.51);

		td.add<C::NONE>(cd, 1.46);
		td.add<C::NONE>(cd, 1.47);
	}

	auto getBad(CentroidData *cd){
		td.add<C::NONE>(cd, 31);
		td.add<C::NONE>(cd, 41);
		td.add<C::NONE>(cd, 51);

		td.add<C::NONE>(cd, 10);
		td.add<C::NONE>(cd, 20);
	}

	auto getExtreme(CentroidData *cd){
		td.add<C::NONE>(cd, 1.5);
		td.add<C::NONE>(cd, 1.5);
		td.add<C::NONE>(cd, 1.5);
		td.add<C::NONE>(cd, 1.5);
		td.add<C::NONE>(cd, 1.5);
	}

} // anonymous namespace

int main(){
	const auto memsize = sizeof(CentroidData) + sizeof(CentroidData::Centroid) * (MyRawTDigest::capacity() - 1);

	CentroidData *cd = reinterpret_cast<CentroidData *>(malloc(memsize));

	printf("None...\n");
	get(cd);
	td.add<C::NONE>(cd, 1.52);
	td.print(cd);
	printf("%10.6f\n", td.percentile_50(cd));



	printf("Standard...\n");
	get(cd);
	td.add<C::STANDARD>(cd, 1.52);
	td.print(cd);
	printf("%10.6f\n", td.percentile_50(cd));



	printf("Aggresive...\n");
	get(cd);
	td.add<C::AGGRESSIVE>(cd, 1.52);
	td.print(cd);
	printf("%10.6f\n", td.percentile_50(cd));



	printf("Bad input + standard...\n");
	getBad(cd);
	td.add<C::STANDARD>(cd, 1.52);
	td.print(cd);
	printf("%10.6f\n", td.percentile_50(cd));


	printf("Bad input + aggresive...\n");
	getBad(cd);
	td.add<C::AGGRESSIVE>(cd, 1.52);
	td.print(cd);
	printf("%10.6f\n", td.percentile_50(cd));



	printf("Extreme input...\n");
	getExtreme(cd);
	td.compress(cd);
	td.print(cd);
	printf("%10.6f\n", td.percentile_50(cd));



	free(cd);
}


