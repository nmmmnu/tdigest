#include <cstdint>
#include <cassert>
#include <cstring>
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



struct Centroid{
	double   mean;
	uint64_t weight;

	// constexpr static Centroid create(double m = 0, uint64_t w = 0){
	// 	return Centroid{ m, w };
	// }

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

static_assert(std::is_trivial_v<Centroid>);



class RawTDigest{
	size_t	capacity_;
	double	delta_;

public:
	constexpr RawTDigest(size_t capacity, double delta) : capacity_(capacity), delta_(delta){
		assert(capacity_ >= 2);
	}

	enum class Compression{
		NONE		,
		STANDARD	,
		AGGRESSIVE
	};

	constexpr size_t capacity() const{
		return capacity_;
	}

	static void clearFast(Centroid *cd){
		cd[0].clear();
	}

	constexpr size_t bytes() const{
		return capacity_ * sizeof(Centroid);
	}

	void print(const Centroid *cd) const{
		printf("Centroids, capacity %zu\n", capacity());

		for(size_t i = 0; i < capacity(); ++i){
			auto const &x = cd[i];
			if (!x)
				break;

			x.print();
		}
	}

public:
	void clear(Centroid *cd) const{
		memset(cd, 0, bytes());
	}

	void load(Centroid *cd, const void *src) const{
		memcpy(cd, src, bytes());
	}

	void store(const Centroid *cd, void *dest) const{
		memcpy(dest, cd, bytes());
	}

private:
	size_t getSize_(const Centroid *cd) const{
		size_t size = 0;

		for(size_t i = 0; i < capacity(); ++i){
			auto const &x = cd[i];
			if (!x)
				break;

			++size;
		}

		return size;
	}

	std::pair<uint64_t, size_t> getWeightAndSize_(const Centroid *cd) const{
		std::pair<uint64_t, size_t> r{ 0, 0 };

		for(size_t i = 0; i < capacity(); ++i){
			auto const &x = cd[i];
			if (!x)
				break;

			r.first += x.weight;
			++r.second;
		}

		return r;
	}

	double percentile_(const Centroid *cd, size_t size, uint64_t weight, double const p) const{
		if (size == 0)
			return 0;

		double const targetRank = p * static_cast<double>(weight);
		double cumulativeWeight = 0;

		for(size_t i = 0; i < size - 1; ++i){
			cumulativeWeight += static_cast<double>(cd[i].weight);
			if (cumulativeWeight >= targetRank)
				return cd[i].mean;
		}

		return cd[size - 1].mean;
	}

public:
	template<Compression C = Compression::AGGRESSIVE>
	void add(Centroid *cd, double value, uint64_t weight = 1) const{
		assert(weight > 0);

		auto size = getSize_(cd);

		auto insert = [&](){
			using namespace t_digest_impl_;

			insertIntoSortedRange(cd, cd + size, Centroid{ value, weight} );

			if (++size < capacity())
				cd[size].clear();
		};

		if (size < capacity())
			return insert();

		if constexpr(C == Compression::NONE)
			return;

		if constexpr(C == Compression::STANDARD)
			size = compressNormal_(cd, size);

		if constexpr(C == Compression::AGGRESSIVE)
			size = compressAggressive_(cd, size);

		if (size < capacity())
			return insert();

		// drop the value
		// should be unreachible if Aggressive,
	}

public:
	size_t compress(Centroid *cd) const{
		size_t const size = getSize_(cd);

		return compressNormal_(cd, size);
	}

public:
	double percentile_50(const Centroid *cd) const{
		return percentile(cd, 0.50);
	}

	double percentile_95(const Centroid *cd) const{
		return percentile(cd, 0.95);
	}

	double percentile(const Centroid *cd, double const p) const{
		assert(p >= 0.00 && p <= 1.00);

		auto [weight, size] = getWeightAndSize_(cd);

		return percentile_(cd, size, weight, p);
	}

	template<typename IT, typename OutIT>
	void percentile(const Centroid *cd, IT first, IT last, OutIT out) const{
		auto [weight, size] = getWeightAndSize_(cd);

		auto f = [&](double p){
			assert(p >= 0.00 && p <= 1.00);
			return percentile_(cd, size, weight, p);
		};

		std::transform(first, last, out, f);
	}

private:
	size_t compressNormal_(Centroid *cd, size_t size) const{
		if (size < 2)
			return size;

		return compressCentroids_<1>(cd, size, delta_);
	}

	size_t compressAggressive_(Centroid *cd, size_t size) const{
		if (size < 2)
			return size;

		auto const distance = findMinDistance__(cd, size);

		if (distance > delta_)
			return compressCentroids_<0>(cd, size, distance);
		else
			return compressCentroids_<1>(cd, size, delta_);
	}

	template<bool UseWeight>
	size_t compressCentroids_(Centroid *cd, size_t size, double delta) const{
		assert(size > 1);

		size_t newSize = 0;
		auto   current = cd[0];

		auto const _ = [](double weight) -> double{
			if constexpr(UseWeight)
				return weight;
			else
				return 1.0;
		};

		for (size_t i = 1; i < size; ++i){
			auto const distance = std::abs(cd[i].mean - current.mean);
			auto const weight_u = current.weight + cd[i].weight;
			auto const weight   = static_cast<double>(weight_u);

			if (_(weight) * distance <= delta) {
				current.mean   = (current.getWeightedMean() + cd[i].getWeightedMean()) / weight;
				current.weight = weight_u;
			}else{
				cd[newSize++] = current;
				current = cd[i];
			}
		}

		cd[newSize++] = current;

		if (newSize < capacity())
			cd[newSize].clear();

		return newSize;
	}

	static double findMinDistance__(const Centroid *cd, size_t const size){
		assert(size > 1);

		double minDistance = std::numeric_limits<double>::max();

		for(auto it = cd; it != cd + size - 1; ++it){
			auto const distance = std::abs(it->mean - std::next(it)->mean);

			if (distance < minDistance)
				minDistance = distance;
		}

		return minDistance;
	}
};



namespace{
	constexpr size_t SIZE  = 5;
	constexpr double DELTA = 0.05;

	constexpr auto CM = RawTDigest::Compression::NONE;

	auto get(Centroid *cd){
		RawTDigest td{ SIZE, DELTA };

		td.clearFast(cd);

		td.add<CM>(cd, 25.0);
		td.add<CM>(cd, 1.50);
		td.add<CM>(cd, 1.51);
		td.add<CM>(cd, 1.46);
		td.add<CM>(cd, 1.47);
	}

	auto getBad(Centroid *cd){
		RawTDigest td{ SIZE, DELTA };

		td.clearFast(cd);

		td.add<CM>(cd, 31);
		td.add<CM>(cd, 41);
		td.add<CM>(cd, 51);
		td.add<CM>(cd, 10);
		td.add<CM>(cd, 20);
	}

	auto getExtreme(Centroid *cd){
		RawTDigest td{ SIZE, DELTA };

		td.clearFast(cd);

		td.add<CM>(cd, 1.5);
		td.add<CM>(cd, 1.5);
		td.add<CM>(cd, 1.5);
		td.add<CM>(cd, 1.5);
		td.add<CM>(cd, 1.5);
	}

} // anonymous namespace

int main(){
	RawTDigest td{ SIZE, DELTA };

	Centroid *cd = reinterpret_cast<Centroid *>(malloc(td.bytes()));

	using C = RawTDigest::Compression;

	constexpr double pp[3]{ 0.05, 0.50, 0.95 };
	double oo[3];



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

	td.percentile(cd, std::cbegin(pp), std::cend(pp), std::begin(oo));
	for(auto const &x : oo)
		printf("-> %10.6f\n", x);



	printf("Extreme input...\n");
	getExtreme(cd);
	td.compress(cd);
	td.print(cd);
	printf("%10.6f\n", td.percentile_50(cd));

	td.percentile(cd, std::cbegin(pp), std::cend(pp), std::begin(oo));
	for(auto const &x : oo)
		printf("-> %10.6f\n", x);



	free(cd);
}


